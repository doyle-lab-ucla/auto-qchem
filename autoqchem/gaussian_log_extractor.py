import re

from autoqchem.descriptor_functions import *

# TODO: develop structure for descriptors and their saving

logger = logging.getLogger(__name__)
float_or_int_regex = "[-+]?[0-9]*\.[0-9]+|[0-9]+"


class gaussian_log_extractor(object):

    def __init__(self, log_file_path):

        with open(log_file_path) as f:
            self.log = f.read()

        # initialize descriptors
        self.descriptors = {}
        self.atom_descriptors = {}

        self.split_parts()  # split parts
        self.get_atom_labels()  # fetch atom labels
        self.get_geometries()  # fetch geometries for each log section
        self.compute_buried_volumes()  # compute buried volumes
        self.get_frequencies_and_moment_vectors()  # fetch vibration table and vectors
        self.get_freq_part_descriptors()  # fetch descriptors from frequency section
        self.get_td_part_descriptors()  # fetch descriptors from TD section

        # concatenate atom_desciptors from various sources
        self.atom_descriptors = pd.concat(self.atom_descriptors.values(), axis=1)
        self.atom_descriptors.index = self.labels

    def get_atom_labels(self):
        """fetch the z-matrix from the log file and save atom labels"""

        # regex logic, fetch part between "Multiplicity =\d\n" and a double line
        # break (empty line may contain spaces)
        z_matrix = re.findall("Multiplicity = \d\n(.*?)\n\s*\n", self.log, re.DOTALL)[0]
        self.labels = [items[0] for items in map(str.split, z_matrix.split('\n'))]

    def split_parts(self):
        """split the log file into parts that correspond to gaussian steps"""

        # regex logic: log parts start with a new line and " # " pattern
        log_parts = re.split("\n\s#\s", self.log)[1:]
        self.parts = {}
        for p in log_parts:
            # regex logic: find first word in the text
            name = re.search("^\w+", p).group(0)
            self.parts[name] = p

    def __process_geom_block(self, geom) -> pd.DataFrame:
        """convert geom text block to dataframe"""

        geom = geom.splitlines()  # split lines
        geom = map(str.strip, geom)  # strip outer spaces
        geom = filter(lambda line: set(line) != {'-'}, geom)  # remove lines that only contain "---"
        geom = map(str.split, geom)  # split each line by space

        # convert to np.array for further manipulation (note: arrays have unique dtype, here it's str)
        geom_arr = np.array(list(geom))
        # create a dataframe
        geom_df = pd.concat([
            pd.DataFrame(geom_arr[:, 1:3].astype(int), columns=['AN', 'Type']),
            pd.DataFrame(geom_arr[:, 3:].astype(float), columns=list('XYZ'))
        ], axis=1)

        # add atom labels
        geom_df.insert(loc=0, column='Atom', value=self.labels)
        return geom_df

    def get_geometries(self):
        """fetch geometries from the log file"""

        logger.info("Extracting geometry.")
        # regex logic: find parts between "Standard orientation.*X Y Z" and "Rotational constants"
        self.geoms = {}
        for part_name, part_text in self.parts.items():
            geoms = re.findall("Standard orientation:.*?X\s+Y\s+Z\n(.*?)\n\s*Rotational constants",
                               part_text, re.DOTALL)
            # fetch all geometry blocks in the text
            geom_dfs = [self.__process_geom_block(geom) for geom in geoms]
            self.geoms[part_name] = geom_dfs[-1] if geom_dfs else None

    def compute_buried_volumes(self, radius=3):
        """get buried volumes for atoms"""

        logger.info(f"Computing buried volumes within radius: {radius} Angstroms.")
        geom = self.geoms['freq']
        self.atom_descriptors['vbur'] = pd.Series(geom.index.map(lambda i: occupied_volume(geom, i, radius)),
                                                  name='VBur')

    def get_frequencies_and_moment_vectors(self):
        """extract the frequency and other features for each mode within a text section"""

        logger.info("Extracting vibrational frequencies and moment vectors.")
        if 'freq' not in self.parts:
            logger.warning("Output file does not have a 'freq' part. Cannot extract frequencies.")
            return

        # regex logic: text between "Harmonic... normal coordinates and Thermochemistry, preceeded by a line of "---"
        freq_part = re.findall("Harmonic frequencies.*normal coordinates:\s*(\n.*?)\n\n\s*-+\n.*Thermochemistry",
                               self.parts['freq'], re.DOTALL)[0]

        # split each section of text with frequencies
        # regex logic, each frequency part ends with a \s\d+\n, note: we do not use DOTALL here!
        freq_sections = re.split("\n.*?\s\d+\n", freq_part)[1:]

        freq_dfs, vector_dfs = [], []
        for freq_section in freq_sections:
            # frequencies
            freqs = re.findall("\n(\s\w+.*?)\n\s+Atom", freq_section, re.DOTALL)[0]
            freqs = [text.split("--") for text in freqs.splitlines()]
            freqs = [[item[0].strip()] + item[1].split() for item in freqs]
            freqs = np.array(freqs).T.tolist()
            freq_dfs.append(pd.DataFrame(freqs[1:], columns=freqs[0]))

            # vectors
            vectors = re.findall("\n(\s+Atom.*)", freq_section, re.DOTALL)[0]
            vectors = [text.split() for text in vectors.splitlines()]
            vector_df = pd.DataFrame(vectors[1:], columns=vectors[0])
            vector_df.drop(["Atom", "AN"], axis=1, inplace=True)
            vector_dfs.append(vector_df)

        # combine into one frame
        frequencies = pd.concat(freq_dfs)
        frequencies['mode_number'] = range(1, len(frequencies) + 1)
        self.modes = frequencies.set_index('mode_number').astype(float)

        vectors = pd.concat(vector_dfs, axis=1)
        vectors.columns = pd.MultiIndex.from_product([list(range(1, len(frequencies) + 1)), ['X', 'Y', 'Z'], ])
        vectors.columns.names = ['mode_number', 'axis']
        vectors.index = self.labels
        self.mode_vectors = vectors.astype(float)

    def get_freq_part_descriptors(self):
        """extract descriptors from frequency part"""

        logger.info("Extracting frequency section descriptors")
        if 'freq' not in self.parts:
            logger.warning("Output file does not have a 'freq' section. Cannot extract descriptors.")
            return

        text = self.parts['freq']

        # single value descriptors
        single_value_desc_list = [
            {"name": "number_of_atoms", "prefix": "NAtoms=\s*", "type": int},
            {"name": "charge", "prefix": "Charge\s=\s*", "type": int},
            {"name": "multiplicity", "prefix": "Multiplicity\s=\s*", "type": int},
            {"name": "dipole", "prefix": "Dipole moment \(field-independent basis, Debye\):.*?Tot=\s*", "type": float},
            {"name": "molar_mass", "prefix": "Molar Mass =\s*", "type": float},
            {"name": "molar_volume", "prefix": "Molar volume =\s*", "type": float},
            {"name": "electronic_spatial_extent", "prefix": "Electronic spatial extent\s+\(au\):\s+<R\*\*2>=\s*",
             "type": float},
            {"name": "E_scf", "prefix": "SCF Done:\s+E.*?=\s*", "type": float},
            {"name": "zero_point_correction", "prefix": "Zero-point correction=\s*", "type": float},
            {"name": "E_thermal_correction", "prefix": "Thermal correction to Energy=\s*", "type": float},
            {"name": "H_thermal_correction", "prefix": "Thermal correction to Enthalpy=\s*", "type": float},
            {"name": "G_thermal_correction", "prefix": "Thermal correction to Gibbs Free Energy=\s*", "type": float},
            {"name": "E_zpe", "prefix": "Sum of electronic and zero-point Energies=\s*", "type": float},
            {"name": "E", "prefix": "Sum of electronic and thermal Energies=\s*", "type": float},
            {"name": "H", "prefix": "Sum of electronic and thermal Enthalpies=\s*", "type": float},
            {"name": "G", "prefix": "Sum of electronic and thermal Free Energies=\s*", "type": float},
        ]

        for desc in single_value_desc_list:
            value = re.search(f"{desc['prefix']}({float_or_int_regex})", text, re.DOTALL).group(1)
            self.descriptors[desc["name"]] = desc['type'](value)

        # stoichiometry
        self.descriptors['stoichiometry'] = re.search("Stoichiometry\s*(\w+)", text).group(1)

        # convergence, regex-logic: last word in each line should be "YES"
        string = re.search("(Maximum Force.*?)\sPredicted change", text, re.DOTALL).group(1)
        self.descriptors['converged'] = all(np.array(re.findall("(\w+)\n", string)) == 'YES')

        # energies, regex-logic: find all floats in energy block, split by occupied, virtual orbitals
        string = re.search("Population.*?SCF density.*?(\sAlph.*?)\n\s*Condensed", text, re.DOTALL).group(1)
        energies = [re.findall(f"({float_or_int_regex})", s_part) for s_part in string.split("Alpha virt.", 1)]
        occupied_energies, unoccupied_energies = [map(float, e) for e in energies]
        homo, lumo = max(occupied_energies), min(unoccupied_energies)
        self.descriptors['homo_energy'] = homo
        self.descriptors['lumo_energy'] = lumo
        self.descriptors['electronegativity'] = -0.5 * (lumo + homo)
        self.descriptors['hardness'] = 0.5 * (lumo - homo)

        # atom_dependent section
        # Mulliken population
        string = re.search("Mulliken charges.*?\n(.*?)\n\s*Sum of Mulliken", text, re.DOTALL).group(1)
        charges = np.array(list(map(str.split, string.splitlines()))[1:])[:, 2]
        mulliken = pd.Series(charges, name='Mulliken_charge')

        # APT charges
        string = re.search("APT charges.*?\n(.*?)\n\s*Sum of APT", text, re.DOTALL).group(1)
        charges = np.array(list(map(str.split, string.splitlines()))[1:])[:, 2]
        apt = pd.Series(charges, name='APT_charge')

        # NPA charges
        string = re.search("Summary of Natural Population Analysis:.*?\n\s-+\n(.*?)\n\s=+\n", text, re.DOTALL).group(1)
        population = np.array(list(map(str.split, string.splitlines())))[:, 2:]
        npa = pd.DataFrame(population, columns=['NPA_charge', 'NPA_core', 'NPA_valence', 'NPA_Rydberg', 'NPA_total'])

        # NMR
        string = re.findall(f"Isotropic\s=\s*({float_or_int_regex})\s*Anisotropy\s=\s*({float_or_int_regex})", text)
        nmr = pd.DataFrame(np.array(string).astype(float), columns=['NMR_shift', 'NMR_anisotropy'])

        self.atom_descriptors['freq'] = pd.concat([mulliken, apt, npa, nmr], axis=1)

    def get_td_part_descriptors(self):
        """extract descriptors from TD part"""

        logger.info("Extracting TD section descriptors")
        if 'TD' not in self.parts:
            logger.warning("Output file does not have a 'TD' section. Cannot extract descriptors.")
            return

        text = self.parts['TD']

        single_value_desc_list = [
            {"name": "ES_root_dipole", "prefix": "Dipole moment \(field-.*?, Debye\):.*?Tot=\s*", "type": float},
            {"name": "ES_root_molar_volume", "prefix": "Molar volume =\s*", "type": float},
            {"name": "ES_root_electronic_spatial_extent",
             "prefix": "Electronic spatial extent\s+\(au\):\s+<R\*\*2>=\s*", "type": float},
        ]

        for desc in single_value_desc_list:
            value = re.search(f"{desc['prefix']}({float_or_int_regex})", text, re.DOTALL).group(1)
            self.descriptors[desc["name"]] = desc['type'](value)

        # excited states
        string = re.findall(f"Excited State.*?({float_or_int_regex})\snm"
                            f".*f=({float_or_int_regex})"
                            f".*<S\*\*2>=({float_or_int_regex})", text)
        self.transitions = pd.DataFrame(np.array(string).astype(float),
                                        columns=['ES_transition', 'ES_osc_strength', 'ES_<S**2>'])

        # atom_dependent section
        # Mulliken population
        string = re.search("Mulliken charges.*?\n(.*?)\n\s*Sum of Mulliken", text, re.DOTALL).group(1)
        charges = np.array(list(map(str.split, string.splitlines()))[1:])[:, 2]
        mulliken = pd.Series(charges, name='ES_root_Mulliken_charge')

        # NPA charges
        string = re.search("Summary of Natural Population Analysis:.*?\n\s-+\n(.*?)\n\s=+\n", text, re.DOTALL).group(1)
        population = np.array(list(map(str.split, string.splitlines())))[:, 2:]
        npa = pd.DataFrame(population, columns=['ES_root_NPA_charge', 'ES_root_NPA_core', 'ES_root_NPA_valence',
                                                'ES_root_NPA_Rydberg', 'ES_root_NPA_total'])

        self.atom_descriptors['TD'] = pd.concat([mulliken, npa], axis=1)
