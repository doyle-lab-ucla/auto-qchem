from config import *

class gaussian_file_generator():
    """generator of gaussian input files"""

    def __init__(self, has_heavy, workflow, directory):
        """configure workflow and basis set to be used"""

        self.has_heavy = has_heavy
        self.workflow = workflow
        self.directory = directory

        if self.has_heavy:
            self.basis_set = config.generic_set
        else:
            self.basis_set = config.light_set

    def generate_file(self, name, resource_block, coords_block,
                      light_elements, heavy_elements, charge, multiplicity):

        if self.has_heavy:
            heavy_block = ""
            heavy_block += "\n" + " ".join(light_elements + ["0"]) + "\n"
            heavy_block += f"{config.light_set}\n****"
            heavy_block += "\n" + " ".join(heavy_elements + ["0"]) + "\n"
            heavy_block += f"{config.heavy_set}\n****\n"
            heavy_block += "\n" + " ".join(heavy_elements + ["0"]) + "\n"
            heavy_block += f"{config.heavy_set}\n"

        output = ""
        output += resource_block
        output += f"%Chk={name}_1.chk\n"

        # geoemtry optimization
        output += f"# {self.workflow.value['opt'](self.basis_set)}\n\n"
        output += f"{name}\n\n"
        output += f"{charge} {multiplicity}\n"
        output += coords_block

        # add heavy block
        if self.has_heavy:
            output += heavy_block

        # frequency calculation
        output += "\n\n\n--Link--1\n"
        output += resource_block
        output += f"%Oldchk={name}_1.chk\n"
        output += f"%Chk={name}_2.chk\n"
        output += f"# {self.workflow.value['freq'](self.basis_set)}\n"

        if self.has_heavy:
            output += heavy_block

        # time_dependent calculation (only for equilibrium workflow)
        if 'td' in self.workflow.value:
            output += "\n\n\n--Link--1\n"
            output += resource_block
            output += f"%Oldchk={name}_2.chk\n"
            output += f"%Chk={name}_3.chk\n"
            output += f"{self.workflow.value['td'](self.basis_set)}\n"

            if self.has_heavy:
                output += heavy_block

        with open(f"{self.directory}/{name}.gjf", "w") as file:
            file.write(output)