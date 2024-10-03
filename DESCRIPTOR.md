This document contains detailed explanations of all descriptor names in Auto-QChem with their units. All calculations are done in Gaussian with Gaussian terminologies. If not specified for excited state, all values represent ground state. 


E: Sum of electronic and thermal Energies (Hartree)

ES_root_dipole: Dipole moment for the first excited state (Debye)

ES_root_electronic_spatial_extent: Electronic spatial extent for the first excited state (a.u.)

ES_root_molar_volume: Molar volume for the first excited state (Bohr^3/mole)

E_scf: 	Energy from Self-Consistent Field method (Hartree)

E_thermal_correction: Thermal correction to Energy (Hartree)

E_zpe: Sum of electronic and zero-point Energies (Hartree)

G: Sum of electronic and thermal Free Energies (Hartree)

G_thermal_correction: Thermal correction to Gibbs Free Energy (Hartree)

H: Sum of electronic and thermal Enthalpies (Hartree)

H_thermal_correction: Thermal correction to Enthalpy (Hartree)

Dipole: Dipole moment (field-independent basis, Debye)

electronegativity: –0.5*(lumo_energy + homo_energy) (Hartree)

electronic_spatial_extent: Electronic spatial extent (a.u.)

hardness: 0.5*(lumo_energy – homo_energy) (Hartree)

homo_energy: HOMO energy (Hartree)

lumo_energy: LUMO energy (Hartree)

molar_mass: Molar Mass (grams/mole)

molar_volume: Molar volume (Bohr^3/mole)

zero_point_correction: Zero-point correction (Hartree)

ES_<S**2>_<number>: Total spin operator for excited state <number>, which corresponds to the same excited state in ES_osc_strength_<number> (unitless)

ES_osc_strength_<number>: The oscillator strength for excited state <number>, ranked from highest to lowest. E.g., ES_osc_strength_0 represent the excited state with the highest oscillator strength out of all excited state calculated (unitless)

ES_transition_<number>: The wavelength for the excited state transition from ground state to excited state <number>, which corresponds to the same excited in ES_osc_strength_<number> (nanometer).

<atom>_APT_charge: <atom>’s Atomic Polar Tensor charge (a.u.)

<atom>_ES_root_Mulliken_charge: <atom>’s Mulliken charge in the first excited state (a.u.)

<atom>_ES_root_NPA_Rydberg: <atom>’s Rydberg orbital charges from NPA for the first excited state (a.u.)

<atom>_ES_root_NPA_charge: <atom>’s natural charge from NPA for the first excited state (a.u.)

<atom>_ES_root_NPA_core: <atom>’s core orbital charges from NPA for the first excited state (a.u.)

<atom>_ES_root_NPA_total: <atom>_ES_root_NPA_core + <atom>_ES_root_NPA_valence + <atom>_ES_root_NPA_Rydberg

<atom>_ES_root_NPA_valence: <atom>’s valence orbital charges from NPA (a.u.)

<atom>_Mulliken_charge: <atom>’s Mulliken charge (a.u.)

<atom>_NMR_anisotropy: <atom>’s NMR shielding anisotropy 

<atom>_NMR_shift: <atom>’s isotropic NMR shielding/chemical shift (ppm)

<atom>_NPA_Rydberg: <atom>’s Rydberg orbital charges from NPA (a.u.)

<atom>_NPA_charge: <atom>’s natural charge from NPA (a.u.)

<atom>_NPA_core: <atom>’s core orbital charges from NPA (a.u.)

<atom>_NPA_total: <atom>_NPA_core + <atom>_NPA_valence + <atom>_NPA_Rydberg

<atom>_NPA_valence: <atom>’s valence orbital charges from NPA (a.u.)	

<atom>_VBur: for a given <atom>, within a 3Å radius, the fraction of the volume that is occupied by other atoms in the molecule (represented by their Van der Waals radii). 
