class config:
    """configuration constants"""

    # framework config block
    data_dir = "C:/Users/AndrzejZuranski/tmp"
    db = "mysql://"  # tbd

    # Gaussian config block
    theory = "B3LYP"
    light_set = "6-31G*"
    generic_set = "genecp"
    heavy_set = "LANL2DZ"
    max_light_atomic_number = 36
    atoms_per_processor = 6
    max_processors = 20
    ram_per_processor = 2