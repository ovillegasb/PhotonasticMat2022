
import pandas as pd
import numpy as np


class MoleculeDefintionError(Exception):
    pass


Elements = {  # g/mol
    'H': {'mass': 1.0080, 'num': 1},
    'C': {'mass': 12.011, 'num': 6},
}


def load_structure(file):

    # Extract information from file
    if file.endswith("xyz"):
        dfatoms = _load_xyz(file)
        dfatoms["mass"] = dfatoms["atsb"].apply(lambda at: Elements[at]["mass"])

    else:
        raise MoleculeDefintionError("Molecule format is not recognized")

    return dfatoms


def _load_xyz(file):
    """Read a file xyz."""
    coord = pd.read_csv(
        file,
        sep=r'\s+',
        skiprows=2,
        header=None,
        names=['atsb', 'x', 'y', 'z'],
        dtype={'x': np.float64, 'y': np.float64, 'z': np.float64}
    )

    return coord
