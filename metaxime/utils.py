from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi
from rdkit.Chem.inchi import MolToInchiKey

from typing import Dict, Tuple, Any, Optional, Iterable, Literal, Set, Union, List
from typing import Callable, Dict, Any, Mapping, Optional

import logging
import os
import tarfile
import pandas as pd
import tempfile

from typing import Dict, Tuple, Optional

def convert_depiction(
    idepic: str,
    itype: Literal["smiles", "inchi"] = "smiles",
    otype: Iterable[str] = ("inchikey",),
) -> Dict[str, str]:
    """Convert a chemical depiction to one or more other formats using RDKit.

    Args:
        idepic: Input depiction string (SMILES or InChI).
        itype: Type of the input depiction, either "smiles" or "inchi".
        otype: Iterable of desired output types. Any of {"smiles", "inchi", "inchikey"}.

    Returns:
        Dict[str, str]: A mapping from requested output types to their string values.

    Raises:
        NotImplementedError: If `itype` or any requested output type is unsupported.
        TypeError: If the input cannot be parsed into an RDKit molecule.
        ValueError: If `otype` is empty.
    """
    logging.debug(f"input: {idepic}")
    logging.debug(f"itype: {itype}")
    otypes: Set[str] = set(otype)
    if not otypes:
        raise ValueError("`otype` must contain at least one output type.")
    # Import
    if itype == "smiles":
        rdmol = MolFromSmiles(idepic, sanitize=True)
    elif itype == "inchi":
        rdmol = MolFromInchi(idepic, sanitize=True)
    else:
        raise NotImplementedError(f'"{itype}" is not a valid input type (use "smiles" or "inchi").')
    if rdmol is None:
        raise TypeError(f'Failed to parse depiction "{idepic}" of type "{itype}".')
    logging.debug("Sanitized the input molecule")
    # Export
    supported = {"smiles", "inchi", "inchikey"}
    unknown = otypes - supported
    if unknown:
        raise NotImplementedError(f"Unsupported output type(s): {sorted(unknown)}. "
                                  f"Supported: {sorted(supported)}")
    out: Dict[str, str] = {}
    if "smiles" in otypes:
        # canonical SMILES by default
        out["smiles"] = MolToSmiles(rdmol)
    if "inchi" in otypes:
        out["inchi"] = MolToInchi(rdmol)
    if "inchikey" in otypes:
        out["inchikey"] = MolToInchiKey(rdmol)
    logging.debug("Exported the requested output depictions")
    return out


def read_compressed_tsv(file_path: str) -> pd.DataFrame:
    """Load a TSV file of reaction recipes, supporting plain or tar.gz formats.

    Args:
        file_path (str): Path to the .tsv or .tar.gz file containing the TSV.

    Returns:
        pd.DataFrame: Loaded DataFrame from the TSV file.
    """
    # Case 1: file is a regular TSV
    if file_path.endswith(".tsv"):
        return pd.read_csv(file_path, comment="#", sep="\t", header=None)

    # Case 2: file is a .tar.gz containing a TSV
    elif file_path.endswith(".tar.gz"):
        with tempfile.TemporaryDirectory() as tmpdir:
            with tarfile.open(file_path, "r:gz") as tar:
                # find the TSV file inside
                tsv_members = [m for m in tar.getmembers() if m.name.endswith(".tsv")]
                if not tsv_members:
                    raise FileNotFoundError("No .tsv file found in tar.gz archive.")
                member = tsv_members[0]
                tar.extract(member, tmpdir)
                extracted_path = os.path.join(tmpdir, member.name)
                return pd.read_csv(extracted_path, comment="#", sep="\t", header=None)

    else:
        raise ValueError(f"Unsupported file type: {file_path}")
