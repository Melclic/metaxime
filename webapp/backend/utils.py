from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi
from rdkit.Chem.inchi import MolToInchiKey

from typing import Dict, Tuple, Any, Optional, Iterable, Literal, Set, Union, List
from typing import Callable, Dict, Any, Mapping, Optional
from typing import Dict, Tuple, Optional

import logging

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
