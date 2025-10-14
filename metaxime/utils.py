from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi
from rdkit.Chem.inchi import MolToInchiKey

from typing import Dict, Tuple, Any, Optional, Iterable, Literal, Set, Union, List
from typing import Callable, Dict, Any, Mapping, Optional
from typing import Dict, Tuple, Optional

import logging
import os
import tarfile
import pandas as pd
import tempfile

from cobra import Model, Reaction, Metabolite


##### Merge

import logging
from cobra import Model, Reaction, Metabolite
from typing import Dict

def merge_models(
    source_model: Model,
    input_target_model: Model,
) -> Model:
    """Merge a COBRApy model into another by matching metabolites via annotation overlap.

    For each metabolite in the source model, this function tries to find a match in the
    target model based on overlapping annotations. If a match is found, the metabolite
    in the source model is mapped to the corresponding target metabolite. Reactions
    from the source model are then copied and added to the target model using this mapping.

    Args:
        source_model (Model): The model whose reactions should be merged into the target.
        target_model (Model): The model to which matching reactions are added.

    Returns:
        Dict[str, str]: Mapping of source metabolite IDs to target metabolite IDs.
    """
    def _annotations_overlap(a: Dict[str, Any], b: Dict[str, Any]) -> bool:
        """
        Return True if there is at least one overlapping value for any key present in both dicts.
        Values may be strings or iterables of strings.
        """
        def _to_norm_set(v: Any) -> Set[str]:
            """Normalize a value (str or iterable) to a lowercase, stripped set of strings."""
            if v is None:
                return set()
            if isinstance(v, (list, tuple, set)):
                items: Iterable[Any] = v
            else:
                items = [v]
            return {str(x).strip().lower() for x in items if x is not None}
        for key in set(a) & set(b):
            if _to_norm_set(a[key]) & _to_norm_set(b[key]):
                return True
        return False

    target_model = input_target_model.copy()
    gen_ori_convert_metabolites: Dict[str, str] = {}
    # Match metabolites based on annotation overlap
    for gen_m in source_model.metabolites:
        for ori_m in target_model.metabolites:
            if _annotations_overlap(gen_m.annotation, ori_m.annotation):
                logging.debug(f"{gen_m.id} matches {ori_m.id}")
                gen_ori_convert_metabolites[gen_m.id] = ori_m.id
                break
    # Copy reactions with mapped metabolites
    new_reactions = []
    for r in source_model.reactions:
        reaction = Reaction(r.id)
        reaction.name = r.name
        reaction.lower_bound = r.lower_bound
        reaction.upper_bound = r.upper_bound
        reaction.annotation = r.annotation
        reac_meta_dict = {}
        for met, coeff in r.metabolites.items():
            mapped_id = gen_ori_convert_metabolites.get(met.id, met.id)
            try:
                meta = target_model.metabolites.get_by_id(mapped_id)
            except KeyError:
                meta = source_model.metabolites.get_by_id(met.id)
            reac_meta_dict[meta] = coeff
        reaction.add_metabolites(reac_meta_dict)
        new_reactions.append(reaction)
    #add the reactions to the model
    target_model.add_reactions(new_reactions)
    logging.info(f"Added {len(new_reactions)} reactions to {target_model.id or 'target model'}")
    return target_model


#COBRA-K connector

#### convert


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


