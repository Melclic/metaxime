

from typing import Any, Dict, Iterable, Set

def annotations_overlap(a: Dict[str, Any], b: Dict[str, Any]) -> bool:
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
    target_model = input_target_model.copy()
    gen_ori_convert_metabolites: Dict[str, str] = {}

    # Match metabolites based on annotation overlap
    for gen_m in source_model.metabolites:
        for ori_m in target_model.metabolites:
            if annotations_overlap(gen_m.annotation, ori_m.annotation):
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

    target_model.add_reactions(new_reactions)
    logging.info(f"Added {len(new_reactions)} reactions to {target_model.id or 'target model'}")
    return target_model


#COBRA-K connector
