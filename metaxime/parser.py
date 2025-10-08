from __future__ import annotations

import copy
import logging
import pandas as pd
import numpy as np
import re
import networkx as nx

from typing import Dict, Tuple, Any, Optional, Iterable, Literal, Set, Union, List
from typing import Callable, Dict, Any, Mapping, Optional

from cobra import Model, Reaction, Metabolite

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity

from metaxime.cache_data import RR_Data
from metaxime.utils import convert_depiction

from biopathopt.utils import merge_annot_dicts

class ParserRP2(RR_Data):
    def __init__(
            self, 
            rp2_scope_path: str,
            rp2_cmp_path: str,
            rp2_paths_path: str,
            use_progressbar=False, 
            low_memory_mode=False,
        ):
        """Class that inherits Data used to build a cobra model
        """
        super().__init__(low_memory_mode=low_memory_mode, use_progressbar=use_progressbar)
        self.rp_strc = self._read_rp2cmp(rp2_cmp_path)
        self.rp_scope = self._read_rp2scope(rp2_scope_path)
        self.rp_paths = self._read_rp2paths(rp2_paths_path)
        self.all_paths = self._extract_all_paths(self.rp_paths)
        self.completed_paths = self._process_all_paths(self.all_paths)


    def _best_inchi_match(
        self,
        query_inchi: str,
        inchi_dict: Dict[str, str],
        radius: int = 2,
        n_bits: int = 2048,
    ) -> Tuple[Optional[str], float]:
        """Find the key of the molecule most similar to the given InChI using RDKit Morgan fingerprints.

        Args:
            query_inchi (str): The query molecule InChI.
            inchi_dict (Dict[str, str]): Dictionary with keys as identifiers and values as InChIs.
            radius (int): Morgan fingerprint radius. Defaults to 2.
            threshold (float): Minimum similarity threshold for a confident match. Defaults to 0.7.

        Returns:
            Tuple[Optional[str], float]: (Best matching key, Tanimoto score). Returns (None, 0.0) if no confident match found.
        """
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=n_bits)
        query_mol = Chem.MolFromInchi(query_inchi)
        if query_mol is None:
            logging.error(f"Invalid query InChI: {query_inchi}")
            return None, 0.0

        query_fp = gen.GetFingerprint(query_mol)
        best_key, best_score = None, -1.0

        for key, inchi in inchi_dict.items():
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                logging.warning(f"Invalid InChI for {key}")
                continue
            sim = TanimotoSimilarity(query_fp, gen.GetFingerprint(mol))
            if sim > best_score:
                best_key, best_score = key, sim
        return best_key, best_score


    #### RP2 output readers


    def _read_rp2paths(self, rp2paths_path: str) -> Union[Dict[int, Dict[int, Dict[int, Dict[str, Any]]]], bool]:
        """Read rp2paths pathways CSV and build nested dictionary of reactions.

        This function parses the output of rp2paths (pathways CSV) and reconstructs
        the nested pathway structure. Each pathway is broken down into steps and
        sub-steps, storing rule, reaction, and compound information.

        The output has the following structure:
            rp_paths[path_id][step][sub_step] = {
                'rule_id': str,
                'rule_ori_reac': str,
                'rule_score': float | None,
                'left': Dict[str, int],
                'right': Dict[str, int],
                'path_id': int,
                'step': int,
                'transformation_id': str
            }

        Args:
            rp2paths_pathways (str): Path to the rp2paths pathway output CSV file.

        Returns:
            Union[Dict[int, Dict[int, Dict[int, Dict[str, Any]]]], bool]:
                - Nested dictionary describing pathways if successful.
                - False if parsing fails due to malformed data or missing file.
        """
        rp_paths = {}
        current_path_id = None
        path_step = 0
        try:
            df = pd.read_csv(rp2paths_path)
            df['Transformation ID'] = [i[:-2] for i in df['Unique ID']]
            required_cols = {"Path ID", "Unique ID", "Rule ID", "Left", "Right"}
            if not required_cols.issubset(df.columns):
                logging.error(f"Missing required columns in {rp2paths_pathways}. Found: {list(df.columns)}")
                return False
            for _, row in df.iterrows():
                # Parse path id and step counter
                try:
                    pid = int(row["Path ID"])
                except (TypeError, ValueError):
                    logging.error(f"Cannot convert Path ID to int ({row['Path ID']})")
                    return False
                if pid != current_path_id:
                    path_step = 1
                    current_path_id = pid
                else:
                    path_step += 1
                # Parse reactants/products
                def parse_side(text: str) -> Dict[str, int]:
                    side: Dict[str, int] = {}
                    txt = str(text).replace("'", "").replace("-", "_").replace("+", "")
                    if not txt or txt == "nan":
                        return side
                    for chunk in txt.split(":"):
                        if not chunk:
                            continue
                        parts = chunk.split(".", 1)
                        if len(parts) != 2:
                            logging.warning(f"Malformed side chunk '{chunk}' in '{txt}'")
                            continue
                        sto_str, name = parts[0].strip(), parts[1].strip()
                        try:
                            sto = int(sto_str)
                        except ValueError:
                            logging.error(f"Cannot convert stoichiometry to int ({sto_str})")
                            return {}
                        cid = self.single_depr_mnxm(name)
                        side[cid] = sto
                    return side
                left = parse_side(row["Left"])
                right = parse_side(row["Right"])
                if (pd.notna(row["Left"]) and not left) or (pd.notna(row["Right"]) and not right):
                    return False
                # For each rule_id, look up reactions from retrorules_prop
                for r_id in str(row["Rule ID"]).split(","):
                    if r_id=='nan' or not r_id:
                        logging.warning(f'The following rule id is empty: {r_id}')
                        continue
                    rr_reacts = self.retrorules_prop.get(r_id.strip(), {})
                    if not rr_reacts:
                        logging.warning(f'Cannot recover the following rule: {r_id}') 
                    for react in rr_reacts:
                        #there can be multiple substrates for each reaction
                        for sub in rr_reacts[react]:
                            rp_paths.setdefault(current_path_id, {}).setdefault(path_step, {}).setdefault(r_id, {}).setdefault(react, {})[sub] = {
                                "rule_id": r_id,
                                "rule_mnxr": react,
                                "rule_mnxm": sub,
                                "rule_score": rr_reacts[react][sub].get("Score", 0.0),
                                "right": dict(right),
                                "left": dict(left),
                                "path_id": pid,
                                "step": path_step,
                                "transformation_id": row.get("Transformation ID"),
                            }
        except FileNotFoundError:
            logging.error(f"Cannot find file: {rp2paths_pathways}")
            return False
        except OSError as e:
            logging.error(f"Error reading {rp2paths_pathways}: {e}")
            return False
        return rp_paths


    def _read_rp2scope(self, scope_path: str) -> Dict[str, Dict[str, Any]]:
        """Parse the RetroPath2 scope CSV file into transformation data.

        This function reads an `out_scope.csv` file from RetroPath2 and aggregates
        information by `Transformation ID`, computing the average score and collecting
        associated reaction SMILES strings.

        Args:
            scope_path (str): Path to the RetroPath2 scope CSV file.

        Returns:
            Dict[str, Dict[str, Any]]: A dictionary where each key is a Transformation ID,
            and the value is another dictionary containing:
                - "score": float, average score of the transformation.
                - "reaction_smiles": str, reaction SMILES (if unique).
        """
        df = pd.read_csv(scope_path)
        # Compute mean scores per transformation
        transf_id__score = (
            df[['Transformation ID', 'Score']]
            .groupby('Transformation ID')['Score']
            .apply(set)
            .to_dict()
        )
        transf_id__score = {tid: np.mean(list(scores)) for tid, scores in transf_id__score.items()}
        # Collect unique reaction SMILES per transformation
        transf_id__reaction_smiles = (
            df[['Transformation ID', 'Reaction SMILES']]
            .groupby('Transformation ID')['Reaction SMILES']
            .apply(set)
            .to_dict()
        )
        for tid, smiles_set in transf_id__reaction_smiles.items():
            if len(smiles_set) == 1:
                transf_id__reaction_smiles[tid] = list(smiles_set)[0]
            else:
                logging.warning(f"Multiple reaction SMILES found for {tid}: {smiles_set}")
        # Combine into structured dictionary
        scope = {
            tid: {
                "score": transf_id__score.get(tid, np.nan),
                "reaction_smiles": transf_id__reaction_smiles.get(tid, None)
            }
            for tid in transf_id__score.keys()
        }
        return scope



    def _read_rp2cmp(self, path: str) -> Dict[str, Dict[str, Any]]:
        """Parse an RP2Paths compounds and enrich with InChI/InChIKey.

        The input file is expected to be tab-separated with a header row.
        This function skips the header and treats the first two columns as:
          - column 0: compound identifier (CID)
          - column 1: SMILES string

        For each row, it attempts to fetch `inchi` and `inchikey` via
        `self.queryCIDstrc(cid)`. If missing, it tries to derive them by
        converting the SMILES using `self._convert_depiction`.

        Args:
            path: Path to the compounds TSV file.

        Returns:
            Dict[str, Dict[str, Any]]: Mapping CID -> {'smiles': str, 'inchi': str?, 'inchikey': str?}

        Raises:
            RuntimeError: If the file cannot be read.
        """
        rp_strc: Dict[str, Dict[str, Any]] = {}
        try:
            # File has a header line; skip it and use positional columns like the original.
            df = pd.read_csv(path)
            for _, row in df.iterrows():
                cid = str(row.iloc[0])
                smiles = str(row.iloc[1])
                rp_strc[cid] = {"SMILES": smiles}
                # InChI
                try:
                    res_conv = convert_depiction(idepic=smiles, itype="smiles", otype={"inchi"})
                    rp_strc[cid]["InChI"] = res_conv["inchi"]
                except NotImplementedError:
                    logging.warning(f"Could not convert SMILES to InChI for CID={cid}: {smiles!r}")
                # InChIKey
                try:
                    res_conv = convert_depiction(idepic=smiles, itype="smiles", otype={"inchikey"})
                    rp_strc[cid]["InChIKey"] = res_conv["inchikey"]
                except NotImplementedError:
                    logging.warning(f"Could not convert SMILES to InChIKey for CID={cid}: {smiles!r}")
                # get the xref
                xref = None
                mnxm = None
                if 'MNXM' in cid:
                    xref, _ = self.mnxm_xref(self.single_depr_mnxm(cid))
                    xref['xref'] = merge_annot_dicts(
                        xref['xref'],
                        rp_strc[cid],
                    )
                    rp_strc[cid] = xref
                else:
                    try:
                        mnxm = self.inchikey_mnxm[rp_strc[cid]["InChIKey"]]
                    except KeyError:
                        try:
                            mnxm = self.inchikey2_mnxm['-'.join(rp_strc[cid]["InChIKey"].split('-')[:2])]
                        except KeyError:
                            pass
                    if mnxm:
                        xref, _ = self.mnxm_xref(mnxm)
                #print(f"xref: {xref['xref']}")
                if xref:
                    xref['xref'] = merge_annot_dicts(
                        xref['xref'],
                        rp_strc[cid],
                    )
                    rp_strc[cid] = xref
        except (FileNotFoundError, OSError, pd.errors.EmptyDataError) as e:
            logging.error(f"Could not read the compounds file ({path}): {e}")
            raise RuntimeError from e
        return rp_strc


    #### Convert Monocomponent Reactions

    def _extract_all_paths(
            self,
            rp_paths: Dict[int, Dict[int, Dict[str, Dict[str, Dict[str, Any]]]]]
        ) -> Dict[int, List[Dict[int, Dict[str, str]]]]:
        """Extract all paths from nested rp_paths by building a directed graph.

        This function creates a directed graph for each pathway, where:
          step → rule → reaction → substrate → next_step

        It then enumerates all simple paths between the first and last step,
        returning a structured representation of each.

        Args:
            rp_paths: Nested dictionary structure:
                rp_paths[path_id][step][rule_id][reaction_id][substrate_id] -> {...}

        Returns:
            Dict mapping path_id -> list of paths, where each path is a dictionary:
                { step: {"rule": str, "reaction": str, "substrate": str}, ... }
        """
        to_ret: Dict[int, List[Dict[int, Dict[str, str]]]] = {}

        for path_id, steps_dict in rp_paths.items():
            G = nx.DiGraph()
            path_steps = sorted(steps_dict.keys())

            if not path_steps:
                to_ret[path_id] = []
                continue

            min_step, max_step = min(path_steps), max(path_steps)

            # Build directed graph
            for step, rule_dict in steps_dict.items():
                G.add_node(step)
                for rule_id, react_dict in rule_dict.items():
                    G.add_node(rule_id)
                    G.add_edge(step, rule_id)

                    for react_id, sub_dict in react_dict.items():
                        G.add_node(react_id)
                        G.add_edge(rule_id, react_id)

                        for sub_id in sub_dict.keys():
                            G.add_node(sub_id)
                            G.add_edge(react_id, sub_id)
                            G.add_edge(sub_id, step + 1)

            # Enumerate all simple paths
            all_paths: List[Dict[int, Dict[str, str]]] = []
            for node_path in nx.all_simple_paths(G, source=min_step, target=max_step + 1):
                current_step = None
                path_data: Dict[int, Dict[str, str]] = {}

                for node in node_path[:-1]:  # skip the last sink node
                    if isinstance(node, int):
                        current_step = node
                        path_data[current_step] = {}
                    elif isinstance(node, str):
                        # Identify type by prefix (no regex)
                        if node.startswith("RR"):
                            path_data[current_step]["rule"] = node
                        elif node.startswith("MNXR"):
                            path_data[current_step]["reaction"] = node
                        elif node.startswith("MNXM"):
                            path_data[current_step]["substrate"] = node

                all_paths.append(path_data)

            to_ret[path_id] = all_paths

        return to_ret


    def _complete_monocomponent_reaction(
        self,
        input_subpath: Dict[str, Any],
        rp_path: Dict,
        match_threshold: float = 0.8,
    ) -> None:
        """
        Complete a single monocomponent RP2 step by reconciling predicted species with
        the original reaction recipe and adding any missing reactants/products.

        Args:
            subpath: Dict with keys 'rule', 'reaction', 'substrate'.
            match_threshold: Minimum similarity to accept the target match without fallback.

        Returns:
            The updated subpath (mutated copy in-place).
        """
        subpath = copy.deepcopy(input_subpath)
        rp_rule = subpath["rule"]
        rp_rule_reac = subpath["reaction"]
        rp_rule_substrate = subpath["substrate"]

        # 1) Recover original recipe reactants/products (prefer rr_recipes, fall back to cache)
        try:
            ori_reactants = self.rr_recipes[rp_rule_reac]["main_reactants"] | self.rr_recipes[rp_rule_reac]["secondary_reactants"]
            ori_products  = self.rr_recipes[rp_rule_reac]["main_products"]  | self.rr_recipes[rp_rule_reac]["secondary_products"]
        except (KeyError, TypeError):
            try:
                base_id = self.single_depr_mnxr(rp_rule_reac)
                ori_reactants = (
                    self.mnxr_prop[base_id]["main_reactants"]
                    | self.mnxr_prop[rp_rule_reac]["secondary_reactants"]
                )
                ori_products = (
                    self.mnxr_prop[base_id]["main_products"]
                    | self.mnxr_prop[rp_rule_reac]["secondary_products"]
                )
            except (KeyError, TypeError):
                raise KeyError(f"Cannot find original reactants/products for {rp_rule_reac}")

        logging.debug(f"\t\tori_reactants: {ori_reactants}")
        logging.debug(f"\t\tori_products: {ori_products}")

        # 2) Build InChI pool for all species mentioned in the original recipe
        ori_strc_inchi: Dict[str, str] = {}
        left_out: list[str] = []
        for mid in (ori_reactants | ori_products):
            inchi = self.rp_strc.get(mid, {}).get("InChI")
            if not inchi:
                try:
                    tmp = self.mnxm_prop[mid]["InChI"]
                    if not pd.isna(tmp):
                        inchi = tmp
                except KeyError:
                    inchi = None
            if inchi:
                ori_strc_inchi[mid] = inchi
            else:
                logging.warning(f"Cannot recover the InChI for {mid}")
                left_out.append(mid)

        # 3) Identify the predicted major product on the RP2 right-hand side
        rp_right = rp_path[rp_rule][rp_rule_reac][rp_rule_substrate]["right"]
        if len(rp_right) != 1:
            raise KeyError(f"Multiple elements on RP2 right-hand side: {rp_right}")

        rp_predict_strc_id = next(iter(rp_right.keys()))
        target_best_mnxm, score = self._best_inchi_match(self.rp_strc[rp_predict_strc_id]["InChI"], ori_strc_inchi)

        if score <= match_threshold and left_out:
            if len(left_out) == 1:
                logging.warning(f"Using unknown as main RHS metabolite: {left_out[0]}")
                target_best_mnxm = left_out[0]
            else:
                raise KeyError("Cannot confidently identify the RHS metabolite")

        logging.debug(f"\t\t{rp_predict_strc_id} -> {target_best_mnxm}")

        # 4) Determine orientation by checking where the target metabolite belongs
        if target_best_mnxm in ori_products:
            reactant_dir, product_dir = "left", "right"
        elif target_best_mnxm in ori_reactants:
            reactant_dir, product_dir = "right", "left"
        else:
            raise KeyError("Cannot recognize reactant/product orientation")

        # 5) Fetch predicted sides for this step
        rp_reactants = rp_path[rp_rule][rp_rule_reac][rp_rule_substrate][reactant_dir]
        rp_products  = rp_path[rp_rule][rp_rule_reac][rp_rule_substrate][product_dir]
        logging.debug(f"\t\trp_reactants: {rp_reactants}")
        logging.debug(f"\t\trp_products: {rp_products}")

        # 6) Map predicted species to original IDs using InChI similarity (reactants)
        reactants_rp2ori: Dict[str, str] = {}
        for mid in rp_reactants:
            if mid not in ori_reactants:
                inchi = self.rp_strc.get(mid, {}).get("InChI")
                if inchi:
                    best_mnxm, _ = self._best_inchi_match(inchi, ori_strc_inchi)
                    reactants_rp2ori[mid] = best_mnxm

        # 7) Map predicted species to original IDs (products)
        products_rp2ori: Dict[str, str] = {rp_predict_strc_id: target_best_mnxm}
        for mid in rp_products:
            if mid not in ori_products and mid != rp_predict_strc_id:
                inchi = self.rp_strc.get(mid, {}).get("InChI")
                if inchi:
                    best_mnxm, _ = self._best_inchi_match(inchi, ori_strc_inchi)
                    products_rp2ori[mid] = best_mnxm

        logging.debug(f"\t\tConverted rp_reactants: {[reactants_rp2ori.get(i, i) for i in rp_reactants.keys()]}")
        logging.debug(f"\t\tConverted rp_products: {[products_rp2ori.get(i, i) for i in rp_products.keys()]}")

        # 8) Determine which species are missing (to add from the recipe)
        to_add_reactants = set(ori_reactants.keys()) - {reactants_rp2ori.get(i, i) for i in rp_reactants.keys()}
        to_add_products  = set(ori_products.keys())  - {products_rp2ori.get(i, i)  for i in rp_products.keys()}
        logging.debug(f"\t\tto_add_reactants: {to_add_reactants}")
        logging.debug(f"\t\tto_add_products: {to_add_products}")

        # 9) Populate the completed subpath
        subpath["reactants"] = copy.deepcopy(rp_reactants)
        subpath["products"]  = copy.deepcopy(rp_products)
        subpath["transformation_id"] = rp_path[rp_rule][rp_rule_reac][rp_rule_substrate]["transformation_id"]

        # Add any missing recipe species with their stoichiometries
        for mid in to_add_products:
            subpath["products"][mid] = ori_products[mid]
        for mid in to_add_reactants:
            subpath["reactants"][mid] = ori_reactants[mid]

        logging.debug(f"\t\tfull reactants: {subpath['reactants']}")
        logging.debug(f"\t\tfull products: {subpath['products']}")
        logging.debug("\t -------")

        # TODO: if you need to reconcile stoichiometries beyond simple union, do it here.

        return subpath


    def _process_all_paths(
        self,
        all_paths: Dict[int, Any],
    ) -> None:
        """Iterate through all pathway structures and process each subpath step.

        This function iterates over all pathway entries, invoking 
        _step_complete_monocomponent_reaction on each step of
        every subpath. It handles missing keys gracefully and logs warnings.

        Args:
            step_function: Function to apply to each path step. It must accept a dict.

        Returns:
            None
        """
        to_ret = {}
        for rp_path_num in all_paths:
            logging.debug(f'------ {rp_path_num} -------')
            to_ret[rp_path_num] = []
            for rp_subpath in all_paths[rp_path_num]:
                to_overwrite = {}
                is_valid = True
                for path_step in rp_subpath:
                    try:
                        logging.debug(f"\t-> {path_step}")
                        to_overwrite[path_step] = self._complete_monocomponent_reaction(
                                rp_subpath[path_step],
                                rp_path=self.rp_paths[rp_path_num][path_step],
                            )
                    except KeyError as e:
                        logging.warning(f"Skipping {rp_path_num} because of the following KeyError: {e}")
                        is_valid = False
                        break
                if is_valid:
                    to_ret[rp_path_num].append(to_overwrite)
        return to_ret


    def return_rp2_models(
        self,
        compartment_id: str = "c",
        reaction_lower_bound: float = 0.0,
        reaction_upper_bound: float = 1000.0,
    ) -> Dict[int, Dict[int, Model]]:
        """Build one COBRA model per (path_num, subpath) from RetroPath-like data.

        Args:
            compartment_id: COBRA compartment ID to assign to created metabolites.
            reaction_lower_bound: Lower bound applied to every reaction.
            reaction_upper_bound: Upper bound applied to every reaction.

        Returns:
            Dict mapping path_num -> {subpath_index -> cobra.Model}.
        """
        to_ret = {}

        for rp_path_num in self.completed_paths:
            logging.debug(f'------ {rp_path_num} -------')
            to_ret[rp_path_num] = {}
            for idx, rp_subpath in enumerate(self.completed_paths[rp_path_num]):
                model = Model(f'rp2_{rp_path_num}_{idx}')
                model_meta = {}
                model_reac = []
                for path_step in rp_subpath: #step in that enumarated path
                    rp_rule = rp_subpath[path_step]['rule']
                    rp_reactants = rp_subpath[path_step]['reactants']
                    rp_products = rp_subpath[path_step]['products']
                    rp_rule_trans_id = rp_subpath[path_step]['transformation_id']
                    #Reaction
                    reaction = Reaction(rp_rule_trans_id)
                    reaction.name = ''
                    reaction.subsystem = ''
                    reaction.lower_bound = reaction_lower_bound  
                    reaction.upper_bound = reaction_upper_bound
                    reaction.annotation.update({
                        'rp_score': self.rp_scope[rp_rule_trans_id]['score'],
                        'rp_step': path_step,
                        'rp_id': rp_rule,
                        'mnxr': self.single_depr_mnxr(rp_subpath[path_step].get('reaction')),
                        #add the ec from out_scope
                    })
                    #Metabolite
                    for cid in rp_reactants|rp_products:
                        if not cid in model_meta:
                            try:
                                xref = self.rp_strc[cid]
                            except KeyError:
                                xref, _ = self.mnxm_xref(cid)
                            model_meta[cid] = Metabolite(
                                    f'{cid}_{compartment_id}',
                                    formula=xref.get('formula', ''),
                                    name=xref.get('name', cid),
                                    charge=xref.get('charge', ''),
                                    compartment=compartment_id,
                                
                            )
                            if 'xref' in xref:
                                model_meta[cid].annotation.update(xref['xref'])
                            else:
                                model_meta[cid].annotation.update(xref)
                    model_reaction_dict = {}
                    for cid in rp_products:
                        model_reaction_dict[model_meta[cid]] = rp_products[cid]
                    for cid in rp_reactants:
                        model_reaction_dict[model_meta[cid]] = -rp_reactants[cid]
                    reaction.add_metabolites(
                        model_reaction_dict
                    )
                    model_reac.append(reaction)
                model.add_reactions(model_reac)
                to_ret[rp_path_num][idx] = model
        return to_ret
