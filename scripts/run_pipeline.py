#!/usr/bin/env python3
import argparse
import logging
from pathlib import Path
import networkx as nx
import tarfile
import json
import tempfile

from metaxime.parser import ParserRP2
from metaxime.utils import merge_models
from biopathopt import ModelBuilder
from cobra.io import write_sbml_model, read_sbml_model

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

def build_cli():
    parser = argparse.ArgumentParser(
        description="Merge RP2 pathway models into a COBRA base model and export as zip"
    )

    parser.add_argument("--scope", required=True, help="RP2 out_scope csv")
    parser.add_argument("--compounds", required=True, help="RP2 out_compounds csv")
    parser.add_argument("--paths", required=True, help="RP2 out_paths csv")
    parser.add_argument("--target_model", required=True, help="Target COBRA model")
    parser.add_argument("--out_tar", required=True, help="Output zip file")
    parser.add_argument("--out_json", required=True, help="Output json summary file for all pathways")
    parser.add_argument("--source_comp", default="c", help="Source compartment id")
    parser.add_argument("--target_comp", default="c", help="Target compartment id")
    parser.add_argument("--use_inchikey2", action="store_true", help="Use InChIKey2 fallback")
    parser.add_argument("--find_all_parentless", action="store_true", help="Do not include models with parentless heterologous molecules")

    return parser


def main():
    args = build_cli().parse_args()

    scope_path = Path(args.scope).resolve()
    compounds_path = Path(args.compounds).resolve()
    paths_path = Path(args.paths).resolve()

    target_model_path = Path(args.target_model).resolve()

    out_tar = Path(args.out_tar).resolve()
    out_tar.parent.mkdir(parents=True, exist_ok=True)
    out_json = Path(args.out_json).resolve()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    tmp_json = {}

    with tempfile.TemporaryDirectory() as tmpdirname:
        tmpdir = Path(tmpdirname)
        logging.info("Temporary directory: %s", tmpdir)
        parser = ParserRP2(
            rp2_scope_path=str(scope_path),
            rp2_cmp_path=str(compounds_path),
            rp2_paths_path=str(paths_path),
        )
        all_models = parser.return_rp2_models(compartment_id=args.source_comp)
        target_builder = ModelBuilder(str(target_model_path))

        for path_id in all_models:
            for sub_path_id in all_models[path_id]:
                model_id = f"rp2_{path_id}_{sub_path_id}"
                logging.info("Processing %s", model_id)
                try:
                    merged = merge_models(
                        all_models[path_id][sub_path_id],
                        target_builder.model,
                        source_compartment=args.source_comp,
                        target_compartment=args.target_comp,
                        find_all_parentless_source=args.find_all_parentless,
                        use_inchikey2=args.use_inchikey2,
                    )

                    # SBML output
                    sbml_file = tmpdir / f"{model_id}.xml"
                    write_sbml_model(merged, str(sbml_file))
                    logging.info("Saved merged SBML: %s", sbml_file)

                    # Graph JSON output
                    G = parser.cobra_model_to_digraph(all_models[path_id][sub_path_id])
                    # find out what are the childless and parentless to 
                    tmp_G = parser.remove_dangling_reactions(G)
                    parentless_nodes = [n for n in tmp_G.nodes if tmp_G.in_degree(n) == 0]
                    childless_nodes = [n for n in tmp_G.nodes if tmp_G.out_degree(n) == 0]
                    tmp_G = None

                    for i in G.nodes:
                        if i in parentless_nodes:
                            G.nodes[i]['topology'] = 'start'
                        elif i in childless_nodes:
                            G.nodes[i]['topology'] = 'end'
                        else:
                            G.nodes[i]['topology'] = 'intermediate'

                    G_json = nx.node_link_data(G)
                    G_json['steps'] = len([i for i in G.nodes if G.nodes[i]['type']=='reaction'])
                    G_json['id'] = model_id
                    tmp_json[model_id] = G_json

                except ValueError as e:
                    logging.warning("Error in %s: %s", model_id, e)

        # Create tar.gz with maximum gzip compression
        with tarfile.open(out_tar, mode="w:gz", compresslevel=9) as tf:
            # SBML models under merged_sbml_model/
            tf.add(tmpdir, arcname=".")

        logging.info("Archive created at: %s", out_tar)
        logging.info("Temporary folder will be removed when context ends")

    # Write the combined graph summary json
    with out_json.open("w", encoding="utf-8") as fh:
        json.dump(tmp_json, fh, ensure_ascii=False, indent=2)
    logging.info("Graph summary JSON written to: %s", out_json)

    logging.info("Done")


if __name__ == "__main__":
    main()