#!/usr/bin/env python3
import argparse
import logging
from pathlib import Path
import zipfile
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

    out_zip = Path(args.out_tar).resolve()
    out_zip.parent.mkdir(parents=True, exist_ok=True)

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
                    out_file = tmpdir / f"{model_id}.xml"
                    write_sbml_model(merged, str(out_file))
                    logging.info("Saved %s", out_file)
                except ValueError as e:
                    logging.warning("Error in %s: %s", model_id, e)
        # Create zip with maximum compression
        with zipfile.ZipFile(
            out_zip,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=9,
        ) as zf:
            for xml_file in tmpdir.glob("*.xml"):
                arcname = f"merged_rp2_models/{xml_file.name}"
                zf.write(xml_file, arcname=arcname)
                logging.debug("Added %s as %s", xml_file, arcname)
        logging.info("Archive created at: %s", out_zip)
        logging.info("Temporary folder will be removed when context ends")
    logging.info("Done")


if __name__ == "__main__":
    main()