#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Merge SBML models
 *
 * Example:
 *   nextflow run merge_model.nf --source_input_model user_model.xml --target_input_model iML1515.xml
 */


/* ----------------------------
 * Parameters
 * ---------------------------- */
params.help = false
params.output_folder = 'results'
params.merge = [ source_input_model: null, target_input_model: null ]
params.parse = [ rp2_scope: null, rp2_compounds: null, rp2_paths: null, compartment_id: 'c' ]

/* ----------------------------
 * Process
 * ---------------------------- */
process merge_models {
    errorStrategy 'terminate'
    container "metaxime:rewrite"

    publishDir (
        path: { "${params.output_folder}/" },
        mode: "copy",
        pattern: "merged_model.xml"
    )

    input:
      path source_model_file
      path target_model_file

    output:
      path "merged_model.xml", emit: merged_model

    script:
    """
#!/usr/bin/env python3
import sys
from cobra.io import write_sbml_model
from metaxime.utils import merge_models
import biopathopt

#source_model = read_sbml_model("${source_model_file}")
source_model =  biopathopt.ModelBuilder("${source_model_file}")
#target_model = read_sbml_model("${target_model_file}")
target_model =  biopathopt.ModelBuilder("${target_model_file}")
merged = merge_models(source_model.model, target_model.model)
write_sbml_model(merged, 'merged_model.xml')
    """
}

process parse_rp2 {
  container "metaxime:rewrite"
  errorStrategy 'terminate'

  publishDir(
    path: { "${params.output_folder}/" },
    mode: 'copy',
    pattern: '*.xml'
  )

  input:
    path rp2_scope_file
    path rp2_compounds_file
    path rp2_paths_file
    val comp_id

  output:
    path "*.xml", emit: models

  script:
  """
#!/usr/bin/env python3
import json
import pathlib
from cobra.io import write_sbml_model
from metaxime.parser import ParserRP2

# Parse and build models
parser = ParserRP2(
    rp2_scope_path="${rp2_scope_file}",
    rp2_cmp_path="${rp2_compounds_file}",
    rp2_paths_path="${rp2_paths_file}",
)
all_models = parser.return_rp2_models(
    compartment_id="${comp_id}"
) 

# Save each model
for path_id in all_models:
    for sub_path_id in all_models[path_id]:
        write_sbml_model(all_models[path_id][sub_path_id], f"{all_models[path_id][sub_path_id].id}.xml")
  """
}

/* ----------------------------
 * Workflow
 * ---------------------------- */

def helpMessageMerge() {
    log.info """
╭─────────────────────────────────────────────────────────────╮
│                       Merge COBRA Models                     │
╰─────────────────────────────────────────────────────────────╯

Usage:
  nextflow run merge_model.nf --input_model <SBML_FILE> [--out_name merged_model.xml]

Options:
  --source_input_model <path>   Path to source input SBML model (required)
  --target_input_model <path>   Path to target input SBML model (required)
  --output_folder <path>    Path to the folder output (default: metaxime)
  --help                 Show this help message

Example:
 nextflow run merge_model.nf --source_input_model user_model.xml --target_input_model iML1515.xml
""".stripIndent()
}

workflow merge {
    def P = (params.merge ?: [:]) as Map

    ch_source_model = Channel.fromPath(P.get('source_input_model')?: '', checkIfExists: true)
    ch_target_model = Channel.fromPath(P.get('target_input_model')?: '', checkIfExists: true)

    if (params.help || !P.source_input_model || !P.target_input_model){
        helpMessageMerge()
        exit 0
    }

    merge_models(ch_source_model, ch_target_model)
}

def helpMessageParse() {
  log.info """
╭────────────────────────────────────────────────────────────╮
│                 RP2 → COBRA Models Exporter                │
╰────────────────────────────────────────────────────────────╯

Options:
    --rp2_scope <path>          Path to RP2 out_scope.csv file (required)
    --rp2_compounds <path>      Path to RP2 out_compounds.csv file (required)
    --rp2_paths <path>          Path to RP2 out_paths.csv file (required)
    --output_folder <path>      Output directory for generated models
                                (default: "metaxime")
    --compartment_id <string>   Compartment ID to assign to metabolites
                                (default: "c")
    --help                      Show this help message

Usage:
    nextflow run parse_rp2_models.nf \\
        --rp2_scope out_scope.csv \\
        --rp2_compounds out_compounds.csv \\
        --rp2_paths out_paths.csv \\
        --output_folder metaxime_models \\
        --compartment_id c
""".stripIndent()
}

workflow parse {
    def P = (params.parse ?: [:]) as Map
    def COMP_ID = (P.get('compartment_id') ?: 'c') as String

    ch_scope     = Channel.fromPath(P.get('rp2_scope')    ?: '', checkIfExists: true)
    ch_compounds = Channel.fromPath(P.get('rp2_compounds')?: '', checkIfExists: true)
    ch_paths     = Channel.fromPath(P.get('rp2_paths')    ?: '', checkIfExists: true)
    ch_compartment = Channel.value(COMP_ID)

    if( params.help || !P.rp2_scope || !P.rp2_compounds || !P.rp2_paths ) {
        helpMessageParse()
        exit 0
    }

    parse_rp2(ch_scope, ch_compounds, ch_paths, ch_compartment)
}
