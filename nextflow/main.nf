nextflow.enable.dsl=2

/*
 * Help message
 */
def helpMessageParse() {
    log.info """
╭────────────────────────────────────────────────────────────╮
│          RP2 to COBRA model merger and zip exporter        │
╰────────────────────────────────────────────────────────────╯

Options:
    --scope <path>              Path to RP2 out_scope csv file (required)
    --compounds <path>          Path to RP2 out_compounds csv file (required)
    --paths <path>              Path to RP2 out_paths csv file (required)
    --target_model <path>       Path to target COBRA SBML model (required)
    --output_folder <path>      Output folder (required)

    --source_comp <string>      Source compartment id for RP2 models
                                (default: "c")
    --target_comp <string>      Target compartment id for COBRA model
                                (default: "c")

    --use_inchikey2             Enable InChIKey2 based fallback matching
    --find_all_parentless       Enable search and rescue of parentless metabolites

    --help                      Show this help message

Usage:
    nextflow run main.nf \\
        --scope out_scope.csv \\
        --compounds out_compounds.csv \\
        --paths out_paths.csv \\
        --target_model iML1515.xml \\
        --out_tar merged_rp2_models.zip \\
        --source_comp c \\
        --target_comp c \\
        --use_inchikey2 \\
        --find_all_parentless
""".stripIndent()
}

params.scope           = null
params.compounds       = null
params.paths           = null
params.target_model    = null

params.out_tar         = "merged_rp2_models.zip"

params.source_comp     = "c"
params.target_comp     = "c"

params.use_inchikey2        = true
params.find_all_parentless  = true

process RUN_MERGE_PIPELINE {

    tag "merge_rp2_models"
    container "melclic/metaxime:latest"

    publishDir (
        path: { "${params.output_folder}/" },
        mode: "copy", 
        pattern: "output.zip"
    )

    input:
        path scope_file
        path compounds_file
        path paths_file
        path target_model_file

    output:
        path "output.zip", emit: merged_zip

    script:
        """
        python /home/MetaXime/run_pipeline.py \
            --scope ${scope_file} \
            --compounds ${compounds_file} \
            --paths ${paths_file} \
            --target_model ${target_model_file} \
            --out_tar output.zip \
            --source_comp ${params.source_comp} \
            --target_comp ${params.target_comp} \
            ${params.use_inchikey2        ? "--use_inchikey2" : ""} \
            ${params.find_all_parentless  ? "--find_all_parentless" : ""}
        """
}


workflow {

    if (params.help) {
        helpMessageParse()
        exit 0
    }

    if (!params.scope || !params.compounds || !params.paths || !params.target_model) {
        log.error "Missing required arguments"
        helpMessageParse()
        exit 1
    }

    scope_ch = channel.fromPath(params.scope, checkIfExists: true)
    compounds_ch = channel.fromPath(params.compounds, checkIfExists: true)
    paths_ch = channel.fromPath(params.paths, checkIfExists: true)
    target_model_ch = channel.fromPath(params.target_model, checkIfExists: true)

    RUN_MERGE_PIPELINE(
        scope_ch,
        compounds_ch,
        paths_ch,
        target_model_ch
    )
}