nextflow.enable.dsl=2


def helpMessage() {
    log.info """
MetaXime retrosynthesis pipeline

Steps
  1. EXTRACT_SINK     extract sink.csv from SBML model
  2. RETROSYNTHESIS   run RP2 retrosynthesis on sink and source InChI
  3. PARSE_RP2        merge RP2 results into target COBRA model and export zip

Required parameters
  --model_file <path>         SBML model xml or json (used for sink and as target model)
  --source_inchi <string>     Source InChI for retrosynthesis
  --output_folder <path>      Output folder

Optional retrosynthesis parameters
  --rules_file <path>         Rules file for RP2 (use NONE.rules to disable, or omit to use RP2 defaults)
                              [${params.rules_file}]
  --std_mode <string>         Standardisation mode
                              [${params.std_mode}]
  --max_steps <int>           Maximum RP2 steps
                              [${params.max_steps}]
  --topx <int>                Top X pathways
                              [${params.topx}]
  --accept_partial_results    Accept partial results (flag)
  --diameters <string>        Comma separated diameters
                              [${params.diameters}]
  --rule_type <string>        Rule type
                              [${params.rule_type}]
  --ram_limit <int>           RAM limit in GB
                              [${params.ram_limit}]

MetaXime merge parameters
  --source_comp <string>      Source compartment id (RP2 models)
                              [${params.source_comp}]
  --target_comp <string>      Target compartment id (COBRA model)
                              [${params.target_comp}]
  --use_inchikey2             Use InChIKey2 fallback matching (flag)
  --find_all_parentless       Enable parentless metabolite search (flag)

Usage example
  nextflow run main.nf \\
    --model_file models/iML1515.xml \\
    --source_inchi "InChI=1S/..." \\
    --rules_file data/rules.json \\
    --out_tar results/rp2_merged.zip \\
    --max_steps 6 \\
    --topx 1000 \\
    --use_inchikey2 \\
    --find_all_parentless
""".stripIndent()
}

// Required
params.model_file = null
params.source_inchi = null
params.out_tar = null

// Optional
params.rules_file = null
params.std_mode   = "H added + Aromatized"
params.max_steps  = 6
params.topx       = 1000
params.accept_partial_results = false

params.diameters  = "2,4,6,8,10,12,14,16"
params.rule_type  = "all"
params.ram_limit  = 18

params.source_comp = "c"
params.target_comp = "c"
params.use_inchikey2       = false
params.find_all_parentless = false

params.help = false

process EXTRACT_SINK {

    tag "sink_extractor"
    container "melclic/biopathopt:latest"

    input:
        path model_file

    output:
        path "sink.csv", emit: sink_csv

    script:
        """
        python /home/extract_sink.py \\
            --model ${model_file} \\
            --out   sink.csv
        """
}

process RETROSYNTHESIS {

    tag "retrosynthesis"
    container "melclic/retrosynthesis:latest"

      input:
        path sink_file
        val source_inchi
        val accept_partial_results
        path rules_file

    output:
        path "out_paths.csv", emit: out_paths
        path "out_compounds.csv", emit: out_compounds
        path "out_scope.csv", emit: out_scope

    script:
        def partial_flag = accept_partial_results ? "--accept-partial-results" : ""
        def rules_flag = rules_file.getName() == "NONE.rules" ? "" : "--rules-file ${rules_file}"
        """
        python /home/rp2/retropipeline.py \
            --sink-file ${sink_file} \
            --source-inchi "${source_inchi}" \
            --std-mode "${params.std_mode}" \
            --max-steps ${params.max_steps} \
            --topx ${params.topx} \
            --diameters "${params.diameters}" \
            --rule-type "${params.rule_type}" \
            --ram-limit "${params.ram_limit}" \
            ${partial_flag} \
            ${rules_flag} \
            --out-paths out_paths.csv \
            --out-compounds out_compounds.csv \
            --out-scope out_scope.csv
        """
}

process PARSE_RP2 {

    tag "parse_rp2"
    container "melclic/metaxime:latest"

    publishDir (
        path: { "${params.output_folder}/" },
        mode: "copy", 
        pattern: "output.zip"
    )
    publishDir (
        path: { "${params.output_folder}/" },
        mode: "copy", 
        pattern: "output.json"
    )

    input:
        path(scope_file)
        path(compounds_file)
        path(paths_file)
        path(target_model_file)

    output:
        path "output.zip", emit: merged_zip
        path "output.json", emit: summary_graphs

    script:
        """
        python /home/MetaXime/run_pipeline.py \\
            --scope        ${scope_file} \\
            --compounds    ${compounds_file} \\
            --paths        ${paths_file} \\
            --target_model ${target_model_file} \\
            --out_tar      output.zip \\
            --out_json     output.json \\
            --source_comp  ${params.source_comp} \\
            --target_comp  ${params.target_comp} \\
            ${params.use_inchikey2        ? "--use_inchikey2" : ""} \\
            ${params.find_all_parentless  ? "--find_all_parentless" : ""}
        """
}

workflow {
    if (params.help) {
        helpMessage()
        System.exit(0)
    }

    // rules_file is no longer required here
    if (!params.model_file || !params.source_inchi) {
        log.error "Missing required parameters: --model_file, --source_inchi"
        helpMessage()
        System.exit(1)
    }

    /*
     * We need the same model file to go to:
     *  - EXTRACT_SINK  (to build sink.csv)
     *  - PARSE_RP2      (as target_model_file)
     */
    // declare channel variables so into{} can assign them
    //def model_for_sink_ch
    model_ch = channel.fromPath(params.model_file, checkIfExists: true)

    // values for source_inchi and accept_partial_results
    def source_inchi_ch = channel.value(params.source_inchi)
    def accept_ch       = channel.value(params.accept_partial_results as boolean)
    def ch_rules = params.rules_file ? channel.fromPath(params.rules_file, checkIfExists: true) : channel.fromPath("NONE.rules")

    model_ch.view { v -> "model_file: $v" }
    source_inchi_ch.view { v -> "source_inchi: $v" }
    accept_ch.view { v -> "accept_partial_results: $v" }
    ch_rules.view { v -> "passed_rules_file: $v" }

    // run EXTRACT_SINK first
    extract_sink_ch = EXTRACT_SINK(model_ch)
    extract_sink_ch.view { v -> "extract_sink_ch: $v" }
    extract_sink_ch.sink_csv.view { v -> "extract_sink_ch.sink_csv: $v" }

    retro_ch = RETROSYNTHESIS(
        extract_sink_ch.sink_csv, 
        source_inchi_ch, 
        accept_ch,
        ch_rules,
    )

    retro_ch.out_scope.view { v -> "retro_ch.out_scope: $v" }
    retro_ch.out_compounds.view { v -> "retro_ch.out_compounds: $v" }
    retro_ch.out_paths.view { v -> "retro_ch.out_paths: $v" }

    PARSE_RP2(
        retro_ch.out_scope,
        retro_ch.out_compounds,
        retro_ch.out_paths,
        model_ch
    )
}