```
nextflow run main.nf -entry merge --merge.source_input_model input/rp2_model.xml --merge.target_input_model input/enriched_model.xml
```

```
nextflow run main.nf -entry parse --parse.rp2_scope input/out_scope.csv --parse.rp2_compounds input/out_compounds.csv --parse.rp2_paths input/out_paths.csv
```
