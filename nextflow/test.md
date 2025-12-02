```
nextflow run main.nf -entry merge --merge.source_input_model input/rp2_model.xml --merge.target_input_model input/enriched_model.xml
```

```
nextflow run main.nf -entry parse --parse.rp2_scope input/out_scope.csv --parse.rp2_compounds input/out_compounds.csv --parse.rp2_paths input/out_paths.csv
```

```
nextflow run complete_workflow.nf --model_file /Users/melchiordulac/workspace/metaxime/notebooks/iML1515.xml --source_inchi "InChI=1S/C6H6O2/c7-5-3-1-2-4-6(5)8/h1-4,7-8H" --output_folder /Users/melchiordulac/workspace/metaxime/nextflow/ --max_steps 3
```
