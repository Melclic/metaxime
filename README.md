# MetaXime

This projects makes it more accessible to study the output of RetroPath2.0 by converting the results to metabolic models for downstream analysis. The project also uses CobraK to perform different 

## Requirements

[BioPathOpt](https://github.com/melclic/bioPathOpt)

## Introduction

This project parses the output of [RetroPath2.0](https://www.myexperiment.org/workflows/4987/versions/16.html) and completes the monocomponent reactions of the tool to full reactions. This includes all the cofactors that makes the reaction complete and extracts them to a cobrapy metabolic model. This includes a sink and source to enable the user to simulate the production of the molecule. The project also has the ability to merge two models together, and works seamlessly with the [retrosynthesis](https://github.com/melclic/retrosynthesis). 

This is part of one of the projects that may be found [here](https://github.com/brsynth) and lives as a [Galaxy platform](https://galaxy-synbiocad.org/). The point of this project is to use an alternative method and enable researchers to use other workflow tools to run the design pipeline.

## How to run 

To predict metabolic pathways, click on "Run Job" on the top banner. 

![Run Job](images/run_job.png)

There you will need to provide the molecule that you would like to produce, the organism in which you would like to produce the molecule of interest, and the maximum number of reactions that the pathways may have.

![Search](images/search_name.png)
![SMILES](images/search_smiles.png)
![Draw](images/search_draw.png)

The results may be inspected once the job is completed. The following error 

## Build

To build the project, navigate to the project folder and run the following commands:

```
docker build -t metaxime -f images/Dockerfile .
docker run -it -v $(pwd)/mx-results:/mx-results -p 80:80 -p 8888:8888 metaxime
```

The results should be accessible at: `http://localhost:80`

## Structure

This project may be separated into three different parts, the retrosynthesis part to predict new pathways, the pathway analysis and the front end of the project

### RetroSynthesis

These include [RetroPath2](https://myexperiment.org/workflows/4987.html), [rp2paths](https://github.com/brsynth/rp2paths) and [RetroRules](https://retrorules.org/). (TODO: finish separating that part into its own docker and API service: https://github.com/Melclic/retrosynthesis). 

RetroPath2 is a [KNIME](https://www.knime.com/) workflow and requires the installation of the workflow execution environment to enable its execution. Rp2paths is a python package and RetroRules is a database of reaction rules. The script in this project parses database and extracts only information that is to be used in a particular execution (see reaction rule diameters. More info [here](https://www.jfaulon.com/retropath2-0-a-retrosynthesis-workflow-for-metabolic-engineers-biorxiv/)).   

### Pathway Analysis

The pathway analysis part of the project may be found under the metaxime/ folder and is the central part of this project. It contains all the scripts that parse the output of RetroPath2 and rp2paths and turn them into [SBML](https://en.wikipedia.org/wiki/SBML) metabolic model files. This project supports the [MIRIAM](https://en.wikipedia.org/wiki/Minimum_information_required_in_the_annotation_of_models) standard for SBML files, that contains the annotations for species, reactions and pathways to link these properties of a metabolic model to various databases.

We call "rpSBML" SBML files that are enriched to contain extra information related to the heterologous pathway. This includes the thermodynamic properties of the chemical species, reaction and pathways. It also holds the FBA run information and any extra annotation related to chemoinformatics properties of the predicted reactions.

![Enriched Annotation](images/rp_annot.png)

More detail on each of these files may be found in the various github projects described in the top of this README.

### Front End and Architecture

The architecture of the project as follows:

![Structure](images/strct.jpg)

The output of the project is a compressed file named `rpcollection.tar.xz`. 

```
├── log.json
├── model_json
│   └── rp_*_rpsbml.json
├── models
│   └─── rp_*_rpsbml.xml
└── networks
    └── rp_*_rpsbml.json
```

The web interface (in the static/ folder) contains all the JS and HTML. The `draw_network.js` renders the network JSON of the parsed SBML by [networkX](https://networkx.org/) and renders it on the page using the [dagre-d3](https://github.com/dagrejs/dagre-d3) project. 

Rendering the molecules in javascript is done by a branched project called [smiledDrawer](https://github.com/Melclic/smilesDrawer).
