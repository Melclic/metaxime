<<<<<<< HEAD
FROM brsynth/rpbase:newrules

WORKDIR /home/
=======
FROM brsynth/rpbase:dev
>>>>>>> dev

RUN apt-get install --quiet --yes --no-install-recommends \
			libxext6  \
    	libxrender-dev \
	 && conda install -y -c rdkit rdkit

RUN mkdir /home/input_cache/
COPY newrules_compounds.tsv /home/input_cache/rr_compounds.tsv
COPY newrules_rules.tsv /home/input_cache/rules_rall.tsv
COPY newrules_rxn_recipes.tsv /home/input_cache/rxn_recipes.tsv

COPY rpCache.py /home/

RUN mkdir /home/input_cache/

COPY input_cache/chem_prop.tsv /home/input_cache/
COPY input_cache/chem_xref.tsv /home/input_cache/
COPY input_cache/comp_xref.tsv /home/input_cache/
COPY input_cache/reac_xref.tsv /home/input_cache/
COPY input_cache/rr_compounds.tsv /home/input_cache/
COPY input_cache/rules_rall.tsv /home/input_cache/
COPY input_cache/rxn_recipes.tsv /home/input_cache/

RUN python rpCache.py
