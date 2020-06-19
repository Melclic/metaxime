FROM brsynth/rpbase:dev

WORKDIR /home/

RUN apt-get install --quiet --yes --no-install-recommends \
			libxext6  \
    	libxrender-dev \
	 && conda install -y -c rdkit rdkit

RUN mkdir /home/input_cache/

### MNXref Version 2019/02/13 ###
COPY input_cache/chem_prop.tsv /home/input_cache/
COPY input_cache/chem_xref.tsv /home/input_cache/
COPY input_cache/comp_xref.tsv /home/input_cache/
COPY input_cache/reac_xref.tsv /home/input_cache/
#### new reaction rules ####
COPY input_cache/rr_compounds.tsv /home/input_cache/
COPY input_cache/rules_rall.tsv /home/input_cache/
COPY input_cache/rxn_recipes.tsv /home/input_cache/

COPY rpCache.py /home/
RUN python rpCache.py

RUN rm -r input_cache/
