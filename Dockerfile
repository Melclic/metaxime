FROM brsynth/rpbase:dev

WORKDIR /home/

RUN apt-get install --quiet --yes --no-install-recommends \
			libxext6  \
    	libxrender-dev \
	 && conda install -y -c rdkit rdkit

RUN mkdir /home/input_cache/
RUN mkdir /home/cache/

### MNXref Version 2019/02/13 ###
COPY input_cache/chem_prop.tsv /home/input_cache/
COPY input_cache/chem_xref.tsv /home/input_cache/
COPY input_cache/comp_xref.tsv /home/input_cache/
COPY input_cache/reac_xref.tsv /home/input_cache/
#### new reaction rules ####
COPY input_cache/rr_compounds.tsv /home/input_cache/
COPY input_cache/rules_rall.tsv /home/input_cache/
COPY input_cache/rxn_recipes.tsv /home/input_cache/
#### thermo ##
COPY thermo_input_cache.tar.xz /
RUN tar xf /thermo_input_cache.tar.xz -C /
RUN mv /thermo_input_cache/* /home/input_cache/
RUN mv /home/input_cache/cc_preprocess.npz /home/cache/cc_preprocess.npz

COPY rpCache.py /home/
RUN python rpCache.py

RUN rm -r input_cache/
