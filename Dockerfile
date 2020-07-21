FROM brsynth/rpbase:v2

WORKDIR /home/

RUN apt-get install --quiet --yes --no-install-recommends \
			libxext6  \
    	libxrender-dev \
	 && conda install -y -c rdkit rdkit

RUN mkdir /home/input_cache/
RUN mkdir /home/cache/

### MNXref Version 2019/02/13 ###
#COPY input_cache/chem_prop.tsv /home/input_cache/
#COPY input_cache/chem_xref.tsv /home/input_cache/
#COPY input_cache/comp_xref.tsv /home/input_cache/
#COPY input_cache/reac_xref.tsv /home/input_cache/
#### new reaction rules ####
#COPY input_cache/rr_compounds.tsv /home/input_cache/
#COPY input_cache/rules_rall.tsv /home/input_cache/
#COPY input_cache/rxn_recipes.tsv /home/input_cache/
#### thermo ##
COPY input_cache_thermo.tar.xz /
RUN tar xf /input_cache_thermo.tar.xz -C /
RUN mv /input_cache_thermo/* /home/input_cache/
RUN mv /home/input_cache/cc_preprocess.npz /home/cache/cc_preprocess.npz
RUN rm /input_cache_thermo.tar.xz
RUN rm -r /input_cache_thermo/

COPY rpCache.py /home/
RUN python rpCache.py

RUN rm -r input_cache/
