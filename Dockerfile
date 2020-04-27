FROM brsynth/rpbase:newrules

WORKDIR /home/

RUN apt-get install --quiet --yes --no-install-recommends \
			libxext6  \
    	libxrender-dev \
	 && conda install -y -c rdkit rdkit

RUN mkdir /home/input_cache/
COPY newrules_compounds.tsv /home/input_cache/rr_compounds.tsv
COPY newrules_rules.tsv /home/input_cache/rules_rall.tsv
COPY newrules_rxn_recipes.tsv /home/input_cache/rxn_recipes.tsv

COPY rpCache.py /home/
RUN python rpCache.py
