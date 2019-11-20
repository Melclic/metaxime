FROM brsynth/rpbase

RUN apt-get install --quiet --yes --no-install-recommends \
			libxext6  \
    	libxrender-dev \
	 && conda install -y -c rdkit rdkit

COPY rpCache.py rpToolCache.p[y] /home/
