FROM brsynth/rpbase

RUN apt-get install --quiet --yes --no-install-recommends \
			libxext6  \
    	libxrender-dev \
	 && conda install -y -c rdkit rdkit

COPY rpCache.py /home/
RUN python rpCache.py

#ONBUILD COPY Dockerfile rpToolCache.p[y] /home/
