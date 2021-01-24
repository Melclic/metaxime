FROM debian:latest

RUN apt-get update \
 && apt-get install -yq --no-install-recommends \
    ca-certificates \
    build-essential \
    cmake \
    wget \
    libboost-dev \
    libboost-iostreams-dev \
    libboost-python-dev \
    libboost-regex-dev \
    libboost-serialization-dev \
    libboost-system-dev \
    libboost-thread-dev \
    libcairo2-dev \
    libeigen3-dev \
    python3-dev \
    python3-pip \
    default-jdk \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN pip3 install numpy

#### install RDKIT #####
RUN cd /
ARG RDKIT_VERSION=Release_2020_09_3
RUN wget --quiet https://github.com/rdkit/rdkit/archive/${RDKIT_VERSION}.tar.gz \
	&& tar -xzf ${RDKIT_VERSION}.tar.gz \
	&& mv rdkit-${RDKIT_VERSION} /rdkit \
	&& rm ${RDKIT_VERSION}.tar.gz
RUN cd rdkit/External/INCHI-API && \
	./download-inchi.sh

WORKDIR /rdkit/build/

RUN cmake -D RDK_BUILD_INCHI_SUPPORT=ON \ 
          -D PYTHON_EXECUTABLE=/usr/bin/python3.7 \
	.. && \
	make && \
	make install 

RUN make -j $(nproc) \
	&& make install

ENV RDBASE /rdkit
ENV LD_LIBRARY_PATH $RDBASE/lib
ENV PYTHONPATH $PYTHONPATH:$RDBASE

###### install the pip packages ####

RUN rm -vf /var/lib/apt/lists/*
RUN apt-get update

RUN apt-get install -y git \
	software-properties-common \
	xz-utils \
	python3-dev \
	python3-setuptools \
	libxml2-dev \
	swig	

RUN pip3 install wheel
###### openbabel #####
WORKDIR /
RUN wget https://github.com/openbabel/openbabel/archive/openbabel-3-1-1.tar.gz
RUN tar -xf openbabel-3-1-1.tar.gz
RUN cd /openbabel-openbabel-3-1-1/
RUN mkdir build
WORKDIR /openbabel-openbabel-3-1-1/build/
RUN cmake -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON  -DPYTHON_EXECUTABLE=/usr/bin/python3.7 ..
RUN make -j $(nproc)
RUN make install

RUN pip3 install python-libsbml \
	networkx==2.3 \
	pandas \
	#openbabel \
	timeout-decorator \
	cython \
	biopython==1.77 \
	equilibrator-api \
	equilibrator-pathway \
	objsize \
	#shared_memory_dict \
	graphviz \
	pydotplus \
	lxml

RUN rm -rf $(dirname  $(which python))/../lib/python3.7/site-packages/ruamel*
RUN pip3 install cobra==0.16

WORKDIR /home/

###### MARVIN ####
WORKDIR /home/extra_packages/
COPY marvin_linux_20.9.deb /home/
COPY license.cxl /home/extra_packages/
ENV CHEMAXON_LICENSE_URL /home/extra_packages/license.cxl
RUN dpkg -i /home/marvin_linux_20.9.deb
RUN rm /home/marvin_linux_20.9.deb

#### extra packages #####
RUN apt-get install -y t-coffee emboss
WORKDIR /home/extra_packages 
RUN git clone -b mel --single-branch https://gitlab.com/Melclic/equilibrator-assets.git
RUN cd equilibrator-assets && pip3 install -e . && cd /home/extra_packages/

RUN git clone https://gitlab.irstea.fr/jacques.fize/GMatch4py.git
RUN cd GMatch4py && pip3 install -e . && cd /home/extra_packages/

RUN apt-get clean
RUN apt-get autoremove -y

WORKDIR /home/

#### generate the cache
COPY init_cache.py /home/
ADD metaxime /home/metaxime/
ADD selenzy /home/selenzy/
RUN mkdir /home/metaxime/input_cache/
COPY metaxime/input_cache/rpselenzyme_data.tar.xz /home/metaxime/input_cache/
RUN python3 init_cache.py
