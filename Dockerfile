FROM continuumio/miniconda3:4.9.2

RUN apt-get update && apt-get install build-essential -y
RUN conda install -c conda-forge rdkit openbabel networkx==2.3 cython python-libsbml 
RUN conda install -c conda-forge cobra
RUN pip3 install objsize

###### MARVIN ####

WORKDIR /home/extra_packages/

ENV MARVIN_VERSION=20.9

COPY marvin_linux_$MARVIN_VERSION.deb /home/extra_packages/
COPY license.cxl /home/extra_packages/
ENV CHEMAXON_LICENSE_URL /home/extra_packages/license.cxl
RUN dpkg -i /home/extra_packages/marvin_linux_$MARVIN_VERSION.deb
RUN rm /home/extra_packages/marvin_linux_$MARVIN_VERSION.deb

##### GMatch4py #####

#graph calculations
RUN conda install -c conda-forge scipy threadpoolctl joblib scikit-learn gensim psutil cython
RUN pip3 install smart-open
RUN git clone https://gitlab.irstea.fr/jacques.fize/GMatch4py.git
RUN cd GMatch4py && pip3 install -e . && cd /home/extra_packages/

##### Equilibrator ####

#install the develop version of equilibrator
#RUN pip install equilibrator-api equilibrator-cache equilibrator-pathway
#RUN git clone --single-branch --branch develop https://gitlab.com/equilibrator/equilibrator-api.git
#RUN cd equilibrator-api && pip3 install -e . && cd ..

#RUN git clone -b mel --single-branch https://gitlab.com/Melclic/equilibrator-assets.git
#RUN git clone https://gitlab.com/equilibrator/equilibrator-assets.git
#RUN cd equilibrator-assets && pip3 install -e . && cd ..

#equilibrator-pathway
#RUN pip install equilibrator-pathway==0.3.1
#RUN git clone --single-branch --branch develop https://gitlab.com/equilibrator/equilibrator-pathway.git
#RUN cd equilibrator-pathway && pip3 install -e . && cd ..

RUN conda install -c conda-forge equilibrator-api equilibrator-cache==0.4.2 equilibrator-pathway
RUN git clone https://gitlab.com/equilibrator/equilibrator-assets.git
RUN cd equilibrator-assets && pip3 install -e . && cd ..

WORKDIR /home/

RUN pip3 install timeout_decorator

#### generate the cache
ADD metaxime /home/metaxime/
ADD selenzy /home/selenzy/
RUN mkdir /home/metaxime/input_cache/
COPY selenzy/rpselenzyme_data.tar.xz /home/metaxime/input_cache/
COPY init_eq.py /home/
RUN python3 init_eq.py
COPY init_cache.py /home/
RUN python3 init_cache.py
