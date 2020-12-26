FROM continuumio/conda-ci-linux-64-python3.8:latest

USER root
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /home/

### install all the requirements ###
COPY requirements.txt /home/
RUN sh -c 'echo "deb http://ftp.us.debian.org/debian sid main" >> /etc/apt/sources.list'
#fix because of debian update-alternatives limitations of not considering anything outside of /usr/share/man
RUN mkdir -p /usr/share/man/man1
RUN apt-get update
RUN sed 's/#.*//' /home/requirements.txt | xargs apt-get install -y
RUN apt-get clean
RUN apt-get autoremove -y
#RUN rm -rf /var/lib/apt/lists/*

##### conda install ###
RUN conda update -n base -c defaults conda
RUN conda install -y -c conda-forge python-libsbml rdkit networkx==2.3 numpy pandas openbabel timeout-decorator cython
RUN conda install -y -c anaconda biopython==1.77
RUN conda install -y -c biobuilds t-coffee
RUN conda install -y -c bioconda emboss

RUN pip install equilibrator-pathway==0.3.1 timeout-decorator objsize shared_memory_dict graphviz pydotplus lxml redis rq flask-restful flask-cors

RUN rm -rf $(dirname  $(which python))/../lib/python3.8/site-packages/ruamel*
RUN pip install cobra==0.16

###### MARVIN ####

RUN mkdir /home/extra_packages/
COPY marvin_linux_20.9.deb /home/
COPY license.cxl /home/extra_packages/
ENV CHEMAXON_LICENSE_URL /home/extra_packages/license.cxl
RUN dpkg -i /home/marvin_linux_20.9.deb
RUN rm /home/marvin_linux_20.9.deb

#### extra install from source or from git ####
WORKDIR /home/extra_packages/
RUN git clone https://gitlab.irstea.fr/jacques.fize/GMatch4py.git
RUN cd GMatch4py && pip install . && cd /home/extra_packages/

RUN git clone --single-branch --branch develop https://gitlab.com/equilibrator/equilibrator-api.git
RUN cd equilibrator-api && pip install -e . && cd /home/extra_packages/

RUN git clone https://gitlab.com/equilibrator/equilibrator-assets.git
RUN cd equilibrator-assets && pip install -e . && cd /home/extra_packages/
RUN cd /home/

#############################################
########## Equilibrator #####################
#############################################

COPY init_equilibrator.py /home/extra_packages/
RUN chmod +x /home/extra_packages/init_equilibrator.py
RUN python /home/extra_packages/init_equilibrator.py

