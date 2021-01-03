FROM informaticsmatters/rdkit-python3-debian:latest

USER root

RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

ENV PATH /opt/conda/bin:$PATH

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

####### equilibrator ########

RUN pip install equilibrator-api
RUN python -c "from equilibrator_api import ComponentContribution; cc = ComponentContribution()"

##### conda install ###
RUN conda update -n base -c defaults conda
RUN conda install -y -c conda-forge python-libsbml rdkit networkx==2.3 numpy pandas openbabel timeout-decorator cython
RUN conda install -y -c anaconda biopython==1.77
RUN conda install -y -c biobuilds t-coffee
RUN conda install -y -c bioconda emboss

RUN rm -rf $(dirname  $(which python))/../lib/python3.8/site-packages/ruamel*
RUN pip install cobra==0.16

###### MARVIN ####

RUN mkdir /home/extra_packages/
COPY marvin_linux_20.9.deb /home/
COPY license.cxl /home/extra_packages/
ENV CHEMAXON_LICENSE_URL /home/extra_packages/license.cxl
RUN dpkg -i /home/marvin_linux_20.9.deb
RUN rm /home/marvin_linux_20.9.deb

#RUN pip install equilibrator-pathway==0.3.1 objsize shared_memory_dict graphviz pydotplus lxml
RUN pip install equilibrator-pathway objsize shared_memory_dict graphviz pydotplus lxml

#### extra install from source or from git ####
WORKDIR /home/extra_packages/

RUN mkdir equilibrator-assets
COPY extra_packages/equilibrator-assets /home/extra_packages/equilibrator-assets/

RUN cd equilibrator-assets && pip install -e . && cd /home/extra_packages/

RUN git clone https://gitlab.irstea.fr/jacques.fize/GMatch4py.git
RUN cd GMatch4py && pip install -e . && cd /home/extra_packages/

#############################################
########## Equilibrator #####################
#############################################

WORKDIR /home/metaxime/

COPY metaxime /home/metaxime/

WORKDIR /home/

COPY selenzy/ /home/selenzy/

COPY init_cache.py /home/
RUN chmod +x /home/init_cache.py
#RUN python /home/init_cache.py

USER rdkit

