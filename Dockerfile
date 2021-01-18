FROM conda/miniconda3

##### conda install ###
RUN conda install -y --no-update-dependencies -c conda-forge python-libsbml rdkit networkx==2.3 numpy pandas openbabel timeout-decorator cython
RUN conda install -y --no-update-dependencies -c anaconda biopython==1.77
RUN conda install -y --no-update-dependencies -c biobuilds t-coffee
RUN conda install -y --no-update-dependencies -c bioconda emboss

RUN rm -rf $(dirname  $(which python))/../lib/python3.7/site-packages/ruamel*
RUN pip install cobra==0.16

###### MARVIN ####
RUN mkdir /home/extra_packages/
COPY marvin_linux_20.9.deb /home/
COPY license.cxl /home/extra_packages/
ENV CHEMAXON_LICENSE_URL /home/extra_packages/license.cxl
RUN dpkg -i /home/marvin_linux_20.9.deb
RUN rm /home/marvin_linux_20.9.deb

#RUN pip install equilibrator-pathway==0.3.1 objsize shared_memory_dict graphviz pydotplus lxml
#RUN pip install equilibrator-pathway objsize shared_memory_dict graphviz pydotplus lxml
RUN pip install equilibrator-pathway objsize graphviz pydotplus lxml

####### equilibrator ########
RUN pip install equilibrator-api
RUN python -c "from equilibrator_api import ComponentContribution; cc = ComponentContribution()"

#### extra install from source or from git ####
WORKDIR /home/extra_packages/
RUN git clone -b mel --single-branch https://gitlab.com/Melclic/equilibrator-assets.git
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
RUN python /home/init_cache.py
