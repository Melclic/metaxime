FROM continuumio/conda-ci-linux-64-python3.8:latest

USER root
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /home/

### install all the requirements ###
COPY docker_files/apt_requirements.txt /home/
RUN sh -c 'echo "deb http://ftp.us.debian.org/debian sid main" >> /etc/apt/sources.list'
#fix because of debian update-alternatives limitations of not considering anything outside of /usr/share/man
RUN mkdir -p /usr/share/man/man1
RUN apt-get update
RUN sed 's/#.*//' /home/apt_requirements.txt | xargs apt-get install -y
RUN apt-get clean
RUN apt-get autoremove -y
#RUN rm -rf /var/lib/apt/lists/*

##### conda install ###
RUN conda update -n base -c defaults conda
RUN conda install -y -c conda-forge python-libsbml rdkit networkx==2.3 numpy pandas openbabel timeout-decorator cython
RUN conda install -y -c anaconda biopython==1.77
RUN conda install -y -c biobuilds t-coffee
RUN conda install -y -c bioconda emboss

RUN pip install equilibrator-pathway==0.3.1
#RUN rm -rf /usr/local/lib/python3/site-packages/ruamel*
RUN rm -rf $(dirname  $(which python))/../lib/python3.8/site-packages/ruamel*
RUN pip install cobra

###### MARVIN ####

COPY docker_files/rp2/marvin_linux_20.9.deb /home/
COPY docker_files/rp2/license.cxl /home/
ENV CHEMAXON_LICENSE_URL /home/license.cxl
RUN dpkg -i /home/marvin_linux_20.9.deb

#### extra install from source or from git ####

RUN git clone https://gitlab.irstea.fr/jacques.fize/GMatch4py.git
RUN cd GMatch4py && pip install . && cd ..

RUN git clone --single-branch --branch develop https://gitlab.com/equilibrator/equilibrator-api.git
RUN cd equilibrator-api && pip install -e . && cd ..

RUN git clone https://gitlab.com/equilibrator/equilibrator-assets.git
RUN cd equilibrator-assets && pip install -e . && cd ..

###############################
########## RETROPATH 2 ########
###############################

WORKDIR /home/

ENV DOWNLOAD_URL_RP2 https://download.knime.org/analytics-platform/linux/knime_4.2.2.linux.gtk.x86_64.tar.gz
ENV INSTALLATION_DIR_RP2 /usr/local
ENV KNIME_DIR $INSTALLATION_DIR_RP2/knime
ENV HOME_KNIME_DIR /home/knime

 # Download KNIME
RUN curl -L "$DOWNLOAD_URL_RP2" | tar vxz -C $INSTALLATION_DIR_RP2 \
    && mv $INSTALLATION_DIR_RP2/knime_* $INSTALLATION_DIR_RP2/knime

# Install Rserver so KNIME can communicate with R
RUN R -e 'install.packages(c("Rserve"), repos="http://cran.rstudio.com/")'

# Build argument for the workflow directory
ONBUILD ARG WORKFLOW_DIR="workflow/"
# Build argument for additional update sites
ONBUILD ARG UPDATE_SITES

# Create workflow directory and copy from host
ONBUILD RUN mkdir -p /payload
ONBUILD COPY $WORKFLOW_DIR /payload/workflow

# Create metadata directory
ONBUILD RUN mkdir -p /payload/meta

# Copy necessary scripts onto the image
COPY docker_files/rp2/docker_conf/getversion.py /scripts/getversion.py
COPY docker_files/rp2/docker_conf/listvariables.py /scripts/listvariables.py
COPY docker_files/rp2/docker_conf/listplugins.py /scripts/listplugins.py
COPY docker_files/rp2/docker_conf/run.sh /scripts/run.sh

# Let anyone run the workflow
RUN chmod +x /scripts/run.sh

# Add KNIME update site and trusted community update site that fit the version the workflow was created with
ONBUILD RUN full_version=$(python /scripts/getversion.py /payload/workflow/) \
&& version=$(python /scripts/getversion.py /payload/workflow/ | awk '{split($0,a,"."); print a[1]"."a[2]}') \
&& echo "http://update.knime.org/analytics-platform/$version" >> /payload/meta/updatesites \
&& echo "http://update.knime.org/community-contributions/trusted/$version" >> /payload/meta/updatesites \
# Add user provided update sites
&& echo $UPDATE_SITES | tr ',' '\n' >> /payload/meta/updatesites

# Save the workflow's variables in a file
ONBUILD RUN find /payload/workflow -name settings.xml -exec python /scripts/listplugins.py {} \; | sort -u | awk '!a[$0]++' > /payload/meta/features

ONBUILD RUN python /scripts/listvariables.py /payload/workflow

# Install required features
ONBUILD RUN "$KNIME_DIR/knime" -application org.eclipse.equinox.p2.director \
-r "$(cat /payload/meta/updatesites | tr '\n' ',' | sed 's/,*$//' | sed 's/^,*//')" \
-p2.arch x86_64 \
-profileProperties org.eclipse.update.install.features=true \
-i "$(cat /payload/meta/features | tr '\n' ',' | sed 's/,*$//' | sed 's/^,*//')" \
-p KNIMEProfile \
-nosplash

# Cleanup
ONBUILD RUN rm /scripts/getversion.py && rm /scripts/listvariables.py && rm /scripts/listplugins.py

############################### Workflow ##############################

ENV RETROPATH_VERSION 9
ENV RETROPATH_URL https://myexperiment.org/workflows/4987/download/RetroPath2.0_-_a_retrosynthesis_workflow_with_tutorial_and_example_data-v${RETROPATH_VERSION}.zip
ENV RETROPATH_SHA256 79069d042df728a4c159828c8f4630efe1b6bb1d0f254962e5f40298be56a7c4

# Download RetroPath2.0
WORKDIR /home/
RUN echo "$RETROPATH_SHA256 RetroPath2_0.zip" > RetroPath2_0.zip.sha256
RUN cat RetroPath2_0.zip.sha256
RUN echo Downloading $RETROPATH_URL
RUN curl -v -L -o RetroPath2_0.zip $RETROPATH_URL && sha256sum RetroPath2_0.zip && sha256sum -c RetroPath2_0.zip.sha256
RUN unzip RetroPath2_0.zip && mv RetroPath2.0/* /home/
RUN rm RetroPath2_0.zip

#install the additional packages required for running retropath KNIME workflow
RUN /usr/local/knime/knime -application org.eclipse.equinox.p2.director -nosplash -consolelog \
-r http://update.knime.org/community-contributions/trunk,\
http://update.knime.com/analytics-platform/4.2,\
http://update.knime.com/community-contributions/trusted/4.2 \
-i org.knime.features.chem.types.feature.group,\
org.knime.features.datageneration.feature.group,\
jp.co.infocom.cheminfo.marvin.feature.feature.group,\
org.knime.features.python.feature.group,\
org.rdkit.knime.feature.feature.group \
-bundlepool /usr/local/knime/ -d /usr/local/knime/

############################# Files and Tests #############################

COPY docker_files/rp2/callRP2.py /home/
COPY docker_files/rp2/rp2_sanity_test.tar.xz /home/

#test
ENV RP2_RESULTS_SHA256 7428ebc0c25d464fbfdd6eb789440ddc88011fb6fc14f4ce7beb57a6d1fbaec2
RUN tar xf /home/rp2_sanity_test.tar.xz -C /home/ 
RUN chmod +x /home/callRP2.py
RUN /home/callRP2.py -sinkfile /home/test/sink.csv -sourcefile /home/test/source.csv -rulesfile /home/test/rules.tar -rulesfile_format tar -max_steps 3 -scope_csv test_scope.csv
RUN echo "$RP2_RESULTS_SHA256 test_scope.csv" | sha256sum --check

############################################
############ RP2paths ######################
############################################

WORKDIR /home/

# Download and "install" rp2paths release
# Check for new versions from 
# https://github.com/brsynth/rp2paths/releases
ENV RP2PATHS_VERSION 1.0.2
ENV RP2PATHS_URL https://github.com/brsynth/rp2paths/archive/v${RP2PATHS_VERSION}.tar.gz
# NOTE: Update sha256sum for each release
ENV RP2PATHS_SHA256 3813460dea8bb02df48e1f1dfb60751983297520f09cdfcc62aceda316400e66
RUN echo "$RP2PATHS_SHA256  rp2paths.tar.gz" > rp2paths.tar.gz.sha256
RUN cat rp2paths.tar.gz.sha256
RUN echo Downloading $RP2PATHS_URL
RUN curl -v -L -o rp2paths.tar.gz $RP2PATHS_URL && sha256sum rp2paths.tar.gz && sha256sum -c rp2paths.tar.gz.sha256
RUN tar xfv rp2paths.tar.gz && mv rp2paths*/* /home/
RUN grep -q '^#!/' RP2paths.py || sed -i '1i #!/usr/bin/env python3' RP2paths.py

COPY docker_files/callRP2paths.py /home/

#############################################
######### RetroRules ########################
#############################################

WORKDIR home/

RUN wget https://retrorules.org/dl/preparsed/rr02/rp2/hs -O /home/rules_rall_rp2.tar.gz && \
    tar xf /home/rules_rall_rp2.tar.gz -C /home/ && \
    mv /home/retrorules_rr02_rp2_hs/retrorules_rr02_rp2_flat_forward.csv /home/rules_rall_rp2_forward.csv && \
    mv /home/retrorules_rr02_rp2_hs/retrorules_rr02_rp2_flat_retro.csv /home/rules_rall_rp2_retro.csv && \
    mv /home/retrorules_rr02_rp2_hs/retrorules_rr02_rp2_flat_all.csv /home/rules_rall_rp2.csv && \
    rm -r /home/retrorules_rr02_rp2_hs && \
    rm /home/rules_rall_rp2.tar.gz

#############################################
########## Equilibrator #####################
#############################################

COPY docker_files/init_equilibrator.py /home/
RUN chmod +x /home/init_equilibrator.py
RUN python /home/init_equilibrator.py