FROM ubuntu:18.04

WORKDIR /home/

###### copy all the necessary files #####

COPY license.cxl /home/
ENV CHEMAXON_LICENSE_URL /home/license.cxl
COPY init_equilibrator.py /home/

### install all the requirements ###

RUN apt-get update
RUN sed 's/#.*//' apt_requirements.txt | xargs apt-get install
RUN apt-get clean
RUN apt-get autoremove
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /usr/local/lib/python3.8/site-packages/ruamel*
RUN pip install -r pip_requirements.txt

#### extra install from source or from git ####

RUN git clone https://gitlab.irstea.fr/jacques.fize/GMatch4py.git
RUN cd GMatch4py && pip install . && cd ..

RUN git clone --single-branch --branch develop https://gitlab.com/equilibrator/equilibrator-api.git
RUN cd equilibrator-api && pip install -e . && cd ..

#equilibrator-assets
RUN git clone https://gitlab.com/equilibrator/equilibrator-assets.git
RUN cd equilibrator-assets && pip install -e . && cd ..

##############################
######### RDkit ##############
##############################

ARG RDKIT_VERSION=Release_2020_03_2
RUN wget --quiet https://github.com/rdkit/rdkit/archive/${RDKIT_VERSION}.tar.gz \
 && tar -xzf ${RDKIT_VERSION}.tar.gz \
 && mv rdkit-${RDKIT_VERSION} rdkit \
 && rm ${RDKIT_VERSION}.tar.gz

RUN mkdir /rdkit/build
WORKDIR /rdkit/build

RUN cmake -Wno-dev \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_INSTALL_PREFIX=/usr \
    -D Boost_NO_BOOST_CMAKE=ON \
    -D PYTHON_EXECUTABLE=/usr/bin/python3 \
    -D RDK_BUILD_AVALON_SUPPORT=ON \
    -D RDK_BUILD_CAIRO_SUPPORT=ON \
    -D RDK_BUILD_CPP_TESTS=OFF \
    -D RDK_BUILD_INCHI_SUPPORT=ON \
    -D RDK_BUILD_FREESASA_SUPPORT=ON \
    -D RDK_INSTALL_INTREE=OFF \
    -D RDK_INSTALL_STATIC_LIBS=OFF \
    ..

# Copy rdkit installation from rdkit-build-env
COPY --from=rdkit-build-env /usr/lib/libRDKit* /usr/lib/
COPY --from=rdkit-build-env /usr/lib/cmake/rdkit/* /usr/lib/cmake/rdkit/
COPY --from=rdkit-build-env /usr/share/RDKit /usr/share/RDKit
COPY --from=rdkit-build-env /usr/include/rdkit /usr/include/rdkit
COPY --from=rdkit-build-env /usr/lib/python3/dist-packages/rdkit /usr/lib/python3/dist-packages/rdkit

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
COPY rp2_docker_conf/getversion.py /scripts/getversion.py
COPY rp2_docker_conf/listvariables.py /scripts/listvariables.py
COPY rp2_docker_conf/listplugins.py /scripts/listplugins.py
COPY rp2_docker_conf/run.sh /scripts/run.sh

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

COPY rpTool.py /home/
COPY galaxy/code/tool_RetroPath2.py /home/
COPY rp2_sanity_test.tar.xz /home/

#test
ENV RP2_RESULTS_SHA256 7428ebc0c25d464fbfdd6eb789440ddc88011fb6fc14f4ce7beb57a6d1fbaec2
RUN tar xf /home/rp2_sanity_test.tar.xz -C /home/ 
RUN /home/tool_RetroPath2.py -sinkfile /home/test/sink.csv -sourcefile /home/test/source.csv -rulesfile /home/test/rules.tar -rulesfile_format tar -max_steps 3 -scope_csv test_scope.csv
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

COPY rpTool.py /home/
COPY galaxy/code/tool_rp2paths.py /home/

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

