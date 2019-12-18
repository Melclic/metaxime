FROM conda/miniconda3

WORKDIR /home/

RUN apt-get --quiet update && \
    #apt-get --quiet --yes dist-upgrade && \
    apt-get install --quiet --yes --no-install-recommends \
    ca-certificates \
    build-essential \
    cmake \
    git \
    wget \
    xz-utils && \
    pip install --upgrade pip && \ 
    conda update -n base -c defaults conda && \
    conda install -y -c SBMLTeam python-libsbml

COPY rpSBML.py /home/
