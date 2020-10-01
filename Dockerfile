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

RUN conda install -c anaconda cython
RUN pip install networkx==2.3 numpy pandas

RUN git clone https://gitlab.irstea.fr/jacques.fize/GMatch4py.git
RUN cd GMatch4py && pip install . && cd ..

COPY rpSBML.py /home/
COPY rpGraph.py /home/
COPY rpMerge.py /home/
COPY rpDraw.py /home/
COPY data/ /home/

ENV PYTHONPATH="/home"
