FROM --platform=linux/amd64 melclic/biopathopt:latest

RUN conda install -n biopathopt -c conda-forge rdkit networkx

RUN mkdir /home/MetaXime/
ADD metaxime/ /home/MetaXime/metaxime/
COPY setup.py /home/MetaXime/
COPY requirements.txt /home/MetaXime/
COPY README.md /home/MetaXime/
COPY scripts/run_pipeline.py /home/MetaXime/

WORKDIR /home/MetaXime/

RUN conda run -n biopathopt pip install -e .