FROM brsynth/rpbase

RUN apt-get install --quiet --yes --no-install-recommends \
			libxext6  \
    	libxrender-dev \
	 && conda install -y -c rdkit rdkit

COPY rpCache.py /home/

ENTRYPOINT ["python"]
CMD ["/home/rpCache.py"]

# Open server port
EXPOSE 8997
