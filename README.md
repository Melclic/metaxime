# MetaXime

Metabolic engineering project that provides different analysis of heterologous pathways


docker build -t metaxime -f Dockerfile .
docker run -it -p 80:80 --entrypoint /bin/bash metaxime
docker run -it -p 80:80 metaxime
