# MetaXime

Metabolic engineering project that provides different analysis of heterologous pathways

TODO: consider seperating the single into threee dockers: 1) RP2+rp2paths+RetroRules 2) MetaXime 3) MetaXime-frontend  

docker build -t metaxime -f Dockerfile .
docker run -it -p 80:80 --entrypoint /bin/bash metaxime
docker run -it -p 80:80 metaxime
docker run -it -p 80:80 -p 8888:8888 metaxime
docker run -it -v $(pwd)/mx-results:/home/mx-results -p 80:80 -p 8888:8888 metaxime



docker build -t metaxime --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) --build-arg USER_NAME=$(id -u -n) -f Dockerfile .
docker run -it --rm -v $(pwd)/mx-results:/home/mx-results --user "$(id -u):$(id -g)" -p 80:80 -p 8888:8888 metaxime

docker build -t metaxime --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) -f Dockerfile .
docker run -it --rm -v $(pwd)/mx-results:/home/mx-results -p 80:80 -p 8888:8888 metaxime


docker run -it -v $(pwd)/mx-results:/mx-results -p 80:80 -p 8888:8888 metaxime
