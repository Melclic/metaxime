mkdir -p $(pwd)/mx-results
docker build -t metaxime -f Dockerfile .
docker run -it -v $(pwd)/mx-results:/home/mx-results -p 80:80 -p 8888:8888 metaxime

