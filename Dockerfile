#FROM ubuntu:22.04 
FROM python:3.12 

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    make \
    liblapack-dev \
    liblapack3 \
    libopenblas-dev \
    mpich \
    && rm -rf /var/lib/apt/lists/*

#RUN apt install python3.12

RUN pip install --upgrade pip
RUN pip install numpy


WORKDIR /boost
#COPY boost_1_88_0.tar.gz .

RUN wget https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.gz
RUN tar xzvf boost_1_88_0.tar.gz

WORKDIR /boost/boost_1_88_0

RUN ./bootstrap.sh --prefix=/usr/
RUN ./b2 -q
RUN ./b2 install

#RUN apt-get install python3-dev

#RUN apt-get install -y build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev libboost-all-dev 

#RUN apt-get install libboost-all-dev


# RUN chmod -R 777 .

RUN pip install dwave-neal

WORKDIR /solver
COPY . . 



# USER 1001:1001

ENV PYBOOST="-I/usr/local/include/python3.12 -I/usr/include -L/usr/lib -lboost_python312 -lpython3.12 -lboost_numpy312"

RUN make clean
RUN make
