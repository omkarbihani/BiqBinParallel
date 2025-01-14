FROM ubuntu:24.04 

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    make \
    liblapack-dev \
    liblapack3 \
    libopenblas-dev \
    mpich \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /parallel_biqbin_maxcut

COPY . /parallel_biqbin_maxcut

CMD ["bash"]