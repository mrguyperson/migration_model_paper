FROM ubuntu:jammy

SHELL ["/bin/bash", "-c"]
WORKDIR /project
COPY environment.yml ./

# RUN mkdir /usr/project/data

# update and install dependencies
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends git \
                                               gh \
                                               unzip \
                                               wget \
    && apt-get clean

# install mamba/conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash ./Mambaforge-Linux-x86_64.sh -bf -p /usr/mambaforge/ && \
    rm ./Mambaforge-Linux-x86_64.sh && \
    ln -s /usr/mambaforge/etc/profile.d/mamba.sh /etc/profile.d/mamba.sh && \
    echo ". /usr/mambaforge/etc/profile.d/mamba.sh" >> ~/.bashrc && \
    echo "mamba activate base" >> ~/.bashrc
ENV PATH="/usr/mambaforge/bin:${PATH}"

RUN mamba env update -n base -f environment.yml && \
    mamba init

