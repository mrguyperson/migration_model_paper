FROM fedora:36

# WORKDIR /usr/project
# SHELL ["/bin/bash", "-c"]
# COPY code/ /usr/project/code
# COPY envs/ /usr/project/envs
# COPY Snakefile /usr/project/.


# RUN mkdir /usr/project/data

# update and install dependencies
RUN dnf update -y && \
    dnf install -y wget unzip which git

# install mamba/conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh  && \
    bash ./Mambaforge-Linux-x86_64.sh -bf -p /usr/mambaforge/ && \
    rm ./Mambaforge-Linux-x86_64.sh
ENV PATH /usr/mambaforge/bin:$PATH
# RUN mamba env create -f environment.yml
WORKDIR /usr/project
COPY environment.yml .
RUN mamba env create -f environment.yml

RUN echo "source activate analysis" > ~/.bashrc



# RUN mamba create -c conda-forge -c bioconda -n snakemake snakemake
# RUN echo "source activate snakemake" > ~/.bashrc
# ENV PATH /opt/mambaforge/envs/env/bin:$PATH
# RUN curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O && \
#     bash Miniconda3-4.7.12.1-Linux-x86_64.sh -bf -p /usr/miniconda3/ && \
#     rm Miniconda3-4.7.12.1-Linux-x86_64.sh && \
#     /usr/miniconda3/bin/conda clean -tipsy && \
#     ln -s /usr/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
#     echo ". /usr/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
#     echo "conda activate base" >> ~/.bashrc

#     # Add conda to PATH and set locale
# ENV PATH="/usr/miniconda3/bin:${PATH}"
# ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# # Configure Conda channels and install Mamba
# RUN conda config --add channels bioconda \
#     && conda config --add channels conda-forge \
#     # && conda config --set channel_priority strict \
#     && conda install mamba \
#     && mamba clean --all

# CMD /bin/bash
# FROM fedora:36

# ENV CONDA_DIR=/opt/conda
# ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
# ENV PATH=${CONDA_DIR}/bin:${PATH}

# # 1. Install just enough for conda to work
# # 2. Keep $HOME clean (no .wget-hsts file), since HSTS isn't useful in this context
# # 3. Install miniforge from GitHub releases
# # 4. Apply some cleanup tips from https://jcrist.github.io/conda-docker-tips.html
# #    Particularly, we remove pyc and a files. The default install has no js, we can skip that
# # 5. Activate base by default when running as any *non-root* user as well
# #    Good security practice requires running most workloads as non-root
# #    This makes sure any non-root users created also have base activated
# #    for their interactive shells.
# # 6. Activate base by default when running as root as well
# #    The root user is already created, so won't pick up changes to /etc/skel
# RUN dnf update -y && \
#     dnf install -y \
#         wget bzip2 ca-certificates \
#         git \
#         tini \
#         which && \
#     # dnf clean && \
#     rm -rf /var/lib/apt/lists/* && \
#     wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -O /tmp/mambaforge.sh && \
#     /bin/bash /tmp/mambaforge.sh -b -p ${CONDA_DIR} && \
#     rm /tmp/mambaforge.sh && \
#     conda clean --tarballs --index-cache --packages --yes && \
#     find ${CONDA_DIR} -follow -type f -name '*.a' -delete && \
#     find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete && \
#     conda clean --force-pkgs-dirs --all --yes  && \
#     echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc && \
#     echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

# ENTRYPOINT ["tini", "--"]
# CMD [ "/bin/bash" ]