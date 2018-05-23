FROM neurodebian:xenial-non-free
MAINTAINER manjebrinkhuis@gmail.com

USER root

# Create user and setup directories
RUN useradd -m -s /bin/bash neuro

RUN mkdir /data \
  && mkdir /opt/conda \
  && mkdir /opt/freesurfer \
  && chown -R neuro:neuro /opt/conda

# Dependencies
RUN apt-get update && apt-get -yq dist-upgrade \
  && apt-get install -yq --no-install-recommends \
    wget \
    gnupg \
    bzip2 \
    ca-certificates \
    locales \
    fonts-liberation \
    libgl1-mesa-dri \
    libgl1-mesa-glx \
    libglu1-mesa \
    libxi-dev \
    libxmu-dev \
    libglu1-mesa-dev \
    libjpeg62-dev \
    fsl-core \
    fslview \
    graphviz \
    libgraphviz-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Lanuage
RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen \
  && locale-gen

ENV LC_ALL=en_US.UTF-8 \
  LANG=en_US.UTF-8 \
  LANGUAGE=en_US.UTF-8

EXPOSE 8888

WORKDIR /home/neuro
USER neuro

ENV FSLDIR=/usr/share/fsl/5.0 \
  FREESURFER_HOME=/opt/freesurfer \
  CONDA_PATH=/opt/conda \
  CONDA_FNAME=Miniconda3-latest-Linux-x86_64.sh \
  CONDA_URL=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  SHELL=/bin/bash \
  USER=neuro

RUN bash ${FSLDIR}/etc/fslconf/fsl.sh
ENV PATH=${FSLDIR}/bin:${PATH}

RUN wget ${CONDA_URL} -P /tmp/ \
  && bash /tmp/${CONDA_FNAME} -b -u -p ${CONDA_PATH}

ENV PATH=${CONDA_PATH}/bin:${PATH}
ENV CONDA_ENVS_PATH=/home/neuro/.conda

# Jupyterlab + environment
RUN conda install -y jupyterlab \
  && conda config --add channels conda-forge

# R
RUN conda install -y -c r r-base r-irkernel

# Neuro
RUN conda create -y -n neuro \
    python=3 \
    pandas \
    numpy \
    scipy \
    rpy2 \
    nipype \
    nibabel \
    ipykernel
RUN /home/neuro/.conda/neuro/bin/python -m ipykernel \
    install --user --name neuro

# Clean up
RUN conda clean -tipsy

RUN mkdir /home/neuro/.jupyter
ADD jupyter_notebook_config.py /home/neuro/.jupyter/jupyter_notebook_config.py

ADD setup.sh /home/neuro/setup.sh
ENTRYPOINT [ "/home/neuro/setup.sh" ]
