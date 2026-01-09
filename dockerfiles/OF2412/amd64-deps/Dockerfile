# Start from the ubuntu Openfoam 2106 image
FROM opencfd/openfoam-dev:2412
USER root
ARG PYTHON_VERSION=3.7
ENV PATH="/root/miniconda3/bin:${PATH}"

RUN rm /etc/apt/sources.list.d/openfoam.list && \
    cp /etc/apt/sources.list /etc/apt/sources.list.backup && \
    grep -v -e "openfoam" /etc/apt/sources.list.backup > /etc/apt/sources.list && \
    echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections && \
    apt-get update && \
    apt-get install -yy -q pwgen npm nodejs cmake git wget bzip2 unzip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Anaconda installing
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
bash Miniconda3-latest-Linux-x86_64.sh -b && \
rm Miniconda3-latest-Linux-x86_64.sh && \
. /root/miniconda3/etc/profile.d/conda.sh && \  
export PATH=/root/miniconda3/bin:$PATH && \
wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.7.1%2Bcpu.zip && \
unzip libtorch-cxx11-abi-shared-with-deps-2.7.1+cpu.zip && \
rm libtorch-cxx11-abi-shared-with-deps-2.7.1+cpu.zip && \
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
conda tos accept --override-channels --channel  https://repo.anaconda.com/pkgs/r && \
conda install -y -c conda-forge muq cmake pybind11 && \
conda clean -y --all
ENV TORCH_LIBRARIES=/libtorch
ENV MUQ_LIBRARIES=/root/miniconda3

# Export OpenFOAM environment variables for use in subsequent builds
RUN /bin/bash -c "source /usr/lib/openfoam/openfoam2412/etc/bashrc && \
    echo 'export FOAM_ETC=${FOAM_ETC}' >> /etc/environment && \
    echo 'export FOAM_APP=${FOAM_APP}' >> /etc/environment && \
    echo 'export FOAM_APPBIN=${FOAM_APPBIN}' >> /etc/environment && \
    echo 'export FOAM_LIBBIN=${FOAM_LIBBIN}' >> /etc/environment && \
    echo 'export FOAM_SRC=${FOAM_SRC}' >> /etc/environment && \
    echo 'export WM_PROJECT_DIR=${WM_PROJECT_DIR}' >> /etc/environment"

ENV FOAM_ETC=/usr/lib/openfoam/openfoam2412/etc
RUN echo 'source /usr/lib/openfoam/openfoam2412/etc/bashrc' >> ~/.bashrc