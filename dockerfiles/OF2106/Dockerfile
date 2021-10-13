# Start from the ubuntu Openfoam 2106 image
FROM ithacafv/openfoam2106

USER root
ARG PYTHON_VERSION=3.7

RUN rm /etc/apt/sources.list.d/openfoam.list && \
    cp /etc/apt/sources.list /etc/apt/sources.list.backup && \
    grep -v -e "openfoam" /etc/apt/sources.list.backup > /etc/apt/sources.list && \
    echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections && \
    apt-get update && \
    apt-get install -yy -q pwgen npm nodejs cmake git wget bzip2 && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Anaconda installing
RUN git clone --recursive https://github.com/pytorch/pytorch && \
cd pytorch && git submodule sync && git submodule update --init --recursive --jobs 0 && \
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh && \
bash Miniconda3-py39_4.10.3-Linux-x86_64.sh -b && \
rm Miniconda3-py39_4.10.3-Linux-x86_64.sh && \
. /root/miniconda3/etc/profile.d/conda.sh && \  
export PATH=/root/miniconda3/bin:$PATH && \
conda install -y numpy ninja pyyaml mkl mkl-include setuptools cmake cffi typing && \
conda install -y -c pytorch magma-cuda90 && \
pip install typing_extensions && \
export CMAKE_PREFIX_PATH=${CONDA_PREFIX:-"$(dirname $(which conda))/../"} && \
python setup.py install && \
cd .. && \
rm -r pytorch && \
conda install -y -c conda-forge muq cmake && \
conda clean -y --all
ENV TORCH_LIBRARIES=/root/miniconda3/lib/python3.9/site-packages/torch
ENV MUQ_LIBRARIES=/root/miniconda3
RUN echo 'source /home/foam/.bashrc' >> ~/.bashrc 

