# Multi-architecture Dockerfile for ITHACA-FV
# Automatically selects the correct base image based on target platform
ARG TARGETARCH
FROM ithacafv/ithacafv-dependencies:${TARGETARCH} AS base

LABEL maintainer="giovannistabile@santannapisa.it"

# Install additional packages
RUN apt-get update && \
    apt-get install -y \
        git \
        vim \
        ssh \
        sudo \
        wget \
        software-properties-common && \
    rm -rf /var/lib/apt/lists/*

# Set working directory and clone ITHACA-FV
WORKDIR /usr/dir
RUN git clone https://github.com/ITHACA-FV/ITHACA-FV.git

# Environment variables for bashrc
ARG of_var="source \$FOAM_ETC/bashrc"
ARG ithaca_var="source /opt/ITHACA-FV/etc/bashrc"

# Update bashrc with OpenFOAM and ITHACA-FV sources
RUN echo $of_var >> /etc/bash.bashrc && \
    echo $ithaca_var >> /etc/bash.bashrc

# Update bashrc with OpenFOAM and ITHACA-FV sources
RUN echo $of_var >> /etc/bash.bashrc && \
    echo $ithaca_var >> /etc/bash.bashrc


RUN /bin/bash -c "source \$FOAM_ETC/bashrc && \
    cd ITHACA-FV && git submodule update --init && source etc/bashrc && \
    if [ \"$TARGETARCH\" = \"amd64\" ]; then \
        ./Allwmake -taumq -j 4; \
    else \
        ./Allwmake -au -j 4; \
    fi"


# Source bashrc on container start
RUN /bin/bash -c "source /etc/bash.bashrc"
ENTRYPOINT ["/bin/bash"]
