# Multi-architecture Dockerfile for ITHACA-FV
# Automatically selects the correct base image based on target platform
ARG TARGETARCH
FROM ithacafv/openfoam2412-muq2-pytorch:${TARGETARCH} AS base

LABEL maintainer="moaadkhamlich@gmail.com"

USER root

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

# Create ithacafv user and group
ARG USER=ithacafv
RUN if id -u 1000 >/dev/null 2>&1; then \
        userdel $(id -nu 1000) || true; \
    fi && \
    if getent group 1000 >/dev/null 2>&1; then \
        groupdel $(getent group 1000 | cut -d: -f1) || true; \
    fi && \
    adduser --disabled-password --gecos '' --uid 1000 $USER && \
    adduser $USER sudo && \
    echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Set environment variables
ENV HOME=/home/$USER
ENV USER=$USER

# Set working directory and clone ITHACA-FV
WORKDIR /opt
RUN git clone https://github.com/mathLab/ITHACA-FV.git && \
    chown -R $USER:$USER ITHACA-FV && \
    chown -R $USER:$USER /home/$USER

# Environment variables for bashrc
ARG of_var="source /usr/lib/openfoam/openfoam2412/etc/bashrc"
ARG ithaca_var="source /opt/ITHACA-FV/etc/bashrc"

# Update bashrc with OpenFOAM and ITHACA-FV sources
RUN echo $of_var >> /etc/bash.bashrc && \
    echo $ithaca_var >> /etc/bash.bashrc

# Switch to ithacafv user
USER $USER

# Build ITHACA-FV
RUN /bin/bash -c "source /usr/lib/openfoam/openfoam2412/etc/bashrc && \
    cd ITHACA-FV && \
    git submodule update --init && \
    source etc/bashrc && \
    ./Allwmake -au -j 4"

# Copy binaries and libraries to system paths (as root)
USER root
RUN if [ -d "/home/$USER/OpenFOAM/$USER-openfoam2412/platforms/linux64GccDPInt32Opt/bin" ]; then \
        cp -r /home/$USER/OpenFOAM/$USER-openfoam2412/platforms/linux64GccDPInt32Opt/bin/* /usr/local/bin/ || true; \
        cp -r /home/$USER/OpenFOAM/$USER-openfoam2412/platforms/linux64GccDPInt32Opt/lib/* /usr/local/lib/ || true; \
    fi

# Final setup
USER $USER
WORKDIR $HOME

# Source bashrc on container start
RUN /bin/bash -c "source /etc/bash.bashrc"

ENTRYPOINT ["/bin/bash"]