FROM ithacafv/openfoam2106:manifest-arm64
LABEL maintainer="moaadkhamlich@gmail.com"

RUN apt-get update \
	&& apt-get install -y \
        git \
		vim \
		ssh \
		sudo\
		wget\
		software-properties-common ;\
		rm -rf /var/lib/apt/lists/*

# add and enable the default user
ARG USER=ithacafv
RUN adduser --disabled-password --gecos '' $USER
RUN adduser $USER sudo; echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

RUN mv /root/* /opt/
WORKDIR /opt/
RUN git clone https://github.com/mathLab/ITHACA-FV.git
RUN chown -R $USER:$USER OpenFOAM-v2106 OpenFOAM-v2106 ITHACA-FV
RUN chown -R $USER:$USER /home/$USER
USER $USER
ENV HOME /home/$USER
ENV USER $USER

RUN /bin/bash -c "source OpenFOAM-v2106/etc/bashrc && \
    cd ITHACA-FV && source etc/bashrc && \
    ./Allwmake -au -j 4";


WORKDIR $HOME
COPY .bashrc $HOME/.bashrc
COPY .vimrc $HOME/.
ENTRYPOINT ["/bin/bash"]
