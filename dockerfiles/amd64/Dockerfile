FROM  opencfd/openfoam2106-dev
LABEL maintainer="moaadkhamlich@gmail.com"


# add and enable the default user
ARG USER=ithacafv
ENV USER $USER

# Create the user
RUN adduser --disabled-password --gecos '' $USER && \
    adduser $USER sudo; echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers && \
    usermod -a -G $USER $USER
#make sure everything is in place
ENV HOME /home/$USER
RUN chown -R $USER:$USER /home/$USER
RUN chown -R $USER:$USER /home/openfoam && rm -r /home/openfoam && \
    chown -R $USER:$USER /home/sudofoam && rm -r /home/sudofoam

USER root
RUN apt-get update &&\
    apt-get install -y \
        git \
        vim \
        ssh \
        sudo\
        wget\
        software-properties-common ;\
        sudo rm -rf /var/lib/apt/lists/*

WORKDIR /usr/lib
RUN git clone https://github.com/mathLab/ITHACA-FV.git
RUN chown -R $USER:$USER openfoam ITHACA-FV
RUN chown -R $USER:$USER /home/$USER
USER $USER

RUN /bin/bash -c "source openfoam/openfoam2106/etc/bashrc && \
    cd ITHACA-FV && source etc/bashrc && \
    ./Allwmake -au -j 4";

USER $USER
WORKDIR $HOME
COPY .bashrc $HOME/.bashrc
ENTRYPOINT ["/bin/bash"]
