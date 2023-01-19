FROM  opencfd/openfoam2106-dev
LABEL maintainer="moaadkhamlich@gmail.com"


# add enviromental variables and enable the default user
ARG USER=ithacafv
ARG of_var="source /usr/lib/openfoam/openfoam2106/etc/bashrc"
ARG ithaca_var="source /usr/lib/ITHACA-FV/etc/bashrc"
ENV USER $USER

# Create the user
RUN adduser --disabled-password --gecos '' $USER && \
    adduser $USER sudo; echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers && \
    usermod -a -G $USER $USER && \
    userdel sudofoam && \
    usermod -u 1000 ithacafv  && \
    groupmod -g 1000 ithacafv
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

RUN /bin/bash -c "source /usr/lib/openfoam/openfoam2106/etc/bashrc && \
    cd ITHACA-FV && git submodule update --init && source /etc/bash.bashrc && \
    ./Allwmake -tau -j 4";

USER root
RUN 	cp -r /home/ithacafv/OpenFOAM/ithacafv-v2106/platforms/linux64GccDPInt32Opt/bin/* /bin/
RUN 	cp -r /home/ithacafv/OpenFOAM/ithacafv-v2106/platforms/linux64GccDPInt32Opt/lib/* /lib/

#Update bashrc 
RUN  echo $of_var >> /etc/bash.bashrc
RUN  echo $ithaca_var >> /etc/bash.bashrc

#Source bashrc
USER $USER
RUN /bin/bash -c "source /etc/bash.bashrc"

USER $USER
WORKDIR $HOME
ENTRYPOINT ["/bin/bash"]

