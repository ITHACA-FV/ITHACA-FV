docker pull ithacafv/openfoam5-muq2-pytorch:v0
docker run -ti -d --name foam5 -v "${PWD}":/home/ofuser/app:rw ithacafv/openfoam5-muq2-pytorch /bin/bash
docker exec foam5 /bin/bash -c "cd /home/ofuser/app; ls"
docker exec foam5 /bin/bash -c "export MUQ_LIBRARIES=/home/Installations/MUQ_INSTALL; export MUQ_EXT_LIBRARIES=/home/Installations/MUQ_INSTALL/muq_external; export TORCH_LIBRARIES=~/pytorch/torch; cd /home/ofuser/app; source /opt/openfoam5/etc/bashrc; source etc/bashrc; ./Allwmake -taumq"
