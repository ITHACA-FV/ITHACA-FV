docker pull ithacafv/openfoam1812-muq2-pytorch
docker run -ti -d --name foam1812 -v "${PWD}":/home/ofuser/app:rw ithacafv/openfoam1812-muq2-pytorch /bin/bash
docker exec foam1812 /bin/bash -c "cd /home/ofuser/app; ls"
docker exec foam1812 /bin/bash -c "export MUQ_LIBRARIES=/home/Installations/MUQ_INSTALL; export MUQ_EXT_LIBRARIES=/home/Installations/MUQ_INSTALL/muq_external; export TORCH_LIBRARIES=/pytorch/torch; cd /home/ofuser/app; source /root/OpenFOAM/OpenFOAM-v1812/etc/bashrc; source etc/bashrc; ./Allwmake -taumq"

