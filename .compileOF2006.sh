docker pull ithacafv/openfoam2006-muq2-pytorch
docker run -ti -d --name foam2006 -v "${PWD}":/home/ofuser/app:rw ithacafv/openfoam2006-muq2-pytorch /bin/bash
docker exec foam2006 /bin/bash -c "cd /home/ofuser/app; ls"
docker exec foam2006 /bin/bash -c "export MUQ_LIBRARIES=/home/Installations/MUQ_INSTALL; export MUQ_EXT_LIBRARIES=/home/Installations/MUQ_INSTALL/muq_external; export TORCH_LIBRARIES=/pytorch/torch; cd /home/ofuser/app; source /usr/lib/openfoam/openfoam2006/etc/bashrc; source etc/bashrc; ./Allwmake -taumq"
