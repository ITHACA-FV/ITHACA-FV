#!/bin/bash
docker pull ithacafv/openfoam2012-muq2-pytorch
docker run -ti -d --name foam2012 -v "${PWD}":/home/ofuser/app:rw ithacafv/openfoam2012-muq2-pytorch /bin/bash
docker exec foam2012 /bin/bash -c "cd /home/ofuser/app; ls"
docker exec foam2012 /bin/bash -c "export MUQ_LIBRARIES=/home/Installations/MUQ_INSTALL; export MUQ_EXT_LIBRARIES=/home/Installations/MUQ_INSTALL/muq_external; export TORCH_LIBRARIES=/pytorch/torch; cd /home/ofuser/app; source /usr/lib/openfoam/openfoam2012/etc/bashrc; source etc/bashrc; ./Allwmake -taumq"
