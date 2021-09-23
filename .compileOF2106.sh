#!/bin/bash
docker pull ithacafv/openfoam2106-muq2-pytorch
docker run -ti -d --name foam2106 -v "${PWD}":/home/ofuser/app:rw ithacafv/openfoam2106-muq2-pytorch /bin/bash
docker exec foam2106 /bin/bash -c "source /usr/lib/openfoam/openfoam2106/etc/bashrc; cd /home/ofuser/app; source etc/bashrc; ./Allwmake -taumq"
