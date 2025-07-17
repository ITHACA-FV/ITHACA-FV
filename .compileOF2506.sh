#!/bin/bash
docker pull ithacafv/openfoam2506-muq2-pytorch
docker run -ti -d --name foam2506 -v "${PWD}":/home/ofuser/app:rw ithacafv/openfoam2506-muq2-pytorch /bin/bash
docker exec foam2506 /bin/bash -c "source /usr/lib/openfoam/openfoam2506/etc/bashrc; cd /home/ofuser/app; git config --global --add safe.directory /home/ofuser/app; source etc/bashrc; git submodule update --init; ./Allwclean; ./Allwmake -taumq"
