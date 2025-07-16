#!/bin/bash
docker pull ithacafv/openfoam2412-muq2-pytorch
docker run -ti -d --name foam2412 -v "${PWD}":/home/ofuser/app:rw ithacafv/openfoam2412-muq2-pytorch /bin/bash
docker exec foam2412 /bin/bash -c "source /usr/lib/openfoam/openfoam2412/etc/bashrc; cd /home/ofuser/app; git config --global --add safe.directory /home/ofuser/app; source etc/bashrc; git submodule update --init; ./Allwclean; ./Allwmake -taumq"
