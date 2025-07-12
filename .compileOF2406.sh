#!/bin/bash
docker pull shenhuiruan/ithaca2406
docker run -ti -d --name foam2406 -v "${PWD}":/home/foam:rw shenhuiruan/ithaca2406 /bin/bash
docker exec foam2406 /bin/bash -c "source /root/OpenFOAM/OpenFOAM-v2406/etc/bashrc; cd /home/foam; source etc/bashrc; git submodule update --init; ./Allwmake -taumq"
