#!/bin/bash
docker pull shenhuiruan/ithaca2312
docker run -ti -d --name foam2312 -v "${PWD}":/home/foam:rw shenhuiruan/ithaca2312 /bin/bash
docker exec foam2312 /bin/bash -c "
  git config --global --add safe.directory /home/foam && \
  source /root/OpenFOAM/OpenFOAM-v2312/etc/bashrc && \
  cd /home/foam && \
  source etc/bashrc && \
  git submodule update --init && \
  ./Allwmake -taumq
"