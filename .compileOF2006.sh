docker pull openfoamplus/of_v2006_centos73
docker run -ti -d --name foam2006 -v ${PWD}:/home/ofuser/app:rw openfoamplus/of_v2006_centos73 /bin/sh
docker exec foam2006 /bin/sh -c "cd /home/ofuser/app; ls"
docker exec foam2006 /bin/sh -c "source /opt/OpenFOAM/setImage_v2006.sh; cd /home/ofuser/app; source etc/bashrc; ./Allwmake -tau"
