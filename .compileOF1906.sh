docker pull openfoamplus/of_v1906_centos73
docker run -ti -d --name foam1906 -v ${PWD}:/home/ofuser/app:rw openfoamplus/of_v1906_centos73 /bin/sh
docker exec foam1906 /bin/sh -c "cd /home/ofuser/app; ls"
docker exec foam1906 /bin/sh -c "source /opt/OpenFOAM/setImage_v1906.sh; cd /home/ofuser/app; source etc/bashrc; ./Allwmake -tau"
