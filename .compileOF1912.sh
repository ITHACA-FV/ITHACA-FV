docker pull openfoamplus/of_v1912_centos73
docker run -ti -d --name foam1912 -v ${PWD}:/home/ofuser/app:rw openfoamplus/of_v1912_centos73 /bin/sh
docker exec foam1912 /bin/sh -c "cd /home/ofuser/app; ls"
docker exec foam1912 /bin/sh -c "source /opt/OpenFOAM/setImage_v1912.sh; cd /home/ofuser/app; source etc/bashrc; ./Allwmake -tau"
