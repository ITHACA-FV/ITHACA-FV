docker exec foam1812 /bin/sh -c "source /opt/OpenFOAM/setImage_v1812.sh"
docker exec foam1812 /bin/sh -c "cd /home/ofuser/app; ./Allwmake -j"