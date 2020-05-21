sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 6C0DAC728B29D817
sudo apt-get update --fix-missing
sudo add-apt-repository http://dl.openfoam.org/ubuntu
sudo apt-get update --fix-missing
sudo apt-get -y install -qq openfoam5 --allow-unauthenticated
