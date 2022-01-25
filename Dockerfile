# getting base image ubuntu
from ubuntu 

MAINTAINER Purusharth Saxena <purusharth.saxena@stud.uni-heidelberg.de>

RUN useradd -ms /bin/bash hiflow
USER hiflow

RUN mkdir /home/hiflow3

WORKDIR /home/hiflow3

RUN apt-get update

RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata


RUN apt-get update; \
    apt-get install -y wget build-essential libopenmpi-dev libtinyxml2-dev cmake cmake-curses-gui libboost-all-dev libmetis-dev; \
    wget --no-check-certificate https://emcl-gitlab.iwr.uni-heidelberg.de/hiflow3.org/hiflow3/-/archive/v3_release/hiflow3-v3_release.tar.gz -O - | tar -xz;   \ 
    cd hiflow3-v3_release; \
    mkdir build; \
    cd build; \
    cmake .. -DHIFLOW_CPP_ISO_STANDARD=c++17; \
    make -j$(nproc)





