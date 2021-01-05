FROM ubuntu:14.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y \
		vim \
    sudo \
		wget \
    make \
		software-properties-common ;\
		rm -rf /var/lib/apt/lists/*

#RUN useradd --user-group --create-home --shell /bin/bash root ;\
#	echo "root ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

RUN sh -c "wget -O - http://dl.openfoam.org/gpg.key | apt-key add -" ;\
  add-apt-repository http://dl.openfoam.org/ubuntu ;\
  apt-get update ;\
  apt-get install -y \
    openfoam30 \
    paraviewopenfoam44 ;\
  echo "source /opt/openfoam30/etc/bashrc" >> ~root/.bashrc

USER root