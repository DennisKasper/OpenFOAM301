FROM ubuntu:14.04

RUN : \
  && apt-get update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  wget \
  make \
  software-properties-common \
  gnuplot \
  && rm -rf /var/lib/apt/lists/*

RUN : \
  && useradd --user-group --create-home --shell /bin/bash foam \
  && echo "foam ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

RUN : \
  && sh -c "wget -O - http://dl.openfoam.org/gpg.key | apt-key add -" \
  && add-apt-repository http://dl.openfoam.org/ubuntu \
  && apt-get update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  openfoam30 \
  && rm -rf /var/lib/apt/lists/* \
  && echo "source /opt/openfoam30/etc/bashrc" >> ~foam/.bashrc

COPY --chown=foam ./applications /home/foam/applications

USER foam

RUN : \
  && . /opt/openfoam30/etc/bashrc \
  && wmake /home/foam/applications/solvers/incompressible/adjointLaplace \
  && wmake /home/foam/applications/solvers/incompressible/adjointSpalartAllmaras
