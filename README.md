# Discrete Adjoint Approach to the Spalart-Allmaras Turbulence Model for Incompressible Flow in OpenFOAM

Masterthesis in [Computational Mechanics](https://www.ed.tum.de/en/ed/studies/degree-program/computational-mechanics-m-sc/) at [Technical University Munich](https://www.tum.de/).

All examples have been tested with [OpenFoam-3.0.1](https://openfoam.org/version/3-0-1/) on [Ubuntu 14.04.6 LTS (Trusty Tahr)](https://releases.ubuntu.com/14.04/) . For convenience a docker image and start script is provided.

## Build Instructions

OS Manjaro Linux, Docker version 20.10.1

Clone this repository an open the root location in a terminal.

The following command will build a local Docker image called `openfoam/openfoam301`
```bash
$ docker build -t openfoam/openfoam301 .
```

## Getting Started

Once build, start the container via the startup script
```bash
$ ./openfoam301-ubuntu
```
Which will start an interactive session in your Docker container with your current working directory on your local filesystem mounted as /home/openfoam/ in your container. The container will be removed once you exit it.

Compile the adjoint solver
```
# wmake /home/openfoam/applications/solvers/incompressible/adjointSpalartAllmaras/
```

Test the solver
```
# adjointSpalartAllmaras -help
```

Expected output
```bash
Usage: adjointSpalartAllmaras [OPTIONS]
options:
  -case <dir>       specify alternate case directory, default is the cwd
  -noFunctionObjects
                    do not execute functionObjects
  -parallel         run in parallel
  -roots <(dir1 .. dirN)>
                    slave root directories for distributed running
  -srcDoc           display source code in browser
  -doc              display application documentation in browser
  -help             print the usage

Using: OpenFOAM-3.0.1 (see www.OpenFOAM.org)
Build: 3.0.1-119cac7e8750
```

## Example

Follow [thesis](./thesis/thesis_dennis_kasper.pdf) chapter six (p. 25) for the unit cube example.

In a terminal navigate to `/home/openfoam/run/unitCubeSpalartAllmaras/` and call the solver
```
# adjointSpalartAllmaras
```

## Reference

[1] Docker build [CFDEngine](https://github.com/CFDEngine/simple-openfoam-dockerfile/)  
[2] Startscript [openfoam4-linux](http://dl.openfoam.org/docker/openfoam4-linux) and [Install guide](https://openfoam.org/download/8-linux/)  
[3] [Master Thesis](./thesis/thesis_dennis_kasper.pdf)

## ToDo

- [ ] Add NACA 2412 Airfoil and ONERA M6 Wing examples
- [ ] fix solver/lib output path -> /home/openfoam/ not /root/OpenFOAM/-3.0.1/
- [ ] fix paraFOAM GUI error
- [ ] compile solvers in docker build process

