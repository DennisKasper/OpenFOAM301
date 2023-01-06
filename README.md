# Discrete Adjoint Approach to the Spalart-Allmaras Turbulence Model for Incompressible Flow in OpenFOAM

Masterthesis in [Computational Mechanics](https://www.ed.tum.de/en/ed/studies/degree-program/computational-mechanics-m-sc/) at [Technical University Munich](https://www.tum.de/).

All examples have been tested with [OpenFoam-3.0.1](https://openfoam.org/version/3-0-1/) on [Ubuntu 14.04.6 LTS (Trusty Tahr)](https://releases.ubuntu.com/14.04/) . For convenience a docker image and startup script is provided.

## Build Instructions

OS Manjaro Linux, Docker version 20.10.22

Clone this repository an open the root location in a terminal.

The following command will build a local Docker image called `openfoam/openfoam301`

```bash
docker build -t openfoam/openfoam301 .
```

## Getting Started

Once build, start the container via the startup script

```bash
./openfoam301-ubuntu
```

which will start an interactive session in your Docker container with your current working directory on your local filesystem mounted as `/home/openfoam/` in your container. The container will be removed once you exit it.

Test the solver

```bash
adjointSpalartAllmaras -help
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

## Usage

Follow [thesis](./thesis/thesis_dennis_kasper.pdf) chapter six (p. 25) for the unit cube example.

Run the container with the startup script, navigate to `/home/openfoam/run/unitCubeSpalartAllmaras/` and call the solver

```bash
adjointSpalartAllmaras
```

Expected Output

```bash
/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  3.0.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
Build  : 3.0.1-119cac7e8750
Exec   : adjointSpalartAllmaras
Date   : Jan 06 2023
Time   : 16:06:08
Host   : "0d907f0f3b09"
PID    : 171
Case   : /home/openfoam/run/unitCubeSpalartAllmaras
nProcs : 1
sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
fileModificationChecking : Monitoring run-time modified files using timeStampMaster
allowSystemOperations : Allowing user-supplied system call operations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Create time

Create mesh for time = 0

Reading field p

Reading field U

...

Iteration = 36
DILUPBiCG:  Solving for nuTilda, Initial residual = 1.31045457787e-12, Final residual = 1.36408329078e-23, No Iterations 3

Iteration = 37
DILUPBiCG:  Solving for nuTilda, Initial residual = 9.22875821365e-13, Final residual = 3.76451815031e-24, No Iterations 3

// * * * * * * * *   End Total finite difference method    * * * * * * * * * //

// * * * * * * * * * *     End sensitivity analysis      * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Total executionTime = 1.6 s  ClockTime = 2 s

End
```

## Reference

[1] Docker build process [CFDEngine1](https://github.com/CFDEngine/simple-openfoam-dockerfile/), [CFDEngine2](https://www.cfdengine.com/blog/how-to-install-openfoam-anywhere-with-docker)  
[2] [OpenFoam](https://openfoam.org/) Startscript [openfoam4-linux](http://dl.openfoam.org/docker/openfoam4-linux) and [Install guide](https://openfoam.org/download/8-linux/)  
[3] [Master Thesis](./thesis/thesis_dennis_kasper.pdf)

## ToDo

- [ ] Add NACA 2412 Airfoil and ONERA M6 Wing examples
