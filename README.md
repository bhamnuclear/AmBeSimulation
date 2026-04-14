# This example is for AmBe source simulation
It also includes examples for other sources like PuBe, CmBe, PoBe. See following for details

## Citation
Is using any parts of this work, please cite the relevant paper that you can find here open access: https://www.sciencedirect.com/science/article/pii/S0168900225010356?via%3Dihub

Cite as:<br>
F. Falezza et al., "Simulation and analysis of an AmBe source"
Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, Volume 1085, 2026, 171233, ISSN 0168-9002, doi.org/10.1016/j.nima.2025.171233

## Depen
dencies
- Geant4
- ROOT
- OpenMPI or other MPI interface 

MPI can be disabled in CMake and main. In any case, it can be chosen at runtime whether to run with or without MPI via the `AMBE_MPI` environment flag: set this either `TRUE` or `FALSE`.

## MPI USAGE
Use with [https://github.com/bhamnuclear/Geant4-MPI/](https://github.com/bhamnuclear/Geant4-MPI/), here available in the `Geant4-MPI` folder.

compile Geant4 MPI first using `Compile.sh` or othewrwise, then compile this simulation:
```
$ mkdir build
$ cd build
$ cmake -DG4mpi_DIR=$PWD/../Geant4-MPI/install/lib/G4mpi-11.3.1 ..
$ make
```
Please replace 11.3.1 with your Geant4 version.


Run with
```$ AMBE_MPI=true mpirun --mca accelerator rocm -n 2 ambe -t 8 macros/shortlong.mac```
or 
```$ AMBE_MPI=true mpirun -n 2 ambe -t 8 macros/shortlong.mac```


## Current Status
1. Scoring has been implemented correctly
2. Full Multithreading and MPI support
3. AmBe and 239PuBe sources. Other sources have only been partially implemented
6. World material is `G4_AIR`


## Usage
Various options have been implemented using cli. It is also possible to use macros, but there is no full support for these and they have not been fully tested. Some macro commands can replace cli commands, but some of the macro commands cannot replace cli commands.

Note that cli commands are more stable and do not cause crashing when run once the simulation is initialised.

### Cli commands
Options list (as listed in main.cc)
- `-i` `--interactive` -> enable GUI. Default OFF is options specified.
- `-t` `--threads` -> Specify number of threads to use. If using MPI, this is the number of threads allocated to each node. DEFAULT 1.
- `-ffs` `--fissionFragmentsScore` -> Enable scoring of fission fragments (neutrons and Ions). DEFAULT OFF.
- `-r` `--isotope` -> Specify Neutron Source to use. 0 is AmBe, 1 is 239PuBe. DEFAULT AmBe. Can also be set to monoenergetic passing SingleX.YZ Where X.YZ is the energy in MeV.
- `-in` `--initialNeutrons` -> Score initial neutrons from the Priamry Generator. DEFAULT OFF.
- `-nw` `--nowater` -> Define is world is water or not (or other materials).
- `-nt` `--neutronTracking` -> Selectively enable or disable neutron tracking scoring
- `-cs` `--casing` -> Define Casing. 0 is cylindrical approximation, 1 is X3 casing.
- `-sg` `--scoregamma` -> Add gamma scoring to resulting tree
- `-ass` `--azimuthalSurfaceScoring` -> Score Azimuthal distirbution of outgoing neutrons.

## Extra sources

To use other sources that are not AmBe, like CmBe, PuBe, RaBe and so on, the following modifications are needed:<br>
- create your input isotope card. This needs to contain all the alpha decays (see nndc) and intensities as displayed there
- load this in PrimaryGenerator (lines 165-168) or otherwise from main or however you prefer
- define your Eloss equations. Suggested use is Lise++, which uses atima1.4. This code supports isotopes of mass above natural uranium. Sample various energies from just above to the isotope maximum down to 0. Get both range (z) as a function of energy (E) and energy as a function of range: z(E) and E(z).
- Fit your Eloss equations (quadratic fit recommended). Implement these in PrimaryGeneratorAction::Eloss
- Run the simulation
- ENJOY!

If you publish your work (or oterwise) do not forget to cite this work :)
