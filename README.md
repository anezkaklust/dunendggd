# DUNE Near Detector Geometries in GeGeDe for TOAD

This is a repository for developing any [DUNE (Near Detector) geometries](https://github.com/DUNE/dunendggd) with [gegede](https://github.com/brettviren/gegede). 

Specifically, this fork of the repository is focused on creating a geometry for [The Teststand of an Overpressure Argon Detector (TOAD)](https://indico.fnal.gov/event/54624/contributions/244519/subcontributions/8549/attachments/156742/204687/HPgTPCTestBeamIntro_050722.pdf).

## Getting Started

See the [gegede installation page](https://github.com/brettviren/gegede/blob/master/doc/install.org) for the prerequisites. Similarly to 
[duneggd](https://github.com/DUNE/duneggd), install `virtualenv`:

```bash
wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.tar.gz
tar -xf virtualenv-1.11.tar.gz
python virtualenv-1.11/virtualenv.py dunendggd_toad
```

And clone this repository:
```bash
git clone https://github.com/anezkaklust/dunendggd_toad
```

## Setup

Each time you login on your machine, you need to activate your virtual environment, setup `duneggd_toad`, and set your python path:

```bash
source dunendggd_toad/bin/activate
cd dunendggd_toad/
python setup.py develop --user 
export PATH=`pwd`/dunendggd_toad/bin:${PATH}
```


## TOAD Geometry
To create `.gdml` file of TOAD, run 

```bash
gegede-cli duneggd_toad/Config/TOAD_Concept.cfg  -o toad_concept.gdml
```

## Visualization
To visualize your geometry file you can use ROOT:
```bash
root -l 'geoDisplay.C("toad_concept.gdml",5)''
```

## Contact
* **TOAD geometry:**
   * Anežka Klustová `a.klustova20@imperial.ac.uk`
* **dunendggd:**
  * Guang Yang `guang.yang.1@stonybrook.edu`
  * Jose Palomino`jose.palominogallo@stonybrook.edu`
* **GeGeDe:**
  * Brett Viren `bviren@bnl.gov`
