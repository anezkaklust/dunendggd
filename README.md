# dunendggd_toad - DUNE Near Detector Geometries in GeGeDe for The Teststand of an Overpressure Argon Detector

This is a repository for developing any DUNE (Near Detector) geometries with [[https://github.com/brettviren/gegede][gegede]]. 

Specifically, this fork of the repository is focused on creating a geoemtry for The Teststand of an Overpressure Argon Detector (TOAD).

* Getting Started

See the [[https://github.com/brettviren/gegede/blob/master/doc/install.org][gegede installation page]] for the prerequisites. Similarly to [[https://github.com/DUNE/duneggd][duneggd]], install virtualenv:

#+BEGIN_EXAMPLE
  $ wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.tar.gz
  $ tar -xf virtualenv-1.11.tar.gz
  $ python virtualenv-1.11/virtualenv.py dunendggd_toad
#+END_EXAMPLE

And clone this repository
#+BEGIN_EXAMPLE
  $ git clone https://github.com/anezkaklust/dunendggd_toad
#+END_EXAMPLE

# Setup

Each time you login on your machine, you need to activate your virtual environment, setup `duneggd_toad`, and set your python path:

#+BEGIN_EXAMPLE
  $ source dunendggd_toad/bin/activate
  $ cd dunendggd_tad/
  $ python setup.py develop --user 
  $ export PATH=`pwd`/dunendggd_toad/bin:${PATH}
#+END_EXAMPLE


# TOAD Geometry
To create `.gdml` file of TOAD, run 

```bash
gegede-cli ../duneggd_toad/Config/TOAD_Concept.cfg  -o toad_concept.gdml
```

# Quick Visualization
To do a quick check or your geometry file you can use ROOT:
```bash
root -l 'geoDisplay.C("toad_concept.gdml",5)''
```

# Contact
* **TOAD geometry:**
   * Anežka Klustová `a.klustova20@imperial.ac.uk`
* **dunendggd:**
  * Guang Yang `guang.yang.1@stonybrook.edu`
  * Jose Palomino`jose.palominogallo@stonybrook.edu`
* **GeGeDe:**
  * Brett Viren `bviren@bnl.gov`
