What you need:
====================================
ROOT 5.34/35+ <br>
gcc compiler (general)


Source codes:
====================================
Located in **src** directory...

**TDISMC_EIC.cpp**         :  pion structure function with ep scattering at JLEIC <br>
**TDISMC_EICn.cpp**        :  pion structure function with eD scattering at JLEIC <br>
**TDISMC_EICK.cpp**        :  kaon structure function with ep scattering at JLEIC <br>
**cteq/**                  :  cteqpdf.h and data based call files (c++ wrapper) <br>
**cteq-tbls/**             :  nucleon PDFs table <br>
**structure_functions/**   :  various regularization form for pion SF/FF


How to change inputs:
====================================
Located in **inputs** directory...

**kinematics.inputs** : edit this document to change simulation kinematics (e.g. number of events, x range, Q2 range, pbeam, kbeam)

All other constants are changed in **src/TDISMC_EIC.h**

&#8595; Below you can see the current kinematics inputs &#8595;


```python
!more inputs/kinematics.input
```

    XMIN=0.001
    XMAX=1.00
    Q2MIN=1.0
    Q2MAX=1100.0
    NEVTS=500000
    PBEAM=135.0
    KBEAM=10.0


How to run:
====================================
**./run_batch.sh <final_state\>** : Final states...(pi/p, pi/n, k/lambda)

&#8595; Below you can see an example for a pion and neutron final state simulation &#8595;


```python
!./run_mesonMC.sh pi/n
```

    
    Pion with neutron final state selected
    
    Your kinematics: [xBj_min:xBj_max] = [ 0.001000: 1.000000] 
    Your kinematics: [Q2_min:Q2_max] = [ 1.000000:1100.000000] 
    Incident Ion Mass   0.93827 GeV 
    Incident Electron, Ion Momenta:  10.0000,   135.00 GeV/c | s_0 =  5400.1019 GeV^2 
    Warning in <TTree::Bronch>: TLorentzVector cannot be split, resetting splitlevel to 0
    Warning in <TTree::Bronch>: TLorentzVector cannot be split, resetting splitlevel to 0
    Total of 57310 events out of 500000 Trials ============================] 100 %
    (int) 57310


ROOT and LUND outputs:
====================================
In the **OUTPUTS** directory are the ROOT and LUND outputs for the simulation for further analysis.


Running GEANT4
====================================
Located in **g4e_files/** directory...

**./run_g4e.sh** : Will run the python script for the GEANT4 simulation from the TDIS_lund.dat file <br>
**eic_g4e.py**   : This python script will run GEANT4 simulation for the lund file specified for detector='jleic' and beamline='erhic'

*This code is maintained by Richard Trotta (trotta@cua.edu).*
