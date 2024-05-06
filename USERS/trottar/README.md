What you need:
====================================
ROOT 5.34/35+ <br>
gcc compiler (general)


Source codes:
====================================
Located in **[src](src)** directory...

**[EIC_mesonMC_n.cpp](src/EIC_mesonMC_n.cpp)**        :  Leading neutron, pion structure function with ep scattering <br>
**[EIC_mesonMC_Lambda.cpp](src/EIC_mesonMC_Lambda.cpp)**        :  Leading lambda, kaon structure function with ep scattering <br>
**[EIC_mesonMC_Sigma.cpp](src/EIC_mesonMC_Sigma.cpp)**        :  Leading sigma, kaon structure function with ep scattering <br>
**[cteq/](src/cteq/)**                  :  cteqpdf.h and data based call files (c++ wrapper) <br>
**[cteq-tbls/](src/cteq-tbls/)**             :  nucleon PDFs table <br>
**[functions/](src/functions/)**   :  various functions (e.g., structure and splitting)


How to change inputs:
====================================
Located in **[inputs](inputs)** directory...

**[kinematics.inputs](inputs/kinematics.inputs)** : edit this document to change simulation kinematics (e.g. number of events, x range, Q2 range, pbeam, kbeam)

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
**./run_mesonMC.sh <flags\> <final_state\>** : Final states...(pi/p, pi/n, k/lambda)

&#8595; Below you can see an example for a pion and neutron final state simulation &#8595;


```python
!./run_mesonMC.sh -g pi/n
! -g, event generator flag
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
In the **[OUTPUTS](OUTPUTS)** directory are the ROOT and LUND outputs for the simulation for further analysis.

Physics Documents:
====================================
Located in **[docs/](docs)** directory...

**[Chekanov_LN_e-p_HERA_2002.pdf](docs/Chekanov_LN_e-p_HERA_2002.pdf)** : Leading neutron production in e + p collisions at HERA (2002) S. Chekanov, et al. <br>
**[Aguilar_pi-k_structure_EIC_2019.pdf](docs/Aguilar_pi-k_structure_EIC_2019.pdf)** : Pion and Kaon Structure at the Electron-Ion Collider (2019) A. C. Aguilar, et al.  <br>
**[fixed_collider_physics.pdf](docs/fixed_collider_physics.pdf)** : Comparison of fixed to collider physics by Richard L. Trotta  <br>
**[meson_split_function.pdf](meson_split_function.pdf)** : Summary of splitting functions used in **[barry_split_func.py](src/functions/split_functions/barry_split_func.py)** by Patrick Barry & Richard L. Trotta  <br>
	
	
*This code is maintained by Richard L. Trotta (trottar.iii@gmail.com).*
