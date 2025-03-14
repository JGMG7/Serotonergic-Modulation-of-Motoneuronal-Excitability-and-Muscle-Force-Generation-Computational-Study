# Serotonergic-Modulation-of-Motoneuronal-Excitability-and-Muscle-Force-Generation-Computational-Study
This study investigates the impact of serotonin (5-HT) on motoneuron electrical activity and muscle force generation. Using a computational model, we explore how 5-HT receptors influence motoneuron excitability and muscle function at different stimulation frequencies.


This package is running with the NEURON simulation program version 8.2 as Python 3.12 module, available on internet at: https://www.neuron.yale.edu/neuron/ and www.python.org

The package contains the follow mechanisms (.mod files) and programs (.py files) needed to do the simulations:

Computational model of a motoneuron incorporating the 5-HT1a and 5-HT2a serotonin receptor, presynaptic neuron that releases serotonin (5-HT) upon stimulation and muscle. The model aims to elucidate the role of serotonin in motor function and its impact on muscle contraction.

PROGRAMS
mn_10hz_input.py, mn_40hz_input.py and mn_100hz_input.py  : which in turn contains the following program:

v_e_moto6_export.hoc    :  Anatomical data corresponding to a cat motoneuron (i.e., v_e_moto6), which is available in the public database (www.neuromorph.org)

the difference between programs is the stimuli frecuency, wich is into the programs name (10, 40 or 100 hz)

each program contain 3 models:

a) Motoneuron without 5-HT receptor
b) Motoneuron with 5-HT receptors and serotonine release = 10000 molecules
c) Motoneuron with 5-HT receptors and serotonine release = 10000000 molecules

MECHANISMS
Receptor_5HT1a.mod : 5-HT1a receptor
Receptor_5HT2.mod : 5-HT2a receptor

Ca_conc.mod : Calcium concentration

CaL_pp.mod : L-type Calcium channel

caL.mod : high-threshold calcium current in the presynaptic terminal

CaN.mod : N-type Calcium channel

HH2.mod : Hippocampal HH channels

KCa.mod : Ca activated potassium channel

KDr.mod : Delayed rectifier potassium channel

Naf.mod : Fast Sodium Channel

Nap.mod : Persistent Sodium Channel

rel2.mod : kinetic model for the release of transmitter

HOW TO RUN
To compile the demo, NEURON and INTERVIEWS must be installed and working on the machine you are using. Just type "nrnivmodl" to compile the mechanisms given in the mod files.


MODIFICATIONS

"""""5-HT concentration:"""""

To modify 5-HT concentration release, change in the "insert presynaptic mechanisms" sections of program.

nt_rel2 = 10000 // 10000000 number of transmitter molecule per vesicle <---- 

For more information about how to get NEURON in Python and how to install it, please refer to the following sites: https://www.neuron.yale.edu/neuron/
