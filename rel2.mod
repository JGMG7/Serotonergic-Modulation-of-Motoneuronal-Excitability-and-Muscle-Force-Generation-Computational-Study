TITLE transmitter release

COMMENT

The models used were obtained and then modified from ModelDB (McDougal et al., 2017; accession number 18198). 

- Destexhe, A., Mainen, Z. F., & Sejnowski, T. J. (1994). Synthesis of models for excitable membranes, synaptic transmission and neuromodulation using a common kinetic formalism. Journal of Computational Neuroscience, 1(3). https://doi.org/10.1007/BF00961734

- Destexhe, A., & Sejnowski, T. J. (1995). G protein activation kinetics and spillover of γ-aminobutyric acid may account for differences between inhibitory responses in the hippocampus and thalamus. Proceedings of the National Academy of Sciences of the United States of America, 92(21). https://doi.org/10.1073/pnas.92.21.9515

- McDougal, R. A., Morse, T. M., Carnevale, T., Marenco, L., Wang, R., Migliore, M., Miller, P. L., Shepherd, G. M., & Hines, M. L. (2017). Twenty years of ModelDB and beyond: building essential modeling tools for the future of neuroscience. In Journal of Computational Neuroscience (Vol. 42, Issue 1). https://doi.org/10.1007/s10827-016-0623-7


 Simple (minimal?) model of transmitter release

 - single compartment, need calcium influx and efflux

 - Ca++ binds to a "fusion factor" protein F leading to an activated form FA.
   Assuming a cooperativity factor of 4 (see Augustine & charlton, 
   J Physiol. 381: 619-640, 1986), one obtains:

	F + 4 Cai <-> FA	(kb,ku)

 - FA binds to presynaptic vesicles and activates them according to:

	FA + V <-> VA		(k1,k2)

   VA represents the "activated vesicle" which is able to bind to the
   membrane and release transmitter.  Presynaptic vesicles (V) are 
   considered in excess.

 - VA releases nt transmitter molecules in the synaptic cleft

	VA  ->  nt T		(k3)

   This reaction is the slowest and a constant number of transmitter per 
   vesicule is considered (nt).  

 - Finally, T is hydrolyzed according to a first-order reaction

	T  ->  ...		(kh)


-----------------------------------------------------------------------------
ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS rel2
	USEION ca READ ica, cai WRITE cai
	RANGE T,FA,CA,Fmax,Ves,b,u,k1,k2,k3,nt,kh
: from cad :
	RANGE depth,kt,kd,cainf,taur
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
: from cad:
	(molar) = (1/liter)			: moles do not appear in units
:	(mM)	= (millimolar)
	(um)	= (micron)
:	(mA)	= (milliamp)
	(msM)	= (ms mM)

}
: from cad:

CONSTANT {
	FARADAY = 96489		(coul)		: moles do not appear in units
:	FARADAY = 96.489	(k-coul)	: moles do not appear in units
}

PARAMETER {

	Ves = 0.1 	(mM)		: conc of vesicles
	Fmax = 0.001	(mM)		: conc of fusion factor F
	b = 1e16 	(/mM4-ms)	: ca binding to F
	u = 0.1  	(/ms)		: ca unbinding 
	k1 = 1000   	(/mM-ms)	: F binding to vesicle
	k2 = 0.1	(/ms)		: F unbinding to vesicle
	k3 = 4   	(/ms)		: exocytosis of T
	nt = 10000			: nb of molec of T per vesicle
	kh = 10  	(/ms)		: cst for hydolysis of T
: from cad:
	depth	= .1	(um)		: depth of shell
	taur	= 700	(ms)		: rate of calcium removal
	cainf	= 1e-8	(mM)
	kt	= 1	(mM/ms)		: estimated from k3=.5, tot=.001
	kd	= 5e-4	(mM)		: estimated from k2=250, k1=5e5
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
	drive_pump	(mM/ms)
}

STATE {
	FA	(mM)
	VA	(mM)
	T	(mM)
	cai	(mM) 
}

INITIAL {
	FA = 0
	VA = 0
	T = 0
:	cai = 1e-8
	cai = kd
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

LOCAL bfc , kfv

DERIVATIVE state {

	bfc = b * (Fmax-FA-VA) * cai^4
	kfv = k1 * FA * Ves
:	this is the old equation incorporated into the below:
:	cai'	= - bfc + 4 * u * FA 
	FA'	= bfc - u * FA - kfv + k2 * VA
	VA'	= kfv - (k2+k3) * VA
	T'	= nt * k3 * VA - kh * T
: from cad:

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump inward

:	drive_pump = -tot * k3 * cai / (cai + ((k2+k3)/k1) )	: quasistat
	drive_pump = -kt * cai / (cai + kd )		: Michaelis-Menten

:	this is the eq for cai prime from cad incorporated into below:
:	cai' = drive_channel + drive_pump + (cainf-cai)/taur
	cai'= -bfc+4*u*FA + drive_channel + drive_pump + (cainf-cai)/taur

}	
