COMMENT
Take from Kim (2017). Obtained from ModelDB (McDougal et al., 2017; accession  number 235769)

- Kim, H. (2017). Muscle length-dependent contribution of motoneuron Cav1.3 channels to force production in model slow motor unit. Journal of Applied Physiology, 123(1). https://doi.org/10.1152/japplphysiol.00491.2016

- McDougal, R. A., Morse, T. M., Carnevale, T., Marenco, L., Wang, R., Migliore, M., Miller, P. L., Shepherd, G. M., & Hines, M. L. (2017). Twenty years of ModelDB and beyond: building essential modeling tools for the future of neuroscience. In Journal of Computational Neuroscience (Vol. 42, Issue 1). https://doi.org/10.1007/s10827-016-0623-7

ENDCOMMENT


TITLE Intracellular calcium dynamics

NEURON {
	SUFFIX Ca_conc
	USEION ca READ ica WRITE cai
	RANGE cai, cca
}

UNITS	{
	(mV) 		= (millivolt)
	(mA) 		= (milliamp)
	FARADAY 	= (faraday) (coulombs)
	(molar) 	= (1/liter)
	(mM) 		= (millimolar)
}

PARAMETER {
	f = 0.004		
	kCa = 8			(/ms)	
	alpha = 1    	(mol/C/cm2)
}

ASSIGNED {
	cai			(mM)
	ica			(mA/cm2)
}

STATE {
	cca		(mM)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

INITIAL {
	cca = 0.0001
}

DERIVATIVE state {
	cca' = f * (- alpha * ica - kCa * cca)
	cai = cca
}