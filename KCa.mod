COMMENT
Take from Kim (2017). Obtained from ModelDB (McDougal et al., 2017; accession  number 235769)

- Kim, H. (2017). Muscle length-dependent contribution of motoneuron Cav1.3 channels to force production in model slow motor unit. Journal of Applied Physiology, 123(1). https://doi.org/10.1152/japplphysiol.00491.2016

- McDougal, R. A., Morse, T. M., Carnevale, T., Marenco, L., Wang, R., Migliore, M., Miller, P. L., Shepherd, G. M., & Hines, M. L. (2017). Twenty years of ModelDB and beyond: building essential modeling tools for the future of neuroscience. In Journal of Computational Neuroscience (Vol. 42, Issue 1). https://doi.org/10.1007/s10827-016-0623-7

ENDCOMMENT


TITLE Ca activated potassium channel

NEURON {
	SUFFIX KCa
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gkcabar, ik, i, g
}

UNITS {
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
	(molar)	=	(1/liter)
	(mM)	=	(millimolar)
	(S)	=	(siemens)
}

PARAMETER {
	gkcabar = 0.02		(S/cm2)
	Kd		= 0.0005		(mM)
}

ASSIGNED {
	ik			(mA/cm2)
	v			(mV)
	ek			(mV)
	cai			(mM)
	i			(mA/cm2)
	g			(S/cm2)
}

INITIAL {
	ik = 0
}

BREAKPOINT {
	g = gkcabar * (cai/(cai+Kd))
	i = g * (v - ek)
	ik = i
}


