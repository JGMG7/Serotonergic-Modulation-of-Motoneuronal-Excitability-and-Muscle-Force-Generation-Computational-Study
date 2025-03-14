COMMENT
Take from Kim (2017). Obtained from ModelDB (McDougal et al., 2017; accession  number 235769)

- Kim, H. (2017). Muscle length-dependent contribution of motoneuron Cav1.3 channels to force production in model slow motor unit. Journal of Applied Physiology, 123(1). https://doi.org/10.1152/japplphysiol.00491.2016

- McDougal, R. A., Morse, T. M., Carnevale, T., Marenco, L., Wang, R., Migliore, M., Miller, P. L., Shepherd, G. M., & Hines, M. L. (2017). Twenty years of ModelDB and beyond: building essential modeling tools for the future of neuroscience. In Journal of Computational Neuroscience (Vol. 42, Issue 1). https://doi.org/10.1007/s10827-016-0623-7

ENDCOMMENT


TITLE Delayed rectifier potassium channel

NEURON {
	SUFFIX KDr
	USEION k READ ek WRITE ik
	RANGE gkdrbar, ik, g, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
}

PARAMETER {
	gkdrbar 	= .8 	(mho/cm2)
}

ASSIGNED {
	v 		(mV)
	ek          (mV)
	ik		(mA/cm2)
	i (mA/cm2)
	g (S/cm2)
	ninf ntau
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gkdrbar * n * n * n * n
	i = g * (v - ek)
	ik = i
}

INITIAL {
	:rates(v)
	:n = ninf
	n = 0.1239
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n) / ntau
}

PROCEDURE rates(v(mV)) {
	ntau = 5 / (exp((v + 50)/40) + exp(-(v + 50)/50))
	ninf = 1 / ((exp(-(v + 31)/15)) + 1)
}




