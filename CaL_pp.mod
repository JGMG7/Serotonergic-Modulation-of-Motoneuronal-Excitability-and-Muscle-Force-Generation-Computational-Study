COMMENT
Take from Kim (2017). Obtained from ModelDB (McDougal et al., 2017; accession  number 235769)

- Kim, H. (2017). Muscle length-dependent contribution of motoneuron Cav1.3 channels to force production in model slow motor unit. Journal of Applied Physiology, 123(1). https://doi.org/10.1152/japplphysiol.00491.2016

- McDougal, R. A., Morse, T. M., Carnevale, T., Marenco, L., Wang, R., Migliore, M., Miller, P. L., Shepherd, G. M., & Hines, M. L. (2017). Twenty years of ModelDB and beyond: building essential modeling tools for the future of neuroscience. In Journal of Computational Neuroscience (Vol. 42, Issue 1). https://doi.org/10.1007/s10827-016-0623-7

ENDCOMMENT


TITLE L-type Calcium channel

NEURON {
	POINT_PROCESS CaL_pp
	USEION ca WRITE ica
	RANGE gcalbar, g, eca, i, ica
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	gcalbar  = 7e-07 (uS)
	mtau	 = 60	    (ms)
	eca      = 60
}

ASSIGNED {
	ica		(nA)
	v		(mV)
	g		(uS)
	i		(nA)
	minf
}

STATE {
	m
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gcalbar * m
	i =  g * (v - eca)
	ica=i
}

INITIAL {
	rates(v)
	m = minf
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m) / mtau
}

PROCEDURE rates(v(mV)) {
	minf = 1 / (1 + exp(-(v + 43)/6))
}