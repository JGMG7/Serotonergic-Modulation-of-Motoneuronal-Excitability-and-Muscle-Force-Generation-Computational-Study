COMMENT
Take from Kim (2017). Obtained from ModelDB (McDougal et al., 2017; accession  number 235769)

- Kim, H. (2017). Muscle length-dependent contribution of motoneuron Cav1.3 channels to force production in model slow motor unit. Journal of Applied Physiology, 123(1). https://doi.org/10.1152/japplphysiol.00491.2016

- McDougal, R. A., Morse, T. M., Carnevale, T., Marenco, L., Wang, R., Migliore, M., Miller, P. L., Shepherd, G. M., & Hines, M. L. (2017). Twenty years of ModelDB and beyond: building essential modeling tools for the future of neuroscience. In Journal of Computational Neuroscience (Vol. 42, Issue 1). https://doi.org/10.1007/s10827-016-0623-7

ENDCOMMENT


TITLE N-type Calcium channel

NEURON {
	SUFFIX CaN
	USEION ca READ cai, cao WRITE ica
	RANGE gcanbar, ica, g, i, eca
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
}

PARAMETER {
	gcanbar  = .01 (mho/cm2)
	mtau	= 15	(ms)
	htau	= 50	(ms)
	R = 8.31441 (VC/Mol/K)
	T = 309.15 (k) 
	Z = 2
	F = 96485.309 (/C)
}

ASSIGNED {
	v 		(mV)
	eca 	(mv)
	ica		(mA/cm2)
	i 		(mA/cm2)
	g 		(S/cm2)
	minf hinf
	cai cao	(mM)
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	eca = ((1000*R*T)/(Z*F))*log(cao/cai)
	g = gcanbar * m * m * h
	i = g * (v - (eca-70))
	ica = i
}

INITIAL {
	m = 0.004199
	h = 0.9219
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m) / mtau
	h' = (hinf - h) / htau
}

PROCEDURE rates(v(mV)) {
	minf = 1 / (1 + exp((v+25)/-5))
	hinf = 1 / (1 + exp((v+43)/5))
}



