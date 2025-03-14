COMMENT
Take from Kim (2017). Obtained from ModelDB (McDougal et al., 2017; accession  number 235769)

- Kim, H. (2017). Muscle length-dependent contribution of motoneuron Cav1.3 channels to force production in model slow motor unit. Journal of Applied Physiology, 123(1). https://doi.org/10.1152/japplphysiol.00491.2016

- McDougal, R. A., Morse, T. M., Carnevale, T., Marenco, L., Wang, R., Migliore, M., Miller, P. L., Shepherd, G. M., & Hines, M. L. (2017). Twenty years of ModelDB and beyond: building essential modeling tools for the future of neuroscience. In Journal of Computational Neuroscience (Vol. 42, Issue 1). https://doi.org/10.1007/s10827-016-0623-7

ENDCOMMENT


TITLE Fast Sodium Channel

NEURON {
	SUFFIX Naf
	USEION na READ ena WRITE ina
	RANGE gnafbar, ina, g, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
}

PARAMETER {
	gnafbar	=0.06 	(mho/cm2) <0,1e9>
}

ASSIGNED {
	v (mV)
	ena (mv)
	ina (mA/cm2)
	i (mA/cm2)
	g (S/cm2)
	minf mtau
	hinf htau
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gnafbar * m * m * m * h
	i = g * (v - ena)
	ina = i
}

INITIAL { 
	m = 0.001
	h = 0.5829
}

DERIVATIVE states { 
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

PROCEDURE rates(v(mV)) {LOCAL a, b
	a = (-0.4*(v+49))/(exp(-(v+49)/5)-1)
	b = (0.4*(v+25))/(exp((v+25)/5)-1)
	mtau = 1/(a + b)
	minf = a/(a + b)
	htau = 30/(exp((v+60)/15)+exp(-(v+60)/16))
	hinf = 1/(1+exp((v+58)/7))
}