NEURON {
    POINT_PROCESS Receptor_5HT2
    POINTER C
    RANGE gmax, Erev, i, mod_gna, mod_gk, K1, K2
    NONSPECIFIC_CURRENT i
    GLOBAL K3, K4, KD, d1, d2, gnabar, gkbar
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (umho) = (micromho)
    (mM) = (milli/liter)
}

PARAMETER {
    K1 = 0.00225 (/ms mM)  : Forward binding rate to receptor
    K2 = 0.02 (/ms)       : Backward (unbinding) rate of receptor
    K3 = 0.083 (/ms)         : Rate of G-protein production
    K4 = 0.0079 (/ms)        : Rate of G-protein decay
    d1 = 0.017 (/ms)         : Rate of desensitization
    d2 = 0.0053 (/ms)        : Rate of re-sensitization
    KD = 100                 : Dissociation constant of Ca++ channel
    n = 4                    : Number of binding sites of G-protein on Ca++
    Erev = 0 (mV)            : Reversal potential (E_Ca)
    gmax = 0 (umho)          : Maximum conductance
    gnabar = 0.12 (S/cm2)    : Maximum sodium conductance
    gkbar = 0.036 (S/cm2)    : Maximum potassium conductance
}

ASSIGNED {
    v (mV)    : Postsynaptic voltage
    i (nA)    : Current = g * (v - Erev)
    g (umho)  : Conductance
    mod_gna (S/cm2)
    mod_gk (S/cm2)
    C (mM)    : Pointer to transmitter concentration
    Gn        : Term for G-protein raised to power n
}

STATE {
    R         : Fraction of activated receptor
    D         : Fraction of desensitized receptor
    G         : Fraction of activated G-protein
}

INITIAL {
    R = 0
    D = 0
    G = 0
}

BREAKPOINT {
    SOLVE bindkin METHOD cnexp
    Gn = G^n
    g = gmax * Gn / (Gn + KD)
    i = g * (v - Erev)
}

DERIVATIVE bindkin {
    R' = K1 * C * (1 - R - D) - K2 * R + d2 * D
    D' = d1 * R - d2 * D
    G' = K3 * R - K4 * G
}
