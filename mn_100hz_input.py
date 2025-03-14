# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 13:24:29 2025

@author: sopho
"""
import neuron
from neuron.units import ms, mV, µM
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
from neuron import h


# =================================================
# ####### * whitout raphe nuclei input * #######
# =================================================

temperature = 37  # float(input("temperature: "))


stim_interv = 10  # 100 Hz

h.celsius = temperature


# ####### *motoneuron* #######
# =====================================

h.load_file("stdrun.hoc")
soma = h.load_file("v_e_moto6_export.hoc")


soma = h.soma
soma.Ra = 100
soma.cm = 1

soma.insert("pas")
g_pas = 1/225


soma.insert("Naf")
soma.gnafbar_Naf = 0.71
soma.ena = 50

soma.insert("KDr")
soma.gkdrbar_KDr = 0.23
soma.ek = -80

soma.insert("CaN")
soma.gcanbar_CaN = 0.013

soma.insert("KCa")
soma.gkcabar_KCa = 0.0258
soma.ek = -80

soma.insert("Ca_conc")

iCaL = h.CaL_pp(soma(0.5))


for i in range(311):
    dend = h.dend[i]
    dend.insert("pas")
    dend.g_pas = 1/11000  # [S/cm^2]
    dend.e_pas = -70        # [mV]
    dend.Ra = 70  # [Ohm*cm]
    dend.cm = 1             # [microF/cm^2]


cals = []

for i in range(311):
    dend = h.dend[i]

    cal = h.CaL_pp(dend(0.5))
    cals.append(cal)


# Axon Hillock
axon_hillock = h.Section(name="axon_hillock")
axon_hillock.L = 20
axon_hillock.diam = 13
axon_hillock.nseg = 11

axon_hillock.pt3dclear()  # to define as a cone
axon_hillock.pt3dadd(0, 0, 0, 20.65)
axon_hillock.pt3dadd(10, 0, 0, 3.83)

axon_hillock.insert("pas")
axon_hillock.g_pas = 1/11000
axon_hillock.e_pas = -70
axon_hillock.Ra = 70
axon_hillock.cm = 1


axon_hillock.insert("Naf")
axon_hillock.gnafbar_Naf = 2.7
axon_hillock.ena = 50

axon_hillock.insert("Nap")
axon_hillock.gnapbar_Nap = 0.000033
axon_hillock.ena = 50

axon_hillock.insert("KDr")
axon_hillock.gkdrbar_KDr = 0.17
axon_hillock.ek = -80


# Initial Segment
initial_segment = h.Section(name="initial_segment")
initial_segment.L = 30
initial_segment.nseg = 3
initial_segment.diam = 3.3

initial_segment.insert("pas")
initial_segment.g_pas = 1/11000
initial_segment.e_pas = -70
initial_segment.Ra = 70
initial_segment.cm = 1

initial_segment.insert("Naf")
initial_segment.gnafbar_Naf = 2.7
initial_segment.ena = 50

initial_segment.insert("Nap")
initial_segment.gnapbar_Nap = 0.000033
initial_segment.ena = 50

initial_segment.insert("KDr")
initial_segment.gkdrbar_KDr = 0.17
initial_segment.ek = -80


axon_hillock.connect(soma(1))
initial_segment.connect(axon_hillock(1))


# presyn and soma synapsis

def create_stimulation(section, tau, start, number, interval, weight, repeat):
    synapses = []
    stimuli = []
    netcons = []

    for i in range(repeat):
        syn = h.ExpSyn(section(0.5))
        syn.tau = tau

        stim = h.NetStim()
        stim.start = start
        stim.number = number
        stim.interval = stim_interv

        nc = h.NetCon(stim, syn)
        nc.weight[0] = weight
        synapses.append(syn)
        stimuli.append(stim)
        netcons.append(nc)

    return synapses, stimuli, netcons


syn_soma, stim_soma, nc_soma = create_stimulation(
    soma, 0.5, 0, 100, stim_interv, 1.5, repeat=1)


# dendrites stimulation

stims = []
netcons = []

for i in range(311):
    dend = h.dend[i]

    stim_syn = h.ExpSyn(dend(0.5))
    stim_syn.tau = 0.5
    stims.append(stim_syn)

    stim1 = h.NetStim()
    stim1.start = 0
    stim1.number = 100
    stim1.interval = stim_interv

    nc_da1 = h.NetCon(stim1, stim_syn)
    nc_da1.weight[0] = 1.5
    netcons.append(nc_da1)
# %%

# ####### *muscle model* #######
# =====================================

muscle = h.Section(name='muscle')
muscle.L = 10
muscle.diam = 10
muscle.insert("pas")
muscle.g_pas = 0.002

calciumObject = h.calcium(muscle(0.5))
forceObject = h.force(muscle(0.5))
# %%

# ####### connect neuron to muscle #######
# ======================================
neuromuscularJunction = \
    h.NetCon(initial_segment(0.5)._ref_v, calciumObject, sec=initial_segment)
neuromuscularJunction.threshold = -40 * mV


# Set Pointers
h.setpointer(calciumObject._ref_A, 'aPointer', forceObject)
h.setpointer(calciumObject._ref_xm, 'xmPointer', forceObject)

# %%

v_mn = h.Vector().record(soma(0.5)._ref_v)
inseg = h.Vector().record(initial_segment(0.5)._ref_v)
dend = h.Vector().record(dend(0.8)._ref_v)
ca = h.Vector().record(calciumObject._ref_Ca)
f = h.Vector().record(forceObject._ref_F)
t = h.Vector().record(h._ref_t)
na = h.Vector().record(soma(0.5)._ref_ina)
k = h.Vector().record(soma(0.5)._ref_ik)
na_is = h.Vector().record(initial_segment(0.5)._ref_ina)
k_is = h.Vector().record(initial_segment(0.5)._ref_ik)


# Run simulation

tstop = 1 * (stim_interv * 100) + 500 * ms
h.load_file('stdrun.hoc')
h.finitialize(-65 * mV)
h.continuerun(tstop)


# Plotting

plt.figure(figsize=(10, 10))
'''
ax1 = plt.subplot(6, 1, 1)
ax1.plot(t, pres_v, 'b-')
ax1.axis([0, tstop, -70, 40])
ax1.set_ylabel('presyn neuron (mV)')
ax1.set_xlabel('t (ms)')

ax2 = plt.subplot(6, 1, 2)
ax2.plot(t, ser, 'b-')
ax2.axis([0, tstop, -2, 60])
ax2.set_ylabel('5-HT rel (nM)')  # arreglar unidad
ax2.set_xlabel('t (ms)')
'''

ax3 = plt.subplot(6, 1, 3)
ax3.plot(t, v_mn, 'b-')
ax3.axis([0, tstop, -100, 50])
ax3.set_ylabel('Motoneuron (mV)')
'''
ax4 = plt.subplot(11, 1, 4)
ax4.plot(t, na, 'b-')
ax4.axis([0, tstop, -10, 1])
ax4.set_ylabel('Na soma')
ax4.set_xlabel('t (ms)')

ax5 = plt.subplot(11, 1, 5)
ax5.plot(t, k, 'b-')
ax5.axis([0, tstop, -1, 5])
ax5.set_ylabel('K soma')
ax5.set_xlabel('t (ms)')

ax6 = plt.subplot(11, 1, 6)
ax6.plot(t, inseg, 'b-')
ax6.axis([0, tstop, -70, 50])
ax6.set_ylabel('init seg (mV)')
ax6.set_xlabel('t (ms)')

ax7 = plt.subplot(11, 1, 7)
ax7.plot(t, na_is, 'b-')
ax7.axis([0, tstop, -30, 1])
ax7.set_ylabel('Na init seg')
ax7.set_xlabel('t (ms)')

ax8 = plt.subplot(11, 1, 8)
ax8.plot(t, k_is, 'b-')
ax8.axis([0, tstop, -1, 2.5])
ax8.set_ylabel('K init seg')
ax8.set_xlabel('t (ms)')
'''
ax9 = plt.subplot(6, 1, 4)
ax9.plot(t, ca / µM, 'b-')
ax9.axis([0, tstop, 0, 0.1])
ax9.set_ylabel('Calcium (ÂµM)')

ax10 = plt.subplot(6, 1, 5)
ax10.plot(t, f, 'b-')
ax10.axis([0, tstop, 0, 25])
ax10.set_ylabel('Muscle Force (N)')
ax10.set_xlabel('t (ms)')

ax11 = plt.subplot(6, 1, 6)
ax11.plot(t, dend, 'b-')
ax11.axis([0, tstop, -90, 100])
ax11.set_ylabel('dend[2] (mV)')
ax11.set_xlabel('t (ms)')

plt.show()
# %%
# save data

time = np.array(t)
vm_mn = np.array(v_mn)
init_seg = np.array(inseg)
vm_dend = np.array(dend)
ca_muscle = np.array(ca)
force_muscle = np.array(f)
na_mn = np.array(na)
k_mn = np.array(k)
na_init_seg = np.array(na_is)
k_init_seg = np.array(k_is)


data = {
    't (ms)': time,
    'Motoneuron (mV)': vm_mn,
    'init seg (mV)': init_seg,
    'dend[2] (mV)': vm_dend,
    'Calcium (µM)': ca_muscle,
    'Muscle Force (N)': force_muscle,
    'Na soma': na_mn,
    'K soma': k_mn,
    'Na init seg': na_init_seg,
    'K init seg': k_init_seg
}


df = pd.DataFrame(data)


df.to_csv('mn_100hz_1.csv', index=False)

#
#
#
#
#

# ==================================================
# ####### * physiological raphe nuclei input * #######
# ==================================================

stim_interv = 10  # 100 Hz

h.celsius = temperature


# ####### *presynaptic cell* #######
# =====================================

presyn = h.Section(name="presyn")
presyn.diam = 10
presyn.L = 10
presyn.Ra = 100   # Axial resistance in Ohm * cm
presyn.cm = 1  # Membrane capacitance in micro Farads / cm^2

presyn.insert("pas")
presyn.g_pas = 1/500
presyn.e_pas = -70


presyn.insert("hh2")
presyn.ek = -90
presyn.gnabar_hh2 = 0.1
presyn.gkbar_hh2 = 0.03


presyn.insert("caL")      # HV Ca++ channel for transmitter release


rel = h.rel2(presyn(0.5))
rel.nt = 10000
# %%

# ####### *motoneuron* #######
# =====================================

h.load_file("stdrun.hoc")
soma = h.load_file("v_e_moto6_export.hoc")


soma = h.soma
soma.Ra = 100
soma.cm = 1

soma.insert("pas")
g_pas = 1/225


soma.insert("Naf")
soma.gnafbar_Naf = 0.71
soma.ena = 50

soma.insert("KDr")
soma.gkdrbar_KDr = 0.23
soma.ek = -80

soma.insert("CaN")
soma.gcanbar_CaN = 0.013

soma.insert("KCa")
soma.gkcabar_KCa = 0.0258
soma.ek = -80

soma.insert("Ca_conc")

iCaL = h.CaL_pp(soma(0.5))


for i in range(311):
    dend = h.dend[i]
    dend.insert("pas")
    dend.g_pas = 1/11000  # [S/cm^2]
    dend.e_pas = -70        # [mV]
    dend.Ra = 70  # [Ohm*cm]
    dend.cm = 1             # [microF/cm^2]


cals = []

for i in range(311):
    dend = h.dend[i]

    cal = h.CaL_pp(dend(0.5))
    cals.append(cal)


# Axon Hillock
axon_hillock = h.Section(name="axon_hillock")
axon_hillock.L = 20
axon_hillock.diam = 13
axon_hillock.nseg = 11

axon_hillock.pt3dclear()  # to define as a cone
axon_hillock.pt3dadd(0, 0, 0, 20.65)
axon_hillock.pt3dadd(10, 0, 0, 3.83)

axon_hillock.insert("pas")
axon_hillock.g_pas = 1/11000
axon_hillock.e_pas = -70
axon_hillock.Ra = 70
axon_hillock.cm = 1


axon_hillock.insert("Naf")
axon_hillock.gnafbar_Naf = 2.7
axon_hillock.ena = 50

axon_hillock.insert("Nap")
axon_hillock.gnapbar_Nap = 0.000033
axon_hillock.ena = 50

axon_hillock.insert("KDr")
axon_hillock.gkdrbar_KDr = 0.17
axon_hillock.ek = -80


# Initial Segment
initial_segment = h.Section(name="initial_segment")
initial_segment.L = 30
initial_segment.nseg = 3
initial_segment.diam = 3.3

initial_segment.insert("pas")
initial_segment.g_pas = 1/11000
initial_segment.e_pas = -70
initial_segment.Ra = 70
initial_segment.cm = 1

initial_segment.insert("Naf")
initial_segment.gnafbar_Naf = 2.7
initial_segment.ena = 50

initial_segment.insert("Nap")
initial_segment.gnapbar_Nap = 0.000033
initial_segment.ena = 50

initial_segment.insert("KDr")
initial_segment.gkdrbar_KDr = 0.17
initial_segment.ek = -80


axon_hillock.connect(soma(1))
initial_segment.connect(axon_hillock(1))


# presyn and soma synapsis

def create_stimulation(section, tau, start, number, interval, weight, repeat):
    synapses = []
    stimuli = []
    netcons = []

    for i in range(repeat):
        syn = h.ExpSyn(section(0.5))
        syn.tau = tau

        stim = h.NetStim()
        stim.start = start
        stim.number = number
        stim.interval = stim_interv

        nc = h.NetCon(stim, syn)
        nc.weight[0] = weight
        synapses.append(syn)
        stimuli.append(stim)
        netcons.append(nc)

    return synapses, stimuli, netcons


syn_presyn, stim_presyn, nc_presyn = create_stimulation(
    presyn, 0.5, 0, 100, stim_interv, 1.5, repeat=1)
syn_soma, stim_soma, nc_soma = create_stimulation(
    soma, 0.5, 0, 100, stim_interv, 1.5, repeat=1)


# dendrites stimulation

stims = []
netcons = []

for i in range(311):
    dend = h.dend[i]

    stim_syn = h.ExpSyn(dend(0.5))
    stim_syn.tau = 0.5
    stims.append(stim_syn)

    stim1 = h.NetStim()
    stim1.start = 0
    stim1.number = 100
    stim1.interval = stim_interv

    nc_da1 = h.NetCon(stim1, stim_syn)
    nc_da1.weight[0] = 1.5
    netcons.append(nc_da1)
# %%

# add 5-HT receptors

rec1a = h.Receptor_5HT1a(initial_segment(0.5))
rec1a.gmax = 0.5
rec1a.Erev = -95
rec1a.K1 = 0.000235
rec1a.K2 = 0.002
rec1a.d2 = 0.00453


rec2 = h.Receptor_5HT2(soma(0.5))
rec2.gmax = 0.001
rec2.Erev = 0
rec2.K1 = 0.00225
rec2.K2 = 0.02

receptors = []

for i in range(311):
    dendrita_actual = h.dend[i]

    rec2 = h.Receptor_5HT2(dendrita_actual(0.5))
    rec2.gmax = 0.001
    rec2.Erev = 0
    rec2.K1 = 0.00225
    rec2.K2 = 0.02

    rec2._ref_C = rel._ref_T
    receptors.append(rec2)


rec1a._ref_C = rel._ref_T
rec2._ref_C = rel._ref_T

# %%

# ####### *muscle model* #######
# =====================================

muscle = h.Section(name='muscle')
muscle.L = 10
muscle.diam = 10
muscle.insert("pas")
muscle.g_pas = 0.002

calciumObject = h.calcium(muscle(0.5))
forceObject = h.force(muscle(0.5))
# %%

# ####### connect neuron to muscle #######
# ======================================
neuromuscularJunction = \
    h.NetCon(initial_segment(0.5)._ref_v, calciumObject, sec=initial_segment)
neuromuscularJunction.threshold = -40 * mV


# Set Pointers
h.setpointer(calciumObject._ref_A, 'aPointer', forceObject)
h.setpointer(calciumObject._ref_xm, 'xmPointer', forceObject)

# %%
# Record data for plots
pres_v = h.Vector().record(presyn(0.5)._ref_v)
ser = h.Vector().record(rel._ref_T)
v_mn = h.Vector().record(soma(0.5)._ref_v)
inseg = h.Vector().record(initial_segment(0.5)._ref_v)
dend = h.Vector().record(dendrita_actual(0.8)._ref_v)
ca = h.Vector().record(calciumObject._ref_Ca)
f = h.Vector().record(forceObject._ref_F)
t = h.Vector().record(h._ref_t)
na = h.Vector().record(soma(0.5)._ref_ina)
k = h.Vector().record(soma(0.5)._ref_ik)
na_is = h.Vector().record(initial_segment(0.5)._ref_ina)
k_is = h.Vector().record(initial_segment(0.5)._ref_ik)


# Run simulation

tstop = 1 * (stim_interv * 100) + 500 * ms
h.load_file('stdrun.hoc')
h.finitialize(-65 * mV)
h.continuerun(tstop)


# Plotting

plt.figure(figsize=(10, 10))

ax1 = plt.subplot(6, 1, 1)
ax1.plot(t, pres_v, 'b-')
ax1.axis([0, tstop, -70, 40])
ax1.set_ylabel('presyn neuron (mV)')
ax1.set_xlabel('t (ms)')

ax2 = plt.subplot(6, 1, 2)
ax2.plot(t, ser, 'b-')
ax2.axis([0, tstop, -2, 10])
ax2.set_ylabel('5-HT rel (nM)')
ax2.set_xlabel('t (ms)')


ax3 = plt.subplot(6, 1, 3)
ax3.plot(t, v_mn, 'b-')
ax3.axis([0, tstop, -100, 50])
ax3.set_ylabel('Motoneuron (mV)')
'''
ax4 = plt.subplot(11, 1, 4)
ax4.plot(t, na, 'b-')
ax4.axis([0, tstop, -10, 1])
ax4.set_ylabel('Na soma')
ax4.set_xlabel('t (ms)')

ax5 = plt.subplot(11, 1, 5)
ax5.plot(t, k, 'b-')
ax5.axis([0, tstop, -1, 5])
ax5.set_ylabel('K soma')
ax5.set_xlabel('t (ms)')

ax6 = plt.subplot(11, 1, 6)
ax6.plot(t, inseg, 'b-')
ax6.axis([0, tstop, -70, 50])
ax6.set_ylabel('init seg (mV)')
ax6.set_xlabel('t (ms)')

ax7 = plt.subplot(11, 1, 7)
ax7.plot(t, na_is, 'b-')
ax7.axis([0, tstop, -30, 1])
ax7.set_ylabel('Na init seg')
ax7.set_xlabel('t (ms)')

ax8 = plt.subplot(11, 1, 8)
ax8.plot(t, k_is, 'b-')
ax8.axis([0, tstop, -1, 2.5])
ax8.set_ylabel('K init seg')
ax8.set_xlabel('t (ms)')
'''
ax9 = plt.subplot(6, 1, 4)
ax9.plot(t, ca / µM, 'b-')
ax9.axis([0, tstop, 0, 0.1])
ax9.set_ylabel('Calcium (ÂµM)')

ax10 = plt.subplot(6, 1, 5)
ax10.plot(t, f, 'b-')
ax10.axis([0, tstop, 0, 25])
ax10.set_ylabel('Muscle Force (N)')
ax10.set_xlabel('t (ms)')

ax11 = plt.subplot(6, 1, 6)
ax11.plot(t, dend, 'b-')
ax11.axis([0, tstop, -90, 100])
ax11.set_ylabel('dend[2] (mV)')
ax11.set_xlabel('t (ms)')

plt.show()
# %%
# save data

time = np.array(t)
presynaptic = np.array(pres_v)
ser = np.array(ser)
vm_mn = np.array(v_mn)
init_seg = np.array(inseg)
vm_dend = np.array(dend)
ca_muscle = np.array(ca)
force_muscle = np.array(f)
na_mn = np.array(na)
k_mn = np.array(k)
na_init_seg = np.array(na_is)
k_init_seg = np.array(k_is)


data = {
    't (ms)': time,
    'presyn neuron (mV)': presynaptic,
    '5-HT rel (nM)': ser,
    'Motoneuron (mV)': vm_mn,
    'init seg (mV)': init_seg,
    'dend[2] (mV)': vm_dend,
    'Calcium (µM)': ca_muscle,
    'Muscle Force (N)': force_muscle,
    'Na soma': na_mn,
    'K soma': k_mn,
    'Na init seg': na_init_seg,
    'K init seg': k_init_seg
}


df = pd.DataFrame(data)


df.to_csv('mn_100hz_2.csv', index=False)

# ==================================================
# ####### * supra physiological raphe nuclei input * #######
# ==================================================


stim_interv = 10  # 100 Hz

h.celsius = temperature


# =====================================
# ####### *presynaptic cell* #######
# =====================================

presyn = h.Section(name="presyn")
presyn.diam = 10
presyn.L = 10
presyn.Ra = 100   # Axial resistance in Ohm * cm
presyn.cm = 1  # Membrane capacitance in micro Farads / cm^2

presyn.insert("pas")
presyn.g_pas = 1/500
presyn.e_pas = -70


presyn.insert("hh2")
presyn.ek = -90
presyn.gnabar_hh2 = 0.1
presyn.gkbar_hh2 = 0.03


presyn.insert("caL")      # HV Ca++ channel for transmitter release


rel = h.rel2(presyn(0.5))
rel.nt = 10000000
# %%

# ####### *motoneuron* #######
# =====================================

h.load_file("stdrun.hoc")
soma = h.load_file("v_e_moto6_export.hoc")


soma = h.soma
soma.Ra = 100
soma.cm = 1

soma.insert("pas")
g_pas = 1/225


soma.insert("Naf")
soma.gnafbar_Naf = 0.71
soma.ena = 50

soma.insert("KDr")
soma.gkdrbar_KDr = 0.23
soma.ek = -80

soma.insert("CaN")
soma.gcanbar_CaN = 0.013

soma.insert("KCa")
soma.gkcabar_KCa = 0.0258
soma.ek = -80

soma.insert("Ca_conc")

iCaL = h.CaL_pp(soma(0.5))


for i in range(311):
    dend = h.dend[i]
    dend.insert("pas")
    dend.g_pas = 1/11000  # [S/cm^2]
    dend.e_pas = -70        # [mV]
    dend.Ra = 70  # [Ohm*cm]
    dend.cm = 1             # [microF/cm^2]


cals = []

for i in range(311):
    dend = h.dend[i]

    cal = h.CaL_pp(dend(0.5))
    cals.append(cal)


# Axon Hillock
axon_hillock = h.Section(name="axon_hillock")
axon_hillock.L = 20
axon_hillock.diam = 13
axon_hillock.nseg = 11

axon_hillock.pt3dclear()  # to define as a cone
axon_hillock.pt3dadd(0, 0, 0, 20.65)
axon_hillock.pt3dadd(10, 0, 0, 3.83)

axon_hillock.insert("pas")
axon_hillock.g_pas = 1/11000
axon_hillock.e_pas = -70
axon_hillock.Ra = 70
axon_hillock.cm = 1


axon_hillock.insert("Naf")
axon_hillock.gnafbar_Naf = 2.7
axon_hillock.ena = 50

axon_hillock.insert("Nap")
axon_hillock.gnapbar_Nap = 0.000033
axon_hillock.ena = 50

axon_hillock.insert("KDr")
axon_hillock.gkdrbar_KDr = 0.17
axon_hillock.ek = -80


# Initial Segment
initial_segment = h.Section(name="initial_segment")
initial_segment.L = 30
initial_segment.nseg = 3
initial_segment.diam = 3.3

initial_segment.insert("pas")
initial_segment.g_pas = 1/11000
initial_segment.e_pas = -70
initial_segment.Ra = 70
initial_segment.cm = 1

initial_segment.insert("Naf")
initial_segment.gnafbar_Naf = 2.7
initial_segment.ena = 50

initial_segment.insert("Nap")
initial_segment.gnapbar_Nap = 0.000033
initial_segment.ena = 50

initial_segment.insert("KDr")
initial_segment.gkdrbar_KDr = 0.17
initial_segment.ek = -80


axon_hillock.connect(soma(1))
initial_segment.connect(axon_hillock(1))


# presyn and soma synapsis

def create_stimulation(section, tau, start, number, interval, weight, repeat):
    synapses = []
    stimuli = []
    netcons = []

    for i in range(repeat):
        syn = h.ExpSyn(section(0.5))
        syn.tau = tau

        stim = h.NetStim()
        stim.start = start
        stim.number = number
        stim.interval = stim_interv

        nc = h.NetCon(stim, syn)
        nc.weight[0] = weight
        synapses.append(syn)
        stimuli.append(stim)
        netcons.append(nc)

    return synapses, stimuli, netcons


syn_presyn, stim_presyn, nc_presyn = create_stimulation(
    presyn, 0.5, 0, 100, stim_interv, 1.5, repeat=1)
syn_soma, stim_soma, nc_soma = create_stimulation(
    soma, 0.5, 0, 100, stim_interv, 1.5, repeat=1)


# dendrites stimulation

stims = []
netcons = []

for i in range(311):
    dend = h.dend[i]

    stim_syn = h.ExpSyn(dend(0.5))
    stim_syn.tau = 0.5
    stims.append(stim_syn)

    stim1 = h.NetStim()
    stim1.start = 0
    stim1.number = 100
    stim1.interval = stim_interv

    nc_da1 = h.NetCon(stim1, stim_syn)
    nc_da1.weight[0] = 1.5
    netcons.append(nc_da1)
# %%

# add 5-HT receptors

rec1a = h.Receptor_5HT1a(initial_segment(0.5))
rec1a.gmax = 0.5
rec1a.Erev = -95
rec1a.K1 = 0.000235
rec1a.K2 = 0.002
rec1a.d2 = 0.00453


rec2 = h.Receptor_5HT2(soma(0.5))
rec2.gmax = 0.001
rec2.Erev = 0
rec2.K1 = 0.00225
rec2.K2 = 0.02

receptors = []

for i in range(311):
    dendrita_actual = h.dend[i]

    rec2 = h.Receptor_5HT2(dendrita_actual(0.5))
    rec2.gmax = 0.001
    rec2.Erev = 0
    rec2.K1 = 0.00225
    rec2.K2 = 0.02

    rec2._ref_C = rel._ref_T

    receptors.append(rec2)

rec1a._ref_C = rel._ref_T
rec2._ref_C = rel._ref_T

# %%

# ####### *muscle model* #######
# =====================================

muscle = h.Section(name='muscle')
muscle.L = 10
muscle.diam = 10
muscle.insert("pas")
muscle.g_pas = 0.002

calciumObject = h.calcium(muscle(0.5))
forceObject = h.force(muscle(0.5))
# %%

# ####### connect neuron to muscle #######
# ======================================
neuromuscularJunction = \
    h.NetCon(initial_segment(0.5)._ref_v, calciumObject, sec=initial_segment)
neuromuscularJunction.threshold = -40 * mV


# Set Pointers
h.setpointer(calciumObject._ref_A, 'aPointer', forceObject)
h.setpointer(calciumObject._ref_xm, 'xmPointer', forceObject)

# %%
# Record data for plots
pres_v = h.Vector().record(presyn(0.5)._ref_v)
ser = h.Vector().record(rel._ref_T)
v_mn = h.Vector().record(soma(0.5)._ref_v)
inseg = h.Vector().record(initial_segment(0.5)._ref_v)
dend = h.Vector().record(dendrita_actual(0.8)._ref_v)
ca = h.Vector().record(calciumObject._ref_Ca)
f = h.Vector().record(forceObject._ref_F)
t = h.Vector().record(h._ref_t)
na = h.Vector().record(soma(0.5)._ref_ina)
k = h.Vector().record(soma(0.5)._ref_ik)
na_is = h.Vector().record(initial_segment(0.5)._ref_ina)
k_is = h.Vector().record(initial_segment(0.5)._ref_ik)


# Run simulation

tstop = 1 * (stim_interv * 100) + 500 * ms
h.load_file('stdrun.hoc')
h.finitialize(-65 * mV)
h.continuerun(tstop)


# Plotting

plt.figure(figsize=(10, 10))

ax1 = plt.subplot(6, 1, 1)
ax1.plot(t, pres_v, 'b-')
ax1.axis([0, tstop, -70, 40])
ax1.set_ylabel('presyn neuron (mV)')
ax1.set_xlabel('t (ms)')

ax2 = plt.subplot(6, 1, 2)
ax2.plot(t, ser, 'b-')
ax2.axis([0, tstop, -2, 10])
ax2.set_ylabel('5-HT rel (nM)')
ax2.set_xlabel('t (ms)')


ax3 = plt.subplot(6, 1, 3)
ax3.plot(t, v_mn, 'b-')
ax3.axis([0, tstop, -100, 50])
ax3.set_ylabel('Motoneuron (mV)')
'''
ax4 = plt.subplot(11, 1, 4)
ax4.plot(t, na, 'b-')
ax4.axis([0, tstop, -10, 1])
ax4.set_ylabel('Na soma')
ax4.set_xlabel('t (ms)')

ax5 = plt.subplot(11, 1, 5)
ax5.plot(t, k, 'b-')
ax5.axis([0, tstop, -1, 5])
ax5.set_ylabel('K soma')
ax5.set_xlabel('t (ms)')

ax6 = plt.subplot(11, 1, 6)
ax6.plot(t, inseg, 'b-')
ax6.axis([0, tstop, -70, 50])
ax6.set_ylabel('init seg (mV)')
ax6.set_xlabel('t (ms)')

ax7 = plt.subplot(11, 1, 7)
ax7.plot(t, na_is, 'b-')
ax7.axis([0, tstop, -30, 1])
ax7.set_ylabel('Na init seg')
ax7.set_xlabel('t (ms)')

ax8 = plt.subplot(11, 1, 8)
ax8.plot(t, k_is, 'b-')
ax8.axis([0, tstop, -1, 2.5])
ax8.set_ylabel('K init seg')
ax8.set_xlabel('t (ms)')
'''
ax9 = plt.subplot(6, 1, 4)
ax9.plot(t, ca / µM, 'b-')
ax9.axis([0, tstop, 0, 0.1])
ax9.set_ylabel('Calcium (ÂµM)')

ax10 = plt.subplot(6, 1, 5)
ax10.plot(t, f, 'b-')
ax10.axis([0, tstop, 0, 25])
ax10.set_ylabel('Muscle Force (N)')
ax10.set_xlabel('t (ms)')

ax11 = plt.subplot(6, 1, 6)
ax11.plot(t, dend, 'b-')
ax11.axis([0, tstop, -90, 100])
ax11.set_ylabel('dend[2] (mV)')
ax11.set_xlabel('t (ms)')

plt.show()
# %%
# save data

time = np.array(t)
presynaptic = np.array(pres_v)
ser = np.array(ser)
vm_mn = np.array(v_mn)
init_seg = np.array(inseg)
vm_dend = np.array(dend)
ca_muscle = np.array(ca)
force_muscle = np.array(f)
na_mn = np.array(na)
k_mn = np.array(k)
na_init_seg = np.array(na_is)
k_init_seg = np.array(k_is)


data = {
    't (ms)': time,
    'presyn neuron (mV)': presynaptic,
    '5-HT rel (nM)': ser,
    'Motoneuron (mV)': vm_mn,
    'init seg (mV)': init_seg,
    'dend[2] (mV)': vm_dend,
    'Calcium (µM)': ca_muscle,
    'Muscle Force (N)': force_muscle,
    'Na soma': na_mn,
    'K soma': k_mn,
    'Na init seg': na_init_seg,
    'K init seg': k_init_seg
}


df = pd.DataFrame(data)


df.to_csv('mn_100hz_3.csv', index=False)