from neuron import h
from neuron.units import ms,mV
import matplotlib.pyplot as plt
import numpy as np

# Create a single section
soma = h.Section(name='soma')
soma.L = 10  # length of the section in microns
soma.diam = 10  # diameter of the section in microns
soma.Ra = 100
soma.cm = 1
soma.insert('hh1')  # insert Hodgkin-Huxley channels

soma(0.5).hh1.gnabar = 0.5  # Sodium conductance in S/cm^2
soma(0.5).hh1.gkbar = 0.036  # Potassium conductance in S/cm^2
soma(0.5).hh1.gl = 0.0003  # Leak conductance in S/cm^2

# Set up current injection
I = h.IClamp(soma(0.5))  # Inject current at the center of the soma
I.delay = 2  # Start time of current injection
I.dur = 40  # Duration of current injection in ms
I.amp = 0.5  # Amplitude of current injection in nA


# Record variables
t = h.Vector().record(h._ref_t)  # Time
v = h.Vector().record(soma(0.5)._ref_v)  # Membrane potential
ina = h.Vector().record(soma(0.5)._ref_ina)  # Sodium current
ik = h.Vector().record(soma(0.5)._ref_ik)  # Potassium current
m = h.Vector().record(soma(0.5).hh1._ref_m)  # m gate variable
h_gate = h.Vector().record(soma(0.5).hh1._ref_h)  # h gate variable
n = h.Vector().record(soma(0.5).hh1._ref_n)  # n gate variable
gNa = h.Vector().record(soma(0.5).hh1._ref_gnabar)  # Sodium conductance
gK = h.Vector().record(soma(0.5).hh1._ref_gkbar)  # Potassium conductance

# Run the simulation
h.load_file('stdrun.hoc')
h.finitialize(-65*mV)
h.continuerun(40*ms)

plt.figure(figsize=(12, 5))

# Plot voltage
plt.subplot(2, 2, 1)
plt.plot(t, v, label='Membrane Potential (mV)')
plt.title('Membrane Potential')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.grid()

# Plot currents
plt.subplot(2, 2, 2)
plt.plot(t, ina, label='INa (nA)')
plt.plot(t, ik, label='IK (nA)')
plt.title('Ionic Currents')
plt.xlabel('Time (ms)')
plt.ylabel('Current (nA)')
plt.legend()
plt.grid()

plt.subplot(2, 2, 3)
plt.plot(t, m, label='m')
plt.plot(t, h_gate, label='h')
plt.plot(t, n, label='n')
plt.title('State Variables')
plt.xlabel('Time (ms)')
plt.ylabel('Gating Variables')
plt.legend()
plt.grid()

# Plot conductance
plt.subplot(2, 2, 4)
plt.plot(t, [soma(0.5).hh1.gnabar * m[i]**3 * h_gate[i] for i in range(len(t))], label='gNa (S/cm2)')
plt.plot(t, [soma(0.5).hh1.gkbar * n[i]**4 for i in range(len(t))], label='gK (S/cm2)')
plt.title('Conductance')
plt.xlabel('Time (ms)')
plt.ylabel('Conductance (S/cm^2)')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()