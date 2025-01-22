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
soma.insert('hh1') 
#soma.insert('cat1i')
#soma.insert('cat1g')
soma.insert('cat1h')

soma(0.5).hh1.gnabar = 0.12  # Sodium conductance in S/cm^2
soma(0.5).hh1.gkbar = 0.036  # Potassium conductance in S/cm^2
soma(0.5).hh1.gl = 0.0003  # Leak conductance in S/cm^2
#soma(0.5).hh1.el = -74

# Set up current injection
I = h.IClamp(soma(0.5))  # Inject current at the center of the soma
I.delay = 2  # Start time of current injection
I.dur = 40  # Duration of current injection in ms
I.amp = 0.1  # Amplitude of current injection in nA

t = h.Vector().record(h._ref_t)  # Time
v = h.Vector().record(soma(0.5)._ref_v)  # Membrane potential
ina = h.Vector().record(soma(0.5)._ref_ina)  # Sodium current
ik = h.Vector().record(soma(0.5)._ref_ik)  # Potassium current
ica = h.Vector().record(soma(0.5)._ref_ica) 
m = h.Vector().record(soma(0.5).hh1._ref_m)  # Na activation
h_gate = h.Vector().record(soma(0.5).hh1._ref_h)  # Na inactivation
n = h.Vector().record(soma(0.5).hh1._ref_n)  # K activation

# Run the simulation
h.load_file('stdrun.hoc')
h.finitialize(-65*mV)
h.continuerun(40*ms)

plt.figure(figsize=(12, 8))

# Plot voltage response
plt.subplot(3, 1, 1)
plt.plot(t, v, label='Membrane Potential (mV)')
plt.title('Membrane Potential with T-type Calcium Channels')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.grid()

# Plot currents
plt.subplot(3, 1, 2)
plt.plot(t, ina, label='INa (nA)', color='red')
plt.plot(t, ik, label='IK (nA)', color='blue')
plt.plot(t, ica, label='ICaT (nA)', color='green')
plt.title('Ionic Currents')
plt.xlabel('Time (ms)')
plt.ylabel('Current (nA)')
plt.legend()
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(t, m, label='m', color='orange')
plt.plot(t, h_gate, label='h', color='purple')
plt.plot(t, n, label='n', color='green')
plt.title('State Variables')
plt.xlabel('Time (ms)')
plt.ylabel('Gating Variables')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()