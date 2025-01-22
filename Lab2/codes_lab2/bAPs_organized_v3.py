import os
import sys
from datetime import datetime

from neuron import h, gui
import numpy as np
from matplotlib import pyplot as plt
import math
from cell import *


path = './'
spine_density = 0.05


seed = 100 
MAX_JITTER = 0 #synchronous activation
cluster_type = "sync_clusters"


np.random.seed(seed)

class config_params():
    pass

config = config_params()

config.CLUSTER_TYPE = None
config.TAU_1_AMPA = 0.3
config.TAU_2_AMPA = 1.8
config.TAU_1_NMDA = 8.019 
config.TAU_2_NMDA = 34.9884
config.N_NMDA = 0.28011
config.GAMMA_NMDA = 0.0765685
config.AMPA_W = 0.00073027
config.NMDA_W = 0.00131038
config.NMDA_W_BLOCKED = 0
config.E_SYN = 0

config.Spike_time = 200
config.SPINE_HEAD_X = 1
config.CLUSTER_L = 20
"""---------- PLEASE PLAY WITH THIS PATAMETER ----------"""
config.CLUSTER_SIZE = 0
"""-----------------------------------------------------"""


#%%

class SimulationNeuron:
    def __init__(self, config, *args):
        self.timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        self.save_path = f'./{self.timestamp}'
        os.makedirs(self.save_path, exist_ok=True)
        self.config = config
        self.create_neuron(*args)
        self.setup_simulation()
        self.create_recorders()
        self.surfix = ''

    def create_neuron(self, path, spine_density):
        self.cell = HPC(path, spine_density)
        self.cell.add_full_spine(self.cell.HCell, 0.25, 1.35, 2.8, self.cell.HCell.soma[0].Ra)

    def visualize_cell(self):
        # use h.plotshape() to visualize the cell
        h.PlotShape(False).plot(plt)
        plt.savefig(f'cell_structure_{self.timestamp}.png')


    def setup_simulation(self):
        h.steps_per_ms = 25
        h.dt = 1.0 / h.steps_per_ms
        h.celsius = 37
        h.v_init = -86
        h.tstop = 500

        self.rd = h.Random(seed)

    def create_recorders(self):
        self.rec_t = h.Vector().record(h._ref_t)
        self.rec_vsoma = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_v)
        # record the voltage of apic 50, 0, 1, 42, 80
        self.rec_vapic_50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_v)
        self.rec_vapic_0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_v)
        self.rec_vapic_1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_v)
        self.rec_vapic_42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_v)
        self.rec_vapic_80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_v)

        #TODO: record the state var, conductance and currrent of soma and apic 50, 0, 1, 42, 80 (如果觉得多就只看Na和Ca，K可以不看)
        #currents
        self.ina_soma = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_ina)
        self.ik_soma = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_ik)
        self.ica_soma = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_ica)

        self.ina_apic50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_ina)
        self.ik_apic50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_ik)
        self.ica_apic50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_ica)

        self.ina_apic0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_ina)
        self.ik_apic0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_ik)
        self.ica_apic0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_ica)

        self.ina_apic1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_ina)
        self.ik_apic1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_ik)
        self.ica_apic1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_ica)

        self.ina_apic42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_ina)
        self.ik_apic42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_ik)
        self.ica_apic42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_ica)

        self.ina_apic80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_ina)
        self.ik_apic80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_ik)
        self.ica_apic80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_ica)

        #state var
        self.soma_na_m = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_m_NaTa_t)#Na+ activation
        self.soma_na_h = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_h_NaTa_t)#Na+ inactivation
        self.soma_nap_m = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_m_Nap_Et2)
        self.soma_nap_h = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_h_Nap_Et2)
        self.soma_k_z = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_z_SK_E2)#SK_E2,calcium-activated potassium channel 
        self.soma_k_m = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_m_SKv3_1)
        self.soma_kp_m = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_m_K_Pst)
        self.soma_kp_h = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_h_K_Pst)
        self.soma_kt_m = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_m_K_Tst)
        self.soma_kt_h = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_h_K_Tst)
        self.soma_cah_m = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_m_Ca_HVA)
        self.soma_cah_h = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_h_Ca_HVA)
        self.soma_cal_m = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_m_Ca_LVAst)
        self.soma_cal_h = h.Vector().record(self.cell.HCell.soma[0](0.5)._ref_h_Ca_LVAst)

        self.apic_na_m50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_m_NaTa_t)#Na+ activation
        self.apic_na_h50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_h_NaTa_t)#Na+ inactivation
        self.apic_k_z50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_z_SK_E2)#SK_E2,calcium-activated potassium channel 
        self.apic_k_m50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_m_SKv3_1)
        self.apic_cah_m50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_m_Ca_HVA)
        self.apic_cah_h50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_h_Ca_HVA)
        self.apic_cal_m50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_m_Ca_LVAst)
        self.apic_cal_h50 = h.Vector().record(self.cell.HCell.apic[50](0.5)._ref_h_Ca_LVAst)

        self.apic_na_m0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_m_NaTa_t)#Na+ activation
        self.apic_na_h0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_h_NaTa_t)#Na+ inactivation
        self.apic_k_z0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_z_SK_E2)#SK_E2,calcium-activated potassium channel 
        self.apic_k_m0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_m_SKv3_1)
        self.apic_cah_m0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_m_Ca_HVA)
        self.apic_cah_h0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_h_Ca_HVA)
        self.apic_cal_m0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_m_Ca_LVAst)
        self.apic_cal_h0 = h.Vector().record(self.cell.HCell.apic[0](0.5)._ref_h_Ca_LVAst)

        self.apic_na_m1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_m_NaTa_t)#Na+ activation
        self.apic_na_h1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_h_NaTa_t)#Na+ inactivation
        self.apic_k_z1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_z_SK_E2)#SK_E2,calcium-activated potassium channel 
        self.apic_k_m1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_m_SKv3_1)
        self.apic_cah_m1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_m_Ca_HVA)
        self.apic_cah_h1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_h_Ca_HVA)
        self.apic_cal_m1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_m_Ca_LVAst)
        self.apic_cal_h1 = h.Vector().record(self.cell.HCell.apic[1](0.5)._ref_h_Ca_LVAst)

        self.apic_na_m42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_m_NaTa_t)#Na+ activation
        self.apic_na_h42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_h_NaTa_t)#Na+ inactivation
        self.apic_k_z42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_z_SK_E2)#SK_E2,calcium-activated potassium channel 
        self.apic_k_m42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_m_SKv3_1)
        self.apic_cah_m42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_m_Ca_HVA)
        self.apic_cah_h42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_h_Ca_HVA)
        self.apic_cal_m42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_m_Ca_LVAst)
        self.apic_cal_h42 = h.Vector().record(self.cell.HCell.apic[42](0.5)._ref_h_Ca_LVAst)

        self.apic_na_m80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_m_NaTa_t)#Na+ activation
        self.apic_na_h80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_h_NaTa_t)#Na+ inactivation
        self.apic_k_z80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_z_SK_E2)#SK_E2,calcium-activated potassium channel 
        self.apic_k_m80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_m_SKv3_1)
        self.apic_cah_m80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_m_Ca_HVA)
        self.apic_cah_h80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_h_Ca_HVA)
        self.apic_cal_m80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_m_Ca_LVAst)
        self.apic_cal_h80 = h.Vector().record(self.cell.HCell.apic[80](0.5)._ref_h_Ca_LVAst)

        #conductance
        self.gbar_soma_Ca_LVAst = self.cell.HCell.soma[0].gCa_LVAstbar_Ca_LVAst
        self.gbar_soma_Ca_HVA = self.cell.HCell.soma[0].gCa_HVAbar_Ca_HVA 
        self.gbar_soma_SKv3_1 =  self.cell.HCell.soma[0].gSKv3_1bar_SKv3_1
        self.gbar_soma_SK_E2 = self.cell.HCell.soma[0].gSK_E2bar_SK_E2
        self.gbar_soma_K_Tst = self.cell.HCell.soma[0].gK_Tstbar_K_Tst
        self.gbar_soma_K_Pst = self.cell.HCell.soma[0].gK_Pstbar_K_Pst 
        self.gbar_soma_Nap_Et2 = self.cell.HCell.soma[0].gNap_Et2bar_Nap_Et2 
        self.gbar_soma_NaTa_t = self.cell.HCell.soma[0].gNaTa_tbar_NaTa_t 

        self.gbar_apic_SK_E2 = self.cell.HCell.apic[50].gSK_E2bar_SK_E2 
        self.gbar_apic_SKv3_1 = self.cell.HCell.apic[50].gSKv3_1bar_SKv3_1 
        self.gbar_apic_NaTa_t = self.cell.HCell.apic[50].gNaTa_tbar_NaTa_t 

        self.gsoma_Ca_LVAst =self.gbar_soma_Ca_LVAst *self.soma_cal_m *self.soma_cal_m *self.soma_cal_h
        self.gsoma_Ca_HVA =self.gbar_soma_Ca_HVA *self.soma_cah_m *self.soma_cah_m *self.soma_cah_h
        self.gsoma_skv3 =self.gbar_soma_SKv3_1 *self.soma_k_m
        self.gsoma_ske2 =self.gbar_soma_SK_E2 *self.soma_k_z
        self.gsoma_kt =self.gbar_soma_K_Tst *self.soma_kt_m *self.soma_kt_m *self.soma_kt_m *self.soma_kt_m *self.soma_kt_h
        self.gsoma_kp =self.gbar_soma_K_Pst *self.soma_kp_m *self.soma_kp_m *self.soma_kp_h
        self.gsoma_nap =self.gbar_soma_Nap_Et2 *self.soma_nap_m *self.soma_nap_m *self.soma_nap_m *self.soma_nap_h
        self.gsoma_nata =self.gbar_soma_NaTa_t *self.soma_na_m *self.soma_na_m *self.soma_na_m *self.soma_na_h

        self.gapic_ske2_50 = self.gbar_apic_SK_E2 * self.apic_k_z50
        self.gapic_skv3_50 = self.gbar_apic_SKv3_1 * self.apic_k_m50
        self.gapic_nata_50 = self.gbar_apic_NaTa_t * self.apic_na_m50 * self.apic_na_m50 * self.apic_na_m50 * self.apic_na_h50

        self.gapic_ske2_0 = self.gbar_apic_SK_E2 * self.apic_k_z0
        self.gapic_skv3_0 = self.gbar_apic_SKv3_1 * self.apic_k_m0
        self.gapic_nata_0 = self.gbar_apic_NaTa_t * self.apic_na_m0 * self.apic_na_m0 * self.apic_na_m0 * self.apic_na_h0

        self.gapic_ske2_1 = self.gbar_apic_SK_E2 * self.apic_k_z1
        self.gapic_skv3_1 = self.gbar_apic_SKv3_1 * self.apic_k_m1
        self.gapic_nata_1 = self.gbar_apic_NaTa_t * self.apic_na_m1 * self.apic_na_m1 * self.apic_na_m1 * self.apic_na_h1
 
        self.gapic_ske2_42 = self.gbar_apic_SK_E2 * self.apic_k_z42
        self.gapic_skv3_42 = self.gbar_apic_SKv3_1 * self.apic_k_m42
        self.gapic_nata_42 = self.gbar_apic_NaTa_t * self.apic_na_m42 * self.apic_na_m42 * self.apic_na_m42 * self.apic_na_h42

        self.gapic_ske2_80 = self.gbar_apic_SK_E2 * self.apic_k_z80
        self.gapic_skv3_80 = self.gbar_apic_SKv3_1 * self.apic_k_m80
        self.gapic_nata_80 = self.gbar_apic_NaTa_t * self.apic_na_m80 * self.apic_na_m80 * self.apic_na_m80 * self.apic_na_h80

    def remove_proximal_apic_sodium_channels(self):
        for sec in h.allsec():
            if 'apic' in sec.name():
                for seg in sec:
                    if h.distance(self.cell.HCell.soma[0](0.5), seg) < 600:
                        seg.gNaTa_tbar_NaTa_t = 0
        self.surfix += 'remove_proximal_apic_sodium_channels_'

    def display_apic_distance(self):   
        # for i in range(0, number of apic(read from file, not a constant)):
        apic_num = 108
        for i in range(0, apic_num+1):
            dist_0 = h.distance(self.cell.HCell.soma[0](0.5), self.cell.HCell.apic[i](0))
            dist_1 = h.distance(self.cell.HCell.soma[0](0.5), self.cell.HCell.apic[i](1))
            print(f"Distance to apic[{i}]: {dist_0}, {dist_1}")

    def create_stim_soma(self, dur, amp):

        self.stim_soma = h.IClamp(self.cell.HCell.soma[0](0.5))
        self.stim_soma.delay = self.config.Spike_time

        self.stim_soma.dur = dur
        self.stim_soma.amp = amp
        
        self.surfix += f'soma_dur_{dur}_'
        self.surfix += f'soma_amp_{amp}_'

        

    def create_stim_apic(self, start):
        self.Stim_1 = h.NetStim()
        self.Stim_1.interval = 1e9
        self.Stim_1.start = start
        self.Stim_1.noise = 0
        self.Stim_1.number = 1
        self.config.stim = self.Stim_1

        self.surfix += f'apic_start_{start}_'

    def add_single_synapse(self, seg):
        synaptic_segments = self.cell.fill_synapse_list_with_spine([seg], self.config)
        self.cell.add_synapses_on_list_of_segments(synaptic_segments, self.cell.synlist, self.cell.conlist, self.config)
        configure_synaptic_delayes(MAX_JITTER, self.cell.conlist, self.rd, self.config, cluster_type=None)

        self.surfix += f'single_synapse_{seg.sec.name()}_'

    def add_clustered_synapses(self, seg, cluster_size):
        self.config.CLUSTER_SIZE = cluster_size
        synaptic_segments = self.cell.fill_clustered_synapses_list_with_spine([seg], self.rd, self.config)
        self.cell.add_synapses_on_list_of_segments(synaptic_segments, self.cell.synlist, self.cell.conlist, self.config)
        configure_synaptic_delayes(MAX_JITTER, self.cell.conlist, self.rd, self.config, cluster_type=None)

        self.surfix += f'clustered_synapses_{seg.sec.name()}_'


    def add_syn_on_shaft(self, seg, start, delay=0, tau1=0.5, tau2=5, e=0, i=0.5,weight=0.00073027):

        self.syn_shaft = h.Exp2Syn(seg)
        self.syn_shaft.tau1 = tau1
        self.syn_shaft.tau2 = tau2
        self.syn_shaft.e = e
        self.syn_shaft.i = i

        self.Stim_2 = h.NetStim()
        self.Stim_2.number = 1
        self.Stim_2.start = start

        self.nc2 = h.NetCon(self.Stim_2, self.syn_shaft)
        self.nc2.delay = delay
        self.nc2.weight[0] = weight
        
        # add name of seg to surfix 
        self.surfix += f'syn_on_shaft_{seg.sec.name()}_'
        self.surfix += f'syn_on_shaft_start_{start}_'
        self.surfix += f'syn_on_shaft_i_{i}_'


    def run_simulation(self):
        h.finitialize(h.v_init)
        h.continuerun(h.tstop)

    def plot_voltage(self):
        plt.figure(figsize=(15, 8))
        plt.plot(self.rec_t, self.rec_vsoma, label='soma')
        plt.plot(self.rec_t, self.rec_vapic_50, label='apic 50')
        plt.plot(self.rec_t, self.rec_vapic_0, label='apic 0')
        plt.plot(self.rec_t, self.rec_vapic_1, label='apic 1')
        plt.plot(self.rec_t, self.rec_vapic_42, label='apic 42')
        plt.plot(self.rec_t, self.rec_vapic_80, label='apic 80')
        plt.xlabel('Time (ms)')
        plt.ylim(-90, 40)
        plt.ylabel('Voltage (mV)')
        plt.legend()
        # plt.savefig(f'voltage_response_{}.png')
        # plt.savefig(f'{self.save_path}/voltage_response.png')
        plt.savefig(f'{self.save_path}/voltage_response_.png')
        # create a txt to store surfix
        with open(f'{self.save_path}/surfix.txt', 'w') as f:
            f.write(self.surfix)


     # TODO: plot state var, conductance and current of soma and apic 50, 0, 1, 42, 80   

    def plot_state_var(self):
        #pass
        plt.figure(figsize=(15, 8))

        plt.subplot(3,1,1)
        plt.plot(self.rec_t, self.soma_na_m, label='m (NaTa_t activation)')
        plt.plot(self.rec_t, self.soma_na_h, label='h (NaTa_t inactivation)')
        plt.plot(self.rec_t, self.soma_nap_m, label='m (Nap_Et2 activation)')
        plt.plot(self.rec_t, self.soma_nap_h, label='h (Nap_Et2 inactivation)')
        plt.xlabel('time (ms)')
        plt.ylabel('State Variable of Na')
        plt.legend()

        plt.subplot(3,1,2)
        plt.plot(self.rec_t, self.soma_k_z, label='z(SK_E2)')
        plt.plot(self.rec_t, self.soma_k_m, label='m(SKv3_1)')
        plt.plot(self.rec_t, self.soma_kp_m, label='m (K_Pst activation)')
        plt.plot(self.rec_t, self.soma_kp_h, label='h (K_Pst inactivation)')
        plt.plot(self.rec_t, self.soma_kt_m, label='m (K_Tst activation)')
        plt.plot(self.rec_t, self.soma_kt_h, label='h (K_Tst inactivation)')
        plt.xlabel('time (ms)')
        plt.ylabel('State Variable of K')
        plt.legend()

        plt.subplot(3,1,3)
        plt.plot(self.rec_t, self.soma_cah_m, label='m (CaH activation)')
        plt.plot(self.rec_t, self.soma_cah_h, label='h (CaH inactivation)')
        plt.plot(self.rec_t, self.soma_cal_m, label='m (CaL activation)')
        plt.plot(self.rec_t, self.soma_cal_h, label='h (CaL inactivation)')
        plt.xlabel('time (ms)')
        plt.ylabel('State Variable of Ca')

        plt.title(f'State Variables Over Time In Soma')
        plt.legend()
        plt.savefig(f'{self.save_path}/state_var_soma_.png')

        plt.figure(figsize=(15, 8))

        plt.subplot(3,1,1)
        plt.plot(self.rec_t, self.apic_na_m50, label='m (NaTa_t_50 activation)')
        plt.plot(self.rec_t, self.apic_na_h50, label='h (NaTa_t_50 inactivation)')
        plt.plot(self.rec_t, self.apic_na_m0, label='m (NaTa_t_0 activation)')
        plt.plot(self.rec_t, self.apic_na_h0, label='h (NaTa_t_0 inactivation)')
        plt.plot(self.rec_t, self.apic_na_m1, label='m (NaTa_t_1 activation)')
        plt.plot(self.rec_t, self.apic_na_h1, label='h (NaTa_t_1 inactivation)')
        plt.plot(self.rec_t, self.apic_na_m42, label='m (NaTa_t_42 activation)')
        plt.plot(self.rec_t, self.apic_na_h42, label='h (NaTa_t_42 inactivation)')
        plt.plot(self.rec_t, self.apic_na_m80, label='m (NaTa_t_80 activation)')
        plt.plot(self.rec_t, self.apic_na_h80, label='h (NaTa_t_80 inactivation)')


        plt.xlabel('time (ms)')
        plt.ylabel('State Variable of Na')
        plt.legend()

        plt.subplot(3,1,2)
        plt.plot(self.rec_t, self.apic_k_z50, label='z(SK_E2_50)')
        plt.plot(self.rec_t, self.apic_k_m50, label='m(SKv3_1_50)')
        plt.plot(self.rec_t, self.apic_k_z0, label='z(SK_E2_0)')
        plt.plot(self.rec_t, self.apic_k_m0, label='m(SKv3_1_0)')
        plt.plot(self.rec_t, self.apic_k_z1, label='z(SK_E2_1)')
        plt.plot(self.rec_t, self.apic_k_m1, label='m(SKv3_1_1)')
        plt.plot(self.rec_t, self.apic_k_z42, label='z(SK_E2_42)')
        plt.plot(self.rec_t, self.apic_k_m42, label='m(SKv3_1_42)')
        plt.plot(self.rec_t, self.apic_k_z80, label='z(SK_E2_80)')
        plt.plot(self.rec_t, self.apic_k_m80, label='m(SKv3_1_80)')
        plt.xlabel('time (ms)')
        plt.ylabel('State Variable of K')
        plt.legend()

        plt.subplot(3,1,3)
        plt.plot(self.rec_t, self.apic_cah_m50, label='m (CaH_50 activation)')
        plt.plot(self.rec_t, self.apic_cah_h50, label='h (CaH_50 inactivation)')
        plt.plot(self.rec_t, self.apic_cal_m50, label='m (CaL_50 activation)')
        plt.plot(self.rec_t, self.apic_cal_h50, label='h (CaL_50 inactivation)')
        plt.plot(self.rec_t, self.apic_cah_m0, label='m (CaH_0 activation)')
        plt.plot(self.rec_t, self.apic_cah_h0, label='h (CaH_0 inactivation)')
        plt.plot(self.rec_t, self.apic_cal_m0, label='m (CaL_0 activation)')
        plt.plot(self.rec_t, self.apic_cal_h0, label='h (CaL_0 inactivation)')
        plt.plot(self.rec_t, self.apic_cah_m1, label='m (CaH_1 activation)')
        plt.plot(self.rec_t, self.apic_cah_h1, label='h (CaH_1 inactivation)')
        plt.plot(self.rec_t, self.apic_cal_m1, label='m (CaL_1 activation)')
        plt.plot(self.rec_t, self.apic_cal_h1, label='h (CaL_1 inactivation)')
        plt.plot(self.rec_t, self.apic_cah_m42, label='m (CaH_42 activation)')
        plt.plot(self.rec_t, self.apic_cah_h42, label='h (CaH_42 inactivation)')
        plt.plot(self.rec_t, self.apic_cal_m42, label='m (CaL_42 activation)')
        plt.plot(self.rec_t, self.apic_cal_h42, label='h (CaL_42 inactivation)')
        plt.plot(self.rec_t, self.apic_cah_m80, label='m (CaH_80 activation)')
        plt.plot(self.rec_t, self.apic_cah_h80, label='h (CaH_80 inactivation)')
        plt.plot(self.rec_t, self.apic_cal_m80, label='m (CaL_80 activation)')
        plt.plot(self.rec_t, self.apic_cal_h80, label='h (CaL_80 inactivation)')
        plt.xlabel('time (ms)')
        plt.ylabel('State Variable of Ca')

        plt.title(f'State Variables Over Time In Apic')
        plt.legend()
        plt.savefig(f'{self.save_path}/state_var_apic_.png')


    def plot_conductance(self):
        #pass
        plt.figure(figsize=(15, 8))

        plt.subplot(3,1,1)
        plt.plot(self.rec_t, self.gbar_soma_NaTa_t *self.soma_na_m *self.soma_na_m *self.soma_na_m *self.soma_na_h, label='gNaTa_t')
        plt.plot(self.rec_t, self.gbar_soma_Nap_Et2 *self.soma_nap_m *self.soma_nap_m *self.soma_nap_m *self.soma_nap_h, label='gNap_Et2')
        plt.xlabel('time (ms)')
        plt.ylabel('Conductance of Na')
        plt.legend()

        plt.subplot(3,1,2)
        plt.plot(self.rec_t, self.gbar_soma_SK_E2 *self.soma_k_z, label='gSK_E2')
        plt.plot(self.rec_t, self.gbar_soma_SKv3_1 *self.soma_k_m, label='gSKv3_1')
        plt.plot(self.rec_t, self.gbar_soma_K_Pst *self.soma_kp_m *self.soma_kp_m *self.soma_kp_h, label='gK_Pst')
        plt.plot(self.rec_t, self.gbar_soma_K_Tst *self.soma_kt_m *self.soma_kt_m *self.soma_kt_m *self.soma_kt_m *self.soma_kt_h, label='gK_Tst')
        plt.xlabel('time (ms)')
        plt.ylabel('Conductance of K')
        plt.legend()

        plt.subplot(3,1,3)
        plt.plot(self.rec_t, self.gbar_soma_Ca_HVA *self.soma_cah_m *self.soma_cah_m *self.soma_cah_h, label='gCa_HVA')
        plt.plot(self.rec_t, self.gbar_soma_Ca_LVAst *self.soma_cal_m *self.soma_cal_m *self.soma_cal_h, label='gCa_LVAst')
        plt.xlabel('time (ms)')
        plt.ylabel('Conductance of Ca')
        plt.legend()

        plt.title(f'Conductance Over Time In Soma')

        plt.legend()
        plt.savefig(f'{self.save_path}/conductance_soma_.png')

        plt.figure(figsize=(15, 8))

        plt.plot(self.rec_t, self.gbar_apic_NaTa_t * self.apic_na_m50 * self.apic_na_m50 * self.apic_na_m50 * self.apic_na_h50, label='gNaTa_t_50')
        plt.plot(self.rec_t, self.gbar_apic_SK_E2 * self.apic_k_z50, label='gSK_E2_50')
        plt.plot(self.rec_t, self.gbar_apic_SKv3_1 * self.apic_k_m50, label='gSKv3_1_50')

        plt.plot(self.rec_t, self.gbar_apic_NaTa_t * self.apic_na_m0 * self.apic_na_m0 * self.apic_na_m0 * self.apic_na_h0, label='gNaTa_t_0')
        plt.plot(self.rec_t, self.gbar_apic_SK_E2 * self.apic_k_z0, label='gSK_E2_0')
        plt.plot(self.rec_t, self.gbar_apic_SKv3_1 * self.apic_k_m0, label='gSKv3_1_0')

        plt.plot(self.rec_t, self.gbar_apic_NaTa_t * self.apic_na_m1 * self.apic_na_m1 * self.apic_na_m1 * self.apic_na_h1, label='gNaTa_t_1')
        plt.plot(self.rec_t, self.gbar_apic_SK_E2 * self.apic_k_z1, label='gSK_E2_1')
        plt.plot(self.rec_t, self.gbar_apic_SKv3_1 * self.apic_k_m1, label='gSKv3_1_1')

        plt.plot(self.rec_t, self.gbar_apic_NaTa_t * self.apic_na_m42 * self.apic_na_m42 * self.apic_na_m42 * self.apic_na_h42, label='gNaTa_t_42')
        plt.plot(self.rec_t, self.gbar_apic_SK_E2 * self.apic_k_z42, label='gSK_E2_42')
        plt.plot(self.rec_t, self.gbar_apic_SKv3_1 * self.apic_k_m42, label='gSKv3_1_42')

        plt.plot(self.rec_t, self.gbar_apic_NaTa_t * self.apic_na_m80 * self.apic_na_m80 * self.apic_na_m80 * self.apic_na_h80, label='gNaTa_t_80')
        plt.plot(self.rec_t, self.gbar_apic_SK_E2 * self.apic_k_z80, label='gSK_E2_80')
        plt.plot(self.rec_t, self.gbar_apic_SKv3_1 * self.apic_k_m80, label='gSKv3_1_80')
        plt.xlabel('time (ms)')
        plt.ylabel('Conductance')

        plt.title(f'Conductance Over Time In Apic')

        plt.legend()
        plt.savefig(f'{self.save_path}/conductance_apic_.png')



    def plot_current(self):
        plt.figure(figsize=(15, 8))
        plt.plot(self.rec_t, self.ina_soma, label='INa (nA)')
        plt.plot(self.rec_t, self.ik_soma, label='IK (nA)')
        plt.plot(self.rec_t, self.ica_soma, label='ICa (nA)')
        plt.xlabel('time (ms)')
        plt.ylabel('Current (nA)')
        plt.title(f'INa, IK and ICaT Over Time in Soma')
        plt.legend()
        plt.savefig(f'{self.save_path}/current_soma_.png')

        plt.figure(figsize=(15, 8))
        plt.plot(self.rec_t, self.ina_apic50, label='INa_50 (nA)')
        plt.plot(self.rec_t, self.ik_apic50, label='IK_50 (nA)')
        plt.plot(self.rec_t, self.ica_apic50, label='ICa_50 (nA)')
        plt.plot(self.rec_t, self.ina_apic0, label='INa_0 (nA)')
        plt.plot(self.rec_t, self.ik_apic0, label='IK_0 (nA)')
        plt.plot(self.rec_t, self.ica_apic0, label='ICa_0 (nA)')
        plt.plot(self.rec_t, self.ina_apic1, label='INa_1 (nA)')
        plt.plot(self.rec_t, self.ik_apic1, label='IK_1 (nA)')
        plt.plot(self.rec_t, self.ica_apic1, label='ICa_1 (nA)')
        plt.plot(self.rec_t, self.ina_apic42, label='INa_42 (nA)')
        plt.plot(self.rec_t, self.ik_apic42, label='IK_42 (nA)')
        plt.plot(self.rec_t, self.ica_apic42, label='ICa_42 (nA)')
        plt.plot(self.rec_t, self.ina_apic80, label='INa_80 (nA)')
        plt.plot(self.rec_t, self.ik_apic80, label='IK_80 (nA)')
        plt.plot(self.rec_t, self.ica_apic80, label='ICa_80 (nA)')
        plt.xlabel('time (ms)')
        plt.ylabel('Current (nA)')
        plt.title(f'INa, IK and ICaT Over Time in Apic')
        plt.legend()
        plt.savefig(f'{self.save_path}/current_apic_.png')
        #pass

    def display_distance(self):
        dist_0 = h.distance(self.cell.HCell.soma[0](0.5), self.cell.HCell.apic[50](0))
        dist_1 = h.distance(self.cell.HCell.soma[0](0.5), self.cell.HCell.apic[50](1))
        print(f"Distance to dend]: {dist_0}, {dist_1}")



#%%

pyr_neuron = SimulationNeuron(config, path, spine_density)

pyr_neuron.display_distance()
# pyr_neuron.remove_proximal_apic_sodium_channels()
#pyr_neuron.display_apic_distance()
## change parameters
#pyr_neuron.create_stim_soma(dur=5, amp=1.9)
#pyr_neuron.create_stim_apic(start=200)
# pyr_neuron.add_single_synapse(seg=pyr_neuron.cell.HCell.apic[50](0.5))
#pyr_neuron.add_clustered_synapses(seg=pyr_neuron.cell.HCell.dend[50](0.5), cluster_size=20)
# pyr_neuron.add_syn_on_shaft(seg=pyr_neuron.cell.HCell.dend[50](0.5), start=250, i=0.5)


#pyr_neuron.run_simulation()
#pyr_neuron.plot_conductance()
#pyr_neuron.plot_voltage()
#pyr_neuron.plot_state_var()
#pyr_neuron.plot_current()




