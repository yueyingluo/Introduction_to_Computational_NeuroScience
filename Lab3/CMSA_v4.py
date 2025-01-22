import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import numpy as np

torch.set_default_dtype(torch.float16)  # 由于LIF神经元相对比较简单，可以使用float16以节省显存
if torch.cuda.is_available():  # 优先使用gpu进行模拟
    device = torch.device('cuda')
else:
    device = torch.device('cpu')

print(f'use device:{device}')
torch.set_grad_enabled(False) # 不使用自动梯度
torch.cuda.empty_cache() # 清空显存
dt = 0.5  # 模拟时间步长 ms
__t__ = 0  # 全局时间变量，用于计算不应期

class LIFlayer:
    def __init__(self, n:int, threshold=-50.0, reset_value=-70.0, membrane_capacitance=1.0, gL=0.05, refractory_period=5):
        self.shape = (n,n)  # 网络形态
        self.threshold = threshold  # 发放阈值
        self.reset_value = reset_value  # 回复和静息电位，单位mV
        self.membrane_capacitance = membrane_capacitance  # 膜电容，单位nF
        self.gL = gL  # 漏电导，单位μS
        self.refractory = refractory_period  # 不应期，单位ms
        self.potential = torch.ones(self.shape).to(device) * (self.reset_value)  # 膜电位，单位mV
        self.spike = torch.zeros(self.shape).to(device)
        self.spike_time = torch.ones(self.shape).to(device) * (-self.refractory*5)  # 上一次的发放时间

    def update(self, input:torch.Tensor):
        assert input.shape == self.shape
        #TODO 请你完成膜电位、神经元发放和不应期的更新
        # 计算电压变化率 dV/dt
        dV_dt = (-self.gL * (self.potential - self.reset_value) + input) / self.membrane_capacitance
        # dt = 1.0  # 时间步长，可以根据需要调整

        # 更新膜电位 V(t) = V(t-1) + dV_dt * dt
        self.potential += dV_dt * dt

        # 判断是否满足发放条件
        is_refractory = (__t__ - self.spike_time) < self.refractory  # 在不应期内
        is_threshold_crossed = (self.potential >= self.threshold) & ~is_refractory  # 满足发放阈值且不在不应期内

        # 更新发放矩阵
        self.spike = is_threshold_crossed.float()

        # 更新发放时间
        self.spike_time[is_threshold_crossed] = __t__

        # 在发放后将膜电位重置
        self.potential[is_threshold_crossed] = self.reset_value

        # 在不应期内将膜电位保持在静息电位
        self.potential[is_refractory] = self.reset_value
        
        return self.potential, self.spike

class Synapseslayer:
    def __init__(self, in_neuron:int, out_neuron:int, m_synapses:int, W=0.02, sigma=18, time_constant=3.0, Vrest = 0.0):
        '''m_synapse must be odd. because in the unfold function, if m is even, the index will shift'''
        self.in_neurons = in_neuron
        self.out_neurons = out_neuron
        assert out_neuron/in_neuron % 1 == 0 or in_neuron/out_neuron % 1 == 0  # 确保E_neurons和I_neurons的数量比是整数，以便于后续进行缩放

        self.shape = (out_neuron, out_neuron, m_synapses, m_synapses)
        self.time_constant = time_constant
        self.weight = self.gaussian(m_synapses, W, sigma)
        self.Vrest = Vrest

        self.i =  torch.zeros(self.shape).to(device)  # 突触电流，单位nA
        self.g = torch.zeros(self.shape).to(device)  # 突触电导，单位μS

    def gaussian(self, n, W, sigma):
        #TODO 请你完成高斯波包函数，返回一个n*n矩阵，其中最大值位于正中间（n为奇数）

        # 创建坐标网格
        x = torch.linspace(-(n-1)/2, (n-1)/2, n)
        y = torch.linspace(-(n-1)/2, (n-1)/2, n)
        X, Y = torch.meshgrid(x, y)
        
        # 计算高斯函数
        gaussian = W * torch.exp(-(X**2 + Y**2)/sigma)
        gaussian = gaussian.to(device)

        return gaussian
    
    def update(self, input: torch.Tensor, potential:torch.Tensor):
        assert input.shape == (self.in_neurons, self.in_neurons)

        if self.in_neurons<self.out_neurons:
            input = self.scale_up(input,self.out_neurons//self.in_neurons)
        else:
            input = self.scale_down(input,self.in_neurons//self.out_neurons)
        
        #TODO 请你完成突触电导和电流的更新（提示：使用torch.einsum计算input和weight的逐项乘积，使用repeat将单个神经元的膜电位展开成和突触张量同阶数的张量）
        # print ("g: "+str(self.g.shape))
        # print ("w: "+str(self.weight.shape))
        # print ("in: "+str(input.shape))
        # print ("Vr: "+str(self.Vrest.shape))
        # 更新电导 g (使用指数衰减模型)
        self.g += (-self.g + torch.einsum('ij,abij->abij', self.weight, input)) / self.time_constant * dt

        # 计算膜电位的重复张量
        # print("Before: "+str(potential.shape))
        potential = potential.reshape(self.out_neurons, self.out_neurons, 1, 1)
        potential_expanded = potential.repeat(1, 1, self.weight.shape[0], self.weight.shape[1])
        # print("V: "+str(potential_expanded.shape))

        # 计算电流 i
        self.i = torch.einsum('ijkl,ijkl->ijkl', self.g, (self.Vrest - potential_expanded))

        return self.i
    
    def scale_up(self, input:torch.Tensor, zoom_rate:int):
        a1 = self.shape[3]//2
        a2 = (self.shape[3]-1)//2

        # 以下四个表达式完成了二维矩阵的扩展，以便后面进行平移展开操作
        input = torch.cat((input, input[:,:a1]), dim=1)
        input = torch.cat((input[:,-a1-a2:-a1], input), dim=1)
        input = torch.cat((input, input[:a1,:]), dim=0)
        input = torch.cat((input[-a2-a1:-a1,:], input), dim=0)

        # 平移展开，得到每个突触对应神经元的spike情况
        input = input.unfold(0, self.shape[2], 1).unfold(1, self.shape[3], 1)

        # 将较小的synapselayer的相邻元素重复，得到较大的layer的输入
        input = input.repeat_interleave(zoom_rate,dim=0).repeat_interleave(zoom_rate,dim=1)
        return input
    
    def scale_down(self, input:torch.Tensor, zoom_rate:int): 
        # 和上面的同理
        a1 = self.shape[3]//2
        a2 = (self.shape[3]-1)//2
        input = torch.cat((input, input[:,:a1]), dim=1)
        input = torch.cat((input[:,-a1-a2:-a1], input), dim=1)
        input = torch.cat((input, input[:a1,:]), dim=0)
        input = torch.cat((input[-a2-a1:-a1,:], input), dim=0)

        input = input.unfold(0, self.shape[2], zoom_rate).unfold(1, self.shape[3], zoom_rate)
        return input
    
class Network:
    def __init__(self, En, In, rp, We, Wi) -> None:
        self.E_neurons = LIFlayer(n=En,refractory_period=rp)
        self.synapsesEE = Synapseslayer(En,En,101,W=We)
        self.synapsesEI = Synapseslayer(En,In,101,W=We)

        self.I_neurons = LIFlayer(n=In, refractory_period=rp)
        self.synapsesIE = Synapseslayer(In,En,101,Vrest=-80,sigma=400,W=Wi)
        self.synapsesII = Synapseslayer(In,In,101,Vrest=-80,sigma=400,W=Wi)

    def update(self, inputE:torch.Tensor, inputI:torch.Tensor):
        E_potential, E_spike = self.E_neurons.update(inputE+self.synapsesEE.i.sum(dim=(2,3))+self.synapsesIE.i.sum(dim=(2,3)))
        I_potential, I_spike = self.I_neurons.update(inputI+self.synapsesII.i.sum(dim=(2,3))+self.synapsesEI.i.sum(dim=(2,3)))
        self.synapsesEE.update(E_spike, E_potential)
        self.synapsesEI.update(E_spike, I_potential)
        self.synapsesIE.update(I_spike, E_potential)
        self.synapsesII.update(I_spike, I_potential)

        return E_potential,E_spike
    
En = 150   # 兴奋性神经元网格边长数
In = 75   # 抑制性神经元网格边长数，En与In之比需要为整数
runtime1 = 50  # 有外界输入的时间
runtime = 100  # 自发动力学的时间
rp = 7.25  # 不应期时间
We = 0.08  # 兴奋性连接权重
Wi = 0.04  # 抑制性连接权重
net = Network(En,In,rp=rp,We=We,Wi=Wi)
voltage_list = []

for i in range(int(runtime1/dt)):
    __t__+=dt
    E_potential, E_spike= net.update(torch.rand((En,En)).to(device)*5, torch.rand((In,In)).to(device)*5)  # 平均外界输入电流2.5nA
    voltage_list.append(E_potential.clone().cpu())

for i in range(int(runtime/dt)):
    __t__+=dt
    E_potential, E_spike= net.update(torch.rand((En,En)).to(device)*0, torch.rand((In,In)).to(device)*0)
    voltage_list.append(E_potential.clone().cpu())

plt.figure(1)
plt.imshow(E_potential.cpu(), cmap='coolwarm')
plt.title(f'We={We}, Wi={Wi}, rp={rp}ms')
plt.colorbar()
plt.show()

# 可以保存数据，使用Watch_Pattern.py观看动画
V = np.array([v.numpy() for v in voltage_list])

plt.figure(2)
for neuron in [(50,55), (55, 55), (60, 55), (65, 55)]:
    plt.plot(V[:, neuron[0], neuron[1]], label = f'{neuron}')
plt.legend()
plt.title(f'Membrane Potential \n We={We}, Wi={Wi}, rp={rp}ms')

np.save(f'E={We}-I={Wi}-R={rp}',V)


