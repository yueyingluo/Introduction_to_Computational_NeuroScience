'''code modified from Lab 4 Discussion, special thanks go to my teammate Jiarui Sun and Yichen Xu, and Chaoming Wang and Tianhao Chu who offer generous help to us.'''
'''This script is used to simulate the spatial CANN model with 1D feature space. 
    The model consists of a place cell network (P_CANN1D) and three grid cell networks (G_CANN1D). 
    The place cell network is responsible for encoding the position of the agent in the environment, 
    while the grid cell networks are responsible for encoding the agent's position in the environment in a grid-like manner. 
    The place cell network and the grid cell networks are coupled through the synaptic connections between them. 
    The model is used to decode the position of the agent in the environment from the activity of the place cell network.
'''

import brainpy as bp
import brainpy.math as bm
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

# Set random seed
bp.math.random.seed(0)

CUE = True

class P_CANN1D(bp.dyn.NeuDyn):
  def __init__(self, num, rho, tau=1., k=8.1, a=0.5,  J0=4., Jpg=0.5, z_min=-bm.pi, z_max=bm.pi,noise=0, **kwargs):
    super().__init__(size=num, **kwargs)

    self.tau = tau
    self.k = k
    self.a = a
    self.J0 = J0
    self.Jpg = Jpg
    self.noise_level=noise

    self.z_min = z_min
    self.z_max = z_max
    self.z_range = z_max - z_min
    self.x = bm.linspace(z_min, z_max, num)
    self.rho = rho
    self.dx = self.z_range / num

    self.A = 1 / (4 * bm.sqrt(bm.pi) * self.a * self.rho * self.k) * (self.rho * self.J0 + bm.sqrt(bm.square(self.rho * self.J0) - 8 * bm.sqrt(2 * bm.pi) * self.a * self.rho * self.k))

    self.u = bm.Variable(bm.zeros(num))
    self.r = bm.Variable(bm.zeros(num))
    self.e_input = bm.Variable(bm.zeros(num)) # external input
    self.g_input = bm.Variable(bm.zeros(num))
    
    self.conn_mat = self.make_conn(self.x)  # 连接矩阵

    # 定义积分函数
    self.integral = bp.odeint(self.derivative)

  def derivative(self, u, t, Ipp, Igp, Iext):
    du = (-u+Ipp+Igp+Iext)/self.tau
    return du

  def dist(self, d):
    return d

  # 计算连接矩阵
  def make_conn(self, x):
    assert bm.ndim(x) == 1
    d = self.dist(x - x[:, None])  # 距离矩阵
    Jxx = self.J0 * bm.exp(-0.5 * bm.square(d / self.a)) / (bm.sqrt(2 * bm.pi) * self.a) 
    return Jxx

  # 获取各个神经元到pos处神经元的输入
  def get_stimulus_by_pos(self, pos):
    # pos = bm.remainder(pos/self.L, 1)*self.z_range
    # gaussian_noise = bm.random.normal(self.x.shape) * self.noise_level
    #gaussian_noise = bm.random.normal(0., self.noise_level, self.x.shape)
    return self.A * bm.exp(-0.25 * bm.square(self.dist(self.x - pos) / self.a)) #+gaussian_noise

  # 网络更新函数
  def update(self, x=None):
    _t = bp.share['t']
    u2 = bm.square(self.u)
    r = u2 / (1.0 + self.k * bm.sum(u2))
    self.r = r
    Ipp = self.rho*bm.dot(self.conn_mat, r)
    self.u[:] = self.integral(self.u, _t,Ipp, self.g_input, self.e_input)
    self.g_input[:] = 0
    self.e_input[:] = 0.  # 重置外部电流

  # 寻找峰值
  def find_peaks_u(self):
    return bm.find_peaks(self.u)
  
  def find_peaks_input(self):
    return bm.find_peaks(self.input)
  
  def decode(self):
    # z = \frac{\sum_i R_p(x_i) x_i}{\sum_i R_p(x_i)},
    zp = bm.sum(self.r * self.x) / bm.sum(self.r)
    return zp


class G_CANN1D(bp.dyn.NeuDyn):
  def __init__(self, num, rho, lmd, tau=1., k=8.1, a=0.5, J0=4.,Jpg=0.5, z_min=-bm.pi, z_max=bm.pi, noise=0,**kwargs):
    super().__init__(size=num, **kwargs)

    self.tau = tau
    self.k = k
    self.a = a
    # self.A = A
    self.J0 = J0
    self.Jpg = Jpg
    self.lmd = lmd # lambda
    self.noise_level=noise

    self.z_min = z_min
    self.z_max = z_max
    self.z_range = z_max - z_min
    self.x = bm.linspace(z_min, z_max, num)
    self.rho = rho
    self.dx = self.z_range / num

    self.A = 1 / (4 * bm.sqrt(bm.pi) * self.a * self.rho * self.k) * (self.rho * self.J0 + bm.sqrt(bm.square(self.rho * self.J0) - 8 * bm.sqrt(2 * bm.pi) * self.a * self.rho * self.k))    

    self.u = bm.Variable(bm.zeros(num))
    self.r = bm.Variable(bm.zeros(num))
    self.e_input = bm.Variable(bm.zeros(num))
    self.p_input = bm.Variable(bm.zeros(num))
    self.conn_mat = self.make_conn(self.x)  # 连接矩阵

    # 定义积分函数
    self.integral = bp.odeint(self.derivative)

  # 微分方程
  def derivative(self, u, t, Igg, Ipg, Iext):
    du = (-u + Igg + Ipg + Iext) / self.tau 
    return du

  # 将距离转换到[-z_range/2, z_range/2)之间
  def dist(self, d):
    d = bm.remainder(d, self.z_range)
    d = bm.where(d > 0.5 * self.z_range, d - self.z_range, d)
    return d

  # 计算连接矩阵
  def make_conn(self, x):
    assert bm.ndim(x) == 1
    d = self.dist(x - x[:, None])  # 距离矩阵
    
    Jxx = self.J0 * bm.exp(-0.5 * bm.square(d / self.a)) / (bm.sqrt(2 * bm.pi) * self.a) 
    return Jxx

  # 获取各个神经元到pos处神经元的输入
  def get_stimulus_by_pos(self, pos):
    pos = bm.remainder(pos/self.lmd, 1) * 2 * bm.pi 
    pos = bm.where(pos > bm.pi, pos - 2 * bm.pi, pos)
    #gaussian_noise = bm.random.normal(self.x.shape) * self.noise_level
    # pos = bm.where(pos < -bm.pi, pos + 2 * bm.pi, pos)
    return self.A * bm.exp(-0.25 * bm.square(self.dist(self.x - pos) / self.a)) #+ gaussian_noise

  # 网络更新函数
  def update(self, x=None):
    _t = bp.share['t']
    u2 = bm.square(self.u)
    r = u2 / (1.0 + self.k * bm.sum(u2))
    # A_p = \frac{\rho_p J_p}{\sqrt{2}} \hat{R}_p + \frac{1}{\sqrt{2}} \sum_i \rho_g \hat{R}_g J_{g,p},
    self.r = r
    
    Igg = self.rho*bm.dot(self.conn_mat, r)
    self.u[:] = self.integral(self.u, _t,Igg, self.p_input, self.e_input)
    self.p_input[:] = 0
    self.e_input[:] = 0.  # 重置外部电流

  # 寻找峰值
  def find_peaks_u(self):
    return bm.find_peaks(self.u)
  
  def find_peaks_input(self):
    return bm.find_peaks(self.input)
  
  
params = {
  'P_CANN1D': {
    'num': 200,
    'tau': 1.,
    'k': 20,
    'a': 0.3,
    'J': 20.,
    'rho': 3.3,
    'Jpg': 0.5,
    'noise':0
  },
  'G_CANN1D': {
    'num': 20,
    'tau1': 2.09,
    'tau2': 1.57,
    'tau3': 1.26,
    'k1': 9.55,
    'k2': 12.73,
    'k3': 15.91,
    'a1': 0.63,
    'a2': 0.47,
    'a3': 0.38,
    'J': 20,
    'rho': 3.3,
    'Jpg': 0.5,
    'lambda1': 3,
    'lambda2': 4,
    'lambda3': 5,
    'noise':0
  }
}


class SpatialCANN(bp.DynamicalSystem):
  def __init__(self, params):
    super().__init__()

    assert params['P_CANN1D']['Jpg'] == params['G_CANN1D']['Jpg']
    self.Jpg = params['P_CANN1D']['Jpg']
    self.z_max = bm.pi
    self.z_min = -bm.pi
    self.z_range = self.z_max - self.z_min
    self.L=30
    self.ag1p = params['G_CANN1D']['a1']
    self.ag2p = params['G_CANN1D']['a2']
    self.ag3p = params['G_CANN1D']['a3']

    # Variables for monitoring
    self.P_u = bm.Variable(np.zeros(params['P_CANN1D']['num']))
    self.G1_u = bm.Variable(np.zeros(params['G_CANN1D']['num']))
    self.G2_u = bm.Variable(np.zeros(params['G_CANN1D']['num']))
    self.G3_u = bm.Variable(np.zeros(params['G_CANN1D']['num']))

    self.P = P_CANN1D(num=params['P_CANN1D']['num'], tau=params['P_CANN1D']['tau'], k=params['P_CANN1D']['k'], a=params['P_CANN1D']['a'], J0=params['P_CANN1D']['J'], rho=params['P_CANN1D']['rho'], Jpg=params['P_CANN1D']['Jpg'], z_max=self.L, z_min=-self.L)
    self.G1 = G_CANN1D(num=params['G_CANN1D']['num'], tau=params['G_CANN1D']['tau1'], k=params['G_CANN1D']['k1'], a=params['G_CANN1D']['a1'], J0=params['G_CANN1D']['J'], rho=params['G_CANN1D']['rho'], Jpg=params['G_CANN1D']['Jpg'], lmd=params['G_CANN1D']['lambda1'])
    self.G2 = G_CANN1D(num=params['G_CANN1D']['num'], tau=params['G_CANN1D']['tau2'], k=params['G_CANN1D']['k2'], a=params['G_CANN1D']['a2'], J0=params['G_CANN1D']['J'], rho=params['G_CANN1D']['rho'], Jpg=params['G_CANN1D']['Jpg'], lmd=params['G_CANN1D']['lambda2'])
    self.G3 = G_CANN1D(num=params['G_CANN1D']['num'], tau=params['G_CANN1D']['tau3'], k=params['G_CANN1D']['k3'], a=params['G_CANN1D']['a3'], J0=params['G_CANN1D']['J'], rho=params['G_CANN1D']['rho'], Jpg=params['G_CANN1D']['Jpg'], lmd=params['G_CANN1D']['lambda3'])
    self.P.noise_level = params['P_CANN1D']['noise']
    self.G1.noise_level = params['G_CANN1D']['noise']
    self.G2.noise_level = params['G_CANN1D']['noise']
    self.G3.noise_level = params['G_CANN1D']['noise']

  def dist(self, d):
    d = bm.remainder(d, self.z_range)
    d = bm.where(d > 0.5 * self.z_range, d - self.z_range, d)
    return d

  def make_conn_gp(self, pnet, gnet, agp):
    p = pnet.x
    g = gnet.x
    phi = bm.remainder(p / gnet.lmd, 1) * 2 * bm.pi
    phi = bm.where(phi > bm.pi, phi - 2 * bm.pi, phi)
    phi = bm.where(phi < -bm.pi, phi + 2 * bm.pi, phi)
    d = self.dist(phi - g[:, None])
    Wpg = self.Jpg * bm.exp(-0.5 * bm.square(d / agp)) / (bm.sqrt(2 * bm.pi) * agp)
    #Wpg += bm.random.normal(0., 0.9, Wpg.shape)
    return Wpg
    
  def make_conn_gp1(self,):
    return self.make_conn_gp(self.P, self.G1, self.ag1p)

  def make_conn_gp2(self,):
    return self.make_conn_gp(self.P, self.G2, self.ag2p)
  
  def make_conn_gp3(self,):
    return self.make_conn_gp(self.P, self.G3, self.ag3p)
  
  def couple_update(self):
    Wpg1 = self.make_conn_gp1()
    Wpg2 = self.make_conn_gp2()
    Wpg3 = self.make_conn_gp3()

    Ipg = bm.dot(Wpg1.T, self.G1.r)*self.G1.rho + bm.dot(Wpg2.T, self.G2.r)*self.G2.rho + bm.dot(Wpg3.T, self.G3.r)*self.G3.rho
    Igp1 = bm.dot(Wpg1, self.P.r) * self.P.rho
    Igp2 = bm.dot(Wpg2, self.P.r) * self.P.rho
    Igp3 = bm.dot(Wpg3, self.P.r) * self.P.rho

    self.P.g_input = Ipg
    # if False:
    self.G1.p_input = Igp1
    self.G2.p_input = Igp2
    self.G3.p_input = Igp3

  def update(self,pos=None):
    _t = bp.share['t']
    dt = bm.get_dt()
    Ip=self.P.get_stimulus_by_pos(pos)
    
    Ig1=self.G1.get_stimulus_by_pos(pos)
    Ig2=self.G2.get_stimulus_by_pos(pos)
    Ig3=self.G3.get_stimulus_by_pos(pos)

    if CUE:
      self.P.e_input = Ip
    self.G1.e_input = Ig1
    self.G2.e_input = Ig2
    self.G3.e_input = Ig3

    self.couple_update()

    self.P.update()
    self.G1.update()
    self.G2.update()
    self.G3.update()

    self.P_u[:] = self.P.u
    self.G1_u[:] = self.G1.u
    self.G2_u[:] = self.G2.u
    self.G3_u[:] = self.G3.u

    return self.P_u, self.G1_u, self.G2_u, self.G3_u
      
PGcann = SpatialCANN(params)  

def run(t,pos):
  bp.share.save(t=t)
  return PGcann(pos=pos)

with bm.environment(dt=0.1):
  v_ext = 0.05
  #dur1, dur2, dur3 = 10, 100, 100
  dur1, dur2, dur3 = 10, 100, 0
  num_1 = int(dur1 / bm.get_dt())
  num_2 = int(dur2 / bm.get_dt())
  num_3 = int(dur3 / bm.get_dt())
  indices = np.arange(num_1+num_2+num_3)
  pos=0
  position = np.zeros(num_1+num_2+num_3)
  position[:num_1] = pos
  for i in range(num_2):
    #pos = position[i+num_1-1] + v_ext * bm.dt
    #assert pos < PGcann.L
    position[i+num_1] = 0.05
  for i in range(num_3):
    #pos = position[i+num_1+num_2-1] - v_ext * bm.dt
    #assert pos > -PGcann.L
    position[i+num_1+num_2] = 0
  position = position.reshape(-1,1)
  
  uP, uG1, uG2, uG3 = bm.for_loop(run, (indices, position))



# Decode the position from P_CANN1D and compare with the true position
  decoded_positions = []
  true_positions = position.flatten()

  for pos in true_positions:
    PGcann.update(pos=pos)
    decoded_positions.append(PGcann.P.decode())

  decoded_positions = np.array(decoded_positions)

# Plot the true positions vs decoded positions
  plt.figure(figsize=(10, 6))
  plt.plot(true_positions, label='True Position')
  plt.plot(decoded_positions, label='Decoded Position', linestyle='dashed')
  plt.xlabel('Time Step')
  plt.ylabel('Position')
  plt.legend()
  plt.title('True Position vs Decoded Position')
  plt.show()

# Calculate the Pearson correlation coefficient
  correlation, _ = pearsonr(true_positions, decoded_positions)
  print(f'Pearson correlation coefficient: {correlation:.3f}')

