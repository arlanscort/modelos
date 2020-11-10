'''
--------------------------------------------------------------------------------
Modelos Hidrologicos de Propagacao (Hydrologic Routing)
--------------------------------------------------------------------------------
Implementacao - Arlan Scortegagna, dez/2019
Ultima atualizacao - Arlan Scortegagna, out/2020
Revisao - ???
--------------------------------------------------------------------------------
'''

import numpy as np
from scipy import special, integrate


def nash_sol_analitica(Qmon, k, n, dt=1):
    # k deve estar com a mesma dimensao do passo de tempo
    # 1 - Obtencao da Impulse Response Function - u(t)
    # (equivale ao HUI/IUH - Hidrograma Unitario Instantaneo)
    u = lambda t: 1/(k*special.gamma(n))*(t/k)**(n-1)*np.exp(-t/k)
    # 2 - Obtencao da Step Response Function - g(t)
    g = [0.0]
    i = 1
    while g[-1] <= 0.999:
        g.append(integrate.quad(u, 0, i*dt)[0])
        i+=1
    # 3 - Obtencao da Pulse Response Function - h(t)
    # (corresponde as ordenadas do Hidrograma Unitario)
    h = np.diff(g)
    # 4 - Inicializacao dos vetores
    HU    = np.zeros(len(h))
    Qprop = []
    # 5 - Convolucao
    for qmon in Qmon:
        HU += qmon*h
        Qprop.append(HU[0])
        HU = np.roll(HU, -1)
        HU[-1] = 0
    # Retorna as vazoes propagadas
    return np.asarray(Qprop)


def nash_sol_diferencial(Qmon, k, n, qini=0, dt=1):
    # k deve estar com a mesma dimensao do passo de tempo
    # 1 - Inicializacao dos vetores
    Qprop = [qini]
    Q_ant = np.zeros(n+1)
    Q_pos = np.zeros(n+1)
    # 2 - Iterar na dimensao tempo
    for i in range(1, len(Qmon)):
        Q_ant[0] = Qmon[i-1]
        Q_pos[0] = Qmon[i]
        # 3 - Iterar nos reservatorios
        for rsv in range(1, n+1):
            Q_pos[rsv] = (2*k-dt)/(2*k+dt)*Q_ant[rsv] + dt/(2*k+dt)*(Q_pos[rsv-1] + Q_ant[rsv-1])
        Qprop.append(Q_pos[-1])
        Q_ant[1:] = Q_pos[1:]
    # Retorna as vazoes propagadas
    return np.asarray(Qprop)


def muskingum(Qmon, k, x, qini=0, dt=1):
    # dt=1 se a dimensao de k for a mesma do passo de tempo
    # k eh a velocidade da onda, ou tempo de residencia (S = k*Q)
    # ex1: se o passo de tempo for de 1 hora e k estiver em horas, dt = 1
    # ex2: se o passo de tempo for de 6 horas e k estiver em horas, dt = 6
    C1 = (dt - 2*k*x)/(2*k*(1-x) + dt)
    C2 = (dt + 2*k*x)/(2*k*(1-x) + dt)
    C3 = (2*k*(1-x) - dt)/(2*k*(1-x) + dt)
    Qprop = np.zeros(len(Qmon))
    Qprop[0] = qini
    for i in range(1, len(Qmon)):
        Qprop[i] = C1*Qmon[i] + C2*Qmon[i-1] + C3*Qprop[i-1]

    return Qprop


# ### Teste nash
# # Exercicio de propagacao hidrologica
# # Ven te Chow, pg. 262
# t = [0, 3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 60]
# chuva = [0, 100, 300, 200, 100, 0, 0, 0, 0, 0, 0, 0]
# vazao = [0, 0, 10, 70, 165, 180, 142, 79, 38, 13, 3, 0]
# dt = 6   # horas
# k = 3.12 # horas
# n = 5.31
# Qprop_analitica = nash_sol_analitica(chuva, k, n, dt=6)
# Qprop_diferencial = nash_sol_diferencial(chuva, k, int(n), dt=6)
# import matplotlib.pyplot as plt
# plt.plot(t, chuva, label='Chuva (m3/s)', color='black')
# plt.plot(t, vazao, label='Vazão observada (m3/s)', color='black')
# plt.plot(t, Qprop_analitica, label='Vazão simulada - Nash solução analítica (m3/s)')
# plt.plot(t, Qprop_diferencial, label='Vazão simulada - Nash solução diferencial (m3/s)')
# plt.legend()
# plt.show()


# ### Teste muskingum
# # Exercicio de propagacao hidrologica
# # Ven te Chow, pg. 262
# t = [i for i in range(1,21)]
# Qmon = [93,137,208,320,442,546,630,678,691,675,634,571,477,390,329,247,184,134,108,90]
# Qresp = [85,91,114,159,233,324,420,509,578,623,642,635,603,546,479,413,341,274,215,170]
# k = 2.3
# x = 0.15
# Qmeu = muskingum(Qmon, k, x, qini=85, dt=1)
# import matplotlib.pyplot as plt
# plt.plot(t, Qmon, label='Inflow (m3/s)', color='blue')
# plt.plot(t, Qresp, label='Outflow - exercicio (m3/s)', color='black')
# plt.plot(t, Qmeu, label='Outflow - obtido (m3/s)', color='red', linestyle=':')
# plt.legend()
# plt.show()


# ### Testes Calibracao
# # Teste para verificar se, ajustando os parametros, o metodo de Nash com a solu-
# # cao analitica pode ser melhor que a solucao diferencial
# from spotpy.analyser import *
# from spotpy.algorithms import sceua
# from spotpy.parameter import Uniform
# # Dados
# t = [0, 3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 60]
# Qmon = [0, 100, 300, 200, 100, 0, 0, 0, 0, 0, 0, 0]
# Qjus = [0, 0, 10, 70, 165, 180, 142, 79, 38, 13, 3, 0]
# # Calibracao 1 - Nash Analitico
# class setup_nash_analitico(object):
#     k = Uniform(low = 0.01, high = 10)
#     n = Uniform(low = 0.01, high = 10)
#     def __init__(self, Qmon, Qjus, dt=1):
#         self.Qmon = Qmon
#         self.Qjus = Qjus
#         self.dt = dt
#     def simulation(self, x):
#         Qprop = nash_sol_analitica(self.Qmon, x[0], x[1], self.dt)
#         return Qprop
#     def evaluation(self):
#         Qobs = self.Qjus
#         return Qobs
#     def objectivefunction(self, simulation, evaluation):
#         rmse = np.sqrt(np.mean((simulation-evaluation)**2))
#         fmin = rmse
#         return fmin
# spot_setup = setup_nash_analitico(Qmon, Qjus, dt=6)
# sampler = sceua(spot_setup)
# sampler.sample(5000, ngs=5, kstop=3, peps=0.1, pcento=0.1)
# results = sampler.getdata()
# params1 = get_best_parameterset(results,maximize=False)
# Qsim1 = nash_sol_analitica(Qmon, params1[0][0], params1[0][1], dt=6)
# # Calibracao 2 - Nash Diferencial
# class setup_nash_diferencial(object):
#     k = Uniform(low =  3, high = 10)
#     n = Uniform(low =  5, high = 10)
#     def __init__(self, Qmon, Qjus, dt=1):
#         self.Qmon = Qmon
#         self.Qjus = Qjus
#         self.dt = dt
#     def simulation(self, x):
#         Qprop = nash_sol_diferencial(self.Qmon, x[0], int(x[1]), self.dt)
#         return Qprop
#     def evaluation(self):
#         Qobs = self.Qjus
#         return Qobs
#     def objectivefunction(self, simulation, evaluation):
#         rmse = np.sqrt(np.mean((simulation-evaluation)**2))
#         fmin = rmse
#         return fmin
# spot_setup = setup_nash_diferencial(Qmon, Qjus, dt=6)
# sampler = sceua(spot_setup)
# sampler.sample(5000, ngs=5, kstop=3, peps=0.1, pcento=0.1)
# results = sampler.getdata()
# params2 = get_best_parameterset(results,maximize=False)
# Qsim2 = nash_sol_diferencial(Qmon, params2[0][0], int(params2[0][1]), dt=6)
# # Plotar
# import matplotlib.pyplot as plt
# plt.plot(t, Qjus, label='Qobs', color='black')
# plt.plot(t, Qsim1, label='Qsim1 - Nash Analitico', color='red')
# plt.plot(t, Qsim2, label='Qsim2 - Nash Diferencial', color='blue')
# plt.legend()
# plt.show()
# # Resultados calibracao - muito mais facil encontrar o minimo global com a fun-
# # cao analitica!
