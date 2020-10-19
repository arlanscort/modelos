'''
Implementacao - Arlan Scortegagna, set/2020
Revisao -
'''
# Modelo proposto - serie de n reservatorio lineares com constante k
# Para cada reservatorio, S = k*Q
# Qcmon - (vazao) condicao de contorno de montante
# Qcjus - (vazao) condicao de contorno de jusante
import numpy as np
from scipy import integrate, special


def rsvs_lin_sol_analitica(Qcmon, dt, k, n):
    # 1 - Impulse Response Function - u(t)
    #   u(t) equivale ao HUI/IUH - Hidrograma Unitario Instantaneo
    u = lambda t: 1/(k*special.gamma(n))*(t/k)**(n-1)*np.exp(-t/k)
    # 2 - Step Response Function - g(t)
    g = [0.0]
    i = 1
    while g[-1] <= 0.999:
        g.append(integrate.quad(u, 0, i*dt)[0])
        i+=1
    # 3 - Pulse Response Function - h(t)
    #   h(t) corresponde as ordenadas do Hidrograma Unitario
    h = np.diff(g)
    ## Convolucao
    # Inicializacao do vetor de vazoes propagadas (saida)
    Qprop = pd.Series([], dtype='float64')
    HU = [0]*len(h)
    for t, q in Qcmon.items():
        HU = [ HU[i] + q*h[i] for i in range(len(h))]
        Qprop.loc[t] = HU[0]
        HU = HU[1:]
        HU.append(0)
    return Qprop


def rsvs_lin_sol_diferencial(Qcmon, dt, k, n):
    idx = Qcmon.index
    Qcmon = Qcmon.to_numpy()
    for rsv in range(n):
        Qprop = [0] # condicao inicial
        for i in range(1, len(Qcmon)):
            qprop = (2*k-dt)/(2*k+dt)*Qprop[i-1] + dt/(2*k+dt)*(Qcmon[i] + Qcmon[i-1])
            Qprop.append(qprop)
        Qcmon = Qprop
    Qprop = pd.Series(index=idx, data=Qprop)
    return Qprop


### Teste com exemplo do Ven te Chow (1988), pg. 262
import pandas as pd
idx = [6*i for i in range(0,11)]
Qcmon = pd.Series(index=idx, data=[0]*len(idx))
Qcmon.loc[6]  = 100
Qcmon.loc[12] = 300
Qcmon.loc[18] = 200
Qcmon.loc[24] = 100
dt = 6   # hr
k = 3.14 # hr (mesma unidade de dt)
n = 5    # hr

Qprop_analitica = rsvs_lin_sol_analitica(Qcmon, dt, k, n)
Qprop_diferencial = rsvs_lin_sol_diferencial(Qcmon, dt, k, n)

import matplotlib.pyplot as plt
Qcmon.plot(label='Montante')
Qprop_analitica.plot(label='Propagacano - solucao analitica')
Qprop_diferencial.plot(label='Propagacao - solucao diferencial')
plt.legend()
plt.show()
