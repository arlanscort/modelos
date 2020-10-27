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

def nash(Qmon, k, n, huini=0):
    # 1 - Obtencao da Impulse Response Function - u(t)
    # (equivale ao HUI/IUH - Hidrograma Unitario Instantaneo)
    u = lambda t: 1/(k*special.gamma(n))*(t/k)**(n-1)*np.exp(-t/k)

    # 2 - Obtencao da Step Response Function - g(t)
    g = [0.0]
    i = 1
    while g[-1] <= 0.999:
        g.append(integrate.quad(u, 0, i)[0])
        i+=1

    # 3 - Obtencao da Pulse Response Function - h(t)
    # (corresponde as ordenadas do Hidrograma Unitario)
    h = np.diff(g)

    # 4 - Inicializacao dos vetores
    HU = np.full(len(h), hu_ini, float)
    Qprop = np.array([], float)

    # 5 - Convolucao
    for qmon in np.nditer([Qmon]):
        HU += qmon*h
        Qprop = np.append(Qprop, HU[0])
        HU = np.roll(HU, -1)
        HU[-1] = 0

    # Retorna as vazoes propagadas
    return Qprop


def muskingum(Qmon, k, x, qini=0):
    # dt=1 se a dimensao de k for a mesma do passo de tempo
    # k eh a velocidade da onda, ou tempo de residencia (S = k*Q)
    # ex1: se o passo de tempo for de 1 hora e k estiver em horas, dt = 1
    # ex2: se o passo de tempo for de 6 horas e k estiver em horas, dt = 6

    C1 = (1 - 2*k*x)/(2*k*(1-x) + 1)
    C2 = (1 + 2*k*x)/(2*k*(1-x) + 1)
    C3 = (2*k*(1-x) - 1)/(2*k*(1-x) + 1)

    Qprop = [qini]
    for i,_ in enumerate(Qmon[1:]):
        qprop = C1*Qmon[i] + C2*Qmon[i-1] + C3*Qmon[i-1]
        Qprop.append(qprop)

    Qprop = np.array(Qprop)
    return Qprop
