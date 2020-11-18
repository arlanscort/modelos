'''
--------------------------------------------------------------------------------
Modelos Hidrologicos de Propagacao (Hydrologic Routing)
--------------------------------------------------------------------------------
Implementacao - Arlan Scortegagna, dez/2019
Ultima atualizacao - Arlan Scortegagna, nov/2020
Revisao - ???
--------------------------------------------------------------------------------
'''

import numpy as np
import math
from scipy import integrate, stats


def rsvs_lineares_nash_sol_analitica(Qmon, k, n, dt=1):
    '''
    ----------------------------------------------------------------------------
    Reservatorios Lineares de Nash - Solucao Analitica
    ----------------------------------------------------------------------------
    Descricao
        Realiza a propagacao por meio de uma cascata de n reservatorios lineares
        (S = kQ). Trata-se de um metodo que serve tanto para propagacao, quanto
        para transferencia de escoamento por meio de um Hidrograma Unitario (HU)
    ----------------------------------------------------------------------------
    Entradas
        Qmon - Hidrograma de vazoes de montante [array_like]
        k  - coeficiente de armazenamento; equivale ao tempo de residencia
        n  - numero de reservatorios da cascata
        dt - passo de tempo
    Saida
        Qprop - Hidrograma de vazoes propagadas [numpy.ndarray]
    ----------------------------------------------------------------------------
    Observacoes
        k e dt devem estar na mesma unidade de tempo (ex: dias ou horas)
    ----------------------------------------------------------------------------
    Referencias
        Nash (1957)
        Ven Te Chow (1988) - p.261
    ----------------------------------------------------------------------------
    '''
    # # Metodo 1 - Integracao numerica do Hidrograma Unitario Instantaneo (HUI)
    # # 1.1 - Definicao de u(t) - "Impulse Response Function" (equivale ao HUI)
    # u = lambda t: 1/(k*math.gamma(n))*(t/k)**(n-1)*np.exp(-t/k)
    # # 1.2 - Obtencao de g(t) - "Step Response Function"
    # int_u = 0
    # i = 0
    # g = []
    # while int_u <= 0.999:
    #     int_u = integrate.quad(u, 0, i*dt)[0]
    #     i += 1
    #     g.append(int_u)
    # # 1.3 - Obtencao de h(t) - "Pulse Response Function" (equivale ao HU)
    # h = np.diff(g)

    # Metodo 2 - Obtencao das ordenadas do HU por meio da CDF Gama
    # 2.1 - Definicao da pdf gama "frozen"
    u = stats.gamma(n, loc=0, scale=k)
    # 2.2 - Obtencao da CDF
    cdf = 0
    i = 0
    g = []
    while cdf <= 0.995:
        cdf = u.cdf(i*dt)
        i += 1
        g.append(cdf)
    # 2.3 - Obtencao das ordenadas h(t)
    h = np.diff(g)

    # Convolucao
    HU = np.zeros(len(h))
    Qprop = []
    for qmon in Qmon:
        HU += qmon*h
        Qprop.append(HU[0])
        HU = np.roll(HU, -1)
        HU[-1] = 0

    return np.asarray(Qprop)


def rsvs_lineares_nash_sol_diferencial(Qmon, k, n, qini=0, dt=1):
    '''
    ----------------------------------------------------------------------------
    Reservatorios Lineares de Nash - Solucao Diferencial
    ----------------------------------------------------------------------------
    Descricao
        Realiza a propagacao por meio de uma cascata de n reservatorios lineares
        (S = kQ)
    ----------------------------------------------------------------------------
    Entradas
        Qmon - Hidrograma de vazoes de montante [array_like]
        k  - coeficiente de armazenamento; equivale ao tempo de residencia
        n  - numero de reservatorios da cascata (inteiro)
        dt - passo de tempo
    Saida
        Qprop - Hidrograma de vazoes propagadas [numpy.ndarray]
    ----------------------------------------------------------------------------
    Observacoes
        k e dt devem estar na mesma unidade de tempo (ex: dias ou horas)
        Para garantir estabilidade, k >= dt/2
    ----------------------------------------------------------------------------
    Referencias
        Todini (1996) - Eq.31, p.354
    ----------------------------------------------------------------------------
    '''
    # Requisito
    n = int(n)

    # 1 - Inicializacao dos vetores
    Qprop = [qini]
    Q_ant = np.zeros(n+1)
    Q_pos = np.zeros(n+1)

    # 2 - Iterar no tempo
    for i in range(1, len(Qmon)):
        Q_ant[0] = Qmon[i-1]
        Q_pos[0] = Qmon[i]

        # 3 - Iterar nos reservatorios
        for rsv in range(1, n+1):
            Q_pos[rsv] = (2*k-dt)/(2*k+dt)*Q_ant[rsv] + dt/(2*k+dt)*(Q_pos[rsv-1] + Q_ant[rsv-1])
        Qprop.append(Q_pos[-1])
        Q_ant[1:] = Q_pos[1:]

    return np.asarray(Qprop)


def muskingum(Qmon, k, x, qini=0, dt=1):
    '''
    ----------------------------------------------------------------------------
    Metodo de Muskingum
    ----------------------------------------------------------------------------
    Descricao
        Realiza a propagacao por meio de um elemento de rio contendo um prisma,
        com volume de armazenamento S = kQ, e uma cunha, com volume de
        mento S = (I-Q)*k
    ----------------------------------------------------------------------------
    Entradas
        Qmon - Hidrograma de vazoes de montante [array_like]
        k  - coeficiente de armazenamento; equivale ao tempo de viagem
        x  - peso
        dt - passo de tempo
    Saida
        Qprop - Hidrograma de vazoes propagadas [numpy.ndarray]
    ----------------------------------------------------------------------------
    Observacoes
        k e dt devem estar na mesma unidade de tempo (ex: dias ou horas)
        x varia entre 0 (level pool) e 0.5 (cunha total), na media x = 0.2
        Para garantir estabilidade, verificar se x <= dt/(2.k) <= (1-x)
    ----------------------------------------------------------------------------
    Referencias
        Ven Te Chow - p.257
    ----------------------------------------------------------------------------
    '''
    C1 = (dt - 2*k*x)/(2*k*(1-x) + dt)
    C2 = (dt + 2*k*x)/(2*k*(1-x) + dt)
    C3 = (2*k*(1-x) - dt)/(2*k*(1-x) + dt)
    Qprop = np.zeros(len(Qmon))
    Qprop[0] = qini
    for i in range(1, len(Qmon)):
        Qprop[i] = C1*Qmon[i] + C2*Qmon[i-1] + C3*Qprop[i-1]
    return Qprop
