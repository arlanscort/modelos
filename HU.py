'''
--------------------------------------------------------------------------------
Hidrogramas Unitarios
--------------------------------------------------------------------------------
Implementacao - Arlan Scortegagna, out/2020
Ultima atualizacao - Arlan Scortegagna, nov/2020
Revisao - ???
--------------------------------------------------------------------------------
'''

import numpy as np
from scipy import integrate, stats


def HUI_nash(ERH, n, k, dt=1):
    '''
    ----------------------------------------------------------------------------
    Hidrograma Unitario - Reservatorios Lineares em Cascata de Nash
    ----------------------------------------------------------------------------
    Descricao
        Realiza a transferencia de vazoes por meio de uma sequencia de reserva-
        torios lineares (S = kQ). Apenas o primeiro reservatorio eh alimentado
        com precipitacao efetiva. Trata-se de um metodo que serve tanto para
        propagacao, quanto para transferencia de escoamento por meio da
        convolucao de um Hidrograma Unitario Instantaneo (HUI)
    ----------------------------------------------------------------------------
    Entradas
        ERH - Effective Rainfall Hyetograph [array_like]
        k - coeficiente de armazenamento; equivale ao tempo de residencia
        n - numero de reservatorios da cascata
    Saida
        DRH - Direct Runoff Hydrograph [numpy.ndarray]
    ----------------------------------------------------------------------------
    Observacoes
        k e dt devem estar na mesma unidade de tempo (ex: dias ou horas)
        Como a solucao adotada eh analitica, por meio de uma funcao gama, nao
        ha necessidade de n ser inteiro
    ----------------------------------------------------------------------------
    Referencias
        Ven Te Chow (1988) - p.261
        Jeng & Coon (2003)
    ----------------------------------------------------------------------------
    '''
    # 1 - Definicao do HUI (neste caso, uma pdf gama "frozen")
    u = stats.gamma(n, loc=0, scale=k)

    # 2 - Ordenadas do HU
    h = []
    i = 0
    while sum(h) <= 0.995:
        ord = u.cdf((i+1)*dt) - u.cdf(i*dt)
        h.append(ord)
        i += 1
    h = np.array(h)

    # 3 - Convolucao
    HU = np.zeros(len(h))
    DRH = []
    for p_ef in ERH:
        HU += p_ef*h
        DRH.append(HU[0])
        HU = np.roll(HU, -1)
        HU[-1] = 0

    return np.array(DRH)


def HUI_3rsv_p_dist(ERH, k, C1, C2, dt=1):
    '''
    ----------------------------------------------------------------------------
    Hidrograma Unitario - Tres Reservatorios Lineares e Precipitacao Distribuida
    ----------------------------------------------------------------------------
    Descricao
        Realiza a transferencia de vazoes por meio de 3 reservatorios alimenta-
        dos simultaneamente com precipitacao. Cada reservatorio representa uma
        fracao de area incremental da bacia, definidas por C1, C2 e C3.
        O Hidrograma Unitario Instantaneo (HUI) dessa configuracao foi formulado
        por Jen e Coon (2003) e sua primeira ordenada eh diferente de zero
    ----------------------------------------------------------------------------
    Entradas
        ERH - Effective Rainfall Hyetograph [array_like]
        k  - coeficiente de armazenamento; equivale ao tempo de residencia
        C1 - fracao de area da bacia 1 (mais proxima do reservatorio)
        C2 - fracao de area da bacia 2 (descarrega na bacia 1)
        dt - passo de tempo
    Saida
        DRH - Direct Runoff Hydrograph [numpy.ndarray]
    ----------------------------------------------------------------------------
    Observacoes
        k e dt devem estar na mesma unidade de tempo (ex: dias ou horas)
        A condicao C1 + C2 + C3 = 1 deve ser satisfeita na calibracao
    ----------------------------------------------------------------------------
    Referencias
        Jeng & Coon (2003)
    ----------------------------------------------------------------------------
    '''
    # Requisito : C1 + C2 + C3 = 1
    C3 = 1 - C1 - C2
    if C3 < 0:
        print('ERRO!')
        return 0

    # 1 - Definicao do HUI
    u = lambda t: C1*(1/k)*np.exp(-t/k) + C2*(1/k)*(t/k)*np.exp(-t/k) + \
            (1 - C1 - C2)*(1/(2*k))*(t/k)**2*np.exp(-t/k)

    # 3 - Ordenadas do HU
    h = []
    i = 0
    while sum(h) <= 0.995:
        ord = integrate.quad(u, i*dt, (i+1)*dt)[0]
        h.append(ord)
        i += 1
    h = np.array(h)

    # 4 - Convolucao
    HU = np.zeros(len(h))
    DRH = []
    for p_ef in ERH:
        HU += p_ef*h
        DRH.append(HU[0])
        HU = np.roll(HU, -1)
        HU[-1] = 0

    return np.array(DRH)
