'''
--------------------------------------------------------------------------------
Modelo GR4J
--------------------------------------------------------------------------------
Implementacao - Arlan Scortegagna, dez/2019
Ultima atualizacao - Arlan Scortegagna, out/2020
Revisao - ???
--------------------------------------------------------------------------------
Forcantes:
    PME - numpy.array 1D contendo a precipitacao medial espacial (mm)
    ETP - numpy.array 1D contendo a evapotranspiracao potencial (mm)
--------------------------------------------------------------------------------
Parametros:
    x1 - capacidade do reservatorio de producao (mm)
    x2 - qtde maxima de agua que pode ser trocada com o aquifero (mm)
    x3 - capacidade maxima do reservatorio de propagacao (mm)
    x4 - parametro relacionado ao tempo de base dos HUs (dias, > 0.5d)
--------------------------------------------------------------------------------
Variaveis de Estado (inseridas no dicionario "Estados"):
    S - armazenamento do reservatorio de producao (mm)
    R - armazenamento do reservatorio de propagacao (mm)
    HU1 - numpy.array 1D contendo os estados iniciais do HU 1 (mm)
    HU2 - numpy.array 1D contendo os estados iniciais do HU 1 (mm)
--------------------------------------------------------------------------------
Outros:
    area - area da bacia em km2 para conversao mm->m3/s (parametro constante)
--------------------------------------------------------------------------------
'''

# from modelos.propagacao import muskingum
import numpy as np

def ordenadas_HU1(x4, D):
    n = int(np.ceil(x4))
    SH1 = np.zeros(n+1)
    for t in range(0, n+1):
        if (t<=0):
            SH1[t] = 0
        elif (t>0) & (t<x4):
            SH1[t] = (t/x4)**D
        else:
            SH1[t] = 1
    OrdHU1 = np.diff(SH1)
    return OrdHU1, n


def ordenadas_HU2(x4, D):
    m = int(np.ceil(2*x4))
    SH2 = np.zeros(m+1)
    for t in range(0, m+1):
        if (t<=0):
            SH2[t] = 0
        elif (t>0) & (t<=x4):
            SH2[t] = (1/2)*(t/x4)**D
        elif (t>x4) & (t<2*x4):
            SH2[t] = 1 - (1/2)*(2-t/x4)**D
        else:
            SH2[t] = 1
    OrdHU2 = np.diff(SH2)
    return OrdHU2, m


def sim_muskingum(area, PME, ETP, x1, x2, x3, x4, k, x, Qmon=None, Estados=None):
    '''
    Variaveis internas
        P1 - altura de precipitacao do passo de tempo
        E  - altura de evapotranspiracao potencial do passo de tempo
        PN - precipitacao liquida
        EN - evapotranspiracao potencial liquida
        PS - montante de precipitacao que entra no reservatorio de SMA
        ES - montante que sai por evapotranspiracao do reservatorio de SMA
        PERC - montante percolado
        PR - 'precipitacao efetiva' (na verdade, considera tb o PERC)
    '''
    # Constantes (passiveis de analise e estudo de caso)
    power = 4
    split = 0.9
    D     = 2.5 # p/ modelos horarios, D = 1.25 (ver Ficchi, 2017, p. 51)
    beta  = 2.25 # p/ modelos horarios, beta = 5.25 (Ficchi, 2017, p. 266)

    # Calcula as ordenadas do HUs
    OrdHU1, n = ordenadas_HU1(x4, D)
    OrdHU2, m = ordenadas_HU2(x4, D)

    # Atribui os estados iniciais
    if Estados is None:
        Estados = {}
    S = Estados.get('S', 0.6*x1)
    R = Estados.get('R', 0.7*x3)
    HU1 = Estados.get('HU1', np.zeros(n))
    HU2 = Estados.get('HU2', np.zeros(m))

    # Executa o processo iterativo
    Q = np.array([], float)
    for P1, E in np.nditer([PME,ETP]):

        # Executa interceptacao e balanco hidrico no reservatorio de SMA/prod.
        if (P1-E) <= 0:
            EN = E - P1
            # !!! acelera calculo de ES !!!
            TWS = 1 if EN/x1 > 13 else np.tanh(EN/x1)
            ES  = S*(2 - S/x1)*TWS / (1 + (1 - S/x1)*TWS)
            # !!! acelera calculo de ES !!!
            S = S - ES
            PR = 0 # (manter pq depois vai somar com o PERC)
        else:
            # enchimento
            PN = P1 - E
            # !!! acelera o calculo de PS !!!
            TWS = 1 if PN/x1 > 13 else np.tanh(PN/x1)
            PS  = x1*(1 - (S/x1)**2)*TWS / (1 + S/x1*TWS)
            # !!! acelera o calculo de PS !!!
            S = S + PS
            PR = PN - PS

        # Percolacao
        PERC = S*(1 - (1 + (S/(beta*x1))**power)**(-1/4))
        S = S - PERC

        # 'Precipitacao efetiva'
        PR += PERC

        # Convolucao do HU1
        HU1 += OrdHU1*(PR*split)
        Q9 = HU1[0]
        HU1 = np.roll(HU1, -1)
        HU1[-1] = 0

        # Convolucao do HU2
        HU2 += OrdHU2*(PR*(1-split))
        Q1 = HU2[0]
        HU2 = np.roll(HU2, -1)
        HU2[-1] = 0

        # Montante de que PODE ser transferido para o aquifero
        EXCH = x2*(R/x3)**(7/2)

        # Escoamento do reservatorio de propagacao
        R = max(0, R+Q9+EXCH)
        QR = R*(1 - (1 + (R/x3)**power)**(-1/4))
        R = R - QR

        # Escoamento direto
        QD = max(0, Q1+EXCH)

        # Escoamento total
        Q = np.append(Q, QR + QD)

    # Consolida as vazoes em m3/s (mm -> m3/s)
    Q = Q*(area/86.4)

    # Propagacao das vazoes de montante
    if Qmon is not None:
        from modelos.propagacao import muskingum
        Qprop = muskingum(Qmon, k, x)
        Q += Qprop

    return Q


def sim_detalhada(PME, ETP, x1, x2, x3, x4, Estados=None):

#     ### TEM QUE TERMINAR, FORAM FEITAS MODIFICAÇ˜OES
#     # Calcula as ordenadas do HUs com base no parametro x4
#     OrdHU1 = f_gera_OrdHU1(x4)
#     OrdHU2 = f_gera_OrdHU2(x4)
#
#     # Atribui os estados iniciais, se nao foram passados
#     S = Estados.get('S', x1*0.3)
#     R = Estados.get('R', x3*0.5)
#     HU1 = Estados.get('HU1', np.zeros(len(OrdHU1)))
#     HU2 = Estados.get('HU2', np.zeros(len(OrdHU2)))
#
#     # Inicializa o DataFrame das saidas
#     HU1cols = ['EstHU1[{}]'.format(i) for i in range(len(OrdHU1))]
#     HU2cols = ['EstHU2[{}]'.format(i) for i in range(len(OrdHU2))]
#     cols = ['Pn','Ps','AE','Perc','PR','PExch', 'AExch1','AExch2', 'AExch', \
#             'QR','QD','Qsim','S','R'] + HU1cols + HU2cols
#     DF = pd.DataFrame(columns = cols)
#
#     # Inicio do processo iterativo
#     Forcantes = pd.concat([PME, ETP], axis=1)
#     Forcantes.columns = ['Precip','PotEvap']
#     for row in Forcantes.itertuples():
#         P1 = row[1]
#         E  = row[2]
#
#         # Interceptacao e balanco no reservatorio de SMA
#         if (P1-E) <= 0:
#             # esvaziamento
#             PN = 0
#             PS = 0
#             EN = E - P1
#             # !!! acelera calculo de PS !!!
#             TWS = 1 if EN/x1 > 13 else np.tanh(EN/x1)
#             ES  = S*(2 - S/x1)*TWS / (1 + (1 - S/x1)*TWS)
#             # !!! acelera calculo de PS !!!
#             S = S - ES
#             AE = ES + P1
#         else:
#             # enchimento
#             EN = 0
#             ES = 0
#             PN = P1 - E
#             # !!! acelera o calculo de PS !!!
#             TWS = 1 if PN/x1 > 13 else np.tanh(PN/x1)
#             PS  = x1*(1 - (S/x1)**2)*TWS / (1 + S/x1*TWS)
#             # !!! acelera o calculo de PS !!!
#             S = S + PS
#             AE = E
#
#         # Percolacao
#         PERC = S*(1 - (1 + (S/(5.25*x1))**4)**(-0.25))
#         # (B = 5.25 para os modelos horarios, ver Tese do Ficchi (2017) pg. 266)
#         S = S - PERC
#
#         # Separacao da precipitacao efetiva
#         PR = PERC + (PN - PS)
#         PRHU1 = PR*0.9
#         PRHU2 = PR*0.1
#
#         # Convolucao do Hidrograma Unitario 1
#         for i in range(len(OrdHU1)-1):
#             HU1[i] = HU1[i+1] + OrdHU1[i]*PRHU1
#         HU1[-1] = OrdHU1[-1]*PRHU1
#
#         # Convolucao do Hidrograma Unitario 2
#         for i in range(len(OrdHU2)-1):
#             HU2[i] = HU2[i+1] + OrdHU2[i]*PRHU2
#         HU2[-1] = OrdHU2[-1]*PRHU2
#
#         # Montante que PODE ser transferido para o aquifero
#         EXCH = x2*(R/x3)**(7/2)
#
#         # Escoamento do reservatorio de propagacao
#         R = max(0, R + HU1[0] + EXCH)
#         QR = R*(1 - (1 + (R/x3)**4)**(-0.25))
#         R = R - QR
#         # (apenas para calcular o AEXCH1)
#         if (R + HU1[0] + EXCH) < 0:
#             AEXCH1 = -(R + HU1[0])
#         else:
#             AEXCH1 = EXCH
#
#         # Escoamento direto
#         QD = max(0, HU2[0] + EXCH)
#         if (HU2[0]+EXCH) < 0:
#             AEXCH2 = - HU2[0]
#         else:
#             AEXCH2 = EXCH
#
#         # Escoamento total
#         Q = QR + QD
#
#         # Insere as saidas em DF
#         AEXCH = AEXCH1 + AEXCH2
#
#         data = [PN, PS, AE, PERC, PR, EXCH, AEXCH1, AEXCH2, AEXCH, QR, QD, Q, S, R]
#         DF = DF.append(pd.Series(index=DF.columns, \
#                         data=data+HU1.tolist()+HU2.tolist(), name=row[0]))
#
#     DF = pd.concat([Forcantes, DF], axis=1)
#     print('Simulacao concluida!')
    return 0
