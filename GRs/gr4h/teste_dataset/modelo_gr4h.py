'''
--------------------------------------------------------------------------------
Parametros

Forcantes  - DataFrame
VarsEstado - dicionario contento np.array para os estados S e R

P1 - altura de precipitacao do passo de tempo
E  - altura de evapotranspiracao potencial do passo de tempo
PN - precipitacao liquida
EN - evapotranspiracao potencial liquida
ES - montante que sai por evapotranspiracao do reservatorio de SMA
AE - evapotranspiracao real ('actual evapotranspiration')
--------------------------------------------------------------------------------
'''

import numpy as np
import pandas as pd

def f_gera_OrdHU1(x4, D=1.25):
    # D = 2.5 para os modelos diarios
    # D = 1.25 para os modelos horarios - ver Tese do Ficchi (2017) pg. 51
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
    return OrdHU1


def f_gera_OrdHU2(x4, D=1.25):
    # D=2.5 para os modelos diarios
    # D=1.25 para os modelos horarios - ver Tese do Ficchi (2017) pg. 51
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
    return OrdHU2


def f_simula_detalhado(x1, x2, x3, x4, PME, ETP, S=None, R=None, HU1=None, HU2=None):
    print("Executando GR4H em modo de simulacao detalhada...")

    # Calcula as ordenadas do HUs com base no parametro x4
    OrdHU1 = f_gera_OrdHU1(x4)
    OrdHU2 = f_gera_OrdHU2(x4)

    # Atribui os estados iniciais, se nao foram passados
    if S is None:
        S = x1*0.3
    if R is None:
        R = x3*0.5
    if HU1 is None:
        HU1 = np.zeros(len(OrdHU1))
    if HU2 is None:
        HU2 = np.zeros(len(OrdHU2))

    # Inicializa o DataFrame das saidas
    HU1cols = ['EstHU1[{}]'.format(i) for i in range(len(OrdHU1))]
    HU2cols = ['EstHU2[{}]'.format(i) for i in range(len(OrdHU2))]
    cols = ['Pn','Ps','AE','Perc','PR','PExch', 'AExch1','AExch2', 'AExch', \
            'QR','QD','Qsim','S','R'] + HU1cols + HU2cols
    DF = pd.DataFrame(columns = cols)

    # Inicio do processo iterativo
    Forcantes = pd.concat([PME, ETP], axis=1)
    Forcantes.columns = ['Precip','PotEvap']
    for row in Forcantes.itertuples():
        P1 = row[1]
        E  = row[2]

        # Interceptacao e balanco no reservatorio de SMA
        if (P1-E) <= 0:
            # esvaziamento
            PN = 0
            PS = 0
            EN = E - P1
            # !!! acelera calculo de PS !!!
            TWS = 1 if EN/x1 > 13 else np.tanh(EN/x1)
            ES  = S*(2 - S/x1)*TWS / (1 + (1 - S/x1)*TWS)
            # !!! acelera calculo de PS !!!
            S = S - ES
            AE = ES + P1
        else:
            # enchimento
            EN = 0
            ES = 0
            PN = P1 - E
            # !!! acelera o calculo de PS !!!
            TWS = 1 if PN/x1 > 13 else np.tanh(PN/x1)
            PS  = x1*(1 - (S/x1)**2)*TWS / (1 + S/x1*TWS)
            # !!! acelera o calculo de PS !!!
            S = S + PS
            AE = E

        # Percolacao
        PERC = S*(1 - (1 + (S/(5.25*x1))**4)**(-0.25))
        # (B = 5.25 para os modelos horarios, ver Tese do Ficchi (2017) pg. 266)
        S = S - PERC

        # Separacao da precipitacao efetiva
        PR = PERC + (PN - PS)
        PRHU1 = PR*0.9
        PRHU2 = PR*0.1

        # Convolucao do Hidrograma Unitario 1
        for i in range(len(OrdHU1)-1):
            HU1[i] = HU1[i+1] + OrdHU1[i]*PRHU1
        HU1[-1] = OrdHU1[-1]*PRHU1

        # Convolucao do Hidrograma Unitario 2
        for i in range(len(OrdHU2)-1):
            HU2[i] = HU2[i+1] + OrdHU2[i]*PRHU2
        HU2[-1] = OrdHU2[-1]*PRHU2

        # Montante que PODE ser transferido para o aquifero
        EXCH = x2*(R/x3)**(7/2)

        # Escoamento do reservatorio de propagacao
        R = max(0, R + HU1[0] + EXCH)
        QR = R*(1 - (1 + (R/x3)**4)**(-0.25))
        R = R - QR
        # (apenas para calcular o AEXCH1)
        if (R + HU1[0] + EXCH) < 0:
            AEXCH1 = -(R + HU1[0])
        else:
            AEXCH1 = EXCH

        # Escoamento direto
        QD = max(0, HU2[0] + EXCH)
        if (HU2[0]+EXCH) < 0:
            AEXCH2 = - HU2[0]
        else:
            AEXCH2 = EXCH

        # Escoamento total
        Q = QR + QD

        # Insere as saidas em DF
        AEXCH = AEXCH1 + AEXCH2

        data = [PN, PS, AE, PERC, PR, EXCH, AEXCH1, AEXCH2, AEXCH, QR, QD, Q, S, R]
        DF = DF.append(pd.Series(index=DF.columns, \
                        data=data+HU1.tolist()+HU2.tolist(), name=row[0]))

    DF = pd.concat([Forcantes, DF], axis=1)
    print('Simulacao concluida!')
    return DF


def f_simula_rapido(x1, x2, x3, x4, PME, ETP, S=None, R=None, HU1=None, HU2=None):
    print("Executando GR4H em modo de simulacao rapida...")

    # Calcula as ordenadas do HUs com base no parametro x4
    OrdHU1 = f_gera_OrdHU1(x4)
    OrdHU2 = f_gera_OrdHU2(x4)

    # Atribui os estados iniciais, se nao foram passados
    if S is None:
        S = x1*0.3
    if R is None:
        R = x3*0.5
    if HU1 is None:
        HU1 = np.zeros(len(OrdHU1))
    if HU2 is None:
        HU2 = np.zeros(len(OrdHU2))

    Qsim = []
    Forcantes = pd.concat([PME, ETP], axis=1)
    for row in Forcantes.itertuples():
        P1 = row[1]
        E  = row[2]

        # Interceptacao e balanco no reservatorio de SMA
        if (P1-E) <= 0:
            # esvaziamento
            EN = E - P1
            # !!! acelera calculo de PS !!!
            TWS = 1 if EN/x1 > 13 else np.tanh(EN/x1)
            ES  = S*(2 - S/x1)*TWS / (1 + (1 - S/x1)*TWS)
            # !!! acelera calculo de PS !!!
            S = S - ES
            PR = 0 # (depois vai somar com PERC)
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
        PERC = S*(1 - (1 + (S/(5.25*x1))**4)**(-0.25))
        # (B = 5.25 para os modelos horarios, ver Tese do Ficchi (2017) pg. 266)
        S = S - PERC

        # Separacao da precipitacao efetiva
        PR += PERC

        # Convolucao dos Hidrogramas Unitarios
        n = len(OrdHU1)-1
        for i in range(len(OrdHU2)-1):
            if i < n:
                HU1[i] = HU1[i+1] + OrdHU1[i]*(PR*0.9)
            HU2[i] = HU2[i+1] + OrdHU2[i]*(PR*0.1)
        HU1[-1] = OrdHU1[-1]*(PR*0.9)
        HU2[-1] = OrdHU2[-1]*(PR*0.1)

        # Montante de que PODE ser transferido para o aquifero
        EXCH = x2*(R/x3)**(7/2)

        # Escoamento do reservatorio de propagacao
        R = max(0, R + HU1[0] + EXCH)
        QR = R*(1 - (1 + (R/x3)**4)**(-0.25))
        R = R - QR

        # Escoamento direto
        QD = max(0, HU2[0] + EXCH)

        # Escoamento total
        QT = QR + QD

        # Consolida serie de saida
        Qsim.append(QT)

    Qsim = pd.Series(data=Qsim, index=Forcantes.index, name='Qsim_mm')
    return Qsim


def f_calibra_NSE(X, *args):
    print('Chamada de simulacao!')
    # Argumentos
    PME     = args[0]
    ETP     = args[1]
    Qobs    = args[2]
    DFcal   = args[3]
    S0_prop = args[4]
    R0_prop = args[5]

    # Parametros
    x1 = X[0]
    x2 = X[1]
    x3 = X[2]
    x4 = X[3]

    # Estados Iniciais
    S = x1*S0_prop
    R = x3*R0_prop
    OrdHU1 = f_gera_OrdHU1(x4)
    OrdHU2 = f_gera_OrdHU2(x4)
    HU1 = np.zeros(len(OrdHU1))
    HU2 = np.zeros(len(OrdHU2))

    # Simulacao
    Qobs = Qobs.to_numpy()
    Qsim = []
    Forcantes = pd.concat([PME, ETP], axis=1)
    for row in Forcantes.itertuples():
        P1 = row[1]
        E  = row[2]

        # Interceptacao e balanco no reservatorio de SMA
        if (P1-E) <= 0:
            # esvaziamento
            EN = E - P1
            # !!! acelera calculo de PS !!!
            TWS = 1 if EN/x1 > 13 else np.tanh(EN/x1)
            ES  = S*(2 - S/x1)*TWS / (1 + (1 - S/x1)*TWS)
            # !!! acelera calculo de PS !!!
            S = S - ES
            PR = 0 # (depois vai somar com PERC)
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
        PERC = S*(1 - (1 + (S/(5.25*x1))**4)**(-0.25))
        # (B = 5.25 para os modelos horarios, ver Tese do Ficchi (2017) pg. 266)
        S = S - PERC

        # Separacao da precipitacao efetiva
        PR += PERC

        # Convolucao dos Hidrogramas Unitarios
        n = len(OrdHU1)-1
        for i in range(len(OrdHU2)-1):
            if i < n:
                HU1[i] = HU1[i+1] + OrdHU1[i]*(PR*0.9)
            HU2[i] = HU2[i+1] + OrdHU2[i]*(PR*0.1)
        HU1[-1] = OrdHU1[-1]*(PR*0.9)
        HU2[-1] = OrdHU2[-1]*(PR*0.1)

        # Montante de que PODE ser transferido para o aquifero
        EXCH = x2*(R/x3)**(7/2)

        # Escoamento do reservatorio de propagacao
        R = max(0, R + HU1[0] + EXCH)
        QR = R*(1 - (1 + (R/x3)**4)**(-0.25))
        R = R - QR

        # Escoamento direto
        QD = max(0, HU2[0] + EXCH)

        # Escoamento total
        QT = QR + QD

        # Consolida serie de saida
        Qsim.append(QT)
    Qsim = np.asarray(Qsim)

    # Nash-Sutcliffe Efficiency
    num = np.sum((Qobs[1440:]-Qsim[1440:])**2)
    den = np.sum((Qobs[1440:]-np.mean(Qobs[1440:]))**2)
    NSE = 1 - num/den
    return 1-NSE


# def f_calibra_KGE(X, *args):
