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

def f_simula_GR4H(Parametros, Forcantes, VarsEstado, OrdHU1, OrdHU2):

    x1  = Parametros['x1']
    x2  = Parametros['x2']
    x3  = Parametros['x3']
    x4  = Parametros['x4']
    PME = Forcantes['pme']
    ETP = Forcantes['etp']
    S   = VarsEstado['S']
    R   = VarsEstado['R']
    HU1 = VarsEstado['HU1']
    HU2 = VarsEstado['HU2']

    # Inicializa o DataFrame que ira conter todas as saidas relevantes
    DF = pd.DataFrame()

    for t in PME.index:
        print(t)
        P1 = PME.loc[t]
        E = ETP.loc[t]

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

        # Convolucao dos Hidrograma Unitario 1
        for i in range(len(OrdHU1)-1):
            HU1[i] = HU1[i+1] + OrdHU1[i]*PRHU1
        HU1[-1] = OrdHU1[-1]*PRHU1

        # Convolucao do Hidrograma Unitario 2
        for i in range(len(OrdHU2)-1):
            HU2[i] = HU2[i+1] + OrdHU2[i]*PRHU2
        HU2[-1] = OrdHU2[-1]*PRHU2

        # Montante de que PODE ser transferido para o aquifero
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
        DF.loc[t,'Precip'] = P1
        DF.loc[t,'PotEvap'] = E
        DF.loc[t,'PN'] = PN
        DF.loc[t,'PS'] = PS
        DF.loc[t,'AE'] = AE
        DF.loc[t,'PERC'] = PERC
        DF.loc[t,'PR'] = PR
        DF.loc[t,'QD'] = QD
        DF.loc[t,'QR'] = QR
        DF.loc[t,'AExch'] = AEXCH1 + AEXCH2
        DF.loc[t,'Qsim'] = Q
        DF.loc[t,'Prod'] = S
        DF.loc[t,'Rout'] = R

    return DF

def f_calibra_GR4H(teste):
    return teste+1
