# '''
# --------------------------------------------------------------------------------
# Parametros
#
# Forcantes  - DataFrame
# VarsEstado - dicionario contento np.array para os estados S e R
#
# P1 - altura de precipitacao do passo de tempo
# E  - altura de evapotranspiracao potencial do passo de tempo
# PN - precipitacao liquida
# EN - evapotranspiracao potencial liquida
# ES - montante que sai por evapotranspiracao do reservatorio de SMA
# AE - evapotranspiracao real ('actual evapotranspiration')
# --------------------------------------------------------------------------------
# '''

import numpy as np
import pandas as pd

# def f_gera_OrdHU1(x4, D=1.25):
#     # D = 2.5 para os modelos diarios
#     # D = 1.25 para os modelos horarios - ver Tese do Ficchi (2017) pg. 51
#     n = int(np.ceil(x4))
#     SH1 = np.zeros(n+1)
#     for t in range(0, n+1):
#         if (t<=0):
#             SH1[t] = 0
#         elif (t>0) & (t<x4):
#             SH1[t] = (t/x4)**D
#         else:
#             SH1[t] = 1
#     OrdHU1 = np.diff(SH1)
#     return OrdHU1


def sacsma_detalhado(X, PME, ETP, St=None):
    print("Executando SAC-SMA em modo de simulacao detalhada...")

    UZTWM = X.get("UZTWM")
    UZFWM = X.get("UZFWM")
    LZTWM = X.get("LZTWM")
    LZFSM = X.get("LZFSM")
    LZFPM = X.get("LZFPM")
    RSERV = X.get("RSERV", 0.3) # SAVED? em mm??
    UZK   = X.get("UZK")
    LZSK  = X.get("LZSK")
    LZPK  = X.get("LZPK")
    PFREE = X.get("PFREE")
    ZPERC = X.get("ZPERC")
    REXP  = X.get("REXP")
    PCTIM = X.get("PCTIM")
    ADIMP = X.get("ADIMP")
    RIVA  = X.get("RIVA", 0.0)
    SIDE  = X.get("SIDE", 0.0)

    if St is None:
        UZTWC = UZTWM * 0.5
        UZFWC = UZFWM * 0.5
        LZTWC = LZTWM * 0.5
        LZFPC = LZFPM * 0.5
        LZFSC = LZFSM * 0.5
        ADIMC = UZTWC + LZTWC
    else:
        UZTWC = St.get("UZTWM")
        UZFWC = St.get("UZFWM")
        LZTWC = St.get("LZTWM")
        LZFPC = St.get("LZFPM")
        LZFSC = St.get("LZFSM")
        ADIMC = St.get("ADIMC")

    for PXV, EP in zip(PME, ETP):
################################################################################
# EVAPOTRANSPIRACAO
################################################################################
        # EDMND - evapotranspiracao potencial
        # E1    - evapotranspiracao real da Zona Superior (UZ) - (livre+tensao)
        # E2    - evapotranspiracao real de agua livre da Zona Superior (UZ)
        # E3    - evapotanspiracao real da Zona Inferior - (livre+tensao)
        # E4
        # E5    - evapotranspiracao da area impearmeavel variavel (ADIMP)
        # RED   - demanda residual, que passa da UZ para a LZ
        EDMND = EP
        ### ZONA SUPERIOR (UZ)
        E1 = EDMND*(UZTWC/UZTWM)
        if UTZWC >= E1: # Todo E1 eh atendido pela agua de tensao da UZ
            # E1 = E1
            UZTWC = UZWTC - E1
            RED = EDMND - E1
            E2 = 0
            # Equilibra
            if ((UZTWC/UZTWM) < (UZTWC/UZFWM)):
                UZRAT = (UZTWC + UZFWC) / (UZTWM + UZFWM)
                UZTWC = UZTWM*UZRAT
                UZFWC = UZFWM*UZRAT
        else: # Seca a agua de tensao da UZ e passa para agua livre da UZ
            E1 = UTZWC
            UZTWC = 0
            RED = EDMND - E1
            if UZFWC >= RED: # Toda residual eh atendida pela agua livre
                E2 = RED
                UZFWC = UZFWC - E2
                RED = 0
                # Equilibra
                if ((UZTWC/UZTWM) < (UZTWC/UZFWM)):
                    UZRAT = (UZTWC + UZFWC) / (UZTWM + UZFWM)
                    UZTWC = UZTWM*UZRAT
                    UZFWC = UZFWM*UZRAT
            else: # Seca a agua livre da UZ e o residual segue pra LZ
                E2 = UZFWC
                UZFWC = 0
                RED = RED - E2
        ### ZONA INFERIOR (LZ)
        E3 = RED*LZTWC/(UZTWM + LZTWM)
        if LZTWC >= E3: # Todo E3 eh atendido pela agua de tensao da LZ
            # E3 = E3
            LZTWC = LZTWC - E3
        else: # Seca a agua de tensao da LZ
            E3 = LZWTC
            LZWTC = 0
        SAVED = RSERV*(LZFPM + LZFSM) # isso nao tinha o no original? Oo
        # Transferencia de agua livre para o reservatorio de tensao da LZ
        RATLZT = LZTWC / LZTWM
        RATLZ  = (LZTWC + LZFPC + LZFSC - SAVED)/(LZTWM + LZFPM + LZFSM - SAVED)
        if (RATLZT < RATLZ):
            DEL = (RATLZ - RATLZT)*LZTWM
            LZTWC = LZTWC + DEL
            # Demonstrei que DEL<(LZFSC+LZFPC), portanto os rsvs sempre suprem
            if LZFSC >= DEL: # Tenta primeiro suprir DEL com agua do suplementar
                LZFSC = LZFSC - DEL
            else: # Se nao der, seca o suplementar e supre com agua do primario
                DEL = DEL - LZFSC
                LZFSC = 0
                LZFPC = LZFPC - DEL
        ### AREA ADIMP
        E5 = E1 + (RED + E2)*((ADIMC - E1 - UZTWC) / (UZTWM + LZTWM)) # ???
        if ADIMC >= E5:
            E5 = E5*ADIMP
            ADIMC = ADIMC - E5
        else:
            E5 = ADIMC
            ADIMC = 0

################################################################################
# PERCOLACAO E ESCOAMENTO SUPERFICIAL
################################################################################
        ### LAMINA NA ZONA SUPERIOR (TWX)
        if PXV >= (UZTWM - UZTWC) # A lamina TWX vai infiltrar
            TWX = PXV - (UZTWM - UZTWC)
            UZTWC = UZTWM
        else: # A lamina TWX fica toda retida no reservatorio UZTW
            TWX = 0
            UZTWC = UZTWC + PXV
        ### LAMINA NA AREA ADIMP
        ADIMC = ADIMC + PXV - TWX
            # Minha interpretacao da linha acima:
            # ADIMC eh uma altura pluviometrica que representa a area saturada
            # da bacia, que eh variavel, limitada e proprocional a umidade rela-
            # tiva dos reservatorios de tensao.
            # Quando o UZTW esta na capacidade maxima, ou seja, UZTWC = UZTWM,
            # presume-se que a area sataurada tb esteja em seu limite maximo.
            # Considerando que ADIMC = ADIMC + (PXV - TWX), temos:
            # - se o UZTW estava cheio, TWX = PXV e (PXV - TWX) eh minimo (=0),
            # logo nao ocorre aumento de ADIMC;
            # - se o UZTW estava vazio, TWX = 0 e (PXV-TWX) eh maximo (=PXV),
            # logo toda a precipitacao vai contribuir para o aumento de ADIMC.
        ### COMPUTA ESCOAMENTO DA AREA IMPERMEAVEL
        ROIMP = PXV * PCTIM
        # SIMPVT = SIMPVT + ROIMP (necessario?)
        # INICIALIZANDO OS ESCOAMENTOS
        SBF   = 0
        SSUR  = 0
        SIF   = 0
        SPREC = 0
        SDRO  = 0
        SPBF  = 0

        NINC  = 1 + 0.2*(UZFWC + TWX)
        DTINC = 1/NINC
        PINC  = TWX/NINC





# def f_simula_rapido(Params, PME, ETP, Qmon=None, Estados=None):
