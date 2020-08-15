'''
Implementacao: A. Scortegagna
Data: agosto de 2020
Verificacao:
'''

import numpy as np
import pandas as pd

def sac_sma_detalhado(X, PME, ETP, Est=None):
    print("Executando SAC-SMA em modo de simulacao detalhada...")

    ## X - parametros
    # 13 obrigatorios
    UZTWM = X.get("UZTWM") #
    UZFWM = X.get("UZFWM") #
    LZTWM = X.get("LZTWM") #
    LZFSM = X.get("LZFSM") #
    LZFPM = X.get("LZFPM") #
    UZK   = X.get("UZK")   #
    LZSK  = X.get("LZSK")  #
    LZPK  = X.get("LZPK")  #
    PFREE = X.get("PFREE") #
    ZPERC = X.get("ZPERC") #
    REXP  = X.get("REXP")  #
    PCTIM = X.get("PCTIM") #
    ADIMP = X.get("ADIMP") #
    # 3 opcionais
    RIVA  = X.get("RIVA", 0.0)
    SIDE  = X.get("SIDE", 0.0)
    RSERV = X.get("RSERV", 0.3)

    # Est - variaveis de estado
    if Est is None:
        UZTWC = UZTWM * 0.5
        UZFWC = UZFWM * 0.5
        LZTWC = LZTWM * 0.5
        LZFPC = LZFPM * 0.5
        LZFSC = LZFSM * 0.5
        ADIMC = UZTWC + LZTWC
    else:
        UZTWC = Est.get("UZTWM")
        UZFWC = Est.get("UZFWM")
        LZTWC = Est.get("LZTWM")
        LZFPC = Est.get("LZFPM")
        LZFSC = Est.get("LZFSM")
        ADIMC = Est.get("ADIMC")

    ############################################################################
    # INICIO DO LOOP EXTERNO
    ############################################################################
    for PXV, EP in zip(PME, ETP):

    # Siglas:
    # UZ - Zona Superior (Upper Zone)
    # LZ - Zona Inferior (Lower Zone)
    # UTZW - Reservatorio de agua de tensao da Zona Superior
    # UZFW - Reservatorio de agua livre da Zona Superior
    # LZTW - Reservatorio de agua de tensao da Zona Inferior
    # LZPW - Reservatorio primario de agua livre da Zona Inferior
    # LZSW - Reservatorio suplementar de agua livre da Zona Inferior

        ########################################################################
        # EVAPOTRANSPIRACAO
        ########################################################################
        # Existem algumas premissas relativas a evapotranspiracao:
        #   1 - a evapotranspiracao potencial nos reservatorios de agua de
        #   tensao eh proporcional a disponibilidade relativa de umidade nos
        #   mesmos, o que nao ocorre nos reservatorios de agua livre;
        #   2 - a evapotranspiracao consome prioritariamente agua de tensao,
        #   consumindo agua livre somente de forma indireta, por meio das
        #   demandas residuais ou de transferencia de agua livre para fins de
        #   equilibrio dos armazenamentos);
        #   3 - em decorrencia da 2a premissa, se houver consumo parcial do UZTW
        #   a demanda restante vai para o LZTW, sem passar pelo UZFW.
        # Significados das variaveis:
        #   EP     - evapotranspiracao potencial global
        #   EP1    - evapotranspiracao potencial no UZTW
        #   E1     - evapotranspiracao real do UZTW
        #   EP2    - evapotranspiracao potencial no UZFW (remanescente)
        #   E2     - evapotranspiracao real do UZFW
        #   RED    - evapotranspiracao potencial residual, passa p/ a LZ
        #   UZRAT  - disponibilidade relativa de agua da UZ, para o equilibrio
        #   EP3    - evapotranspiracao potencial no LZTW
        #   E3     - evapotranspiracao real do LZTW
        #   SAVED  - agua livre reservada, nao sujeita a evapotranspiracao
        #   RATLZT - Relacao conteudo/capacidade de agua de tensao na LZ
        #   RATLZ  - Relacao conteudo/capacidade de umidade total na LZ
        EP1 = EP*(UZTWC/UZTWM)
        if UTZWC >= EPUZ:
            E1 = EP1
            E2 = 0
        else:
            E1 = UTZWC
            EP2 = EP - E1
            if UZFWC >= EP2:
                E2 = EP2
            else:
                E2 = UZFWC
        UZTWC = UTZWC - E1
        UZFWC = UZFWC - E2
        if UZTWC <= 0.00001 : UZTWC = 0
        if UZFWC <= 0.00001 : UZFWC = 0
        RED = EP - E1 - E2
        if (UZTWC/UZTWM) < (UZFWC/UZFWM):
            UZRAT = (UZTWC + UZFWC)/(UZTWM + UZFWM)
            UZTWC = UZTWM*UZRAT
            UZFWC = UZFWM*UZRAT
        EP3 = RED*(LZTWC/(UZTWM + LZTWM))
        if LZTWC >= EP3:
            E3 = EP3
        else:
            E3 = LZWTC
        LZTWC = LZTWC - E3
        if LZTWC <= 0.00001 : LZTWC = 0
        SAVED = RSERV*(LZFPM + LZFSM)
        RATLZT = LZTWC / LZTWM
        RATLZ  = (LZTWC + LZFPC + LZFSC - SAVED)/(LZTWM + LZFPM + LZFSM - SAVED)
        if (RATLZT < RATLZ):
            DEL = (RATLZ - RATLZT)*LZTWM
            LZTWC = LZTWC + DEL # DEL < LZFSC+LZFPC, sempre (ver anexos)
            if LZFSC >= DEL:
                LZFSC = LZFSC - DEL
            else:
                DEL = DEL - LZFSC
                LZFSC = 0
                LZFPC = LZFPC - DEL
        EP5 = E1 + (RED + E2)*((ADIMC - E1 - UZTWC)/(UZTWM + LZTWM)) # ???
        if ADIMC >= E5:
            ADIMC = ADIMC - E5
            E5 = E5*ADIMP
        else:
            E5 = ADIMC
            ADIMC = 0

        ################################################################################
        # PERCOLACAO
        ###############################################################################
        # A percolacao eh a passagem de agua livre da UZ para a LZ
        # Primeiramente, calcula-se a altura pluviometrica superavitaria da UZ,
        # ou seja, a quantidade de precipitacao que supera o deficit de agua de
        # tensao e que eh passivel de percolar, juntamente com a agua do UZFW.
        # Significados das variaveis:
        #   PXV - altura pluviometrica global
        #   TWX - altura pluviometrica superavitaria
        if PXV >= (UZTWM - UZTWC)
            TWX = PXV - (UZTWM - UZTWC)
            UZTWC = UZTWM
        else:
            TWX = 0
            UZTWC = UZTWC + PXV

        ########################################################################
        # AREA IMPERMEAVEL ADICIONAL
        ########################################################################
        # ADIMC eh uma altura pluviometrica que representa a area saturada
        # da bacia. Essa area eh variavel, limitada e proprocional a umidade
        # relativa dos UZTW, ou seja, quanto mais saturado estiver de agua de
        # tensao, maior eh essa area impermeavel adicional.
        # Quando o UZTW esta na capacidade maxima, ou seja, UZTWC = UZTWM,
        # presume-se que a area sataurada tb esteja em seu limite maximo.
        # Considerando que ADIMC = ADIMC + (PXV - TWX), temos:
        #   - se o UZTW ja estava cheio, UZTWC = UZTWM acima, TWX = PXV e
        #   PXV - TWX = 0 (max), logo nao ocorre aumento de ADIMC
        #   - se o UZTW estava vazio, UZTWC = 0, TWX = 0 e PXV - TWX =PXV (max),
        #   logo toda a precipitacao vai contribuir para o aumento de ADIMC.
        ADIMC = ADIMC + PXV - TWX

        ########################################################################
        # ESCOAMENTO GERADO NA AREA IMPERMEAVEL
        ########################################################################
        ROIMP = PXV * PCTIM
        SIMPVT = SIMPVT + ROIMP # ???

        # Acumuladoes do loop interno
        SSUR  = 0 # Escoamento superficial
        SIF   = 0 # Escoamento subsuperficial (interflow)
        SPERC = 0 # Percolacao
        SDRO  = 0 # Escoamento superficial direto (direct runoff)
        SBF   = 0 # Escoamento de base (baseflow - primario + suplementar)
        SPBF  = 0 # Escoamento de base primario

        # Parametros do loop interno
        NINC = int(1 + 0.2*(UZFWC + TWX))
        DINC = 1/NINC
        PINC = TWX/NINC
        DUZ  = 1 - ((1-UZK)**DINC)
        DLZP = 1 - ((1-LZPK)**DINC)
        DLZS = 1 - ((1-LSSK)**DINC)
        # Observacao - o DINC nao considera o passo de tempo, uma vez que presu-
        # m

        # NINC - Numero de incrementos/passos de tempo do loop interno
        # DINC - Intervalo de tempo de cada passo
        # PINC - Umidade incremental de cada passo (max 5mm, ver pasta anexos)
        # DUZ  - Fracao de deplecionamento do reservatorio UZFW
        # DLZP - Fracao de deplecionamento do reservatorio da LZ primario
        # DLZS - Fracao de deplecionamento do reservatorio da LZ suplementar
        ### Numero de incrementos computados no loop interno

            # Frac = (-1)*(S[t+dt]-S[t])/S[t] = 1-exp(1-k.dt) ~= 1-(1-k)^dt
            # Servem para calcular as perdas dos reservatorios...
            # (ver anexos para entender a logica!)
            # Observacao: enquanto as taxas UZK, LZPK e LZSK e ZPERC estiverem com a mesma
            # unidade do passo de tempo basico do modelo (hrˆ-1, diaˆ-1 ou ateh
            # 6hrˆ-1), o modelo nao tera dependencia temporal, pois UZK*DINC,
            # por exemplo, sera adimensional.


        ########################################################################
        # INICIO DO LOOP INTERNO
        ########################################################################
        for i in range(NINC):

            ADSUR = 0 # ??? VAI USAR LA EMBAIXO

            # Escoamento direto gerado na area ADIMP
            RATIO = (ADIMC - UZTWC) / LZTWM
            if RATIO < 0 : RATIO = 0
            ADDRO = PINC*(RATIO**2)
                # Soh quem fez o modelo pra entender as 3 linhas acima...

            # 1o. Retira a agua livre que escoa dos reservatorios inferiores
            DEL = LZFPC * DLZP
            if LZFPC > DEL:
                # DEL = DEL
                LZFPC = LZFPC - DEL
            else:
                DEL = LZFPC
                LZFPC = 0
            SBF  = SBF + DEL
            SPBF = SPBF + DEL

            DEL = LZFSC * DLZS
            if LZFSC > DEL:
                # DEL = DEL
                LZFSC = LZFSC - DEL
            else:
                DEL = LZFSC
                LZFSC = 0
            SBF  = SBF + DEL
            SSBF = SSBF + DEL

            # 2o. Percolacao de agua livre (UZFW + PINC) para os rsvs inferiores
            if (PINC + UZFWC) > 0.01: # Tem agua disponivel para percolar
                # DEFR  - Deficit relativo de umidade da zona inferior
                # PERCM - Limite minimo de percolacao na condicao em que os rsvs
                # DEFR  - Deficit relativo de umidade da zona inferior (DEWET)
                # PERCT - Montante percolado que vai para o rsv de agua de tensao
                # PERCF - Montante de percolacao que vai para os rsvs de agua livre
                PERCM   = LZFPM*DLZP + LZFSM*DLZS
                SLZ_DEF = LZTWM + LZFPM + LZFSM - LZTWC - LZFPC -LZFSC
                DEFR    = SLZ_DEF/(LZTWM + LZFPM + LZFSM)
                PERC    = PERCM*(1 + ZPERC*(DEFR**REXP))*(UZFWC/UZFWM)

                # Primeiro retira-se a agua, depois adiciona-se TWX ??
                if PERC > UZFWC:
                    PERC = UZFWC
                if PERC > SLZ_DEF:
                    PERC = SLZ_DEF
                UZFWC = UZFWC - PERC
                SPERC = SPERC - PERC

                # Escoamento subsuperficial (interflow)
                DEL = UZFWC * DUZ
                UZFSC = UZFSC - DEL
                SIF = SIF + DEL

                # Distribuicao da agua percolada na zona inferior

                # Primeiro de agua de tensao (PERCT)
                PERCT = PERC*(1 - PFREE)
                if (PERCT + LZTWC) > LZTWM:
                    # Agua em excesso eh adicionada ao PERCF
                    PERCT_EXC = PERCT + LZTWC - LZTWM
                    LZTWC = LZTWM
                else:
                    # Nao vai para agua livre
                    PERCT_EXC = 0
                    LZTWC = LZTWC + PERCT


                # Depois agua livre (PERCF)
                PERCF = PERCT_EXC + PERC*PFREE
                # Depois de agua livre
                if PERCF > 0:
                    HPL   = LZFPM / (LZFPM + LZFSM)
                    RATLP = LZFPC/LZFPM
                    RATLS = LZFSC/LZFSM
                    FRACP = min(1, HPL*2*(1-RATLP)/((1-RATLP)+(1-RATLS)))

                    PERCP = PERCF*FRACP

                    # Primeiro coloca no suplementar (tentativo)
                    PERCS = PERCF - PERCP
                    if (LZFSC + PERCS) <= LZFSM:
                        # PERCS = PERCS
                        LZFSC = LZFSC + PERCS
                    else:
                        PERCS = LZFSM - LZFSC
                        LZFSC = LZFSM

                    LZFPC = LZFPC + (PERCF - PERCS)
                    if (LZFPC > LZFPM):
                        EXCESS = LZFPC - LZFPM
                        LZTWC = LZTWC + EXCESS
                        LZFPC = LZFPM

                # Adiciona PINC no UZFWC E NO RSV SUPERFICIAL
                if PINC > 0:

                    if (PINC + UZFWC) > UZFWM:
                        # Parte que nao eh do escoamento direto
                        SUR = (PINC + UZFWC) - UZFWM
                        UZFWC = UZFWM
                        SSUR = SSUR + SUR*PAREA
                        # Do escoament direto
                        ADSUR = SUR*(1 - ADDRO/PINC)
                        SSUR = SSUR + ADSUR*ADIMP

                    else:
                        UZFWC = UZFWC + PINC

            ADIMC = ADIMC + PINC - ADDRO - ADSUR
            if ADIMC > UZTWM + LZTWM:
                ADDRO = ADDRO + ADIMC - (UZTWM + LZTWM)
                ADIMC = UZTWM + LZTWM
            SDOR = SDOR + ADDRO*ADIMP
            if ADIMC < 0.00001 : ADIMC = 0
        ########################################################################
        # FIM DO LOOP INTERNO
        ########################################################################

    ############################################################################
    # FIM DO LOOP EXTERNO
    ############################################################################

    TBF = SBR*PAREA
    BFCC = TBF*(1/(1+SIDE))
    BFP = SPBF*PAREA/(1+SIDE)
    BFS = BFCC - BFP
    BFNCC = TBF - BFCC

    TCI = ROIMP + SDRO + SSUR + SIF + BFCC

    EUSED = E1 + E2 + E3
    E4 = (EDMND - EUSED)*RIVA

    TCI = TCI - E4
    if (TCI < 0):
        E4 = E4 + TCI
        TCI = 0

    SROT = SROT + TCI
