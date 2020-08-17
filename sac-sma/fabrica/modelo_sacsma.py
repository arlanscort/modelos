'''
Implementacao: Arlan Scortegagna, agosto de 2020
Revisao:
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

    # AREA PEARMEAVEL PERMANENTE (PAREA)
    # O Sac-SMA considera 3 areas: uma sempre permeavel (PAREA); outra sempre
    # impermeavel (PCTIM); e outra permeavel variavel limitada a ADIMP.
    PAREA = 1 - ADIMP - PCTIM
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

        # EVAPOTRANSPIRACAO
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

        # PERCOLACAO (a)
        # A percolacao representa a transferencia de agua livre do UZFW p/ a LZ.
        # Primeiramente, calcula-se a altura pluviometrica superavitaria da UZ,
        # ou seja, a quantidade de precipitacao que supera o deficit de agua de
        # tensao e que eh passivel de percolar junto com a agua livre do UZFW.
        # Significados das variaveis:
        #   PXV - altura pluviometrica global
        #   TWX - altura pluviometrica superavitaria
        if PXV >= (UZTWM - UZTWC)
            TWX = PXV - (UZTWM - UZTWC)
            UZTWC = UZTWM
        else:
            TWX = 0
            UZTWC = UZTWC + PXV

        # AREA IMPERMEAVEL ADICIONAL
        # ADIMC eh uma altura pluviometrica que representa a area saturada
        # da bacia. Essa area eh variavel, limitada e proprocional a umidade
        # relativa dos UZTW, ou seja, quanto mais saturado estiver de agua de
        # tensao, maior eh essa area impermeavel adicional.
        # Quando o UZTW esta na capacidade maxima, ou seja, UZTWC = UZTWM,
        # presume-se que a area sataurada tb esteja em seu limite maximo.
        # Considerando que ADIMC = ADIMC + (PXV - TWX), temos:
        #   1 - se o UZTW ja estava cheio, UZTWC = UZTWM acima, TWX = PXV e
        #   PXV - TWX = 0 (max), logo nao ocorre aumento de ADIMC;
        #   2 - se o UZTW estava vazio, UZTWC = 0, TWX = 0 e PXV - TWX = PXV
        # (max), logo toda a precipitacao vai incrementar ADIMC.
        ADIMC = ADIMC + PXV - TWX

        # AREA IMPERMEAVEL PERMANENTE
        ROIMP = PXV * PCTIM
        SIMPVT = SIMPVT + ROIMP # ???

        # Variaveis do loop interno
        SSUR  = 0 # Somatorio do escoamento superficial
        SIF   = 0 # Somatorio do escoamento subsuperficial (interflow)
        SPERC = 0 # Somatorio da percolacao
        SDRO  = 0 # Somatorio do escoamento superficial direto (direct runoff)
        SPBF  = 0 # Somatorio do escoamento de base primario
        SSBF  = 0 # Somatorio do escoamento de base suplementar
        NINC = int(1 + 0.2*(UZFWC + TWX)) # Numero de incrementos/iteracoes
        DINC = 1/NINC                     # Duracao dos incrementos
        PINC = TWX/NINC                   # Altura incremental
        DUZ  = 1 - ((1-UZK)**DINC)  # Fracao de deplecionamento do UZFW
        DLZP = 1 - ((1-LZPK)**DINC) # Fracao de deplecionamento do LZFP
        DLZS = 1 - ((1-LSSK)**DINC) # Fracao de deplecionamento do LZFS
        # Observacoes:
        #   1 - Em DINC foi omitido o passo de tempo do modelo (DT), assim,
        #   cada incremento tem a mesma unidade de tempo do passo basico do
        #   loop externo (1hr, 6hrs, 1dia, etc.);
        #   2 - PINC maximo = 5 mm (ver anexos).
        #   3 - As fracoes de deplecionamento servem para calcular as retiradas
        #   dos reservatorios a cada incremento e equivalem a
        #   (-1)*(S[t+dt]-S[t])/S[t] = 1-exp(1-k*dt) ~= 1-(1-k)^dt (ver anexos).
        #   4 - IMPORTANTE: enquanto os parametros de recessao (UZK, LZPK e
        #   LZSK estiverem na mesma unidade de tempo do loop externo (hrs, dia,
        #   etc.), o modelo nao tera dependencia temporal, uma vez que
        #   UZK*DINC, por exemplo, sera adimensional na exponencial.

        ########################################################################
        # INICIO DO LOOP INTERNO
        ########################################################################
        for i in range(NINC):

            ADSUR = 0 # Escoamento superficial gerado na area ADIMP (0 a priori)

            # ESCOAMENTO DIRETO DA AREA ADIMP
            RATIO = (ADIMC - UZTWC) / LZTWM
            if RATIO < 0 : RATIO = 0
            ADDRO = PINC*(RATIO**2) # ???

            # PERCOLACAO (b)
            # Retira agua livre do LZFP e LZFS, liberando volume desses reserva-
            # torios para receber a percolacao.
            # Significado das variaveis:
            #   DEL_PBF - agua livre que escoa do rsv primario
            #   DEL_SBF - agua livre que escoa do rsv suplementar
            DEL_PBF = LZFPC * DLZP
            if LZFPC < DEL_PBF:
                DEL_PBF = LZFPC
            DEL_SBF = LZFSC * DLZS
            if LZFSC < DEL_SBF:
                DEL_SBF = LZFSC
            LZFPC = LZFPC - DEL_PBF
            if LZFPC < 0.0001 : LZFPC = 0
            LZFSC = LZFSC - DEL_SBF
            if LZFSC < 0.0001 : LZFSC = 0
            SPBF = SPBF + DEL_PBF
            SSBF = SSBF + DEL_SBF


            if (PINC + UZFWC) > 0.01:
                # HA agua livre disponivel na UZ; a percolacao eh efetivada

                # PERCOLACAO (c)
                # Transfere agua livre da UZ para a LZ, antes de adicionar o
                # PINC no UZFW.
                # Significado das variaveis:
                #   PERCM - Limite min de percolacao com rsvs da LZ saturados
                #   DEFR  - Deficit relativo combinado de umidade da LZ
                #   PERC  - Volume de percolacao potencial
                #   DEFA  - Deficit absoluto combinado de umidade da LZ
                # Observacao: no livro do Singh, PERCM = PBASE e DEFR = DEWET
                PERCM = LZFPM*DLZP + LZFSM*DLZS
                DEFR  = 1 - (LZTWC + LZFPC + LZFSC)/(LZTWM + LZFPM + LZFSM)
                PERC  = PERCM*(1 + ZPERC*(DEFR**REXP))*(UZFWC/UZFWM)

                # PERCOLACAO (d) - Calcula a percolacao real de acordo com as
                # disponibilidades em UZFW e nos rsvs da LZ.
                if PERC > UZFWC:
                    PERC = UZFWC
                DEFA = LZTWM + LZFPM + LZFSM - (LZTWC + LZFPC + LZFSC)
                if PERC > DEFA:
                    PERC = DEFA
                UZFWC = UZFWC - PERC
                SPERC = SPERC + PERC

                # ESCOAMENTO SUBSUPERFICIAL (INTERFLOW)
                # Calcula SIF aqui, mas soh vai utilizar a atualizacao de UZFWC
                # no final deste condicional, na etapa de infiltracao
                DEL = UZFWC * DUZ
                SIF = SIF + DEL
                UZFWC = UZFWC - DEL

                # PERCOLACAO (e)
                # Distribui o volume percolado entre os rsvs da LZ.
                # Significado das variaveis:
                #   PERCT - Volume de percolacao que vai para o LZTWM
                #   EXC   - Volume que excede a capacidade do reservatorio
                #   PERCF - Volume de percolacao que vai para o LZFP e LZFS
                #   HPL   - Razao entre a capacidade maxima do rsv primario e a
                #   capacidade maxima de armazenamento de agua livre da LZ
                #   RATLP - Razao conteudo/capacidade do rsv primario
                #   RATLS - Razao conteudo/capacidade do rsv suplementar
                #   FRACP - Fracao de agua livre indo pro rsv primario
                #   PERCP - Volume de percolacao que vai para o rsv primario
                #   PERCS - Volume de percolacao que vai para o rsv suplementar
                PERCT = PERC*(1 - PFREE)
                if (PERCT + LZTWC) > LZTWM:
                    EXC = PERCT - (LZTWM - LZTWC)
                    LZTWC = LZTWM
                else:
                    EXC = 0
                    LZTWC = LZTWC + PERCT
                PERCF = EXC + PERC*PFREE
                if PERCF > 0:
                    HPL   = LZFPM/(LZFPM + LZFSM)
                    RATLP = LZFPC/LZFPM
                    RATLS = LZFSC/LZFSM
                    FRACP = HPL*2*(1-RATLP)/((1-RATLP)+(1-RATLS))
                    if FRACP > 1 : FRACP= 1
                    PERCP = PERCF*FRACP
                    PERCS = PERCF - PERCP
                    if (LZFSC + PERCS) <= LZFSM:
                        LZFSC = LZFSC + PERCS
                    else:
                        PERCS = LZFSM - LZFSC
                        LZFSC = LZFSM
                    LZFPC = LZFPC + (PERCF - PERCS)
                    if (LZFPC > LZFPM):
                        EXC   = LZFPC - LZFPM
                        LZTWC = LZTWC + EXC
                        LZFPC = LZFPM

                # INFILTRACAO
                # Passagem de agua excedente da superficie para o solo,
                # ou seja, adiciona PINC ao UZFW
                if PINC > 0:

                    # ESCOAMENTO SUPERFICIAL
                    if (PINC + UZFWC) > UZFWM
                        SUR = PINC - (UZFWM - UZFWC)
                        UZFWC = UZFWM
                        SSUR = SSUR + SUR*PAREA
                        ADSUR = SUR*(1 - AADRO/PINC)
                        SSUR = SSUR + ADSUR*ADIMP
                    else:
                        # PINC = 0
                        UZFWC = UZFWC + PINC
            else:
                # (PINC+UZFWC) <= 0.01
                # NAO HA agua livre disponivel na UZ; nao ocorre percolacao
                UZFWC = UZFWC + PINC

            # ESCOAMENTO DIRETO DA AREA IMPERMEAVEL ADICIONAL
            ADIMC = ADIMC + PINC - ADDRO - ADSUR
            if ADIMC > UZTWM + LZTWM:
                ADDRO = ADDRO + ADIMC - (UZTWM + LZTWM)
                ADIMC = UZTWM + LZTWM
            SDOR = SDOR + ADDRO*ADIMP
            if ADIMC < 0.00001 : ADIMC = 0
        ########################################################################
        # FIM DO LOOP INTERNO
        ########################################################################


        ### AQUI PRECISA MELHORAR A DESCRICAO
        # Computa os somatorios dos escoamentos e ajusata-os as areas nas quais
        # sao gerados...
        SBF = SPBF + SSBF
        TBF = SBF*PAREA
        BFCC = TBF*(1/(1+SIDE))
        BFP = SPBF*PAREA/(1+SIDE)
        BFS = BFCC - BFP
        if BFS < 0 : BFS = 0.0
        BFNCC = TBF - BFCC

        TCI = ROIMP + SDRO + SSUR + SIF + BFCC

        EUSED = E1 + E2 + E3
        E4 = (EDMND - EUSED)*RIVA

        TCI = TCI - E4
        if (TCI < 0):
            E4 = E4 + TCI
            TCI = 0

        SROT = SROT + TCI

        if ADIMC < UZTWC : ADIMC = UZTWC
        ### ATE AQUI


    ############################################################################
    # FIM DO LOOP EXTERNO
    ############################################################################
