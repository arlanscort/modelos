from numpy import array
from math import log, exp
from random import random
import numpy as np
import pandas as pd

def SAC_SMA(etp, cmb, params, qin=None, states_ini=None, get_states=None):

    """ Modelo chuva-vazão Sacramento Soil Moisture Accouting (SAC-SMA)
        Modelo de propagação em canal (escoamento superficial) com três reservatórios lineares em cascata

        ### Comentario de Mino Sorribas, agosto de 2019
            ASPECTOS REVISADOS
            - opcao de retornar serie de groundwater, argumento 'flowsep'
            - controle de DRY na evaporacao da zona superior
            - 'else' errado no calculo de  SDRO
            - readaptacao de SIDE para formulacao original
            - fator de distribucao de interflow 'fSIF' entre SURF e GRND, (default=1., 100% para surf)
            POSSIVEIS MELHORIAS
            - verificar balanco hidrico
            - simplificar propagacao de reservatorio linear, com controle de fluxo
            - melhorar padrao de argumentos f(cmb,etp,params,states=None,qin=None,**kwargs)
            - remover q0 dos argumentos
            - melhorar __docs__ com padrao PEP ou similar
            - inicializar parametros com variaveis locais
            - utilizar .get para inicializar parametros opcionais

        ### Comentario de Arlan ALmeida, agosto de 2019
            MODIFICACOES
            - nos comentarios, melhorei (no meu entendimento) o esclarecimento de algumas linhas de codigo, deixando as definicoes mais precisas em relacao a literatura
            - tirei cedilhas, acentos, etc. dos comentarios...
            - simplifiquei o algoritmo na parte das perdas por evapotranspiracao, sem modificar a mecanica e o resultado final
            - o codigo esta de acordo com as seguintes versoes de fland1:

            REFERENCIAS
            FORTRAN - Dan Bronman USBR
                https://github.com/danbroman/NWS_SacSMA_source/blob/master/sacsma/src/model_execution/fland1.f
            Em C - Projeto Hydromad
                https://github.com/floybix/hydromad/blob/master/src/fland1.c
            Em R
                https://github.com/tanerumit/sacsmaR/blob/master/R/sacsma.R

                A SEQUENCIA EH
                EVAPOTRANSPIRACAO
                PERCOLACAO
                INFILTRACAO
                ESCOAMENTO
                Ocorre simultanemanete...


def SAC_SMA(etp, cmb, params, qin=None, states_ini=None, get_states=None):

    ENTRADAS
        etp - serie de dados de evapotranspiração potencial[mm]
        cmb - serie de dados de chuva media na bacia [mm]
        qin - serie de dados de vazao afluente a montante da bacia [m3/s]
        params - dicionario de parametros para os quais os limites, as dimensoes e os significados (segundo Vrugt et al., 2006) sao:
            UZTWM [1.0 - 150.0]    (mm)
            UZFWM [1.0 - 150.0]    (mm)
            LZTWM [1.0 - 500.0]    (mm)
            LZFPM [1.0 - 1000.0]   (mm)
            LZFSM [1.0 - 1000.0]   (mm)
            ADIMP [0.0 - 0.40]     (fracao decimal)
            UZK   [0.1 - 0.5 ]     (1/d)
            LZPK  [0.0001 - 0.001] (1/d)
            LZSK  [0.01 - 0.25]    (1/d)
            ZPERC [1.0 - 250.0]    (adimensional)
            REXP  [1.0 - 5.0]      (adimensional)
            PCTIM [0.0 - 0.1]      (fracao decimal)
            PFREE [0.0 - 0.6]      (fracao decimal)
            RIVA  [0.0]            (fracao decimal)
            SIDE  [0.0]            (adimensional)
            RSERV [0.3]            (fracao decimal)
            NCK   [0.01 - ]= tempo de residencia em dias; constante dos reservatorios lineares de Nash em cascata (NC)
            area  = area projetada da bacia em km2 (recomendo utilizar a projecao conica de Albers)
            dt    = passo de tempo do loop principal em unidades de dias; corresponde a resolucao das series de entrada


    area  = área [km2] da sub-bacia simulada (incremental);
    dt    = passo de tempo em dias (e.g. dt = 1.0 para passo diario; dt = 1.0/24.0 para passo horario)
    q0    = (opcional), vazão inicial para reservatórios lineares
    qub   = (opcional, lista com os dados de vazão afluente [m3/s]
    state = lista de condição inicial  [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,Q1,Q2,Q3] #apesar que só opera no solo

    get_states = None, "all" ou "last", controla a saída de variáveis de estado [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,Q1,Q2,Q3]

    Saída:

    output = lista com os dados de vazão calculado pelo modelo [m3/s] (get_states=None)
    output,states = lista com os dados de vazão calculado pelo modelo [m3/s] e variáveis de estado (get_states = "all" ou "last")


    fargs = {
    "feval": sacsma3s,
    "xl":    [  10.0,      5.0,     1.0,     1.0,    1.0,      0.0,   0.1, 0.0001,  0.01,      1.0,   0.0,     0.0,     0.0,  0.,   0.]
    "xu":    [ 150.0,    150.0,   500.0,  1000.0, 1000.0,     0.40,   0.5,  0.025,  0.25,    250.0,   5.0,     0.1,     0.1, 10.,   1.]
    "xname": ['UZTWM', 'UZFWM', 'LZTWM', 'LZFSM', 'LZFPM', 'ADIMP', 'UZK', 'LZPK', 'LZSK', 'ZPERC', 'REXP', 'PCTIM','PFREE', 'k','fSIF']
    "descript":"fland"
    }

    """


    ### Aquisicao de parametros

    # SAC-SMA
    UZTWM = params.get("UZTWM")
    UZFWM = params.get("UZFWM")
    LZTWM = params.get("LZTWM")
    LZFSM = params.get("LZFSM")
    LZFPM = params.get("LZFPM")
    RSERV = params.get("RSERV", 0.3)
    UZK   = params.get("UZK")
    LZSK  = params.get("LZSK")
    LZPK  = params.get("LZPK")
    PFREE = params.get("PFREE")
    ZPERC = params.get("ZPERC")
    REXP  = params.get("REXP")
    PCTIM = params.get("PCTIM")
    ADIMP = params.get("ADIMP")
    RIVA  = params.get("RIVA", 0.0)
    SIDE  = params.get("SIDE", 0.0)
    # Nash cascade
    NCK   = params.get("NCK", 1.0)
    # Dados
    xETP  = params.get("xETP", 1.0)
    xCMB  = params.get("xCMB", 1.0)
    # Modelo
    dt   = params.get("dt")
    area = params.get("area")

    # Inicializando armazenamentos da fase bacia
    if states_ini is None:
        UZTWC = UZTWM * 0.5
        UZFWC = UZFWM * 0.5
        LZTWC = LZTWM * 0.5
        LZFPC = LZFPM * 0.5
        LZFSC = LZFSM * 0.5
        ADIMC = UZTWC + LZTWC
    else: # Ajeitar aqui tb...
        UZTWC = states_ini[0]
        UZFWC = states_ini[1]
        LZTWC = states_ini[2]
        LZFPC = states_ini[3]
        LZFSC = states_ini[4]
        ADIMC = states_ini[5]

    estado_UZFWC = []

    #--------------------------------------------------------------------------
    # MODELO SACRAMENTO SOIL MOISTURE ACCOUNTING (SAC-SMA)
    #--------------------------------------------------------------------------

    # #Lista onde os dados de vazão gerada pela fase bacia serão armazenados
    # qbac = []
    #
    # qZ1 = []
    # qZ2 = []
    # qZ3 = []
    # qTCI = []




    POTIM = PCTIM + ADIMP # Area potencialmente impermeavel (R.J.C. Burnash, 1995)
    # Listas para os 5 escoamentos que compoem o hidrograma
    direto         = []
    superficial    = []
    subsuperficial = []
    suplementar    = []
    primario       = []
    #Numero de passos de tempo
    lenCMB = len(cmb)

    #===========================================================================================================================
    #FASE BACIA
    #===========================================================================================================================
    # Iterando o modelo a cada registro da serie de dados - intervalo de tempo DT
    for i in range(lenCMB):


        #===# PERDAS POR EVAPOTRANSPIRACAO #===#
        EDMND = etp[i] * xETP # Demanda por evapotranspiracao (potencial) do intervalo DT ajustada pelo multiplicador opcional

        # ZONA SUPERIOR
        # E1  - Evapotranspiracao efetiva no reservatorio UZTW
        # E2  - Evapotranspiracao efetiva no reservatorio UZFW
        E1 = EDMND * (UZTWC / UZTWM)
        if E1 < UZTWC:
            E1 = E1 # Coloquei aqui (e nos proximos casos similares) para enfatizar que etr = etp
            UZTWC = UZTWC - E1
            E2 = 0.0
        else:
            E1 = UZTWC # etr1 = toda agua disponivel
            UZTWC = 0.0
            E2 = EDMND - E1
            if E2 > UZFWC:
                E2 = UZFWC # etr2 = toda agua disponivel
                UZFWC = 0.0
            else:
                E2 = E2
                UZFWC = UZFWC - E2
        if (UZFWC/UZFWM) > (UZTWC/UZTWM):
            # Equilibrar o armazenamento relativo da zona superior (UZRAT), transferindo agua livre para o rsv de tensao
            UZRAT = (UZTWC+UZFWC) / (UZTWM+UZFWM)
            UZFWC = UZFWM * UZRAT
            UZTWC = UZTWM * UZRAT

        # ZONA INFERIOR
        # E3     - Evapotranspiracao efetiva no reservatorio LZTW
        # RATLZT - Armazenamento relativo do LZTW
        # RATLZ  - Armazenamento relativo dos reservatorios da zona inferior (LZTW+LZFP+LZFS) considerando SAVED
        E3 = (EDMND - E1 - E2) * LZTWC / (UZTWM + LZTWM)
        if E3 > LZTWC:
            E3 = LZTWC # etr3 = toda agua disponivel
            LZTWC = 0.0
        else:
            E3 = E3
            LZTWC = LZTWC - E3
        # Verificar os armazenamentos relativos...
        SAVED  = RSERV * (LZFPM + LZFSM)
        RATLZT = LZTWC / LZTWM
        RATLZ  = (LZTWC + (LZFPC+LZFSC-SAVED)) / (LZTWM + (LZFPM+LZFSM-SAVED))

        if (RATLZT < RATLZ):
            # Equilibrar os armazenamentos relativos da zona inferior, transferindo agua livre para o rsv de tensao
            DEL = (RATLZ - RATLZT) * LZTWM # DEL significa transferencia de agua livre (ocorre em outras partes do codigo tb)
            LZTWC = LZTWC + DEL
            # Primeiro, tentar suprir a demanda com DEL do reservatorio suplementar
            if DEL <= LZFSC:
                LZFSC = LZFSC - DEL
            else:
                # Suprir a demanda com DEL do reservatorio primario
                DEL = DEL - LZFSC
                LZFSC = 0.0
                LZFPC = LZFPC - DEL

        # AREA IMPERMEAVEL ADICIONAL
        # E5 - Evapotranspiracao efetiva na "Area Adicional Impermeavel" (ADIM)
        # O limite superior de ADIM eh representado por ADIMP, que eh uma fracao decimal que deve ser multiplicada pela area total da bacia
        E5 = E1 + (EDMND - E1) * ((ADIMC - E1 - UZTWC) / (UZTWM + LZTWM))
        if E5 > ADIMC:
            E5 = ADIMC
            ADIMC = 0.0
        else:
            E5 = E5
            ADIMC = ADIMC - E5
        E5 = E5 * ADIMP


        #===# LAMINA EXCEDENTE (TWX) E ESCOAMENTO IMPERMEAVEL #===#
        PXV = cmb[i] * xCMB # Altura de precipitacao acumulada no intervalo DT ajustada pelo multiplicador opcional
        # TWX - Lamina em excesso da zona superior (agua livre), passivel de infiltracao na area permeavel ou de geracao de escoamento direto na area impermeavel
        if (UZTWC + PXV) > UZTWM:
            TWX = UZTWC + PXV - UZTWM
            UZTWC = UZTWM
        else:
            TWX = 0.0
            UZTWC = UZTWC + PXV

        # Atualizando o reservatorio da Area Adiciona Impermeavel (TWX nao entra por motivos obvios)
        ADIMC = ADIMC + PXV - TWX

        # Calculando o escoamento superficial da Area Impermeavel Permanente (PCTIM)
        ROIMP = PXV * PCTIM

        #===========================================================================================
        # INICIO DO LOOP INTERNO DE PERCOLACAO/INFILTRACAO
        #===========================================================================================

        # Inicializando os acumuladores que serao atualizados no loop interno
        SDRO  = 0.0 # SOMATORIO do escoamento direto
        SSUR  = 0.0 # SOMATORIO do escoamento superficial
        SIF   = 0.0 # SOMATORIO do escoamento subsuperficial (interflow)
        SPERC = 0.0 # SOMATORIO da percolacao
        SBF_P = 0.0 # SOMATORIO do escoamento de base primario
        SBF_S = 0.0 # SOMATORIO do escoamento de base suplementar

        # Determinando o numero, a duracao e a lamina dos incrementos que serao processados no loop interno
        NINC = int(round(1.0 + 0.20 * (UZFWC + TWX), 0)) # Numero de incrementos
        DINC = dt / NINC                                 # Duracao de cada incremento [dias]
        PINC = TWX / NINC                                # Lamina processada a cada iteracao do loop interno [mm]
                                                         # Segundo a formulacao, substituindo NINC em PINC e fazendo o limite, temos que PINC max = 5.0 mm

        # Fracoes de deplecionamento: [S(t) - S(t+dt)] / S(t) = 1 - exp(-k.dt) ~= 1 - (1-k)**dt (aproximacao da exponencial por Taylor)
        DUZ  = 1.0 - (1.0 -  UZK)**DINC # Fracao do reservatorio superior de agua livre
        DLZP = 1.0 - (1.0 - LZPK)**DINC # Fracao do reservatorio inferior primaerio
        DLZS = 1.0 - (1.0 - LZSK)**DINC # Fracao do reservatorio inferior suplementar

        for j in range(NINC):

            # Calculando o escoamento direto gerado na Area Adicional Impermeavel (representada pela lamina ADIMC)
            RATIO = (ADIMC - UZTWC) / LZTWM
            if RATIO < 0.0: RATIO = 0.0
            ADDRO = PINC * (RATIO**2)

            # Calculando o escoamento de base gerado nos reservatorios inferiores de agua livre

            # RESERVATORIO PRIMARIO
            BF_P = LZFPC * DLZP
            if LZFPC > BF_P:
                LZFPC = LZFPC - BF_P
            else:
                BF_P = LZFPC
                LZFPC = 0.0
            SBF_P = SBF_P + BF_P*(1 - POTIM)

            # RESERVATORIO SUPLEMENTAR
            BF_S = LZFSC * DLZS
            if LZFSC > BF_S:
                LZFSC = LZFSC - BF_S
            else:
                BF_S = LZFSC
                LZFSC = 0.0
            SBF_S = SBF_S + BF_S*(1 - POTIM)


            #===# PERCOLACAO #===#
            # Transferencia de agua livre da zona superior para a zona inferior e balanco entre os reservatorios....

            if (PINC + UZFWC) > 0.01:
            # Condicao 1: HA disponibilidade de agua livre na zona superior, logo A PERCOLACAO EH REALIZADA
                PBASE  = LZFPM*DLZP + LZFSM*DLZS                        # Limite superior de percolacao (na condicao de saturacao do solo)
                SLZ_DEF = (LZTWM-LZTWC) + (LZFSM-LZFSC) + (LZFPM-LZFPC) # Somatorio dos deficits dos reservatorios da zona inferior
                SLZ_CAP = LZTWM + LZFSM + LZFPM                         # Somatorio das capacidades dos reservatorios da zona inferior
                PERC = PBASE * (1.0 + ZPERC*(SLZ_DEF/SLZ_CAP)**REXP)    # Demanda por percolacao, considerando as deficiencias da zona inferior
                PERC = PERC * (UZFWC/UZFWM)                             # Demanda por percolacao, ponderada pela disponibilidade de agua livre na zona superior
                if PERC >= UZFWC:  # A percolacao efetiva nao pode superar o montante de agua livre da zona superior
                    PERC = UZFWC
                if PERC > SLZ_DEF: # A percolacao efetiva tambem nao pode superar o somatorio dos deficit dos reservatorios da zona inferior
                    PERC = SLZ_DEF
                UZFWC = UZFWC - PERC
                SPERC = SPERC + PERC

                # Calculando o escoamento subsuperficial (interflow), ja descontando a percolacao do reservatorio superior de agua livre
                DEL = UZFWC * DUZ
                SIF = SIF + DEL*(1 - POTIM)
                UZFWC = UZFWC - DEL

                # Distribuindo a agua que percola a partir da zona superior entre os reservatorios da zona inferior
                PERCT = PERC * (1.0 - PFREE)
                PERCF = PERC * PFREE

                if (PERCT + LZTWC) > LZTWM:
                    # Enchendo o reservatorio de tensao e transferindo o excesso para os reservatorios de agua livre
                    DEL = PERCT + LZTWC - LZTWM
                    LZTWC = LZTWM
                    PERCF = PERCF + DEL
                else:
                    LZTWC = LZTWC + PERCT

                # A distribuicao eh feita por ponderacao, tendo o reservatorio suplementar maior prioridade
                # HPL, RATLP e RATLS sao fracoes relacionadas a capacidade ou disponibilidade (usadas na ponderacao)
                # PERCP  = montante de percolacao que vai para o reservatorio primario
                # PERCS  = montante de percolacao que vai para o reservatorio suplementar
                if PERCF > 0.0:
                    # Distribuir percolacao entre os reservatorios suplementar (que tem prioridade) e primario
                    HPL = LZFPM / (LZFPM + LZFSM)
                    RATLP = LZFPC / LZFPM
                    RATLS = LZFSC / LZFSM
                    # Calculando a fracao de agua de percolacao que deve ser transferida para o reservatorio primario (em condicoes ideais)
                    FRACP = min(1.0, (HPL * 2.0 * (1.0-RATLP)) / ((1.0-RATLP) + (1.0-RATLS)))
                    # Calculando o montante de agua que deve percolar para o reservatorio suplementar (em condicoes ideias)
                    PERCS = PERCF * (1.0 - FRACP)
                    if (LZFSC + PERCS) > LZFSM:
                        # Nesse caso, a percolacao efetiva corresponde ao deficit do reservatorio suplementar
                        PERCS = LZFSM - LZFSC
                        LZFSC = LZFSM
                    else:
                        # Nesse caso, a percolacao efetiva eh mantida, i.e. PERCS...
                        LZFSC = LZFSC + PERCS
                    # Calculando a percolacao efetiva que vai para o reservatorio primario
                    PERCP = PERCF - PERCS
                    LZFPC = LZFPC + PERCP
                    # Verificar se o reservatorio primario excedeu sua capacidade e, se necessario, transferir agua pro reservatorio de tensao
                    if LZFPC > LZFPM:
                        EXCESS = LZFPC - LZFPM
                        LZTWC = LZTWC + EXCESS
                        LZFPC = LZFPM


                #===# INFILTRACAO #===#
                # Transferencia de agua de chuva da superficie para dentro do solo (na zona superior)
                ADSUR = 0.0 # Escoamento superficial gerado na area ADIM no loop interno (sera atualizado se houver lamina excedente)
                if PINC > 0.0:
                    if (PINC + UZFWC) > UZFWM:
                    # Condicao em que ocorre saturacao do reservatorio superior de agua livre e lamina excedente
                        UZFWC = UZFWM                      # Saturando o reservatorio superior de agua livre
                        SUR = PINC + UZFWC - UZFWM         # Lamina excedente que gera escoamento superficial
                        ADSUR = SUR * (1.0 - ADDRO/PINC)   # Cuidado! ADSUR eh uma lamina, SSUR eh um volume (jah esta ponderado pelas areas de geracao...)
                        ACTIM = PCTIM + ADIMP*(ADDRO/PINC) # Area ativamente impermeavel
                        SSUR = SSUR + SUR*(1 - ACTIM)      # Em conformidade com R.J.C. Burnash (1995)
                    else:
                    # Condicao em que toda a lamina infiltra
                        UZFWC = UZFWC + PINC

            else:
            # Condicao 2: NAO HA disponibilidade de agua livre na zona superior, logo A PERCOLACAO NAO EH REALIZADA
                UZFWC = UZFWC + PINC
                ADSUR = 0.0

            # Balanco de agua na area ADIM (para encontrar o escoamento direto)
            ADIMC = ADIMC + PINC - ADDRO - ADSUR
            if ADIMC > (UZTWM + LZTWM):
                ADDRO = ADDRO + ADIMC - (UZTWM + LZTWM)
                ADIMC = UZTWM + LZTWM
            SDRO = SDRO + ADDRO*ADIMP

        #===========================================================================================
        # FIM DO LOOP INTERNO DE PERCOLACAO/INFILTRACAO
        #===========================================================================================

        estado_UZFWC.append(UZFWC)
        fconv = (area/dt)*(1000.0/86400.0) # Conversao de lamina [mm] em vazao [m3/s]
        # Escoamento superior - Z1 (escoamento na area impermeavel + escoamento direto + escoamento superficial + escoamento subsuperficial)
        # Todos ja estao ponderados pelas respectivas areas de geracao...

        # Escoamentos gerados pelo SAC-SMA
        direto.append((ROIMP + SDRO) * fconv)
        superficial.append(SSUR * fconv)
        subsuperficial.append(SIF * fconv)
        suplementar.append(SBF_S * fconv)
        primario.append(SBF_P * fconv)
        #
        # # Propagacao do escoamento superior
        # Z1.append(SUP*fconv) # m3/s
        #
        # ##########################
        # # !!!IMPLEMENTAR AQUI!!! #
        # ##########################
        #
        # # Escoamento de base suplementar (Z2)
        # Z2.append(SBF_P*fconv)
        #
        # # Escoamento de base primario (Z3)
        # Z3.append(SBF_S*fconv)



    return direto, superficial, subsuperficial, suplementar, primario, estado_UZFWC



        # # APLICAR O SIDE...
        # # Montantes de escoamento de base do loop principal
        # # TBF   - é o escoamento de base total = observado no canal + nao observado no canal (recarga profunda)
        # # BFCC  - escoamento de base observado no canal (baseflow channel component)
        # # BFNCC - escoamento de base nao observado no canal (baseflow NON-channel component)
        # # SIDE  - escoamento de base nao observado no canal / escoamento de base observado no canal
        # # Portanto: SIDE = (TBF - BFCC) / BFCC
        # TBF  = SBF * PAREA
        # BFCC = TBF / (1.0 + SIDE)
        # BFP  = SPBF * PAREA /(1.0 + SIDE) # Escoamento de base primario
        # BFS  = BFCC - BFP                 # Escoamento de base suplementar
        # if BFS < 0.0: BFS = 0.0
        # BFNCC = TBF - BFCC



    #
    #
    #     # Escoamento da zona inferior
    #
    #     # O montante de escoamento direto ocorre na area ADIM e ja foi multiplicado pela area respectiva
    #     # O montante de escoamento superficial, que ocorre nas areas permeavel e na area ADIM, tb ja foi multiplicado
    #     # Portanto, o escoamento total no canal (area PCTIM + direto + superficial + subsuperficial + de base)
    #     TCI = ROIMP + SDRO + SSUR + SIF + BFCC
    #
    #     #APLICAR O RIVA
    #     # Descontando a evapotranspiracao na vegetacao riparia e afeta o escoamento que vai para o canal
    #     EUSED = E1 + E2 + E3      # Evapotranspiração da area permeavel

# from matplotlib import pyplot as plt
# dados['qobs'].plot(label='qobs', linestyle='--', color='black')
# dados['qsim'].plot(label='qsim')
# dados['direto'].plot(label='direto')
# dados['subsuperficial'].plot(label='subsuperficial')
# dados['superficial'].plot(label='superficial')
# plt.legend()
# plt.grid()
# plt.show()

    #     E4 = (EDMND - EUSED)*RIVA # Evapotranspiracao da area com vegetacao riparia, ponderada pela area
    #
    #     if E4 > TCI:
    #         E4 = TCI
    #         TCI = 0.0
    #     else:
    #         E4 = E4
    #         TCI = TCI - E4
    #     qTCI.append(TCI*conv)
    #
    # return qTCI

        ### AGORA ENTRA A PROPAGACAO PROPOSTA POR VRUGT

        # Propagacao

        # Nao propaga
        # Z2 =
        # Z3 =
        # Z = Z1 + Z2 + Z3


    #     # # GRND = Escoamento subterrâneo
    #     # # SURF = Escoamento superficial
    #
    #     GRND = SIF*(1.-params.get(["fSIF"]) + BFCC   #interflow is distrbuted with fSIF
    #     SURF = TCI - GRND
    #     GRND = GRND - E4
    #     if GRND < 0.0:
    #        SURF = SURF + GRND
    #        GRND = 0.0
    #        if SURF < 0.0: SURF = 0.0
    #
    #     #Calcula a evapotranspiração total que ocorreu efetivamente
    #     # TET = Evapotranspiração total ocorrida (mm)
    #     EUSED = EUSED * PAREA
    #     TET = EUSED + E5 + E4
    #
    #     #Verifica se armazenamento da área impermeável é igual ou maior que da zona de tensão superior
    #     if ADIMC < UZTWC: ADIMC = UZTWC
    #
    #
    #     #Armazena vazão gerada pela bacia, de acordo com a origem (separação)
    #     qbac.append(TCI*conv)    # total
    #     qsurf.append(SURF*conv)  # superficial
    #     qgrnd.append(GRND*conv)  # subterrâneo
    #
    #     #Armazena state variablesBF_P = SPBF * PAREA
    #     #if get_states in ["all","last"]:
    #     #    st = [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,None,None,None]
    #     #    states.append(st)
    #     st = [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,None,None,None]
    #     states.append(st)
    #
    #
    # #===========================================================================================================================
    # #FASE CANAL (propagação do escoamento supercial com 3 reservatórios lineares)
    # #===========================================================================================================================
    #
    # #parametros
    # n_rsv = 3
    # k = params.get(["k"] #adimensional
    #
    # #number of stepsBF_P = SPBF * PAREA
    # nsteps = len(qsurf)
    #
    # #condicao de contorno (Qin)
    # #montante
    # if qub is None:
    #     qub = [0. for _ in range(nsteps)]
    #
    # #soma vazao de montante (qub) com bacia (qsurf)
    # qin = [sum(x) for x in zip(qub,qsurf)]
    #
    # #condicao inicial (para cada um dos 3 reservatorios)
    # if q0 is not None:
    #     qinitial = [q0 for _ in range(n_rsv)]
    #     #teste:
    #     #2 primeiros
    #     #qinc = q0 - qmont[0] - qsurf[0]
    #     #qinitial[:-1] = [qinc,qinc]
    #     #ultimo
    #     #qinitial[-1] = q0
    # else:
    #     qinitial = [0. for _ in range(n_rsv)]
    #
    #
    # #lista para vazao de saida ao final de cada passo/reservatorio
    # qout = [0. for _ in range(nsteps)]
    #
    # for i in range(n_rsv):
    #     j=i+1                  #auxiliar para posicao dq Q1,Q2,Q3 em states
    #     qout[0] = qinitial[i]
    #
    #     for t in range(1, nsteps):
    #         qout[t] = qout[t-1]*(2.*k-1.)/(2.*k+1.) + ( qin[t] + qin[t-1] )/(1.+2.*k)
    #         qout[t] = max(0.,qout[t])
    #
    #         #if get_states in ["all","last"]:
    #         states[t][j+5] = qout[t]
    #
    #     #defluente deste reservatorio é afluente do proximo
    #     qin = [x for x in qout]
    #
    # #Calcula vazao total: canal+subterraneo
    # qsim = [sum(x) for x in zip(qout,qgrnd)]
    #
    #
    # #===========================================================================================================================
    # # AJUSTA SAIDA DO MODELO COM VAZAO (OU VAZAO + VARIAVEIS DE ESTADO)
    # #===========================================================================================================================
    # if get_states is not None:
    #     if get_states == "all":
    #         return qsim,states #_final
    #
    #     if get_states == "last":
    #         return qsim,statdiretoes[-1] #_final[-1]
    #
    # if flowsep == True:
    #     return qsim,qgrnd
    #
    # #retorna somente vazao
    # return qsim

# Execucao




params = {}
with open("parametros_jaguari.txt") as arq:
    for line in arq:
       (key, val) = line.split(";")
       params[key] = float(val)

#dados = pd.read_excel("dados_jaguari_diario.xlsx", index_col='date', parse_dates=True, sep=';')
dados = pd.read_csv("dados_maringa.txt", index_col='date', parse_dates=True, sep=';')
etp = dados['etp']
cmb = dados['cmb']
#qin = dados['qin']

direto, superficial, subsuperficial, suplementar, primario, estado_UZFWC = SAC_SMA(etp, cmb, params)
#
dados['qsim'] = [sum(x) for x in zip(direto, superficial, subsuperficial, suplementar, primario)]
dados['direto']         = direto
dados['superficial']    = superficial
dados['subsuperficial'] = subsuperficial
dados['UZFWC'] = estado_UZFWC
# dados['suplementar']    = suplementar
# dados['primario']       = primario
#
#
#
from matplotlib import pyplot as plt
dados['qobs'].plot(label='qobs', linestyle='--', color='black')
dados['qsim'].plot(label='qsim')
dados['direto'].plot(label='direto')
dados['subsuperficial'].plot(label='subsuperficial')
dados['superficial'].plot(label='superficial')
plt.legend()
plt.grid()
plt.show()
