import numpy as np
import pandas as pd

def SaCSMA_arlan2019(etp, cmb, params, dt, qin=None, state=None, get_states=None):

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

        ### Comentario de Arlan Almeida, agosto de 2019
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
        ###############################################################################################################
        dt    = passo de tempo do loop principal em unidades de dias; corresponde a resolucao das series de entrada
        state = lista de condição inicial  [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,Q1,Q2,Q3] #apesar que só opera no solo
        get_states = None, "all" ou "last", controla a saída de variáveis de estado [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,Q1,Q2,Q3]
        output = lista com os dados de vazão calculado pelo modelo [m3/s] (get_states=None)
    """
    # Aquisicao de parametros
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
    NCK   = params.get("NCK",  1.0)
    xETP  = params.get("xETP", 1.0)
    xCMB  = params.get("xCMB", 1.0)
    area = params.get("area")
    # Inicializacao das variaveis de estado
    if state is None:
        UZTWC = UZTWM * 0.5
        UZFWC = UZFWM * 0.5
        LZTWC = LZTWM * 0.5
        LZFPC = LZFPM * 0.5
        LZFSC = LZFSM * 0.5
        ADIMC = UZTWC + LZTWC
    else:
        UZTWC = state[0]
        UZFWC = state[1]
        LZTWC = state[2]
        LZFPC = state[3]
        LZFSC = state[4]
        ADIMC = state[5]
    # POTIM - percentual referente a Area Potencialmente Impermeavel (R.J.C. Burnash, 1995)
    POTIM = PCTIM + ADIMP
    # Inicializacao do DataFrame que ira conter os escoamentos e variaveis de estado finais do loop externo
    colunas_estados   = ['UZTWC','UZFWC','LZTWC','LZFPC','LZFSC','ADIMC']
    colunas_escoamentos = ['PERC','ROIMP+SDRO','SSUR','SIF','SBF_S','SBF_P']
    df = pd.DataFrame(columns=(colunas_estados + colunas_escoamentos))
    #-----------------------------------------------------------------------------------------------
    ### INICIO DO LOOP EXTERNO
    #-----------------------------------------------------------------------------------------------
    for i in range(len(etp)):
        ### PERDAS POR EVAPOTRANSPIRACAO
        EDMND = etp[i] * xETP # Demanda por evapotranspiracao / evapotranspira potencial no intervalo dt
        ## ZONA SUPERIOR
        # E1 - Evapotranspiracao efetiva no reservatorio superior de agua de tensao (UZTW)
        # E2 - Evapotranspiracao efetiva no reservatorio superior de agua livre (UZFW)s
        E1 = EDMND * (UZTWC / UZTWM)
        if E1 < UZTWC:
            E1 = E1 # Essa sintaxe, e as proximas similares, embora redundante, enfatiza que a ETR = ETP
            UZTWC = UZTWC - E1
            E2 = 0.0
        else:
            E1 = UZTWC
            UZTWC = 0.0
            E2 = EDMND - E1
            if E2 > UZFWC:
                E2 = UZFWC
                UZFWC = 0.0
            else:
                E2 = E2
                UZFWC = UZFWC - E2
        if (UZFWC/UZFWM) > (UZTWC/UZTWM):
        # Equilibrar os armazenamentos da zona superior, transferindo agua livre para o reservatorio de tensao
            UZRAT = (UZTWC + UZFWC) / (UZTWM + UZFWM)
            UZFWC = UZFWM * UZRAT
            UZTWC = UZTWM * UZRAT
        ## ZONA INFERIOR
        # E3     - Evapotranspiracao efetiva no reservatorio inferior de agua de tensao (LZTW)
        # RATLZT - Armazenamento relativo do LZTW
        # RATLZ  - Armazenamento relativo dos reservatorios da zona inferior (LZTW+LZFP+LZFS) descontando SAVED
        E3 = (EDMND - E1 - E2) * LZTWC / (UZTWM + LZTWM)
        if E3 > LZTWC:
            E3 = LZTWC
            LZTWC = 0.0
        else:
            E3 = E3
            LZTWC = LZTWC - E3
        SAVED  = RSERV * (LZFPM + LZFSM)
        RATLZT = LZTWC / LZTWM
        RATLZ  = (LZTWC + (LZFPC+LZFSC-SAVED)) / (LZTWM + (LZFPM+LZFSM-SAVED))
        if (RATLZT < RATLZ):
        # Equilibrar os armazenamentos da zona inferior, transferindo agua livre para o reservatorio de tensao
            DEL = (RATLZ - RATLZT) * LZTWM # DEL significa transferencia de agua livre (ocorre em outras partes do codigo tambem)
            LZTWC = LZTWC + DEL
            if DEL <= LZFSC:
            # Primeiro, tentar suprir a demanda (DEL) com agua livre do reservatorio suplementar
                LZFSC = LZFSC - DEL
            else:
            # Se nao for suficiente, suprir a demanda (DEL) com agua livre do reservatorio primario
                DEL = DEL - LZFSC
                LZFSC = 0.0
                LZFPC = LZFPC - DEL
        ## AREA IMPERMEAVEL ADICIONAL (ADIM)
        # E5 - Evapotranspiracao efetiva na Area Impermeavel Adicional (ADIM)
        E5 = E1 + (EDMND - E1) * ((ADIMC - E1 - UZTWC) / (UZTWM + LZTWM))
        if E5 > ADIMC:
            E5 = ADIMC
            ADIMC = 0.0
        else:
            E5 = E5
            ADIMC = ADIMC - E5
        ### LAMINA EXCEDENTE NA ZONA SUPERIOR (TWX)
        PXV = cmb[i] * xCMB # Altura de precipitacao acumulada no intervalo DT ajustada pelo multiplicador opcional
        # TWX - Lamina em excesso da zona superior (agua livre), passivel de infiltracao na area permeavel ou de geracao de escoamento direto na area impermeavel
        if (UZTWC + PXV) > UZTWM:
            TWX = UZTWC + PXV - UZTWM
            UZTWC = UZTWM
        else:
            TWX = 0.0
            UZTWC = UZTWC + PXV
        ADIMC = ADIMC + PXV - TWX # Desconsidera TWX porque vira escoamento direto na parte impermeavel ou escoamento superficial na parte permeavel
        ### ESCOAMENTO DIRETO - AREA IMPERMEAVEL PERMANENTE
        ROIMP = PXV * PCTIM
        # Inicializacao de constantes e acumuladores para o loop interno
        NINC = int(round(1.0+0.20*(UZFWC+TWX), 0)) # NINC - Numero de incrementos
        DINC = dt / NINC                           # DINC - Duracao de cada incremento [dias]
        PINC = TWX / NINC                          # PINC - Altura de precipitacao processada a cada incremento [mm] (max PINC = 5.0)
        DUZ  = 1.0 - (1.0 -  UZK)**DINC            # DUZ - Fator de deplecao do reservatorio superior de agua livre (escoamento subsuperficial)
        DLZP = 1.0 - (1.0 - LZPK)**DINC            # DLZP - Fator de deplecao do reservatorio inferior primario (escoamento de base lento)
        DLZS = 1.0 - (1.0 - LZSK)**DINC            # DLZS - Fator de deplecao do reservatorio inferior suplementar (escoamento de base rapido)
        SDRO  = 0.0                                # Acumulador dos montantes de escoamento direto (direct runoff)
        SSUR  = 0.0                                # Acumulador dos montantes de escoamento superficial (surface runoff)
        SIF   = 0.0                                # Acumulador dos montantes de escoamento subsuperficial (interflow)
        SPERC = 0.0                                # Acumulador dos montantes de percolacao
        SBF_P = 0.0                                # Acumulador dos montantes do escoamento de base primario
        SBF_S = 0.0                                # Acumulador dos montantes do escoamento de base suplementar
        #-------------------------------------------------------------------------------------------
        ### INICIO DO LOOP INTERNO
        #-------------------------------------------------------------------------------------------
        for j in range(NINC):
            ### ESCOAMENTO DIRETO - AREA IMPERMEAVEL ADICIONAL (ADIM) - PARTE 1
            RATIO = (ADIMC - UZTWC) / LZTWM
            if RATIO < 0.0: RATIO = 0.0
            ADDRO = PINC * (RATIO**2)
            ### ESCOAMENTO DE BASE
            ## RESERVATORIO PRIMARIO
            BF_P = LZFPC * DLZP
            if LZFPC > BF_P:
                LZFPC = LZFPC - BF_P
            else:
                BF_P = LZFPC
                LZFPC = 0.0
            SBF_P = SBF_P + BF_P*(1 - POTIM)
            ## RESERVATORIO SUPLEMENTAR
            BF_S = LZFSC * DLZS
            if LZFSC > BF_S:
                LZFSC = LZFSC - BF_S
            else:
                BF_S = LZFSC
                LZFSC = 0.0
            SBF_S = SBF_S + BF_S*(1 - POTIM)
            ### PERCOLACAO (transferencia de agua livre da zona superior para os reservatorios da zona inferior)
            if (PINC + UZFWC) > 0.01:
            # Como existe disponibilidade de agua livre na zona superior, ocorre a percolacao...
                ## PERCOLACAO - PROCESSAMENTO DA ZONA SUPERIOR
                PBASE  = LZFPM*DLZP + LZFSM*DLZS                        # PBASE   - Percolacao maxima que ocorre na condicao de saturacao do solo
                SLZ_DEF = (LZTWM-LZTWC) + (LZFSM-LZFSC) + (LZFPM-LZFPC) # SLZ_DEF - Somatorio dos deficits dos reservatorios da zona inferior
                SLZ_CAP = LZTWM + LZFSM + LZFPM                         # SLZ_CAP - Somatorio das capacidades dos reservatorios da zona inferior
                PERC = PBASE * (1.0 + ZPERC*(SLZ_DEF/SLZ_CAP)**REXP)    # PERC    - Demanda por percolacao, considerando as deficiencias da zona inferior
                PERC = PERC * (UZFWC/UZFWM)                             # PERC    - Demanda por percolacao, ponderada pela disponibilidade de agua livre na zona superior
                if PERC >= UZFWC:
                # A percolacao efetiva nao pode superar o montante de agua livre da zona superior
                    PERC = UZFWC
                if PERC > SLZ_DEF:
                # A percolacao efetiva tambem nao pode superar o somatorio dos deficit dos reservatorios da zona inferior
                    PERC = SLZ_DEF
                UZFWC = UZFWC - PERC
                SPERC = SPERC + PERC
                ### ESCOAMENTO SUBSUPERFICIAL (INTERFLOW) - Notar que eh calculado apos a percolacao de agua livre da zona superior...
                DEL = UZFWC * DUZ
                SIF = SIF + DEL*(1 - POTIM)
                UZFWC = UZFWC - DEL
                ## PERCOLACAO - PROCESSAMENTO DA ZONA INFERIOR
                PERCT = PERC * (1.0 - PFREE) # PERCT - montante que, a priori, percola para o reservatorio inferior de agua de tensao
                PERCF = PERC * PFREE         # PERCF - montante que, a priori, percola para o reservatorio inferior de agua livre
                if (PERCT + LZTWC) > LZTWM:
                # Enchendo o reservatorio de tensao e transferindo o excesso para os reservatorios de agua livre
                    DEL = PERCT + LZTWC - LZTWM
                    LZTWC = LZTWM
                    PERCF = PERCF + DEL
                else:
                    LZTWC = LZTWC + PERCT
                # A distribuicao da percolacao entre os reservatorios inferiores eh feita por ponderacao, tendo o reservatorio suplementar maior prioridade...
                # HPL, RATLP e RATLS sao fracoes relacionadas a capacidade ou disponibilidade (usadas na ponderacao)
                # PERCP - Montante de percolacao que vai para o reservatorio primario
                # PERCS  = montante de percolacao que vai para o reservatorio suplementar
                if PERCF > 0.0:
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
                    # Nesse caso, a percolacao efetiva eh mantida
                        LZFSC = LZFSC + PERCS
                    # Calculando a percolacao efetiva que vai para o reservatorio primario
                    PERCP = PERCF - PERCS
                    LZFPC = LZFPC + PERCP
                    if LZFPC > LZFPM:
                    # Verificar se o reservatorio primario excedeu sua capacidade e, se necessario, transferir agua pro reservatorio de tensao
                        EXCESS = LZFPC - LZFPM
                        LZTWC = LZTWC + EXCESS
                        LZFPC = LZFPM
                ### INFILTRACAO (transferencia de agua da chuva para o solo - notar que a precipitacao entra nessa etapa)
                ADSUR = 0.0 # Montante de escoamento superficial gerado na ADIM (a priori eh nulo)
                if PINC > 0.0:
                    if (PINC + UZFWC) > UZFWM:
                    # Condicao em que ocorre saturacao, geracao de lamina excedente e consequentemente SUR e ADSUR
                        SUR = PINC + UZFWC - UZFWM         # Lamina excedente que gera escoamento superficial nas areas permeaveis que nao compoem ADIMP
                        UZFWC = UZFWM                      # Saturando o reservatorio superior de agua livre
                        ADSUR = SUR * (1.0 - ADDRO/PINC)   # Lamina excedente que gera escoamento superficial nas areas permeaveis de ADIMP
                        ACTIM = PCTIM + ADIMP*(ADDRO/PINC) # ACTIM - Area Ativamente Impermeavel, dadas das condicoes de umidade do solo
                        SSUR = SSUR + SUR*(1 - ACTIM)      # (1-ACTIM) = "Area Ativamente Permeavel"; calculo em conformidade com R.J.C. Burnash (1995)
                    else:
                    # Condicao em que toda a agua infiltra, nao gera lamina excedente e, portanto, nao gera SUR nem ADSUR
                        UZFWC = UZFWC + PINC
            else:
            # Como nao existe disponibilidade de agua livre na zona superior, nao ocorre percolacao
                UZFWC = UZFWC + PINC
                ADSUR = 0.0
            ### ESCOAMENTO DIRETO - AREA IMPERMEAVEL ADICIONAL (ADIM) - PARTE 2
            ADIMC = ADIMC + PINC - ADDRO - ADSUR
            if ADIMC > (UZTWM + LZTWM):
                ADDRO = ADDRO + ADIMC - (UZTWM + LZTWM)
                ADIMC = UZTWM + LZTWM
            SDRO = SDRO + ADDRO*ADIMP
        #-------------------------------------------------------------------------------------------
        ### FIM DO LOOP INTERNO
        #-------------------------------------------------------------------------------------------
        if ADIMC < UZTWC: ADIMC = UZTWC

        df.loc[i,['UZTWC','UZFWC','LZTWC','LZFPC','LZFSC','ADIMC']] = UZTWC, UZFWC, LZTWC, LZFPC, LZFSC, ADIMC
        df.loc[i,['PERC','ROIMP+SDRO','SSUR','SIF','SBF_S','SBF_P']] = PERC, ROIMP+SDRO, SSUR, SIF, SBF_S, SBF_P
    #-----------------------------------------------------------------------------------------------
    ### FIM DO LOOP EXTERNO
    #-----------------------------------------------------------------------------------------------

    # Conversao das laminas em vazoes...
    fconv = (area/dt)*(1000.0/86400.0)
    df[['PERC','ROIMP+SDRO','SSUR','SIF','SBF_S','SBF_P']] = df[['PERC','ROIMP+SDRO','SSUR','SIF','SBF_S','SBF_P']].applymap(lambda x: x*fconv)

    # Pos processamento das laminas geradas...

        #fconv = (area/dt)*(1000.0/86400.0) # Conversao de lamina [mm] em vazao [m3/s]
        # Escoamento superior - Z1 (escoamento na area impermeavel + escoamento direto + escoamento superficial + escoamento subsuperficial)
        # Todos ja estao ponderados pelas respectivas areas de geracao...

        # Escoamentos gerados pelo SAC-SMA

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



    return df



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
    # E5 = E5 * ADIMP
    #     TET = EUSED + E5 + E4
    #
    # E5 = E5 * ADIMP # O limite superior de ADIM eh por ADIMP, que eh uma fracao decimal a ser multiplicada pela area da bacia
            # O limite superior de ADIM eh representado por ADIMP, que eh uma fracao decimal que deve ser multiplicada pela area total da bacia

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

bacia = 'atibaia_valinhos'
params = {}
with open(bacia + '/parametros.txt') as arq:
    for line in arq:
       (key, val) = line.split()
       params[key] = float(val)

epq = pd.read_csv(bacia + '/epq.txt', sep=";")
etp = epq['etp'].to_numpy()
cmb = epq['cmb'].to_numpy()
qin = epq['qin'].to_numpy()

dt = 1.0
df = SaCSMA_arlan2019(etp, cmb, params, dt)
df['qsim'] = df['ROIMP+SDRO'] + df['SSUR'] + df['SIF'] + df ['SBF_S'] + df['SBF_P']
df['qobs'] = epq['qobs']

df['qsim'].plot(label='qsim')
df['qobs'].plot(label='qobs')
from matplotlib import pyplot as plt
plt.legend()
plt.show()
df.to_excel(bacia+'/resultados_arlan.xlsx')
