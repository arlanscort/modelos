from numpy import zeros, empty, append
import pandas as pd

def SaCSMA_arlan2019(etp, cmb, params, dt, qin=None, state=None):

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
            SIDE  [0.0]            (adimensional)    Razao entre (recarga profunda) / (escoamento de base do hidrograma)
            RSERV [0.3]            (fracao decimal)
            k     [0.01 - ]= tempo de residencia em dias; constante dos reservatorios lineares de Nash em cascata (NC)
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
    k     = params.get("k", 1.0)
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
        q0rsv1 = 0.0 # mudar
        q0rsv2 = 0.0 # mudar
        q0rsv3 = 0.0 # mudar
    else:
        UZTWC = state[0]
        UZFWC = state[1]
        LZTWC = state[2]
        LZFPC = state[3]
        LZFSC = state[4]
        ADIMC = state[5]
        q0rsv1 = state[6]
        q0rsv2 = state[7]
        q0rsv3 = state[8]
    # Atribuicao de qin, que eh oopcional
    if qin is None:
        qin = zeros([len(etp),])
    # POTIM - percentual referente a Area Potencialmente Impermeavel (R.J.C. Burnash, 1995)
    POTIM = PCTIM + ADIMP
    # Inicializacao do DataFrame que ira conter os escoamentos e variaveis de estado finais do loop externo
    c_entradas    = ['EDMND','PXV','qin']
    c_estados     = ['UZTWC','UZFWC','LZTWC','LZFPC','LZFSC','ADIMC']
    c_escoamentos = ['TWX','PERC','ROIMP','SDRO','SSUR','SIF','SBF_S','SBF_P']
    c_resultados  = ['TET','SUP','QRSV1','QRSV2','QRSV3','Z1','Z2','Z3']
    df_sim = pd.DataFrame(columns=(c_entradas + c_estados + c_escoamentos + c_resultados))
    df_sim.index.name = 'passo'
    # Inicilizando os vetores que irao conter os escoamentos finais
    TET   = zeros(len(etp)) ### DESCREVER...
    SUPq  = zeros(len(etp))
    Z1    = zeros(len(etp))
    Z2    = empty(len(etp))
    Z3    = empty(len(etp))
    QRSVs = empty([len(etp),3])
    # Fator de conversao altura de escoamento [mm ]-> vazao [m3/s]
    fconv = (area/dt)*(1000.0/86400.0)
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
        E5 = E5 * ADIMP
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
            ADDRO = PINC * (RATIO**2.0)
            ### ESCOAMENTO DE BASE
            ## RESERVATORIO PRIMARIO
            BF_P = LZFPC * DLZP
            if LZFPC > BF_P:
                LZFPC = LZFPC - BF_P
            else:
                BF_P = LZFPC
                LZFPC = 0.0
            SBF_P = SBF_P + BF_P*(1.0 - POTIM)
            ## RESERVATORIO SUPLEMENTAR
            BF_S = LZFSC * DLZS
            if LZFSC > BF_S:
                LZFSC = LZFSC - BF_S
            else:
                BF_S = LZFSC
                LZFSC = 0.0
            SBF_S = SBF_S + BF_S*(1.0 - POTIM)
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
                SIF = SIF + DEL*(1.0 - POTIM)
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
                        SSUR = SSUR + SUR*(1.0 - ACTIM)      # (1-ACTIM) = "Area Ativamente Permeavel"; calculo em conformidade com R.J.C. Burnash (1995)
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
        ### APLICACAO DE "RIVA"
        # E4 - Evapotranspiracao na vegetacao riparia (mata ciliar + superficie livre do canal)
        EUSED = E1 + E2 + E3
        E4 = (EDMND - EUSED) * RIVA
        SUP = qin[i]/fconv + ROIMP + SDRO + SSUR + SIF - E4
        if SUP < 0.0:
            E4  = E4 + SUP
            SUP = 0.0
        # Notar que todas as alturas que compoem SUP jah haviam sido ponderadas pelas respectivas areas de geracao
        # Essas areas sao PCTIM, 1-ACTIM, 1-POTIM, RIVA...
        ### CONVERSAO DAS ALTURAS EM VAZOES E APLICACAO DE "SIDE"
        TET[i]  = EUSED*(1-POTIM) + E5 + E4 # TET - Evapotranspiracao total
        SUPq[i] = SUP*fconv                 # SUP - Escoamento "superior" (zona superior + superficie); sera propagado
        Z2[i]   = SBF_S*fconv/(1.0 + SIDE)  # Z2  - Escoamento de base suplementar
        Z3[i]   = SBF_P*fconv/(1.0 + SIDE)  # Z3  - Escoamento de base primario
        ### Salvando as informacoes no DataFrame de simulacao...
        df_sim.loc[i,['EDMND','PXV','qin']] = EDMND, PXV, qin[i]
        df_sim.loc[i,['UZTWC','UZFWC','LZTWC','LZFPC','LZFSC','ADIMC']]  =\
        UZTWC, UZFWC, LZTWC, LZFPC, LZFSC, ADIMC
        df_sim.loc[i,['TWX','PERC','ROIMP','SDRO','SSUR','SIF','SBF_S','SBF_P']] =\
        TWX, PERC, ROIMP*fconv, SDRO*fconv, SSUR*fconv, SIF*fconv, SBF_S*fconv, SBF_P*fconv
        # Verificando condicao para ADIMC:
        if ADIMC < UZTWC : ADIMC = UZTWC
    #-----------------------------------------------------------------------------------------------
    ### FIM DO LOOP EXTERNO
    #-----------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------
    ### PROPAGACAO DE "SUPq" COM 3 RESERVATORIOS LINEARES EM CASCATA
    #-----------------------------------------------------------------------------------------------
    # Premissas:
    # Q = k.S
    # dS/dt = qin - Q ---> dQ/dt = k*(qin - Q)
    if (k < 2.0/dt):
        print('Ok!')
    else:
        print('Lascou...')
    QRSVs[0] = [q0rsv1, q0rsv2, q0rsv3]
    for i in range(1,len(SUPq)):
        qin_i_1 = SUPq[i-1]
        qin_i_2 = SUPq[i]
        for rsv in range(0,3):
            QRSVs[i][rsv] = (k*dt)/(2.0+k*dt)*(qin_i_1 + qin_i_2) + (2.0-k*dt)/(2.0+k*dt)*QRSVs[i-1][rsv]
            qin_i_1 = QRSVs[i-1][rsv]
            qin_i_2 = QRSVs[i][rsv]
    df_sim['TET']   = TET
    df_sim['SUP']   = SUPq
    df_sim['QRSV1'] = QRSVs[:,0]
    df_sim['QRSV2'] = QRSVs[:,1]
    df_sim['QRSV3'] = QRSVs[:,2]
    df_sim['Z1']    = df_sim['QRSV3']
    df_sim['Z2']    = Z2
    df_sim['Z3']    = Z3
    df_sim['qsim']  = df_sim['Z1'] + df_sim['Z2'] + df_sim['Z3']

    return df_sim


params = {}
with open('parametros.txt') as arq:
    for line in arq:
       (key, val) = line.split('\t')
       params[key] = float(val)
epq = pd.read_csv('epq.txt', sep=";", parse_dates=True)
etp = epq['etp'].to_numpy()
cmb = epq['cmb'].to_numpy()
dt = 1.0
df_sim = SaCSMA_arlan2019(etp, cmb, params, dt)
df_sim['qobs'] = epq['qobs']

from matplotlib import pyplot as plt
# df_sim['qobs'].plot(label='qobs', color='black')
# df_sim['qsim'].plot(label='qsim', color='red')
# plt.legend()
# plt.show()

df_sim.to_excel('df_sim.xlsx')

df2 = pd.DataFrame()
df2['primario'] = df_sim['SBF_P']
df2['suplementar'] = df_sim['SBF_P'] + df_sim['SBF_S']
df2['subsuperficial'] = df_sim['SBF_P'] + df_sim['SBF_S'] + df_sim['SIF']
df2['superficial'] = df_sim['SBF_P'] + df_sim['SBF_S'] + df_sim['SIF'] + df_sim['SSUR']
df2['direto'] = df_sim['SBF_P'] + df_sim['SBF_S'] + df_sim['SIF'] + df_sim['SSUR'] + df_sim['SDRO']
df2['impermeavel'] = df_sim['SBF_P'] + df_sim['SBF_S'] + df_sim['SIF'] + df_sim['SSUR'] + df_sim['SDRO'] + df_sim['ROIMP']
df2['qobs'] = df_sim['qobs']
df2[['primario','suplementar','subsuperficial','direto','impermeavel']].plot()
df2['qobs'].plot(label='qobs', color='black')
plt.legend()
plt.show()
