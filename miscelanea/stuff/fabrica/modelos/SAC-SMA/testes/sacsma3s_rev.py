from numpy import array
from math import log, exp
from random import random
import numpy as np
import time


### MODELO REVISADO PELO MINO

def sacsma3s_rev(et0, cmb, prm, area, dt, q0=None, qub=None, state=None, get_states=None, flowsep=None):

    """ Modelo chuva-vazão Sacramento Soil Moisture Accouting (SAC-SMA)
        Modelo de propagação em canal (escoamento superficial) com três reservatórios lineares em cascata

        #MS: ASPECTOS REVISADOS
        # - opcao de retornar serie de groundwater, argumento 'flowsep'
        # - controle de DRY na evaporacao da zona superior
        # - 'else' errado no calculo de  SDRO
        # - readaptacao de SIDE para formulacao original
        # - fator de distribucao de interflow 'fSIF' entre SURF e GRND, (default=1., 100% para surf)
        #
        #POSSIVEIS MELHORIAS
        # - verificar balanco hidrico
        # - simplificar propagacao de reservatorio linear, com controle de fluxo
        # - melhorar padrao de argumentos f(cmb,etp,params,states=None,qin=None,**kwargs)
        # - remover q0 dos argumentos
        # - melhorar __docs__ com padrao PEP ou similar
        # - inicializar parametros com variaveis locais
        # - utilizar .get para inicializar parametros opcionais


    Entradas:
    et0 = lista com os dados de evapotranspiração potencial[mm];
    cmb = lista com os dados de chuva média na bacia [mm/dt];
    qub = (opcional),lista com os dados de vazão de montante [m3/s];
    prm = dicionário com os 16-20 parâmetros dos modelos (14-16 da fase bacia <+ 2 multip. dos inputs> + 2 da fase canal);
        {"UZTWM", "UZFWM", "LZTWM", "LZFPM", "LZFSM"} = capacidades máximas dos reservatórios do solo [mm];
        {"UZK", "LZPK", "LZSK"} = taxas de depleção dos reservatórios de água livre [fração/dia]
        {"ZPERC", "REXP"} = coeficiente e expoente da equação de percolação [adim.];
        {"PFREE"} = porção da água percolada que vai direto para os reservatórios de água livre [fração];
        {"PCTIM", "ADIMP"} = porção de área impermeável permanente e de área impermeável adicional [fração];
        {"SIDE"} = porção do escoamento subterrâneo que vai para o canal [fração];
        <{"RSERV", "RIVA"}> = volume na zona inferior não acessível para et. [mm], taxa de evap. da mata ciliar [fração];
        <{"xET0", "xPREC"}> = multiplicadores da evapotranspiração e da chuva média na bacia [adim.];
        {"k"} = fator de depleção dos reservatórios de propagação - tempo de residência [dias];

    area = área [km2] da sub-bacia simulada (incremental);
    dt = passo de tempo em unidades de dias (e.g. passo diário=1., passo horário = 1./24.)
    q0   = (opcional), vazão inicial para reservatórios lineares
    qub  = (opcional, lista com os dados de vazão afluente [m3/s]
    state = lista de condição inicial  [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,Q1,Q2,Q3] #apesar que só opera no solo

    get_states = None, "all" ou "last", controla a saída de variáveis de estado [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,Q1,Q2,Q3]

    Saída:

    output = lista com os dados de vazão calculado pelo modelo [m3/s] (get_states=None)
    output,states = lista com os dados de vazão calculado pelo modelo [m3/s] e variáveis de estado (get_states = "all" ou "last")

    Parametros e limites para calibração:
    Vrugt et al. (2003)
    fargs = {
    "feval": sacsma3s,
    "xl":    [  10.0,      5.0,     1.0,     1.0,    1.0,      0.0,   0.1, 0.0001,  0.01,      1.0,   0.0,     0.0,     0.0,  0.,   0.]
    "xu":    [ 150.0,    150.0,   500.0,  1000.0, 1000.0,     0.40,   0.5,  0.025,  0.25,    250.0,   5.0,     0.1,     0.1, 10.,   1.]
    "xname": ['UZTWM', 'UZFWM', 'LZTWM', 'LZFSM', 'LZFPM', 'ADIMP', 'UZK', 'LZPK', 'LZSK', 'ZPERC', 'REXP', 'PCTIM','PFREE', 'k','fSIF']
    "descript":"fland"
    }

    """

    #--------------------------------------------------------------------------
    # MODELO SAC-SMA + PROPAGAÇÃO DO ESCOAMENTO SUPERFICIAL
    #--------------------------------------------------------------------------

    #lista para armazenar esc. superficial e de base
    qsurf,qgrnd=[],[]
    states=[]

    #Inicializando armazenamentos da fase bacia
    if state is None:
        UZTWC = prm["UZTWM"] * 0.5
        UZFWC = prm["UZFWM"] * 0.5
        LZTWC = prm["LZTWM"] * 0.5
        LZFPC = prm["LZFPM"] * 0.5
        LZFSC = prm["LZFSM"] * 0.5
        ADIMC = UZTWC + LZTWC
    else:
        UZTWC = state[0]
        UZFWC = state[1]
        LZTWC = state[2]
        LZFPC = state[3]
        LZFSC = state[4]
        ADIMC = state[5]



    #Inserindo valores padrão de parâmetros que podem não ser fornecidos
    if "fSIF" not in prm: prm["fSIF"] =1.0 #ms distribuicao de interflow
    if "RSERV" not in prm: prm["RSERV"] = 0.30
    if "RIVA" not in prm: prm["RIVA"] = 0.0
    #if "SIDE" not in prm: prm["SIDE"] = 1.0 #ms
    if "SIDE" not in prm: prm["SIDE"] = 0.0
    if "xET0" not in prm: prm["xET0"] = 1.0
    if "xPREC" not in prm: prm["xPREC"] = 1.0






    #Lista onde os dados de vazão gerada pela fase bacia serão armazenados
    qbac = []

    #Número de segundos no passo de tempo
    secs = 86400 * dt

    #Fator de conversão; X [mm/dt] * conv = Y [m3/s]
    conv = area * 1000 / secs

    #Area permeável
    PAREA = 1.0 - prm["PCTIM"] - prm["ADIMP"]

    #Armazenamento total da zona inferior
    LZMAX = prm["LZTWM"] + prm["LZFPM"] + prm["LZFSM"]

    #Tamanho relativo do armazenamento primário comparado com o armazenamento livre inferior total
    HPL = prm["LZFPM"] / (prm["LZFPM"] + prm["LZFSM"])


    #Numero de passos de tempo
    lenCMB = len(cmb)

    # Intervencao Arlan
    colunas_estados   = ['UZTWC','UZFWC','LZTWC','LZFPC','LZFSC','ADIMC']
    colunas_escoamentos = ['percolacao','direto','superficial','subsuperficial','suplementar','primario']
    df = pd.DataFrame(columns=(colunas_estados + colunas_escoamentos))
    # Intervencao Arlan

    #===========================================================================================================================
    #FASE BACIA
    #===========================================================================================================================
    #Iterando o modelo a cada registro da série de dados
    for i in range(lenCMB):

        #if i>1 and any(np.array(states[i-1]))<0:
        #    print(states[i-1])

        #Demanda potencial de evapotranspiração e chuva média na bacia ajustados
        EDMND = et0[i] * prm["xET0"]
        PXV = cmb[i] * prm["xPREC"]

        #Calculando a perda por evapotranspiração, na zona superior, no intervalo
        # EDMND = Demanda de evapotranspiração | Evapotranspiração de referência | Evapotranspiração Potencial [mm]
        # RED   = Diferença entre a EDMND e evapotranspiração ocorrida no intervalo (mm)
        # E1    = Evaporação ocorrida na zona de tensão superior (mm)
        # E2    = Evaporação ocorrida na zona livre superior (mm)
        # UZRAT = Fração de água em toda a zona superior
        E1 = EDMND * (UZTWC / prm["UZTWM"])
        RED = EDMND - E1
        if RED < 0:
            print(i, RED, EDMND, E1, UZTWC, prm["UZTWM"], (UZTWC / prm["UZTWM"]))
            x=input()
        E2 = 0.0

        DRY = False      #FLAG FOR UZ DRYING
        #Descontando a evaporação da zona de tensão superior, porém não pode ser evaporada mais água do que há nesta camada.
        if UZTWC < E1:
            E1 = UZTWC
            UZTWC = 0.0
            RED = EDMND - E1

            #Descontando o resíduo da PET na zona livre superior
            if UZFWC < RED:
                E2 = UZFWC
                UZFWC = 0.0
                RED = RED - E2
                DRY=True

            else:
                E2 = RED
                UZFWC = UZFWC - E2
                RED = 0.0

        else:
            UZTWC = UZTWC - E1

        if not DRY:
            #Verificando demanda de água pela zona de tensão superior
            if (UZTWC/prm["UZTWM"]) < (UZFWC/prm["UZFWM"]):
                #Fração da água na zona livre superior excedeu a fração na zona de tensão superior, então transfere-se água
                #da zona livre para a de tensão.
                UZRAT = (UZTWC + UZFWC) / (prm["UZTWM"] + prm["UZFWM"])
                UZTWC = UZRAT * prm["UZTWM"]
                UZFWC = UZRAT * prm["UZFWM"]


        #Verificando se os armazenamentos da zona superior secaram
        if UZTWC < 1e-6: UZTWC = 0.0
        if UZFWC < 1e-6: UZFWC = 0.0

        #Calculando a perda por evapotranspiração, na zona inferior, no intervalo
        # E3     = Evaporação ocorrida na zona de tensão inferior (mm)
        # RATLZT = Fração de água na zona de tensão inferior
        # RATLZ  = Fração de água em toda a zona inferior
        # DEL    = Coluna de água transferida da zona livre para a zona de tensão
        E3 = RED * (LZTWC / (prm["UZTWM"] + prm["LZTWM"]))

        if LZTWC < E3:
            E3 = LZTWC
            LZTWC = 0.0

        else:
            LZTWC = LZTWC - E3

        #Verificando demanda de agua pela zona de tensão inferior

        RATLZT = LZTWC / prm["LZTWM"]
        #RATLZ = (LZTWC + LZFPC + LZFSC - prm["RSERV"]) / (LZMAX - prm["RSERV"])

        SAVED = prm["RSERV"]*(prm["LZFPM"]+prm["LZFSM"]) #MS: PINB1 SUBROUTINE
        RATLZ = (LZTWC + LZFPC + LZFSC - SAVED) / (LZMAX - SAVED)



        if RATLZT < RATLZ:

            #Recarregando a zona de tensão inferior com água da zona livre inferior, se houver mais água lá.
            DEL = (RATLZ-RATLZT) * prm["LZTWM"]

            #Transfere água da zona livre inferior suplementar (LZFSC) para a zona de tensão inferior (LZTWC)
            LZTWC = LZTWC + DEL
            LZFSC = LZFSC - DEL

            if LZFSC < 0.0:

                #Se a transferência excedeu LZFSC então o resto vem da zona livre inferior primária (LZFPC)
                LZFPC = LZFPC + LZFSC
                LZFSC = 0.0

        #Verificando se o armazenamento da LZTWC secou
        if LZTWC < 1e-6: LZTWC = 0.0


        #Calculando a perda por evapotranspiração da zona impermeável no intervalo
        # E5 = Evaporação ocorrida na zona impermeavel (mm)
        E5 = E1 + (RED + E2) * ((ADIMC - E1 - UZTWC) / (prm["UZTWM"] + prm["LZTWM"]))

        #Descontando a evaporação do armazenamento da área impermeável
        if ADIMC < E5:
            E5 = ADIMC
            ADIMC = 0.0

        else:
            ADIMC = ADIMC - E5

        #Determinando fração do volume da evapotranspiração na área impermeavel, relativo a toda a evapotranspiração
        #ocorrida na bacia.
        E5 = E5 * prm["ADIMP"]

        #Calculando os escoamentos de percolação e superficial.
        # TWX   = Umidade em excesso na zona de tensão superior, no intervalo (mm)
        # ROIMP = Escoamento superficial da área impermeável
        TWX = PXV + UZTWC - prm["UZTWM"]

        if TWX < 0.0:
            #Se não houve excesso de água na zona de tensão superior...
            UZTWC = UZTWC + PXV
            TWX = 0.0
        else:
            #... do contrário a zona de tensão superior fica saturada
            UZTWC = prm["UZTWM"]

        #Umidade disponível (água que não infiltrou) na zona de tensão superior, vai para o armazenamento da zona impermeável.
        ADIMC = ADIMC + PXV - TWX

        #Calculando o escoamento superficial da área impermeável
        ROIMP = PXV * prm["PCTIM"]

        #arq.write("A %6i %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f" % (i+1, ADIMC, UZTWC, UZFWC, LZTWC, LZFPC, LZFSC))
        #arq.write(" %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n" % (EDMND, E1, E2, E3, E5, PXV, TWX))

        #Inicializando acumuladores do intervalo dt
        SBF   = 0.0    #Escoamento de base
        SSUR  = 0.0    #Escoamento superficial
        SIF   = 0.0    #Escoamento interno (subsuperficial)
        SPERC = 0.0    #Percolação
        SDRO  = 0.0    #Run-off direto
        SPBF  = 0.0    #Escoamento de base da zona livre inferior primária

        #Determinando os incrementos computacionais de tempo para o intervalo básico de tempo.
        #Nenhum incremento irá exceder 5.0 milimetros de UZFWC+TWX.
        # NINC = Número de incrementos de tempo em que o intervalo de tempo será dividido para posterior contabilidade
        #        da umidade do solo.
        # DINC = Comprimento de cada incremento em dias
        # PINC = Quantidade de umidade disponível para cada incremento
        #print(i,'UZFWC;TWX;DINC',UZFWC,TWX)
        NINC = int(round(1.0 + 0.20 * (UZFWC + TWX), 0))
        #NINC = int(1.0 + 0.20 * (UZFWC + TWX))
        DINC = (1.0 / NINC)*dt
        PINC = TWX / NINC

        #Calculando frações de deplecionamento da água para o tempo de incremento, sendo que as taxas de depleção são
        #para um dia.
        # DUZ  = Depleção de água da zona superior, por incremento
        # DLZP  = Depleção de água da zona inferior primária, por incremento
        # DLZS  = Depleção de água da zona inferior suplementar, por incremento
        DUZ  = 1.0 - ((1.0 - prm["UZK"])**DINC)
        DLZP = 1.0 - ((1.0 - prm["LZPK"])**DINC)
        DLZS = 1.0 - ((1.0 - prm["LZSK"])**DINC)

        #if i>0: print(NINC,DINC,PINC,states[i-1][0]/UZTWM,states[i-1][1]/UZFWM,states[i-1][2]/LZTWM,states[i-1][3]/UZTWM,states[i-1][4]/UZTWM)
        #Início do loop para os incrementos do intervalo de tempo
        for j in range(NINC):

            ADSUR = 0.0

            #Calculando escoamento superficial direto da área impermeável adicional
            # ADDRO = Volume (coluna) de run-off direto da área impermeável adicional
            RATIO = (ADIMC - UZTWC) / prm["LZTWM"]
            if RATIO < 0.: RATIO = 0.0
            ADDRO = PINC * (RATIO**2)

            #Calculando o escoamento de base da zona livre inferior primária e o acumulado do intervalo de tempo
            # BF = Escoamento de base
            BF = LZFPC * DLZP

            if LZFPC < BF:
                BF = LZFPC
                LZFPC = 0.0

            else:
                LZFPC = LZFPC - BF

            SBF = SBF + BF
            SPBF = SPBF + BF

            #Calculando o escoamento de base da zona livre inferior suplementar e o acumulado do intervalo de tempo
            BF = LZFSC * DLZS

            if LZFSC < BF:
                BF = LZFSC
                LZFSC = 0.0

            else:
                LZFSC = LZFSC - BF

            SBF = SBF + BF

            #Calculando o volume percolado (se não houver água disponível, pula esta etapa)
            if (PINC + UZFWC) <= 0.010:
                UZFWC = UZFWC + PINC

            else:
                #Há água, calcula percolação.
                # PERC = Volume percolado no incremento de tempo
                # DEFR = Taxa de deficiência de umidade da zona inferior do solo
                PERCM = prm["LZFPM"] * DLZP + prm["LZFSM"] * DLZS
                PERC = PERCM * (UZFWC/prm["UZFWM"])
                DEFR = 1. - ((LZTWC + LZFPC + LZFSC) / LZMAX)

                #print(i,j,'PREC;DEFR;',prm["ZPERC"] * DEFR)
                #if i >=1:print(states[i-1])
                PERC = PERC * (1.0 + prm["ZPERC"] * DEFR**prm["REXP"])

                ##Há água, calcula percolação.
                ## PERC = Volume percolado no incremento de tempo
                ## DEFR = Taxa de deficiência de umidade da zona inferior do solo
                ## A = Parâmetro apresentado em SINGH (1995) para melhorar a calibração de ZPERC
                #PERCM = prm["LZFPM"] * DLZP + prm["LZFSM"] * DLZS
                #DEFR = 1.0 - ((LZTWC + LZFPC + LZFSC) / LZMAX)
                #A = ((prm["UZFWM"] - PERCM) / (PERCM * prm["ZPERC"])) ** (1.0/prm["REXP"])
                #if DEFR <= A:
                #    PERC = (PERCM + (prm["UZFWM"] - PERCM) * (DEFR/A)**prm["REXP"]) * UZFWC/prm["UZFWM"]
                #else:
                #    PERC = UZFWC

                #OBS: A percolação ocorre da zona livre superior antes de PINC ser adicionado

                if PERC > UZFWC:
                    #Percolação não pode exceder o armazenamento da zona livre superior
                    PERC = UZFWC

                UZFWC = UZFWC - PERC

                #Verifica se a percolação excedeu a deficiência da zona inferior
                CHECK = LZTWC + LZFPC + LZFSC + PERC - LZMAX

                if CHECK > 0.0:
                    #Volume dos armazenamentos das zonas inferiores mais percolação não deve exceder a capacidade
                    #máxima da zona inferior.
                    PERC = PERC - CHECK

                    #Devolvendo excesso à zona superior
                    UZFWC = UZFWC + CHECK

                #Acumulando a percolação dos incrementos
                SPERC = SPERC + PERC

                #Calculando o escoamento interno e o acumulado
                # DEL = Escoamento interno
                #OBS: A quantidade PINC ainda não foi adicionada
                DEL = UZFWC * DUZ
                SIF = SIF + DEL
                UZFWC = UZFWC - DEL

                #Distribuir a água percolada entre as zonas inferiores, sendo que a zona de tensão deve ser preenchida antes,
                #com excessão da percolação ocorrida na área PFREE.
                # PERCT = Percolação que vai para a zona de tensão inferior
                # PERCF = Percolação que vai direto para a zona livre inferior
                PERCT = PERC * (1.0 - prm["PFREE"])

                if (PERCT + LZTWC) > prm["LZTWM"]:
                    #Excesso irá para a zona livre inferior
                    PERCF = PERCT + LZTWC - prm["LZTWM"]
                    LZTWC = prm["LZTWM"]

                else:
                    #Zona de tensão inferior recebe água percolada
                    LZTWC = LZTWC + PERCT
                    PERCF = 0.0

                #Distribuir a água percolada para a zona livre inferior entre os armazenamentos de água livre.
                # RATLP = Fração de água na zona livre inferior primária
                # RATLS = Fração de água na zona livre inferior suplementar
                # FRACP = Fração da percolação que vai para zona livre primária
                # PERCP = Quantidade da percolação que vai para a zona livre primária
                # PERCS = Quantidade da percolação que vai para a zona livre suplementar
                # EXCESS = Eventual excesso da capacidade máxima da zona livre inferior primária
                PERCF = PERCF + PERC * prm["PFREE"]

                if PERCF > 0.0:
                    #Distribuindo percolação
                    RATLP = LZFPC / prm["LZFPM"]
                    RATLS = LZFSC / prm["LZFSM"]
                    FRACP = min(1.0, (HPL * 2.0 * (1.0-RATLP)) / ((1.0-RATLP) + (1.0-RATLS)))
                    PERCP = PERCF * FRACP
                    PERCS = PERCF - PERCP

                    #Adicionando o excesso de percolação na zona suplementar
                    LZFSC = LZFSC + PERCS

                    if LZFSC > prm["LZFSM"]:
                        #A adição do excesso da percolação não pode exceder a capacidade máxima da zona livre suplementar
                        PERCS = PERCS - LZFSC + prm["LZFSM"]
                        LZFSC = prm["LZFSM"]

                    #Adicionando o excedente de percolação na zona primária
                    LZFPC = LZFPC + (PERCF - PERCS)

                    #Verificar para ter certeza que o armazenamento livre primário não excedeu a capacidade máxima
                    if LZFPC > prm["LZFPM"]:
                        EXCESS = LZFPC - prm["LZFPM"]
                        LZTWC = LZTWC + EXCESS
                        LZFPC = prm["LZFPM"]

                #Distribuir PINC entre a zona superior livre e escoamento superficial
                if PINC > 0.0:

                    #Verificar se o acréscimo de PINC excede a capacidade máxima da zona livre superior
                    if (PINC+UZFWC) > prm["UZFWM"]:
                        #Calcular o escoamento superficial e a soma dos incrementos
                        # SUR   = Escoamento superficial
                        # ADSUR = Quantidade do escoamento direto que provém da porção ADIMP que não está gerando
                        #         escoamento superficial direto no momento
                        # ADDRO/PINC = Fração da área ADIMP que está gerando escoamento superficial direto no momento
                        SUR = PINC + UZFWC - prm["UZFWM"]
                        print(SUR)
                        UZFWC = prm["UZFWM"]
                        SSUR = SSUR + SUR*PAREA
                        ADSUR = SUR * (1.0 - ADDRO/PINC)
                        SSUR = SSUR + ADSUR * prm["ADIMP"]
                    else:
                        #Não excedeu, ou seja, toda a água infiltra, logo não haverá escoamento superficial
                        UZFWC = UZFWC + PINC

            #Balanço de água da área impermeável
            ADIMC = ADIMC + PINC - ADDRO - ADSUR

            if ADIMC > (prm["UZTWM"] + prm["LZTWM"]):
                ADDRO = ADDRO + ADIMC - (prm["UZTWM"] + prm["LZTWM"])
                ADIMC = prm["UZTWM"] + prm["LZTWM"]

            #else: #<--POSSIVEL ERRO ANGELO
            #Acumulando escoamento superficial direto do incremento <<MS:CUIDADO!
            SDRO = SDRO + ADDRO * prm["ADIMP"]

            if ADIMC < 1e-6: ADIMC = 0.0

            #arq.write("B %6i %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n" % (j+1, SBF, SSUR, SIF, SPERC, SDRO, SPBF))


        #Fim do loop para os incrementos do intervalo de tempo

        #Calcula os acumulados e ajusta a quantidade de run-off pela área em que ele foi gerado
        # EUSED = Evaporação ocorrida na fração de área PAREA, durante o intervalo de tempo
        EUSED = E1 + E2 + E3
        SIF = SIF * PAREA

        #Separando componente do escoamento de base que vai para o canal, da componente que não vai para o canal
        # TBF = é o escoamento de base total
        # BFCC = componente do escoamento de base que vai para o canal
        # BFNCC = componente do escoamento de base que NÃO vai para o canal
        TBF = SBF * PAREA
        BFCC = TBF * (1.0/(1.0+prm["SIDE"]))   #MS
        BFP = SPBF * PAREA /(1.0+ prm["SIDE"]) #MS
        BFS = BFCC - BFP
        if BFS < 0.0: BFS = 0.0
        BFNCC = TBF - BFCC

        #Calculando escoamento afluente da bacia para o canal no intervalo de tempo
        # TCI  = Escoamento afluente total
        # GRND = Escoamento subterrâneo
        # SURF = Escoamento superficial
        TCI = ROIMP + SDRO + SSUR + SIF + BFCC
        GRND = SIF*(1.-prm["fSIF"]) + BFCC   #interflow is distrbuted with fSIF
        SURF = TCI - GRND


        #Calcula da evapotranspiração da vegetação ciliar
        # E4 = Evapotranspiração da mata ciliar (mm)
        E4 = (EDMND - EUSED) * prm["RIVA"]

        #Subtrai a evapotranspiração da mata ciliar do escoamento afluente para o canal
        if E4 > TCI:
            E4 = TCI
            TCI = 0.0

        else:
            TCI = TCI - E4

        GRND = GRND - E4
        if GRND < 0.0:
           SURF = SURF + GRND
           GRND = 0.0
           if SURF < 0.0: SURF = 0.0

        #Calcula a evapotranspiração total que ocorreu efetivamente
        # TET = Evapotranspiração total ocorrida (mm)
        EUSED = EUSED * PAREA
        TET = EUSED + E5 + E4

        #Verifica se armazenamento da área impermeável é igual ou maior que da zona de tensão superior
        if ADIMC < UZTWC: ADIMC = UZTWC


        #Armazena vazão gerada pela bacia, de acordo com a origem (separação)
        qbac.append(TCI*conv)    # total
        qsurf.append(SURF*conv)  # superficial
        qgrnd.append(GRND*conv)  # subterrâneo

        #Armazena state variables
        #if get_states in ["all","last"]:
        #    st = [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,None,None,None]
        #    states.append(st)
        st = [UZTWC,UZFWC,LZTWC,LZFPC,LZFSC,ADIMC,None,None,None]
        states.append(st)

        fconv = (area/dt)*(1000.0/86400.0)
        df.loc[i,['UZTWC','UZFWC','LZTWC','LZFPC','LZFSC','ADIMC']] = UZTWC, UZFWC, LZTWC, LZFPC, LZFSC, ADIMC
        df.loc[i,['percolacao','direto','superficial','subsuperficial','primario']] = PERC*fconv, (ROIMP + SDRO)*fconv, SSUR*fconv, SIF*fconv, SPBF*PAREA*fconv


    #===========================================================================================================================
    #FASE CANAL (propagação do escoamento supercial com 3 reservatórios lineares)
    #===========================================================================================================================

    #parametros
    n_rsv = 3
    k = prm["k"] #adimensional

    #number of steps
    nsteps = len(qsurf)

    #condicao de contorno (Qin)
    #montante
    if qub is None:
        qub = [0. for _ in range(nsteps)]

    #soma vazao de montante (qub) com bacia (qsurf)
    qin = [sum(x) for x in zip(qub,qsurf)]

    #condicao inicial (para cada um dos 3 reservatorios)
    if q0 is not None:
        qinitial = [q0 for _ in range(n_rsv)]
        #teste:
        #2 primeiros
        #qinc = q0 - qmont[0] - qsurf[0]
        #qinitial[:-1] = [qinc,qinc]
        #ultimo
        #qinitial[-1] = q0
    else:
        qinitial = [0. for _ in range(n_rsv)]


    #lista para vazao de saida ao final de cada passo/reservatorio
    qout = [0. for _ in range(nsteps)]

    for i in range(n_rsv):
        j=i+1                  #auxiliar para posicao dq Q1,Q2,Q3 em states
        qout[0] = qinitial[i]

        for t in range(1, nsteps):
            qout[t] = qout[t-1]*(2.*k-1.)/(2.*k+1.) + ( qin[t] + qin[t-1] )/(1.+2.*k)
            qout[t] = max(0.,qout[t])

            #if get_states in ["all","last"]:
            states[t][j+5] = qout[t]

        #defluente deste reservatorio é afluente do proximo
        qin = [x for x in qout]

    #Calcula vazao total: canal+subterraneo
    qsim = [sum(x) for x in zip(qout,qgrnd)]


    #===========================================================================================================================
    # AJUSTA SAIDA DO MODELO COM VAZAO (OU VAZAO + VARIAVEIS DE ESTADO)
    #===========================================================================================================================
    if get_states is not None:
        if get_states == "all":
            return qsim,states #_final

        if get_states == "last":
            return qsim,states[-1] #_final[-1]

    if flowsep == True:
        return qsim,qgrnd

    #retorna somente vazao
    return qsim, df


bacia = 'atibaia_valinhos'

params = {}
with open(bacia + '/parametros.txt') as arq:
    for line in arq:
       (key, val) = line.split()
       params[key] = float(val)

import pandas as pd
epq = pd.read_csv(bacia + '/epq.txt', sep=";")
et0 = epq['etp']
cmb = epq['cmb']
prm = params
area = params["area"]
dt = 1.0

epq['qsim'], df = sacsma3s_rev(et0, cmb, prm, area, dt)

from matplotlib import pyplot as plt
epq['qobs'].plot(label='qobs')
epq['qsim'].plot(label='qsim')
plt.show()
df.to_excel(bacia + '/resultados_revmino.xlsx')
