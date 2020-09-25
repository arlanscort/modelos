def SACSMA(et0, cmb, qin, prm, area, q0=None):
    """ Modelo chuva-vazão Sacramento Soil Moisture Accouting.
        Modelo de propagação em canal Cascata de reservatórios conceituais lineares

    Entradas:
    et0 = lista com os dados de evapotranspiração potencial (= et. de referência) [mm em 1 hora];
    cmb = lista com os dados de chuva média na bacia [mm em 1 hora];
    qin = lista com os dados de vazão de montante [m3/s];
    prm = dicionário com os 16-20 parâmetros dos modelos (14-16 da fase bacia <+ 2 multip. dos inputs> + 2 da fase canal);
        {"UZTWM", "UZFWM", "LZTWM", "LZFPM", "LZFSM"} = capacidades máximas dos reservatórios do solo [mm];
        {"UZK", "LZPK", "LZSK"} = taxas de depleção dos reservatórios de água livre [fração/dia]
        {"ZPERC", "REXP"} = coeficiente e expoente da equação de percolação [adim.];
        {"PFREE"} = porção da água percolada que vai direto para os reservatórios de água livre [fração];
        {"PCTIM", "ADIMP"} = porção de área impermeável permanente e de área impermeável adicional [fração];
        {"SIDE"} = porção do escoamento subterrâneo que vai para o canal [fração];
        <{"RSERV", "RIVA"}> = volume na zona inferior não acessível para et. [mm], taxa de evap. da mata ciliar [fração];
        <{"xET0", "xPREC"}> = multiplicadores da evapotranspiração e da chuva média na bacia [adim.];
        {"Kprop"} = taxa de depleção dos reservatórios de propagação [fração/hora];
        {"lag"} = tempo de advecção (deslocamento temporal) da vazão propagada [horas].
    area = área [km2] da sub-bacia simulada (incremental);
    <q0> = vazão inicial para os reservatórios do modelo de propagação. Se não for fornecido usa qin[0] [m3/s].
    
    Saída:
    qcalc = lista com os dados de vazão calculado pelo modelo [m3/s].
    """
    
    #FASE BACIA
    #===========================================================================================================================
    #Inicializando armazenamentos da fase bacia
    UZTWC = prm["UZTWM"] * 0.5
    UZFWC = prm["UZFWM"] * 0.5
    LZTWC = prm["LZTWM"] * 0.5
    LZFPC = prm["LZFPM"] * 0.5
    LZFSC = prm["LZFSM"] * 0.5
    ADIMC = 0.0
    
    #Inserindo valores padrão de parâmetros que podem não ser fornecidos
    if "RSERV" not in prm: prm["RSERV"] = 0.30
    if "RIVA" not in prm: prm["RIVA"] = 0.0
    if "xET0" not in prm: prm["xET0"] = 1.0
    if "xPREC" not in prm: prm["xPREC"] = 1.0
    
    #Lista onde os dados de vazão gerada pela fase bacia serão armazenados
    qbac = []
    
    #Passo de tempo do intervalo [dias]
    DT = 1.0/24.0
    
    #Fator de conversão; X [mm/DT] * conv = Y [m3/s]
    conv = area / (86.40*DT)
    
    #Area permeável
    PAREA = 1.0 - prm["PCTIM"] - prm["ADIMP"]
    
    #Armazenamento total da zona inferior
    LZMAX = (prm["LZTWM"] + prm["LZFPM"] + prm["LZFSM"])
    
    #Tamanho relativo do armazenamento primário comparado com o armazenamento livre inferior total
    HPL = prm["LZFPM"] / (prm["LZFPM"] + prm["LZFSM"])
    
    
    #Iterando o modelo a cada registro da série de dados
    for i in range(len(cmb)):
        
        #Demanda potencial de evapotranspiração e chuva média na bacia ajustados
        EDMND = et0[i] * prm["xET0"]
        PXV = cmb[i] * prm["xPREC"]
        
        #Calculando a perda por evapotranspiração, na zona superior, no intervalo
        # EDMND = Demanda de evapotranspiração | Evapotranspiração de referência | Evapotranspiração Potencial [mm]
        # RED   = Diferença entre a EDMND e evapotranspiração ocorrida no intervalo (mm)
        # E1    = Evaporação ocorrida na zona de tensão superior (mm)
        # E2    = Evaporação ocorrida na zona livre superior (mm)
        # UZRAT = Fração de água em toda a zona superior
        E1 = EDMND * (UZTWC / prm["UZFWM"])
        RED = EDMND - E1
        E2 = 0.0
        
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
                
            else:
                E2 = RED
                UZFWC = UZFWC - E2
                RED = 0.0
            
        else:
            UZTWC = UZTWC - E1
            
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
        RATLZ = (LZTWC + LZFPC + LZFSC - prm["RSERV"]) / (LZMAX - prm["RSERV"])
        
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
        
        #Inicializando acumuladores do intervalo DT
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
        NINC = int(round(1.0 + 0.20 * (UZFWC + TWX), 0))
        DINC = (1.0 / NINC) * DT
        PINC = TWX / NINC
        
        #Calculando frações de deplecionamento da água para o tempo de incremento, sendo que as taxas de depleção são
        #para um dia.
        # DUZ  = Depleção de água da zona superior, por incremento
        # DLZP  = Depleção de água da zona inferior primária, por incremento
        # DLZS  = Depleção de água da zona inferior suplementar, por incremento
        DUZ  = 1.0 - ((1.0 - prm["UZK"])**DINC)
        DLZP = 1.0 - ((1.0 - prm["LZPK"])**DINC)
        DLZS = 1.0 - ((1.0 - prm["LZSK"])**DINC)
        
        
        #Início do loop para os incrementos do intervalo de tempo
        for j in range(NINC):
            
            ADSUR = 0.0
            
            #Calculando escoamento superficial direto da área impermeável adicional
            # ADDRO = Volume (coluna) de run-off direto da área impermeável adicional
            RATIO = (ADIMC - UZTWC) / prm["LZTWM"]
            if RATIO < 0: RATIO = 0.0
            ADDRO = PINC * RATIO**2
            
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
                DEFR = 1.0 - ((LZTWC + LZFPC + LZFSC) / LZMAX)
                PERC = PERC * (1.0 + prm["ZPERC"] * DEFR**prm["REXP"])
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
                        PERCS = PERCS - LZFSC + LZFSM
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
                
            else:
                #Acumulando escoamento superficial direto do incremento
                SDRO = SDRO + ADDRO * prm["ADIMP"]
                
            if ADIMC < 1e-6: ADIMC = 0.0
        
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
        BFCC = TBF * prm["SIDE"]
        BFP = SPBF * PAREA * prm["SIDE"]
        BFS = BFCC - BFP
        if BFS < 0.0: BFS = 0.0
        BFNCC = TBF - BFCC

        #Calculando escoamento afluente da bacia para o canal no intervalo de tempo
        # TCI  = Escoamento afluente total
        # GRND = Escoamento subterrâneo
        # SURF = Escoamento superficial
        TCI = ROIMP + SDRO + SSUR + SIF + BFCC
        GRND = SIF + BFCC   #interflow is part of ground flow
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
        
        #Apenda vazão gerada pela bacia na respectiva lista
        qbac.append(TCI*conv)
        
    #Concluída a fase bacia
    #===========================================================================================================================
    arq.close()
    
    
    #FASE CANAL
    #===========================================================================================================================
    #Número de reservatórios da fase canal
    NRSV = 2
    
    #Vetor de reservatórios e de volumes de propagação
    RSV, PROP = [None for i in range(NRSV)], [None for i in range(NRSV)]
    
    #Vetor de vazão calculada
    qcalc = []
    
    #Inicializando reservatórios de propagação
    if q0 == None:
        q0 = qin[0]
    
    for i in range(NRSV):
        RSV[i] = (q0 * 86400 * DT) / prm["Kprop"]
    
    #Executando modelo de propagação por reservatórios conceituais lineares
    for i in range(len(qbac)):
        
        RSV[0] = RSV[0] + (qbac[i] + qin[i]) * 3600
        PROP[0] = prm["Kprop"] * RSV[0]
        RSV[0] = RSV[0] - PROP[0]
        
        for j in range(1, NRSV):
            
            RSV[j] = RSV[j] + PROP[j-1]
            PROP[j] = prm["Kprop"] * RSV[j]
            RSV[j] = RSV[j] - PROP[j]
        
        qcalc.append(PROP[-1])
    
    #Deslocando série em prm["lag"] horas e transformando em m3/s
    for i in range(len(qcalc)-1, -1, -1):
        
        j = i - prm["lag"]
        
        if j < 0:
            qcalc[i] = qcalc[i] / 3600
            
        else:
            qcalc[i] = qcalc[j] / 3600
    
    #===========================================================================================================================
    return qcalc    

#1230123012301230123012301230123012301230123012301230123012301230123012301230123012301230123012301230123012301230123012301230123