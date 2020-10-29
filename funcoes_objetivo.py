'''
Arquivo contendo as funcoes objetivo
'''

import numpy as np


def NSE(simulation, evaluation, h_aq):
    # Nash-Sutcliffe efficiency (NSE)
    # 1 - Limpar
    Qsim = np.array([])
    Qobs = np.array([])
    for parQ in np.array(list(zip(simulation, evaluation))[h_aq:]):
        if np.isnan(parQ[0]) or np.isnan(parQ[1]):
            continue
        else:
            Qsim = np.append(Qsim, parQ[0])
            Qobs = np.append(Qobs, parQ[1])
    # 2 - Calcular
    NSE = 1 - np.sum((Qsim-Qobs)**2)/np.sum((Qobs-np.mean(Qobs))**2)
    return NSE


def sqrtNSE(simulation, evaluation):
    # Nash-Sutcliffe efficiency (NSE) aplicada na raiz das vazoes
    # 1 - Limpar
    Qsim = np.array([])
    Qobs = np.array([])
    for parQ in np.array(list(zip(simulation, evaluation))[h_aq:]):
        if np.isnan(parQ[0]) or np.isnan(parQ[1]):
            continue
        else:
            Qsim = np.append(Qsim, parQ[0])
            Qobs = np.append(Qobs, parQ[1])
    # 2 - Calcular
    sqrtQsim = np.sqrt(Qsim)
    sqrtQobs = np.sqrt(Qobs)
    sqrtQmed = np.sqrt(np.mean(Qobs))
    sqrtNSE = 1-np.sum((sqrtQsim-sqrtQobs)**2)/np.sum((sqrtQobs-sqrtQmed)**2)
    return sqrtNSE


def logNSE(simulation, evaluation):
    # Nash-Sutcliffe efficiency (NSE) aplicada no logaritmo das vazoes
    # 1 - Limpar serie observada
    Qsim = np.array([])
    Qobs = np.array([])
    for parQ in np.array(list(zip(simulation, evaluation))[h_aq:]):
        if np.isnan(parQ[0]) or np.isnan(parQ[1]):
            continue
        else:
            Qsim = np.append(Qsim, parQ[0])
            Qobs = np.append(Qobs, parQ[1])
    # 2 - Calcular criterio de eficiencia
    v = min(Qobs[Qobs>0])
    lnQsim = np.log(Qsim + v)
    lnQobs = np.log(Qobs + v)
    lnQmed = np.log(np.mean(Qobs)+v)
    logNSE = 1-np.sum((lnQsim-lnQobs)**2)/np.sum((lnQobs-lnQmed)**2)
    return logNSE


def KGE(simulation, evaluation):
    # Kling-Gupta efficiency (KGE)
    # 1 - Limpar serie observada
    Qsim = np.array([])
    Qobs = np.array([])
    for parQ in np.array(list(zip(simulation, evaluation))[h_aq:]):
        if np.isnan(parQ[0]) or np.isnan(parQ[1]):
            continue
        else:
            Qsim = np.append(Qsim, parQ[0])
            Qobs = np.append(Qobs, parQ[1])
    # 2 - Calcular criterio de eficiencia
    r = np.corrcoef(Qsim, Qobs)[0,1]   # correlacao linear
    alfa = np.std(Qsim)/np.std(Qobs)   # termo de variabilidade
    beta = np.mean(Qsim)/np.mean(Qobs) # termo de vies
    KGE = 1 - ((r-1)**2 + (alfa-1)**2 + (beta-1)**2)**(1/2)
    return KGE
