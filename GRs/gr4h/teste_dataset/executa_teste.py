import numpy as np
import pandas as pd
import datetime as dt
import time
from modelo_gr4h import f_gera_OrdHU1, f_gera_OrdHU2
from modelo_gr4h import f_simula_detalhado, f_simula_rapido
from modelo_gr4h import f_calibra_NSE
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt


# Teste para verificar as funcoes de simulacao
start = time.perf_counter()
dataset = pd.read_csv('L0123003.csv', parse_dates=True, index_col="DatesR")
x1 = 521.113
x2 = -2.918
x3 = 218.009
x4 = 4.124
# DF = f_simula_detalhado(x1, x2, x3, x4, dataset['P'], dataset['E'])
# end = time.perf_counter()
# print('Segundos para concluir simulacao detalhada =', end-start)
start = time.perf_counter()
Qsim = f_simula_rapido(x1, x2, x3, x4, dataset['P'], dataset['E'])
end = time.perf_counter()
print('Segundos para concluir simulacao rapida =', end-start)


# # Teste de calibracao
# # 1 - Leitura dos dados (forcantes e vazoes observadas)
# dataset = pd.read_csv('L0123003.csv', parse_dates=True, index_col="DatesR")
# PME = dataset['P']
# ETP = dataset['E']
# Qobs = dataset['Qmm']
# # 2 - Produzir o DataFrame de calibracao (vazoes observadas e indicativo de uso ou nao)
# aquecimento = [0]*1440
# calibracao  = [1]*(len(Qobs)-1440)
# DFcal = pd.Series(data=aquecimento+calibracao, index=Qobs.index, name='Calibra?')
# DFcal = pd.concat([Qobs, DFcal], axis=1)
# # 3 - Passar os argumentos da funcao e a proporcao dos estados iniciais em tupla
# S0_prop = 0.3
# R0_prop = 0.5
# args = (PME, ETP, Qobs, DFcal, S0_prop, R0_prop)
# # 4 - Passar os limites dos parametros a serem otimizados
# bounds = [(50,850), (-10,10), (100,600), (1,10)]
# # 4 - Rodar algoritmo de calibracao
# resultado = differential_evolution(f_calibra_NSE, bounds, updating='deferred', workers=-1, args=args, disp=True)
