import numpy as np
import pandas as pd
import datetime as dt
import time
from modelo_gr4h import f_gera_OrdHU1, f_gera_OrdHU2
from modelo_gr4h import f_simula_detalhado, f_simula_rapido
from modelo_gr4h import f_calibra_NSE
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt


# Teste para verificar funcao de simulacaoß detalhada
start = time.perf_counter()
dataset = pd.read_csv('L0123003.csv', parse_dates=True, index_col="DatesR")
PME  = dataset['P']
ETP  = dataset['E']
Qobs = dataset['Qmm']
x1 = 521.113
x2 = -2.918
x3 = 218.009
x4 = 4.124
Qsim = f_simula_rapido(x1, x2, x3, x4, PME, ETP, x1*0.3, x3*0.5)
end = time.perf_counter()
minutes = (end-start)/60


# Teste de calibracao


# Teste de convergencia do periodo de aquecimento


# Teste da funcao de simulacao
# start_time = dt.datetime.now()
# for i in range(5):
#     Parametros = {'x1':X[0], 'x2':X[1], 'x3':X[2], 'x4':X[3]}
#     S = 0.3 * X[0]
#     R = 0.5 * X[2]
#     OrdHU1 = f_gera_OrdHU1(X[3])
#     OrdHU2 = f_gera_OrdHU2(X[3])
#     HU1 = np.zeros(len(OrdHU1))
#     HU2 = np.zeros(len(OrdHU2))
#     VarsEstado = {'S': S, 'R': R, 'HU1':HU1, 'HU2':HU2}
#     DF = f_simula_GR4H(Parametros, Forcantes, VarsEstado, OrdHU1, OrdHU2)
# end_time = dt.datetime.now()
# print('Duracao: {}'.format(end_time-start_time))


# Teste da funcao de calibracao
# start_time = dt.datetime.now()
# for i in range(5):
#     Parametros = X_otim
#     Qsim = f_calibra_GR4H(Parametros, Forcantes.pme, Forcantes.etp, Forcantes.qobs, 1440)
# end_time = dt.datetime.now()
# print('Duracao: {}'.format(end_time-start_time))


# # Para efetivar uma calibracao...
# args = (Forcantes.pme, Forcantes.etp, Forcantes.qobs, 1440)
# bounds = [(50,850), (-10,10), (100,600), (1,10)]
# resultado = differential_evolution(f_calibra_GR4H, bounds, updating='deferred', workers=-1, args=args, disp=True)


# Testando os periodos de aquecimento
# dados = pd.read_csv('L0123003.csv', parse_dates=True, index_col="DatesR")
# Forcantes = dados[['P','E','Qmm']].astype('float')
# Forcantes = Forcantes.rename(columns={'P':'pme', 'E':'etp', 'Qmm':'qobs'})
# Forcantes.index.rename('datahora_utc', inplace=True)
# Parametros = {'x1':560, 'x2':-2.5, 'x3':214.2, 'x4':4.11}
# OrdHU1 = f_gera_OrdHU1(4.11)
# OrdHU2 = f_gera_OrdHU2(4.11)
# HU1 = np.zeros(len(OrdHU1))
# HU2 = np.zeros(len(OrdHU2))
# S = np.arange(0,1.25,0.25)
# R = np.arange(0,1.25,0.25)
# for i in S:
#     for j in R:
#         S = i*560.0
#         R = j*214.2
#         VarsEstado = {'S': S, 'R': R, 'HU1':HU1, 'HU2':HU2}
#         print('Simulacao i={}, j={}'.format(i,j))
#         DF = f_simula_GR4H(Parametros, Forcantes, VarsEstado, OrdHU1, OrdHU2)
#         DF.Qsim.plot(label='i={}, j={}'.format(i,j))
#
#
#
#         break
#     break
#
# plt.show()
