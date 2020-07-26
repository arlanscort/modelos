import numpy as np
import pandas as pd
import datetime as dt

from modelo_gr4h import f_gera_OrdHU1, f_gera_OrdHU2
from modelo_gr4h import f_simula_GR4H, f_calibra_GR4H
from scipy.optimize import differential_evolution


dados = pd.read_csv('L0123003.csv', parse_dates=True, index_col="DatesR")
Forcantes = dados[['P','E','Qmm']].astype('float')
Forcantes = Forcantes.rename(columns={'P':'pme', 'E':'etp', 'Qmm':'qobs'})
Forcantes.index.rename('datahora_utc', inplace=True)

X_orig = [521.113, -2.918, 218.009, 4.124]
X_otim = [559.60081023, -2.5, 214.18937429, 4.11512303]

import datetime as dt

X = X_otim

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


# Para efetivar uma calibracao...
bounds = [(50,850), (-10,10), (100,600), (1,10)]
resultado = differential_evolution(f_calibra_GR4H, bounds, updating='deferred', workers=-1, args=args, disp=True)
