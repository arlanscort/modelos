import numpy as np
import pandas as pd
import datetime as dt

from modelo_gr4h import  f_gera_OrdHU1, f_gera_OrdHU2, f_simula_GR4H

dados = pd.read_csv('L0123003.csv', parse_dates=True, index_col="DatesR")
Forcantes = dados[['P','E','Qmm']].astype('float')
Forcantes = Forcantes.rename(columns={'P':'pme', 'E':'etp', 'Qmm':'Qobs'})
Forcantes.index.rename('datahora_utc', inplace=True)
Parametros = {'x1':521.113, 'x2':-2.918, 'x3':218.009, 'x4':4.124}
OrdHU1 = f_gera_OrdHU1(Parametros['x4'])
OrdHU2 = f_gera_OrdHU2(Parametros['x4'])
HU1 = np.zeros(len(OrdHU1))
HU2 = np.zeros(len(OrdHU2))
S = 0.3*Parametros['x1']
R = 0.5*Parametros['x3']
VarsEstado = {'S':S, 'R':R, 'HU1':HU1, 'HU2':HU2}

DF = f_simula_GR4H(Parametros, Forcantes, VarsEstado, OrdHU1, OrdHU2)


teste


# import matplotlib.pyplot as plt


# t0 = dt.datetime(2004,3,1,1)
#
# Qobs = Forcantes.loc[t0:,'Qobs'].to_numpy()
# Qobtido = DF.loc[t0:,'Qsim'].to_numpy()
# Qesperado = pd.read_csv('Qsaidas.csv').to_numpy()
#
# plt.plot(Qobtido, label='obtido', color='red')
# plt.plot(Qesperado, label='esperado', color='blue')
#
# plt.show()
