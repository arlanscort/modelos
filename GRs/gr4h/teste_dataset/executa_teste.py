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

# Calibra
#X = [521.113, -2.918, 218.009, 4.124]
#Qsim, NSE = f_calibra_GR4H(X, Forcantes.pme, Forcantes.etp, Forcantes.qobs, 1440)


# bounds = [(100,850), (-2.5,2.5), (100,600), (1,10)]
# args = (Forcantes.pme, Forcantes.etp, Forcantes.qobs, 1440)
#
# resultado = differential_evolution(f_calibra_GR4H, bounds, args=args, disp=True)
#
#
# from scipy.optimize import rosen, differential_evolution
#
# bounds = [(0,2), (0, 2), (0, 2), (0, 2), (0, 2)]
# result = differential_evolution(rosen, bounds, popsize=150, disp=True)
