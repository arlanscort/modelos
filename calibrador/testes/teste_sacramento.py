import pandas as pd
from modelo_sacsma import sac_sma_d_detalhado



Forcantes = pd.read_excel('forcantes.xlsx', index_col = 'data', parse_dates=True)
Parametros = pd.read_excel('parametros.xlsx', index_col = 'parametro')
area = 1000 # km2
DF = sac_sma_d_detalhado(area, Forcantes, Parametros)

qobs = Forcantes['qobs']
qsim = DF['BFCC'] + DF['QUZ_prop']

import matplotlib.pyplot as plt
qobs.plot(color='black')
qsim.plot(color='red')
plt.show()
