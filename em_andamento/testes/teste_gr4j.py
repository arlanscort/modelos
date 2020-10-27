'''
Teste do GR4J
Dataset disponivel em:
    https://webgr.inrae.fr/en/models/daily-hydrological-model-gr4j/
Arlan Scortegagna - out/2020
'''

import sys
sys.path.append('/Users/arlan/github/modelos/')
import gr4j
import pandas as pd
import spotpy
import numpy as np


# Leitura das forcantes
df = pd.read_excel('GR4J_EN_ORIGINAL.xlsx', sheet_name='GR4J', skiprows=38)
idx = df['Date'].to_numpy()
PME = df['Pluie (mm)'].to_numpy()
ETP = df['ETP (mm)'].to_numpy()
Qjus = df['Débit (m3/s)'].to_numpy()
area = 260


# Teste de verificacao da implementacao do calibrador

import matplotlib.pyplot as plt

def calibra_e_plota(fobj):
    spot_setup = gr4j.spot_setup(area, PME, ETP, Qjus, h_aq=365, fobj=fobj)
    sampler = spotpy.algorithms.sceua(spot_setup, dbname='SCEUA_gr4j_20out2020', dbformat='csv')
    sampler.sample(5000, ngs=6, kstop=3, peps=0.1, pcento=0.1)
    results = sampler.getdata()
    params = spotpy.analyser.get_best_parameterset(results,maximize=False)
    x1, x2, x3, x4 = params[0][0], params[0][1], params[0][2], params[0][3]
    Qsim = gr4j.sim(area, PME, ETP, x1, x2, x3, x4)
    return Qsim

fobj = 'NSE'
Qsim = calibra_e_plota(fobj)
plt.plot(idx, Qsim, label=fobj)

fobj = 'sqrtNSE'
Qsim = calibra_e_plota(fobj)
plt.plot(idx, Qsim, label=fobj)

fobj = 'logNSE'
Qsim = calibra_e_plota(fobj)
plt.plot(idx, Qsim, label=fobj)

fobj = 'KGE'
Qsim = calibra_e_plota(fobj)
plt.plot(idx, Qsim, label=fobj)

plt.plot(idx, Qjus, label='Qobs')

plt.legend()
plt.show()


# # Testar
# 1 - NSE = 0.0922999
# Minimal objective value: 0.0922999
# x1: 281.743
# x2: 2.85714
# x3: 100.75
# x4: 1.72153


# Testar
# 1 - sqrtNSE = 0.120787
# x1: 279.785
# x2: 2.01923
# x3: 66.2981
# x4: 1.40072


# # Testar
# 1 - logNSE = 0.140972
# x1: 251.096
# x2: 1.76681
# x3: 58.0471
# x4: 1.32648


# Testar
# 1 - KGE = 0.0534984
# x1: 201.504
# x2: 3.70369
# x3: 167.493
# x4: 1.1624
