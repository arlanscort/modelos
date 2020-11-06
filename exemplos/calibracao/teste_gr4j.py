'''
Teste do GR4J
Dataset disponivel em:
    https://webgr.inrae.fr/en/models/daily-hydrological-model-gr4j/
Arlan Scortegagna - out/2020
'''

import sys
sys.path.append('/Users/arlan/github/modelos/')
import pandas as pd
import numpy as np
from plotar_hidro import plotar_hidro

# Leitura das forcantes
df = pd.read_excel('GR4J_EN_ORIGINAL.xlsx', sheet_name='GR4J', skiprows=38)
idx = df['Date'].to_numpy()
PME = df['Pluie (mm)'].to_numpy()
ETP = df['ETP (mm)'].to_numpy()
Qjus = df['Débit (m3/s)'].to_numpy()
area = 260

# Calibracao
from spotpy.analyser import *
from spotpy.algorithms import sceua
from spot_setup import setup_gr4j_nash
from plotar_hidro import plotar_hidro
spot_setup = setup_gr4j_nash(area, PME, ETP, Qjus, h_aq=365, fobj='NSE')
sampler = sceua(spot_setup)
sampler.sample(5000, ngs=10, kstop=3, peps=0.1, pcento=0.1)
results = sampler.getdata()
params = get_best_parameterset(results,maximize=False)
bestindex, bestobjf = get_minlikeindex(results)
simulation_fields = get_simulation_fields(results)
Qsim = {}
Qsim['best_NSE'] = list(results[simulation_fields][bestindex][0])

# Plotagem
fig = plotar_hidro(idx, PME, ETP, Qjus, Qmon=None, Qsim=Qsim)
import matplotlib.pyplot as plt
plt.show()
