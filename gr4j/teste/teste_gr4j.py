import sys
sys.path.append('/Users/arlan/github/modelos/GR4J/')

import pandas as pd
import gr4j
df = pd.read_excel('GR4J_EN_ORIGINAL.xlsx', sheet_name='GR4J', skiprows=38)
idx = df['Date']
PME = df['Pluie (mm)'].to_numpy()
ETP = df['ETP (mm)'].to_numpy()
Qobs_mm = df['Débit (mm/j)'].to_numpy()
Qobs = df['Débit (m3/s)'].to_numpy()

area = 260
x1 = 320.11
x2 = 2.42
x3 = 69.63
x4 = 1.39
# fconv
Qsim = gr4j.sim(PME, ETP, area, x1, x2, x3, x4)

import matplotlib.pyplot as plt
plt.plot(idx, Qobs)
plt.plot(idx, Qsim)
plt.show()
