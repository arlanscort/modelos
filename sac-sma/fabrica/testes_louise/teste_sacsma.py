'''

Teste para aplicacao do modelo_sacsma.py
Bacia teste: atibaia_valinhos

'''
import pandas as pd
import modelo_sacsma
from modelo_sacsma import sac_sma_detalhado
from matplotlib import pyplot as plt

############################################################################
# IMPORTACAO DOS PARAMETROS E DAS FORCANTES
############################################################################

parametros = pd.read_excel('parametros.xlsx', index_col = 'parametro')

#print(parametros['valor']['area'])

forcantes = pd.read_excel('forcantes.xlsx')
data = forcantes['data']
PME = forcantes["pme"]
ETP = forcantes["etp"]
QIN1 = forcantes["qin1"]
QOBS = forcantes["qobs"]

############################################################################
# MODELO
############################################################################

df = sac_sma_detalhado(parametros, PME, ETP)
df['qsim'] = df['ROIMP+SDRO'] + df['SSUR'] + df['SIF'] + df ['BFS'] + df['BFP']
df['qobs'] = QOBS
df['data'] = data


plt.figure(figsize=(10, 5))
plt.plot(df['data'], df['qsim'], color='blue', linestyle='-', linewidth=1, label = 'qsim')
plt.plot(df['data'], df['qobs'], 'k--', linewidth=2, label='qobs')
plt.xlabel('data')
plt.ylabel('vazao')
plt.legend(loc='best', fontsize=15)
plt.grid();
plt.savefig('valinhos.png')
plt.show()


#df['qsim'].plot(label='qsim')
#df['qobs'].plot(label='qobs')
#plt.legend()
#plt.show()
#fig.savefig('valinhos.png')

df.to_excel('resultados_valinhos.xlsx')
