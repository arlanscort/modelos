# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 10:25:24 2020
Propagacao do escoamento superficial (referencia: ARNO)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as pyplot

#==============================================================================
# FASE propagação do escoamento supercial com 3 reservatórios lineares
#==============================================================================

Forcantes = pd.read_excel('forcantes2.xlsx')
#Forcantes = pd.read_excel('forcantes.xlsx', index_col = 'data', parse_dates=True)
#Parametros = pd.read_excel('parametros.xlsx', index_col = 'parametro')

# numero de reservatorios
n_rsv = 3
k = 0.75
dt = 1
#k = Parametros['valor'].get(['k'], 0.5)
#I = Forcantes['qin1']

I = np.array([i for i in Forcantes['qin1']])
print(I)
qobs = np.array([i for i in Forcantes['qobs']])
print(qobs)
qin = I

qout = np.zeros(len(I))
qout[0] = qin[0]

for n in range(n_rsv):
    i = 0

    while i+1 < len(I):
        print(i)
        qout[i+1] = (2*k - dt)/(2*k + dt)*qout[i] + (dt)/(2*k + dt) * (qin[i+1] + qin[i])
        print(qout[i+1], qout[i],qin[i+1], qin[i])
        print('tempo')

        i += 1

    qin = [x for x in qout]
    print(qin)
    print(qout)
    print('reservatorio')


pyplot.figure(figsize=(10,5))
pyplot.plot(range(len(I)), I, label= 'I')
pyplot.plot(range(len(I)), qout,label= 'qout')
pyplot.plot(range(len(I)), qobs, color = 'black',label= 'qobs')
pyplot.legend()
pyplot.show()
