'''
Teste para aplicacao do modelo_sacsma.py
Bacia teste: atibaia_valinhos
'''
import pandas as pd
import numpy as np
import modelo_sacsma
from modelo_sacsma import sac_sma_detalhado
from matplotlib import pyplot as plt

############################################################################
# IMPORTACAO DOS PARAMETROS E DAS FORCANTES
############################################################################
Forcantes = pd.read_excel('forcantes.xlsx', index_col = 'data', parse_dates=True)
Parametros = pd.read_excel('parametros.xlsx', index_col = 'parametro')
Estados = None
DF = sac_sma_detalhado(Forcantes, Parametros)
DF.to_excel('teste_valinhos.xlsx')
