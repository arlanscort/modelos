### ARQUIVO PARA A EXECUÇÃO DO OTIMIZADOR MOCOM PARA A CALIBRAÇÃO DOS MODELOS HIDROLÓGICOS DO SIMEPAR ###


from numpy import array
from datetime import datetime
import matplotlib.pyplot as plt
from lib.otimizacao import MOCOM

import csv
import os
import sys


def ReadEPQ(arqname):
    """
        Função de leitura do arquivo EPQ
        Salva no dicionário 'dados' para operar nos modelos de otimização
    """
    ETp, CMB, Qmont, Qexut, peso, data = [], [], [], [], [], []
    dados = { 'data':[], 'ETp':[], 'CMB':[], 'Qmont':[], 'Qexut':[], 'peso':[] }

    with open(arqname,'rt') as f:
        for line in f:
            if line.startswith("#"):
                pass
            else:
                l = line.split(';')
                dados['data'].append( datetime.strptime(l[0],'%Y-%m-%d'))
                dados['ETp'].append( float(l[1]) )
                dados['CMB'].append( float(l[2]) )
                dados['Qmont'].append( float(l[3]) )
                dados['Qexut'].append( float(l[4]) )
                dados['peso'].append( int(float(l[5]) ))

    for p in dados.keys():
        dados[p] = array(dados[p])

    return dados


### SEGUIR OS PASSOS PARA EXECUTAR A CALIBRAÇÃO ###

#1. Definir o nome do arquivo EPQ
nome_epq = "maringa_diario"

#*. Definir a vazão inicial da propagação
Q0 = 25.0 #m³/s

#2. Definir a área da bacia em km²
area = 1238

#3. Definir o passo de tempo EM DIAS (se DT = 1. o passo é de 1 dia; se DT = 1./24. o passo é de 1 hora)
DT = 1.

#4. Definir o diretorio do arquivo de saída que irá conter a frente Pareto
arq_s = "resultado_maringa_diario.txt"

dir_epq = "epq_" + nome_epq + ".txt"
dados = ReadEPQ(dir_epq)
dados['area'] = area
dados['nome'] = 'maringa'
dados['DT'] = DT
dados['Q0'] = Q0

#5. Montar o objeto da classe MOCOM com a chave do modelo, o dicionário de dados, o arquivo de saída, as funções objetivo e o diretório dos arquivos dos ciclos (salvos de 200 em 200)
otimiza = MOCOM('sacsma-3s', dados, arqFinal= arq_s, fObj=['NS','NS_log','erro_vol'], dirCiclos="CiclosMOCOM/")

#6. Execução da otimização
otimiza.Executa()

#7. Plotagem das funções objetivo da frente Pareto
otimiza.PlotaFobj( FilePopInicial = arq_s, funcoes = ['NS','NS_log','erro_vol'], p3D = True)

#8. Simulação com base na solução que deu origem ao melhor resultado de cada função objetivo contida em "efes"
efes = ['KGE']
qsim,fsim,result = otimiza.SimulaMelhor( FilePopInicial = arq_s, funcoes_best = efes, funcoes_calc=efes )

q_sim = result

with open('qsim.csv', 'w') as arq:
    for i in q_sim:
        arq.write(str(i))
        arq.write('\n')
