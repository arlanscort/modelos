'''
Implementacao - Arlan Scortegagna, set/2020
Revisao - Gabriel Bernardini, out/2020
Programa contendo a classe spot_setup para execucao dos modulos do spotpy
'''

import sys
sys.path.append('/Users/arlan/github/modelos/GR4J/')
from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import os
import pandas as pd
import gr4j
import numpy as np

class spot_setup(object):
    x1 = Uniform(low=  1, high=1500, optguess=500)
    x2 = Uniform(low=-10, high=5   , optguess=1  )
    x3 = Uniform(low=1  , high=500 , optguess=250)
    x4 = Uniform(low=0.4, high=4   , optguess=0.5)

    def __init__(self, bacia_nome, bacia_area, obj_func):
        # Definir funcao objetivo
        self.obj_func = obj_func
        # Definir bacia
        self.area = bacia_area
        # Carrega o arquivo PEQ (deve estar na mesma pasta do arquivo de exec)
        self.ETP, self.PME = [], []
        self.data, self.Q_obs = [], []

        path = os.path.dirname(os.path.realpath(__file__))
        PEQ = pd.read_csv('{}.csv'.format(bacia_nome), index_col='data', parse_dates=True)
        self.data  = PEQ.index
        self.pme   = PEQ.pme.to_numpy()
        self.etp   = PEQ.etp.to_numpy()
        self.qin   = PEQ.qin.to_numpy()
        self.qobs  = PEQ.qobs.to_numpy()

    def simulation(self, x):
        # Executa a simulacao com o modelo hidrologico retornando a vazao em mm
        Qsim = gr4j.sim(self.pme, self.etp, self.area, x[0], x[1], x[2], x[3])
        return Qsim

    def evaluation(self):
        Qobs = self.qobs
        return Qobs

    def objectivefunction(self, simulation, evaluation, params=None):
        if params is None:
            params = {}
        wup = Estados.get('wup', 0) # warm-up period

        NSE = (simulation-evaluation)**2 / (evaluation-np.mean(evaluation))**2
        return f
