from spotpy.parameter import Uniform
import modelos.funcoes_objetivo as funcoes_objetivo
import modelos.gr4j as gr4j

class gr4j_muskingum(object):

    x1 = Uniform(low =   1, high = 1500, optguess = 350)
    x2 = Uniform(low = -10, high =    5, optguess =   0)
    x3 = Uniform(low =   1, high =  500, optguess =  90)
    x4 = Uniform(low = 0.5, high =    4, optguess = 1.7)
    k  = Uniform(low = 0.5, high =  10,  optguess = 1.0)
    x  = Uniform(low = 0.0, high = 0.5,  optguess = 0.15)

    def __init__(self, area, PME, ETP, Qjus, Qmon=None, h_aq=0, fobj='NSE'):
        self.area = area
        self.PME  = PME
        self.ETP  = ETP
        self.Qjus = Qjus
        self.Qmon = Qmon
        self.h_aq = h_aq
        self.fobj = fobj

    def simulation(self, x):
        Qsim = gr4j.sim_muskingum(self.area, self.PME, self.ETP, \
                x[0], x[1], x[2], x[3], x[4], x[5], \
                Qmon=self.Qmon)
        return Qsim

    def evaluation(self):
        Qobs = self.Qjus
        return Qobs

    def objectivefunction(self, simulation, evaluation):
        criterio = getattr(funcoes_objetivo, self.fobj)(simulation, evaluation, self.h_aq)
        fmin = 1 - criterio
        return fmin


# class sacramento
# class ...
