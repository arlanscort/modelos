from spotpy.parameter import Uniform
import funcoes_objetivo


class setup_gr4j_nash(object):

    x1 = Uniform(low =   1, high = 1500, optguess =  350)
    x2 = Uniform(low = -10, high =    5, optguess =    0)
    x3 = Uniform(low =   1, high =  500, optguess =   90)
    x4 = Uniform(low = 0.5, high =    4, optguess =  1.7)
    k  = Uniform(low = 0.1, high =    5, optguess =  1.0)
    n  = Uniform(low = 0.1, high =   10, optguess = 0.15)

    def __init__(self, area, PME, ETP, Qjus, Qmon=None, h_aq=0, fobj='NSE'):
        self.area = area
        self.PME  = PME
        self.ETP  = ETP
        self.Qjus = Qjus
        self.Qmon = Qmon
        self.h_aq = h_aq
        self.fobj = fobj

    def simulation(self, x):
        from gr4j import gr4j_nash
        Qsim = gr4j_nash(self.area, self.PME, self.ETP, x[0], x[1], x[2], x[3], x[4], x[5], Qmon=self.Qmon)
        return Qsim

    def evaluation(self):
        Qobs = self.Qjus
        return Qobs

    def objectivefunction(self, simulation, evaluation):
        criterio = getattr(funcoes_objetivo, self.fobj)(simulation, evaluation, self.h_aq)
        fmin = 1 - criterio
        return fmin


class setup_gr4j_muskingum(object):

    x1 = Uniform(low =    1, high = 1500, optguess =  350)
    x2 = Uniform(low =  -10, high =    5, optguess =    0)
    x3 = Uniform(low =    1, high =  500, optguess =   90)
    x4 = Uniform(low =  0.5, high =    4, optguess =  1.7)
    k  = Uniform(low =  0.5, high =    7, optguess =  1.0)
    x  = Uniform(low = 0.01, high =  0.5, optguess = 0.15)

    def __init__(self, area, PME, ETP, Qjus, Qmon=None, h_aq=0, fobj='NSE'):
        self.area = area
        self.PME  = PME
        self.ETP  = ETP
        self.Qjus = Qjus
        self.Qmon = Qmon
        self.h_aq = h_aq
        self.fobj = fobj

    def simulation(self, x):
        from gr4j import gr4j_muskingum
        Qsim = gr4j_muskingum(self.area, self.PME, self.ETP, x[0], x[1], x[2], x[3], x[4], x[5], Qmon=self.Qmon)
        return Qsim

    def evaluation(self):
        Qobs = self.Qjus
        return Qobs

    def objectivefunction(self, simulation, evaluation):
        criterio = getattr(funcoes_objetivo, self.fobj)(simulation, evaluation, self.h_aq)
        fmin = 1 - criterio
        return fmin


class setup_gr4j_prop_interna(object):

    x1 = Uniform(low =   1, high = 1500, optguess =  350)
    x2 = Uniform(low = -10, high =    5, optguess =    0)
    x3 = Uniform(low =   1, high =  500, optguess =   90)
    x4 = Uniform(low = 0.5, high =    4, optguess =  1.7)

    def __init__(self, area, PME, ETP, Qmon, Qjus, h_aq=0, fobj='NSE'):
        self.area = area
        self.PME  = PME
        self.ETP  = ETP
        self.Qjus = Qjus
        self.Qmon = Qmon
        self.h_aq = h_aq
        self.fobj = fobj

    def simulation(self, x):
        from gr4j import gr4j_prop_interna
        Qsim = gr4j_prop_interna(self.area, self.PME, self.ETP, self.Qmon, x[0], x[1], x[2], x[3])
        return Qsim

    def evaluation(self):
        Qobs = self.Qjus
        return Qobs

    def objectivefunction(self, simulation, evaluation):
        criterio = getattr(funcoes_objetivo, self.fobj)(simulation, evaluation, self.h_aq)
        fmin = 1 - criterio
        return fmin


class setup_hymod(object):

    cmax  = Uniform(low=1.0   , high=500)
    bexp  = Uniform(low=0.1   , high=2.0)
    alpha = Uniform(low=0.1   , high=0.99)
    Ks    = Uniform(low=0.001 , high=0.10)
    Kq    = Uniform(low=0.1   , high=0.99)
    k     = Uniform(low=0.5   , high=7)
    x     = Uniform(low=0.01  , high=0.5)

    def __init__(self, area, PME, ETP, Qjus, Qmon=None, h_aq=0, fobj='KGE'):
        self.area = area
        self.PME  = PME
        self.ETP  = ETP
        self.Qjus = Qjus
        self.Qmon = Qmon
        self.h_aq = h_aq
        self.fobj = fobj

    def simulation(self, x):
        from hymod import hymod_nash
        Qsim = hymod_nash(self.area, self.PME, self.ETP, x[0], x[1], x[2], x[3], x[4], x[5], x[6], Qmon=self.Qmon)
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
