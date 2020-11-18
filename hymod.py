# -*- coding: utf-8 -*-
'''
Copyright (c) 2015 by Tobias Houska

This file is part of Statistical Parameter Estimation Tool (SPOTPY).

:author: Tobias Houska and Benjamin Manns

:paper: Houska, T., Kraft, P., Chamorro-Chavez, A. and Breuer, L.:
SPOTting Model Parameters Using a Ready-Made Python Package,
PLoS ONE, 10(12), e0145180, doi:10.1371/journal.pone.0145180, 2015.
'''

# from numba import jit

import propagacao # INTERVENCAO ARLAN
import funcoes_objetivo
from spotpy.parameter import Uniform

# def hymod(Precip, PET, cmax, bexp, alpha, Rs, Rq):
def hymod_nash(area, Precip, PET, cmax, bexp, alpha, Rs, Rq, k, n, Qmon=None): # INTERVENCAO ARLAN

    """
    See https://www.proc-iahs.net/368/180/2015/piahs-368-180-2015.pdf for a scientific paper:

    Quan, Z.; Teng, J.; Sun, W.; Cheng, T. & Zhang, J. (2015): Evaluation of the HYMOD model
    for rainfall–runoff simulation using the GLUE method. Remote Sensing and GIS for Hydrology
    and Water Resources, 180 - 185, IAHS Publ. 368. DOI: 10.5194/piahs-368-180-2015.

    :param cmax:
    :param bexp:
    :param alpha:
    :param Rs:
    :param Rq:
    :return: Dataset of water in hymod (has to be calculated in litres) (???)
    :rtype: list
    """

    # HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
    x_loss = 0.0
    # Initialize slow tank state
    x_slow = 2.3503 / (Rs * 22.5)
    x_slow = 0  # --> works ok if calibration data starts with low discharge
    # Initialize state(s) of quick tank(s)
    x_quick = [0,0,0]
    t = 0
    output = []
    # START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS

    while t <= len(Precip)-1:
        Pval = Precip[t]
        PETval = PET[t]
        # Compute excess precipitation and evaporation
        ER1, ER2, x_loss = excess(x_loss, cmax, bexp, Pval, PETval)
        # Calculate total effective rainfall
        ET = ER1 + ER2
        #  Now partition ER between quick and slow flow reservoirs
        UQ = alpha * ET
        US = (1 - alpha) * ET
        # Route slow flow component with single linear reservoir
        x_slow, QS = linres(x_slow, US, Rs)
        # Route quick flow component with linear reservoirs
        inflow = UQ

        for i in range(3):
            # Linear reservoir
            x_quick[i], outflow = linres(x_quick[i], inflow, Rq)
            inflow = outflow

        # Compute total flow for timestep
        output.append(QS + outflow)
        t = t+1

    # INTERVENCAO ARLAN
    # Propagacao
    import numpy as np
    output = np.array(output)
    if Qmon != None:
        Qprop = propagacao.nash(Qmon, k, n)
        output += Qprop
    # Conversao mm/dia -> m3/s
    output = output * area/86.4

    # INTERVENCAO ARLAN
    return output

# @jit
def power(X,Y):
    X=abs(X) # Needed to capture invalid overflow with netgative values
    return X**Y

# @jit
def linres(x_slow,inflow,Rs):
    # Linear reservoir
    x_slow = (1 - Rs) * x_slow + (1 - Rs) * inflow
    outflow = (Rs / (1 - Rs)) * x_slow
    return x_slow,outflow

# @jit
def excess(x_loss,cmax,bexp,Pval,PETval):
    # this function calculates excess precipitation and evaporation
    xn_prev = x_loss
    ct_prev = cmax * (1 - power((1 - ((bexp + 1) * (xn_prev) / cmax)), (1 / (bexp + 1))))
    # Calculate Effective rainfall 1
    ER1 = max((Pval - cmax + ct_prev), 0.0)
    Pval = Pval - ER1
    dummy = min(((ct_prev + Pval) / cmax), 1)
    xn = (cmax / (bexp + 1)) * (1 - power((1 - dummy), (bexp + 1)))

    # Calculate Effective rainfall 2
    ER2 = max(Pval - (xn - xn_prev), 0)

    # Alternative approach
    evap = (1 - (((cmax / (bexp + 1)) - xn) / (cmax / (bexp + 1)))) * PETval  # actual ET is linearly related to the soil moisture state
    xn = max(xn - evap, 0)  # update state

    return ER1,ER2,xn


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
        Qsim = hymod_nash(self.area, self.PME, self.ETP, x[0], x[1], x[2], x[3], x[4], x[5], x[6], Qmon=self.Qmon)
        return Qsim

    def evaluation(self):
        Qobs = self.Qjus
        return Qobs

    def objectivefunction(self, simulation, evaluation):
        criterio = getattr(funcoes_objetivo, self.fobj)(simulation, evaluation, self.h_aq)
        fmin = 1 - criterio
        return fmin
