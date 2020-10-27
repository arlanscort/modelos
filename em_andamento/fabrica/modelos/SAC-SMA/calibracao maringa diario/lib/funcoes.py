from datetime import datetime
from numpy import zeros, array, ones, mean, std, sqrt, nonzero, log,percentile,arange
from scipy.stats import pearsonr

# Definicão de funcoes estatísticas -------------------------------------------------
def Pearson(s, o):
    return pearsonr(s, o)[0]


def NS(s, o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    obs: dados com peso zero devem ser excluídos antes
    """
    ns = 1 - sum( (s-o)**2 )/sum((o - mean(o) )**2)
    return ns

def NS_log(s, o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    obs: dados com peso zero devem ser excluídos antes
    """
    tiny=0.0001
    s = log(s+tiny)
    o = log(o+tiny)
    ns = 1 - sum( (s-o)**2 )/sum((o - mean(o) )**2)
    return ns

def KGE(s, o,sr=1.,salfa=1.,sbeta=1.):
    """
    Kling-Gupta Efficiency
    input:
        s: simulated
        o: observed
    """
    x1_std, x1_avg = std(o), mean(o)
    x2_std, x2_avg = std(s), mean(s)
  
    alfa = x2_std/x1_std               #variabilidade relativa entre simulado(2) e observado(1)
    #beta = (x2_avg - x1_avg)/x1_std   #bias normalizado pelo desvio padrao da observacao <-termo beta em nse..    
    beta = (x2_avg/x1_avg)             #bias entre vazao simuladao e observada
    
    r = pearsonr(s, o)[0]
    
    #scaling factors    
    sr, salfa, sbeta = sr**2, salfa**2, sbeta**2
    
    #ponderadores
    g1 = ( r - 1 )*( r - 1 ) 
    g2 = ( alfa - 1 )*( alfa - 1 )
    g3 = ( beta - 1 )*( beta - 1 )
    eds  = sqrt( abs(sr*g1) + abs(salfa*g2) + abs(sbeta*g3) ) #formula geral
    kge  = 1 - eds
    
    return kge

def Bias(s, o):
    """
    Retorna valor relativo do Bias (vies)
    """
    
    #bias = abs( sum( s - o )/sum(o) )
    bias = abs(mean(s-o))
    
    return bias
  

def residuo(s, o):
    """
    Adimensional Coefficient
    Entrada
        s: simulated
        o: observed
    Retorna
        f2 = adimensional coeffecient
    """
    p1 = abs( mean(s - o) )/mean(o)
    p2 = std( s - o)/std(o)
    f2 = 0.3*p1 + 0.7*p2
    return f2

def erro_vol(s, o):
    """
    Relação entre os volumes calculados
    Entrada
        s: simulated
        o: observed
    Retorna
        deltaV = relação entre os volumes
    """
    deltaV = abs( ( s.sum() - o.sum() )/ o.sum() )
    
    return deltaV

def rmse(s, o):
    """
    Erro quadrático médio
    Entrada
        s: simulated
        o: observed
    Retorna
        deltaV = relação entre os volumes
    """

    return sqrt( ( (o - s)**2 ).mean() )

def rmsei(s, o):
    """
    Erro quadrático médio inverso
    Entrada
        s: simulated
        o: observed
    Retorna
        deltaV = relação entre os volumes
    """
    tiny = 0.0001
    s = s+tiny
    o = o+tiny
    return sqrt( ( (1./o - 1./s)**2 ).mean() )

def rmseqq(s,o):
    """
    Erro médio quadratico dos quantis
    Entrada
        s: simulated
        o: observed
    Retorna
        deltaV = relação entre os volumes
    """
    lin = arange(5,100,5)
    s = percentile(s,lin)
    o = percentile(o,lin)    

    return sqrt( ( (o - s)**2 ).mean() )    


# Dicionarios auxiliares ---------------------------------------------------------------------------------------------------------------------

Objetivo = {        'NS': [NS, -1],        'NS_log': [NS_log, -1],     'KGE': [ KGE, -1],       'Bias': [Bias,   1], 'Pearson': [Pearson, -1],
              'erro_vol': [erro_vol, 1],     'rmse': [rmse, 1],      'rmsei': [rmsei, 1],    'residuo': [residuo,1], 'rmseqq':[rmseqq,1]}

estatisticas = ['NS', 'NS_log', 'Bias', 'Pearson', 'rmse', 'rmsei', 'erro_vol', 'Erro_pos', 'Erro_neg', 'Erro_total', 'residuo', 'KGE']
