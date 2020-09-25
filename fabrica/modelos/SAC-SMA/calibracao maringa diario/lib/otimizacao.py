from numpy import zeros, array, ones, argsort, mean, std, sqrt, nonzero, all, any, log
from numpy.random import permutation
from random import random
from multiprocessing import Pool
from functools import partial
from math import floor
from shutil import copyfile
import matplotlib.pyplot as plt
import time 
import os
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D


import pandas as pd
import time


#from modelos import Modelos
#from funcoes import estatisticas, Objetivo, Pearson, NS, NS_log, Bias, residuo, rmsei, rmse, erro_vol, KGE


from lib.modelos import Modelos
from lib.funcoes import estatisticas, Objetivo, Pearson, NS, NS_log, Bias, residuo, rmsei, rmse, erro_vol, KGE


# configura graficos
params = {'font.size': 18, 'font.family': 'serif', 'xtick.major.pad': '8'}
mpl.rcParams.update(params)
mpl.rcParams['figure.figsize'] = 16, 12

# num_workers = multiprocessing.cpu_count()  # number of cores . O numero maximo de processos e igual ao numero de cores!

lbl_obs  = 'vazão observada'
lbl_best = 'melhor_'
lbl_vazao = 'Vazão (m³/s)'

lbl_obs  = 'observation'
lbl_best = 'best_'
lbl_vazao = 'Discharge (m³/s)'

class MOCOM(object):
    """
    # Método de calibração Multi-Objective global optimization (MOCOM)
    # Beaseado no artigo Yapo, Gupta, Sorooshian (1998) e na dissertação
    # "A multiobjective global optimization algorithm with application
    #    to calibration of hydrologic models" --- Yapo (1996)
    # Evolui simplex usando processamento paralelo
    
    Parametros:
    
    --> modelo: string 
        Indica o modelo a ser calibrado.
            'sacsma': Sacramento
              'iph2': IPH2 simples
            'iph2rl': IPH2 com propagação por reservatórios lineares
              'smap': modelo SMAP (Soil Moisture Accounting Procedure)
              
    --> entrada: dicionário
        Pode conter quaisquer dados necessários para o modelo a ser calibrado.
        Toda e qualquer informação usada no modelo deverá ser acrescentada neste dicionário, uma vez que o modelo deverá receber
            apenas o dicionário de dados e um array com os parâmetros.
        Obrigatório: vazão observada deve ser um array nomeado 'Qexut' 
                     caso haja pesos (0 ou 1) deve-se incluir o array de pesos com o key 'pesos'. Ex.: entrada['peso'] = array([0,0,0,1,1,1,....1,1])
                     
    --> fObj: lista de strings, opcional
        Lista contendo funções objetivo usadas para calibração (2 ou mais). Default ['NS', 'NS_log', 'Bias']
                'NS': Coeficiente Nash Sutcliffe
            'NS_log': Log do Nash Sutcliffe
               'KGE': Eficiência de Kling-Gupta
              'Bias': Valor absoluto do Bias (vies)
           'Pearson': Coeficiente de correlação de Pearson
          'erro_vol': Relação entre os volumes calculados
              'rmse': Erro quadrático médio
             'rmsei': Erro quadrático médio inverso
           'residuo': média ponderada dos resíduos da média e do desvio-padrão
           
    --> Nciclos: inteiro, opcional
        Número máximo de ciclos caso o calibração não consiga convergir para a fronteira de pareto com todos os pontos não dominados
        O valor default é 20000
        
    --> popTotal: inteiro, opcional
        População total. O valor default é 700 pontos.
        
    --> FilePopInicial: string, opcional
        Nome de arquivo contendo os parâmetros da população inicial. 
        Deve conter uma linha com os nomes (strings) de cada parâmetro. O número de pontos será igual ao número da população total.
        Se não for passado, os pontos iniciais serão gerados por sorteio.
        
    --> parFixos: dicionário, opcional
        Dicionários com os parâmetros que serão mantidos fixos. Podem ser passados quantos parâmetros quiser, só precisa garantir que os  
        nomes dos parametros sejam os mesmos usados pelo modelo. Ex de entrada para SACSMA: parFixos = {'UZK': 0.20, ..., 'lag': 6.0}
        Caso não seja passado nenhum parâmetro, todos os parâmetros do modelo serão calibrados
            
    --> NotaPeso: bool, opcional
        Se True: indica que há um array com pesos (0 ou 1) dentro da estrutura de dados e calculará as estatísticas apenas para dados com peso 1. key do array deverá ser 'peso'.
        Se False: calcula a função objetivo para todos os dados
        
    --> dt: int, opcional
        Passo de tempo das simulação do modelo. Valor default é 3600 (para simulação horária)
        
    """



    def __init__(self, 
                 modelo,                             # string contendo nome do modelo 
                 entrada,                            # estrutura de dados para entrada no modelo escolhido (definido pelo usuario)
                 fObj = ['KGE', 'NS_log','erro_vol'], # lista com funções objetivo 
                 Nciclos = 20000,                    # numero maximo de ciclos 
                 popTotal = 500,                     # população inicial
                 FilePopInicial = None,              # nome do arquivo contendo parametros da pop. inicial (util para continuar calibração interrompida)
                 parFixos = {},                      # se quiser, pode-se entrar com um dicionário contendo parametros fixos ex: {'par1': valor, 'parn': valor}
                 NotaPeso = True,                    # array de pesos Ex.: entrada['peso'] = array([0,0,0,1,1,1,....1,1])
                 dt = 86400,                          # passo de tempo da simulação
                 arqFinal = None,
                 q0 = None,
                 dirCiclos = "CiclosMOCOM/" ):                 # nome do arquivo com resultados finais

        self.DadosEntrada = entrada                                 # dados de entrada do modelo
        self.DadosEntrada['dt'] = dt                                # passo de tempo do modelo
        self.Nciclos = Nciclos                                      # numero de ciclos no MOCOM
        #self.FunObj =       {'fc': [ Objetivo[fc] for fc in fObj ], # dicionário com informações das funções objetivo <- 'Objetivo' é um dicionario!!!
        #                      'nomes': fObj, 'nb': len(fObj)}
        self.FunObj = self.montaFobj(fObj)
        
        self.NotaPeso = NotaPeso
        
        self.dirCiclos = dirCiclos
        self.arqFinal = arqFinal

        # Caso os dados tenha peso, serão separados os índices com peso 1 para cálculo das funções objetivo
        if NotaPeso: 
            self.NotaPeso = NotaPeso
            self.IndicesPeso = nonzero(self.DadosEntrada["peso"])[0]
            self.DadosEntrada["Qexut_peso1"] = array([ self.DadosEntrada["Qexut"][j] for j in self.IndicesPeso ])
        
        self.Modelo = {}
        self.Modelo['nome'] = Modelos[modelo]['nome']   # nome do modelo 
        self.Modelo['fc']   = Modelos[modelo]['fc']     # copia do modelo hidrológico        
        self.DadosEntrada['ParFixos'] = {}              # dicionário de parâmetros fixos
        self.DadosEntrada['calibra'] = []               # lista com parâmetros a serem calibrados
        self.DadosEntrada['OrdemParFixos'] = []         # lista com nomes dos parâmetros fixos (para conservar a ordem)

        linf, lsup, fix, aux = [], [], {}, {}
        
        # Separando parametros fixos e parametros a serem calibrados
        for par in parFixos.keys():            
            self.DadosEntrada['ParFixos'][par] = parFixos[par]
        
        for j,par in enumerate(Modelos[modelo]['params']):
            if par in self.DadosEntrada['ParFixos']:
                self.DadosEntrada['OrdemParFixos'].append(aux[par])
            else:
                linf.append(Modelos[modelo]['Linf'][j])
                lsup.append(Modelos[modelo]['Lsup'][j])
                self.DadosEntrada['calibra'].append(par)

        # Definindo limites de acordo com parametros a serem calibrados
        self.Modelo['Linf'] = array(linf)       # array com limites inferiores dos parametros
        self.Modelo['Lsup'] = array(lsup)       # array com limites superiores dos parametros
        self.ndim = len(self.Modelo['Lsup'])    # dimensao (número de parâmetros a serem calibrados)
        
        # Caso entre com população inicial
        if FilePopInicial is not None: 
            self.Dx = MOCOM.CarregaPop0(self, FilePopInicial)  # Guarda parametros da pop. inicial
            self.popTotal   = len(self.Dx)                     # tamanho da população 
        else: 
            self.ParamsIniciais = None                          
            self.popTotal   = popTotal                         # tamanho da população (valor passado para a funcao)
        
        #~ MOCOM.Executa(self, arqFinal)  # Chama funcao para executar calibração
    

    def montaFobj(self,fObj):
        # dicionário com informações das funções objetivo <- 'Objetivo' é um dicionario!!!
        fObjDic = {'fc': [ Objetivo[fc] for fc in fObj ], 'nomes': fObj, 'nb': len(fObj)}
        return fObjDic


    def PlotaParam(self, FilePopInicial):
        """
         Função para plotar os parâmetros do modelo normalizados pelo valor máximo
           FilePopInicial: arquivo csv com dados da população
         
        """
        params  = MOCOM.CarregaPop0(self, FilePopInicial)
        params = array(params).T
        
        #~ plt.ion()
        fig = plt.figure() 
        axes = fig.add_axes([0.14, 0.15, 0.80, 0.75])    
        
        for i in range(len(params)):
            axes.plot( params[i]/max(params[i]), label = self.DadosEntrada['calibra'][i] )
        axes.grid()
        axes.legend(bbox_to_anchor=(0.5, -0.05), loc='upper center', prop={'size':16}, ncol = 8 )
        axes.set_ylabel("Parametros normalizados")
        plt.show()


    def PlotaFobj(self, FilePopInicial, funcoes = None, p3D = False):
        """
         Função para plotar as estatísticas finais
           funcoes: lista com string das estatísticas (2 ou 3) que se dejesa plotar. 
           p3D: Se true, e forem passadas três estatísticas, será plotado um gráfico em
                3D
        """
        if funcoes is None: funcoes = self.FunObj['nomes']
        STATS = MOCOM.CarregaFObj(self, arqname = FilePopInicial, FobjLista = funcoes)
        
        if len(funcoes) == 1: 
            print("Passar pelo menos duas estatísticas.")
            return 0
        
        if p3D:
            #~ plt.ion()
            if len(funcoes) != 3:
                print("Deve-se passar 3 estatísticas para gráfico 3D.")
                return 0
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            x, y, z = STATS[funcoes[0]], STATS[funcoes[1]], STATS[funcoes[2]]
            ax.scatter(x, y, z)
            ax.set_xlabel(funcoes[0])
            ax.set_ylabel(funcoes[1])
            ax.set_zlabel(funcoes[2])
            plt.show()
        else:
            #~ plt.ion()
            if len(funcoes) > 2:
                print("Foram passadas %d funcoes. Serão plotadas apenas as duas primeiras."%len(funcoes))
            fig = plt.figure()
            ax = fig.add_subplot(111)
            x, y = STATS[funcoes[0]], STATS[funcoes[1]]
            ax.scatter(x, y)
            ax.set_xlabel(funcoes[0])
            ax.set_ylabel(funcoes[1])
            ax.grid()
            plt.show()     


    def SimulaParam(self,params, ns=10, plota = True, funcoes_calc=None, permuta=False):
        """
         Função plota soluções do modelo, a partir de uma lista de parametros fornecida
           FilePopInicial: arquivo csv com dados da população
           ns: número de soluções a serem plotadas
        """
               
        if plota: 
            #~ plt.ion()
            fig = plt.figure()              
            axes = fig.add_axes([0.10, 0.10, 0.80, 0.80])    
            axes.plot(self.DadosEntrada['data'], self.DadosEntrada['Qexut'], color = 'black', label = lbl_obs)


        if funcoes_calc is None:
            ffobj = self.FunObj
        else:
            ffobj = self.montaFobj(funcoes_calc)

        print(ffobj['nomes'])
        
        lpars = range(len(params))
        
        #opcional: sorteia 'ns' pontos aleatoriamente
        if permuta==True:
            lpars = permutation(range(len(params)))
        lpars = lpars[:ns]


        xc,fxc = {},{}
        for i,param in enumerate(params):

            if i not in lpars: continue

            result = self.Modelo['fc']( self.DadosEntrada, param )  #chama modelo
            ofs = self.FObjPontoFunc(param,ffobj)
            
            print([ "  %.2f"%fi for fi in ofs] )

            #salva em dicionario
            xc[i]  = result
            fxc[i] = ofs
            
            if plota:
                axes.plot(self.DadosEntrada['data'], result, label = 'simul_%02d'%i)      
        
        if plota:
            axes.grid()
            axes.legend( )
            axes.set_ylabel(lbl_vazao)
            plt.show()

        return xc,fxc


    def Simula(self,FilePopInicial, ns=10 , plota = True, funcoes_calc=None, permuta=False,inflaction=None):
        """
         Função plota soluções do modelo
           FilePopInicial: arquivo csv com dados da população
           ns: número de soluções a serem plotadas
        """
        #retorna lista de arrays (um para cada conjunto de parametros)
        params  = MOCOM.CarregaPop0(self, FilePopInicial) 
        if inflaction is not None:
            params  = MOCOM.CarregaPop0Gauss(self, FilePopInicial, inflaction=inflaction) #
        
        if plota: 
            #~ plt.ion()
            fig = plt.figure()              
            axes = fig.add_axes([0.10, 0.10, 0.80, 0.80])    
            axes.plot(self.DadosEntrada['data'], self.DadosEntrada['Qexut'], color = 'black', label = lbl_obs)
        

       # print(self.FunObj['nomes'])

        if funcoes_calc is None:
            ffobj = self.FunObj
        else:
            ffobj = self.montaFobj(funcoes_calc)

        print(ffobj['nomes'])

        #sorteia 'ns' pontos aleatoriamente
        if permuta==True:
            lpars = permutation(range(len(params)))
        else:
            lpars = range(len(params))
        lpars = lpars[:ns]

        for i,param in enumerate(params):
            #if i > ns: break
            if i not in lpars: continue
            result = self.Modelo['fc']( self.DadosEntrada, param ) # chama modelo
            #ofs = MOCOM.FObjPonto(self, param )
            ofs = self.FObjPontoFunc(param,ffobj)
            print([ "  %.2f"%fi for fi in ofs] )
            if plota:
                axes.plot(self.DadosEntrada['data'], result, label = 'simul_%02d'%i)
        if plota:
            axes.grid()
            axes.legend( )
            axes.set_ylabel(lbl_vazao)
            plt.show()


    def SimulaMelhor(self, FilePopInicial, funcoes_best = None, funcoes_calc=None, plota = True, fig=None):
        """
         Função plota soluções do modelo
           FilePopInicial: arquivo csv com dados da população
           ns: número de soluções a serem plotadas
           funcoes_best: lista de strings com nome das funcoes objetivo, para as quais queremos as melhores solucoes
           funcoes_calc: lista de strings com nome das fucnoes objetivos, para as quais queremos ver resultados!
        """
        if funcoes_best is None: funcoes_best = self.FunObj['nomes'] #usa nome para buscar no arquivo

        if funcoes_calc is None:
            ffobj = self.FunObj
        else:
            ffobj = self.montaFobj(funcoes_calc)        
        
        #retorna lista de arrays (um para cada conjunto de parametros)
        #print('melhores solucoes salvas no arquivo:'+FilePopInicial)
        params = MOCOM.CarregaFObj(self, arqname = FilePopInicial, FobjLista = funcoes_best, retorna = 'parametros') 
		
        if plota: 
            #~ plt.ion()
            fig = plt.figure()              
            axes = fig.add_axes([0.10, 0.10, 0.80, 0.80])    
            axes.plot(self.DadosEntrada['data'], self.DadosEntrada['Qexut'], color = 'black', label = lbl_obs)
        print(ffobj['nomes'])
        xc,fxc = {},{}
        for fi in funcoes_best:
            
            param = params[fi]
            result = self.Modelo['fc']( self.DadosEntrada, param ) # chama modelo
            #ofs = self.FObjPonto(param) #depende do objeto
            ofs = self.FObjPontoFunc(param,ffobj) #independe do objeto
            print([ "  %.2f"%i for i in ofs] )

            #salva em dicionario
            xc[fi]  = result #lista com duas listas dentros?!
            fxc[fi] = ofs #lista com duas listas dentros?!
            if plota:
                axes.plot(self.DadosEntrada['data'], result, label = lbl_best+'%s'%fi)


        if plota:
            axes.grid()
            axes.legend( )
            axes.set_ylabel(lbl_vazao)
            plt.show()

        return xc,fxc, result #retorna saida e funcao objetivo


    def FObjPontoFunc(self, param, ffobj):
        """
        Calcula função objetivo dado uma lista de parametros e dicionario de funcoes objetivo       
        Param: array com parâmetros a serem avaliados
        ffobj: dicionario de funcoes objetivo para operar
        """        
        modelado = self.Modelo['fc']( self.DadosEntrada, param )

        if self.NotaPeso:
            mod = array([ modelado[j] for j in self.IndicesPeso ])
            result = [ ffobj['fc'][k][0](mod, self.DadosEntrada["Qexut_peso1"])*ffobj['fc'][k][1] for k in range(ffobj['nb']) ]
        else:
            result = [ ffobj['fc'][k][0](modelado, self.DadosEntrada["Qexut"])*ffobj['fc'][k][1] for k in range(ffobj['nb']) ]
        return array(result)  


    #def Executa(self, arqFinal,npro=5):
    def Executa(self, npro=5):
        
        start_time = time.time()

        print("Otimizacao:",start_time)
        
        if not os.path.exists(self.dirCiclos):  # Para escrever arquivos de saida
            os.makedirs(self.dirCiclos)
        pool = Pool( processes = npro)             # inicia processamento paralelo
        
        indices = [i for i in range(0, self.popTotal) ]
        
        # gera populacao inicial (caso nao tenha recebido valores iniciais)        
        #self.Dx = self.CarregaPop0(arqInicial)
        #self.Dx = self.CarregaPop0Gauss(arqInicial,inflaction=1.2)
        print(" Gerando populacao inicial ...")
        if not hasattr(self, 'Dx'):
            self.Dx = [ MOCOM.IniciaPop(self) for i in range(self.popTotal) ]
            
        # Calculando funções objetivo
        self.Df = array(list(pool.imap(partial(MOCOM.FObj, self), indices )))
        
        # Inicia o cálculo da fronteira de pareto, ordena as matrizes e calcula pdf
        MAXR, Nmax, hier_Pareto = MOCOM.ParetoRank(self)        
        hier_Pareto, self.Df, self.Dx = MOCOM.Ordena(self, hier_Pareto)
        dist_prob = MOCOM.pdf(self, hier_Pareto, MAXR)

        print(" Calculou população e funçoes objetivo iniciais. \n Evoluindo ...")
        ciclo, confere = 0, True  # Variaveis auxiliares para checar convergência da calibracao

        at0 = self.dirCiclos+"atual0.csv"
        at1 = self.dirCiclos+"atual1.csv"
        
        while( MAXR > 1 and confere == True):
            
            st = time.time()
            
            # Simplex Generation -- vou gerar Nmax complexos com nsimplex pontos (ndim + 1 ponto cada)
            LNmax = [i for i in range(Nmax)] 
            Sp = [ MOCOM.Simplex(self, Nmax, hier_Pareto, MAXR, dist_prob, l) for l in LNmax ]
            # obs.: Esta funcao demora mais quando em paralelo
            
            # Evoluindo Simplex        
            novo_ponto_fobj = list( pool.imap(partial( MOCOM.MOSIM, self, MAXR), Sp) )
            # MOSIM é a mais demorada -- VERIFICAR O QUE PODE TORNÁ-LA MAIS RÁPIDA (modelos!!!!)
            # obs.: funcao mais rapida em paralelo. 

            # 3 --- Retorno o novo ponto em Dx e Df
            for i in LNmax:
                ind_ruim = int(self.popTotal - Nmax + i)
                self.Dx[ind_ruim][:], self.Df[ind_ruim][:] = novo_ponto_fobj[i][0], novo_ponto_fobj[i][1]
        
            # 4a --- Nova hierarquização de Pareto e distribuição de probabilidade
            MAXR, Nmax, hier_Pareto = MOCOM.ParetoRank(self)
            hier_Pareto, self.Df, self.Dx = MOCOM.Ordena(self, hier_Pareto)
            dist_prob = MOCOM.pdf(self, hier_Pareto, MAXR)
            
            # 4c --- Checando o número de ciclos
            ciclo = ciclo + 1
            
            if ciclo > self.Nciclos: 
                confere = False
                msg = "# PROCESSO INTERROMPIDO POR EXCESSO DE LOOPS "
            else:
                msg = "# CONVERGÊNCIA COMPLETA"

            message = "Calibrando modelo %s na bacia %s. Ciclo = %d. MAXR = %d. Tempo = %.2f minutos."%(self.Modelo['nome'], self.DadosEntrada["nome"], ciclo, MAXR, (time.time() - start_time)/60 )

            if ciclo%50 == 0: 
                print("Calibrando modelo %s. Ciclo = %d. MAXR = %d. Tempo = %.2f minutos.  %.4f"%(self.Modelo['nome'], ciclo, MAXR, (time.time() - start_time)/60, time.time() - st ))
                MOCOM.Escreve(self, message, namefile = at0 ) 
            elif ciclo%100 == 0:
                MOCOM.Escreve(self, message, namefile = at1 ) 
                
            if ciclo%200 == 0:
                copyfile(self.dirCiclos+"atual0.csv", self.dirCiclos+"params_ciclo%05d_%s.csv"%(ciclo, self.Modelo['nome'].upper() ) )
                
        Stats = array(list(pool.imap(partial(MOCOM.EstatFinal, self), indices )))
        pool.close()
        pool.join()
        
        message  = "MOCOM. Status: %s --- %d ciclos. Tempo total = %.2f minutos. "%(msg, ciclo, (time.time() - start_time)/60.)
        
        if self.arqFinal is None:
            MOCOM.Escreve(self, message, FinalStats = Stats, namefile = self.dirCiclos+"FINAL_%s_%s.csv"%(self.Modelo['nome'].upper(), self.DadosEntrada["nome"].upper() ) )
        else: 
            MOCOM.Escreve(self, message, FinalStats = Stats, namefile = self.arqFinal )


    def CarregaFObj(self, arqname, FobjLista, retorna = 'Fobj'):
        """
        arqname: nome do arquivo contendo valores dos parametros inciais
        FobjLista: lista de strings com nomes das estatisticas que deseja retornar
        retorna: matriz de parâmetros (popTotal x ndim)
        
        Observações:
           - arquivo separado por ';'
           - pula primeira linha
           - segunda linha contém strings com nomes dos parâmetros (não importa a ordem)
           - a partir da terceira linha até o fim do arquivo: cada linha contém os parâmetros   
           - irá usar apenas parâmetros a ser calibrados (caso tenha sido passada uma lista com parâmetros
             a serem mantidos fixos)
        """
        
        arq = open(arqname, 'rt')
        aux = []
        for cont,line in enumerate(arq):
            c = line.split(";")
            if cont == 0: continue
            if cont == 1: pi = [ (cj,j) for j,cj in enumerate(c) if cj in self.DadosEntrada['calibra']  ]
            aux.append(array(c[:-1]))
        arq.close()
        
        aux2 = array(aux).T
        
        stats, par = {}, {}
        tam = len(aux2)
        
        aux3 = {}
        
        for i in range(tam):
            st = aux2[tam - i - 1][0] 
            if st in FobjLista and st not in stats:
                stats[st] = aux2[tam - i - 1][1:]
		
        for st in stats.keys():
            aux3[st] = list( map(float, stats[st]) )
            
	
        if retorna == 'Fobj': return aux3
                
        aux4 = {}
        lista_max=['NS', 'NS_log', 'Pearson','KGE'] #Parâmetros que devem ser maximizados
        lista_min=[ 'Bias', 'rmse', 'rmsei', 'erro_vol', 'Erro_pos', 'Erro_neg', 'Erro_total', 'residuo','rmseqq' ] #Parâmetros que devem ser minimizados

        for p in aux3.keys():
            Smax = None
            if p in lista_max:
                Smax = max(aux3[p])
            elif p in lista_min:
                Smax = min(aux3[p]) 
            if Smax is not None:
                x = aux3[p].index(Smax) + 1
                pr = list( map(float, aux[x]))
                if Smax in pr:
                    f = [ pr[ind] for par,ind in pi ]
                    aux4[p] = array(f)
        if retorna == 'parametros':
            return aux4
            

    def CarregaPop0(self, arqname):
        """
        arqname: nome do arquivo contendo valores dos parametros inciais
        retorna: matriz de parâmetros (popTotal x ndim)
        
        Observações:
           - arquivo separado por ';'
           - pula primeira linha
           - segunda linha contém strings com nomes dos parâmetros (não importa a ordem)
           - a partir da terceira linha até o fim do arquivo: cada linha contém os parâmetros   
           - irá usar apenas parâmetros a ser calibrados (caso tenha sido passada uma lista com parâmetros
             a serem mantidos fixos)
        """
        arq = open(arqname, 'rt')
        aux = []
        for cont,line in enumerate(arq):
            c = line.split(";")
            if cont == 0: continue
            elif cont == 1:
                p =   [ (cj,j) for j,cj in enumerate(c) if cj in self.DadosEntrada['calibra']  ]
                continue
            f = [ float(c[ind]) for par,ind in p ]
            aux.append(array(f))
        return aux

    
    def CarregaPop0Gauss(self, arqname, inflaction=1.):
        """
        arqname: nome do arquivo contendo valores dos parametros inciais
        retorna: matriz de parâmetros (popTotal x ndim) supondo uma distrbuicao multinormal
        
        Observações:
           - arquivo separado por ';'
           - pula primeira linha
           - segunda linha contém strings com nomes dos parâmetros (não importa a ordem)
           - a partir da terceira linha até o fim do arquivo: cada linha contém os parâmetros   
           - irá usar apenas parâmetros a ser calibrados (caso tenha sido passada uma lista com parâmetros
             a serem mantidos fixos)
        """
        arq = open(arqname, 'rt')
        aux = []
        for cont,line in enumerate(arq):
            c = line.split(";")
            if cont == 0: continue
            elif cont == 1:
                p =   [ (cj,j) for j,cj in enumerate(c) if cj in self.DadosEntrada['calibra']  ]
                continue
            f = [ float(c[ind]) for par,ind in p ]
            aux.append(array(f))


        from numpy import cov,diag,vstack
        from numpy.random import randn        
        fa = vstack(aux)            #inicializa array
        covmatrix = cov(fa.T)          #assume covariancia igual do conjunto de dados
        sigma = diag(covmatrix)     #variancias
        mu = fa.mean(axis=0)           #medias
        inflaction = 1.               #fator para inflar/relaxar as variancias
        
        naux = len(aux)
        #yaux = mu+inflaction*randn(naux,1)*sigma #monta array
        yaux = [mu+inflaction*randn(1)*sigma for _ in range(naux)] #monta lista

        #verifica limites
        for i,a in enumerate(yaux):
            il,iu=a<self.Modelo['Linf'],a>self.Modelo['Lsup']
            a[il]=self.Modelo['Linf'][il]
            a[iu]=self.Modelo['Lsup'][iu]
            yaux[i]=a
            
        return yaux


    def IniciaPop(self):
        """
           Gera pontos aleatórios como população inicial
           Retorna array com valores para cada parâmetro do modelo (de acordo com os limites permitidos)
        """
        RN = array( [random() for h in range(self.ndim) ] )      
        P = ( ( self.Modelo['Lsup'] - self.Modelo['Linf'] )*RN + self.Modelo['Linf'] )
        return P


    def FObj(self, indice):
        """
        indice corresponde ao ponto avaliado na matriz self.Dx
        Será usada apenas na parte inicial da calibração para permitir o calculo das funções objetivo de toda 
          a matrix self.Dx por processamento paralelo
        Retorna um array contendo as funções objetivo escolhidas para a calibração
        """
        modelado = self.Modelo['fc']( self.DadosEntrada, self.Dx[indice] )          

        if self.NotaPeso:
            mod = array([ modelado[j] for j in self.IndicesPeso ])
            result = [ self.FunObj['fc'][k][0](mod, self.DadosEntrada["Qexut_peso1"])*self.FunObj['fc'][k][1] for k in range(self.FunObj['nb']) ]
        else:
            result = [ self.FunObj['fc'][k][0](modelado, self.DadosEntrada["Qexut"])*self.FunObj['fc'][k][1] for k in range(self.FunObj['nb']) ]

        return array(result)


    def FObjPonto(self, param):
        """
        Calcula função objetivo especificadas em self.FunObj ao rodar o modelo self.Modelo, dada uma lista de parametros        
        Param: array com parâmetros a serem avaliados 
        """        
        modelado = self.Modelo['fc']( self.DadosEntrada, param )
        
        if self.NotaPeso:
            mod = array([ modelado[j] for j in self.IndicesPeso ])
            result = [ self.FunObj['fc'][k][0](mod, self.DadosEntrada["Qexut_peso1"])*self.FunObj['fc'][k][1] for k in range(self.FunObj['nb']) ]
        else:
            result = [ self.FunObj['fc'][k][0](modelado, self.DadosEntrada["Qexut"])*self.FunObj['fc'][k][1] for k in range(self.FunObj['nb']) ]
        return array(result)
     

    def EstatFinal(self, indice):
        """
        Função para cálculo das estatisticas finais
        indice corresponde ao ponto avaliado na matriz self.Dx
        """
        modelado = self.Modelo['fc']( self.DadosEntrada, self.Dx[indice] )
        
        if self.NotaPeso:
            mod = array([ modelado[j] for j in self.IndicesPeso ])
            obs = self.DadosEntrada["Qexut_peso1"]
        else:
            mod = array(modelado)
            obs = array(self.DadosEntrada["Qexut"])
        
        pearson  = Pearson(mod, obs)
        nash     = NS(mod, obs)
        nash_log = NS_log(mod, obs)
        bias     = Bias(mod, obs)
        adimensional = residuo(mod, obs)
        RMSEI    = rmsei(mod, obs)
        RMSE     = rmse(mod, obs)
        er_vol   = erro_vol(mod, obs)
        kge      = KGE(mod, obs)
        
        # Erro positivo, negativo e total
        pos      = mod > obs
        neg      = obs > mod
        
        erro_pos = abs(sum( (mod - obs)*pos )/max(sum(pos), 1) )
        erro_neg = abs(sum( (mod - obs)*neg )/max( sum(neg), 1) )
        erro_total = abs(sum( (mod - obs) )/len(obs))
    
        return array([nash, nash_log, bias, pearson, RMSE, RMSEI, er_vol, erro_pos, erro_neg, erro_total, adimensional, kge])


    def ParetoRank(self):
        """
        Calcula o Rank de Pareto, retornando um vetor com o índice (grau de
        dominância) dado a cada ponto da população
        
        Retorna:
            h: pior índice do rank
            totalMax: Número de pontos com o pior índice
            ranks: vetor (lin) contendo o índice de cada ponto avaliado
        """
        domina = MOCOM.ParetoDomination(self, self.Df)
        tamObj = len(domina)
        ranks = zeros(tamObj)
        Ndomina = sum(domina)
        h = 0
        status = True
        restante = tamObj - Ndomina
        
        while status:
            h = h + 1
            new_objetivo = zeros( shape = (restante, self.FunObj['nb']) )
            idx = zeros( restante )
            j = 0
            
            for i in range( tamObj ):
                if domina[i] == False:
                    new_objetivo[j] = self.Df[i]
                    idx[j] = i
                    j = j + 1
                elif domina[i] == True and ranks[i] == 0:
                    ranks[i] = h
                else:
                    pass
                    
            if j == 0:
                status = False
            else:
                new_domina = MOCOM.ParetoDomination(self, new_objetivo)
            
                for k in range(restante):
                    domina[ int(idx[k]) ] = new_domina[int(k)]
                Ndomina = sum(domina)
                restante = tamObj - Ndomina
        
        valorMax = max(ranks)
        totalMax = (ranks == valorMax).sum()
                
        return h, totalMax, ranks


    def ParetoDomination(self, vetor):
        """
        Entrada:
            vetor: An (n_points, n_costs) array
        Retorna: 
            is_efficient: array do tipo "boolean". Retorna True para os pontos 
            não dominados e False para os pontos dominados
        """
        is_efficient = ones(vetor.shape[0], dtype = bool)
        for i, c in enumerate(vetor):
            is_efficient[i] = all(any(vetor >= c, axis = 1))
        return is_efficient


    def Ordena(self, a):
        """
        Reordenada a matriz de parâmetros e a matriz de funções objetivo de 
        acordo com os índices do rank de Pareto (vetor a)
        """
        # Cria lista com indices em ordem crescente
        a1 = zeros( self.popTotal )
        b1, c1 = zeros( shape=(self.popTotal, self.FunObj['nb']) ), zeros( shape = (self.popTotal, self.ndim) )
    
        rankind = argsort(a)
        Y = len(rankind)
        
        for i in range(Y):
            a1[i] = a[rankind[i]]
            b1[i] = self.Df[rankind[i]]
            c1[i] = self.Dx[rankind[i]]
        return a1, b1, c1
        

    def pdf(self, valores, MAXR):
        """
        Gera uma função discreta de densidade de probabilidade
        de acordo com os índices do rank de Pareto, de forma que
        índices menores terão maior probabilidade de ser escolhidos.
        Eq. número 2 em Yapo (1998)
        
        Entrada:
            valores: array com os índices de cada ponto de acordo
                     com o rank de Pareto
            MAXR: pior índice verificado no rank de Pareto
        
        Retorna:
            prob: array com as probabilidades de cada índice ser selecionado
        """
        denominador = sum(MAXR - valores + 1)
        tam = len(valores)
        prob1 = MAXR - valores + 1
        prob = prob1/denominador
        return prob


    def aleatorio(self, probabilities):
        """
        Baseado na distribuição de probabilidades de selecionar cada um dos
        índices do rank de Pareto, retorna o índice aleatório
        
        probabilities: array com probabilidades de cada índice do
            rank ser selecionado
        """    
        x = random()
        cumulative_probability = 0.0
        indice = 0
        for item_probability in probabilities:
            cumulative_probability += item_probability        
            if x < cumulative_probability: break
            else: indice = indice + 1
        return indice


    def gera_simplex(self, valores, MAXR, prob):
        """
        Gera simplex para etapa de evolução. Escolhe n números aleatórios (com reposição)
        entre os pontos que obtiveram índice < MAXR. 
        
        Entrada:
            valores: array com índices (rank de Pareto) de cada ponto
            MAXR: pior índice do rank
            prob: distribuição de probabilidade
            
        Retorna:
            Bx, Bf: array com conjunto de parâmetros e funções objetivo dos
                    pontos aleatórios
        """
        Bx, Bf = zeros( shape=(self.ndim, self.ndim) ), zeros( shape = (self.ndim, self.FunObj['nb']) )
        L = ones(self.ndim)*-99
        h = 0
        while h < self.ndim:
            x = MOCOM.aleatorio(self, prob)
            #~ if x not in L and valores[x] != MAXR: #(sem reposição)
            if valores[x] != MAXR:
                L[h] = x
                h = h + 1
        L = sorted(L)    
        L = array([int(h) for h in L])
        
        for i in range(0, self.ndim):
            Bx[i, :] = self.Dx[L[i], :]
            Bf[i, :] = self.Df[L[i], :]
            
        return Bx, Bf


    def Simplex(self, Nmax, hier_Pareto, MAXR, dist_prob, m):
        """
        Segunda parte da função simplex
        
        Recebe os números aleatórios escolhidos entre os índices < Nmax e inclui um dos piores
        pontos na posição final para evolução.
        
        Entrada:
            Nmax, popTotal, nsimplex, ndim, nPar: parâmetros e dimensões
            hier_pareto: índices de pareto
            MAXR: maior índice dos rank
            dist_prob: distribuuição de probabilidade
            Dx, Df: matriz de pontos e funções objetivo, respectivamente
            m: índice do ponto com pior rank de pareto
        
        Retorna:
            Bx, Bf: simplex pronto para etapa de evolução
        """
        nsimplex = self.ndim + 1
        Bx, Bf = zeros( shape=(nsimplex, self.ndim) ), zeros( shape = (nsimplex, self.FunObj['nb']) )
        Bx[:-1,:], Bf[:-1,:] = MOCOM.gera_simplex(self, hier_Pareto, MAXR, dist_prob)
        
        ind_ruim = int(self.popTotal - Nmax + m)
        Bx[-1, :], Bf[-1,:]= self.Dx[ind_ruim, :], self.Df[ind_ruim, :]
        return [Bx, Bf]


    def MOSIM(self, MAXR, Sp):
        """
        MOSIM: Multiobjective Downhill Simplex Method. Método de evolução do pior
        ponto do conjunto
        
        Entrada:
            Bf: array sorteado do simplex (ndim + 1, nPar) com as funções objetivo
            Bx: Conjunto de parâmetros sorteado (ndim + 1, ndim) no simplex
            ndim: número de parâmetros/dimensão
            Linf, Lsup: limites dos parâmetros
            ETp, CMB, Qmont, Ainc, Qexut, peso: parâmetros para o cálculo da vazão
            MAXR: pior índice do rank de pareto
            
        Retorna:
            Sref ou Scont: escolhe entre o ponto de reflexão ou contração e retorna o melhor
            OBJref: novo vetor de funções objetivo calculado para o novo ponto 
        """
        #[Bx, Bf]
        Bx = Sp[0]
        Bf = Sp[1]
        
        # 1 - Separo o pior ponto para evoluir
        pior_ponto = Bx[-1,:]
        
        # 2 - Cálculo do centróide do simplex sem o pior ponto
        Sg =  sum(Bx[0:-1, :])/self.ndim
        
        # 3 - Calcula o ponto de reflexão
        Sref = 2*Sg - pior_ponto
        menor = sum(Sref < self.Modelo['Linf'])
        maior = sum(Sref > self.Modelo['Lsup'])
        
        # 3a - Verifica se o ponto de reflexão está fora dos limites
        if (menor + maior ) > 0:
            Sref = MOCOM.IniciaPop(self)
        
        # 3b - Calcula o vetor de funções objetivo e verifica a dominância
        OBJref = MOCOM.FObjPonto(self, Sref)    
        Bf[-1, :] = OBJref    
        dominados = MOCOM.ParetoDomination(self, Bf)
        
        # 3c - Verifica se o ponto é não-dominado pelos demais: se não for (rank = 1),
          # Sref substitui o pior ponto - caso contrário, o ponto de contração é aceito imediatamente
          
        if dominados[-1] == True:
            return [Sref, OBJref]
            # 4 - Se o ponto de reflexão não for aceito, o ponto de contração será aceito imediatamente        
        else:
            Scont = 0.5*(Sg + pior_ponto)
            menor = sum(Scont < self.Modelo['Linf'])
            maior = sum(Scont > self.Modelo['Lsup'])
            # 4a - Verifica se está fora dos limites - se estiver, gera ponto aleatório
            if (menor + maior ) > 0:
                Scont = MOCOM.IniciaPop(self)
            
            OBJref = MOCOM.FObjPonto(self, Scont)
            
            return [Scont, OBJref]


    def Escreve(self, message, FinalStats = None, namefile = None):
        """
        Funcao auxiliar para escrever resultados temporários em um arquivo
        """
        saida = open("%s"%namefile, "wt" )
        
        saida.write("%s \n"%message)
        
        for nm in self.FunObj['nomes']:
            saida.write("%s;"%nm)
            
        if FinalStats is not None:
            for fc in estatisticas:
                saida.write("%s;"%fc)
        
        for pr in self.DadosEntrada['calibra']:
            saida.write("%s;"%pr)
            
        for pr in self.DadosEntrada['OrdemParFixos']:
            saida.write("%s"%pr)
                #~ self.DadosEntrada['OrdemParFixos'] = []

        saida.write('\n')
        
        for j in range(self.popTotal):
            if self.FunObj['nb'] <=1: 
                saida.write("%.8f;"% (self.Df[j] ) )
            else:
                for k in range(self.FunObj['nb']):
                    saida.write("%.8f;"% (self.Df[j][k]) )
                
            if FinalStats is not None:
                for l in range(len(FinalStats[j][:])):
                    saida.write("%.8f;"% FinalStats[j][l] )
            
            for i in range(self.ndim):
                saida.write( "%.8f;" % self.Dx[j][i] )
            
            for par in self.DadosEntrada['OrdemParFixos']:
                saida.write( "%.8f" % self.DadosEntrada['ParFixos'][par] )
            saida.write("\n")
        saida.close()
