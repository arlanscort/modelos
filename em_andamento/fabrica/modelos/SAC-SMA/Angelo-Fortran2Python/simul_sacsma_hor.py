from datetime import datetime
from modelos_hidrologicos import SACSMA

params = "1.000   7.673 0.748994 115.498 5.277121 144.552  32.997  112.301 0.349996 0.034978 0.199997 0.712128 0.000091 0.115548 0.287332   7.9".split()
arqIN = "ecv_simul_zdr_pc.txt"
arqOUT = "result_python.txt"
area = 541.82

#Parâmetros do modelo
prm = {"UZTWM": float(params[0])}
prm["UZFWM"] = float(params[1])
prm["UZK"]   = float(params[2])
prm["ZPERC"] = float(params[3])
prm["REXP"]  = float(params[4])
prm["LZTWM"] = float(params[5])
prm["LZFSM"] = float(params[6])
prm["LZFPM"] = float(params[7])
prm["LZSK"]  = float(params[8])
prm["LZPK"]  = float(params[9])
prm["PFREE"] = float(params[10])
prm["SIDE"]  = float(params[11])
prm["PCTIM"] = float(params[12])
prm["ADIMP"] = float(params[13])
prm["Kprop"] = float(params[14])
prm["lag"]   = int(float(params[15]))

#Dados de entrada
datas, evap, chuva, qmont, qexut, peso = [], [], [], [], [], []
arq = open(arqIN, 'r')

for l in arq:
    
    l = l.split()
    datas.append(datetime(int(l[0]), int(l[1]), int(l[2]), int(l[3]), 0, 0))
    evap.append(float(l[4]))
    chuva.append(float(l[5]))
    qmont.append(float(l[6]))
    qexut.append(float(l[7]))
    peso.append(float(l[8]))
    
arq.close()


#Série de vazão calculada pelo modelo
qmod = SACSMA(evap, chuva, qmont, prm, area, qexut[0])


#Estatísticas de desempenho

#Calculando valor média da vazão na exutória, vazão modelada e dos resíduos
media = [0.0, 0.0, 0.0]

for i in range(len(qmod)):
    
    media[0] = media[0] + qexut[i] * peso[i]              #acumulando valores de vazão observada
    media[1] = media[1] + qmod[i] * peso[i]               #acumulando valores de vazão modelada
    media[2] = media[2] + (qmod[i] - qexut[i]) * peso[i]  #acumulando valores do resíduo

sp = sum(peso)
media = [x/sp for x in media] #média ponderada

#Calculando desvio-padrão das séries de vazão e de resíduo
#Faz o somatório das séries de erro e variância
desvp = [0.0, 0.0, 0.0]
func  = [0.0 for i in range(9)]

for i in range(len(qmod)):
    
    res = qmod[i] - qexut[i]
    desvp[0] = desvp[0] + ((qexut[i] - media[0])**2) * peso[i]
    desvp[1] = desvp[1] + ((qmod[i] - media[1])**2) * peso[i]
    desvp[2] = desvp[2] + ((res - media[2])**2) * peso[i]
    
    func[0] = func[0] + abs(res) * peso[i]
    func[1] = func[1] + ((qmod[i] - qexut[i])**2) * peso[i]
    
    if res > 0:
        func[2] = func[2] + res * peso[i]
        func[3] = func[3] + peso[i]
    
    if res < 0:
        func[4] = func[4] + res * peso[i]
        func[5] = func[5] + peso[i]
    
    func[6] = func[6] + (qexut[i] - media[0]) * (qmod[i] - media[1]) * peso[i]
    func[7] = func[7] + ((qexut[i] - media[0])**2) * peso[i]

desvp = [(x/sp)**0.5 for x in desvp]

#Contabilizando estatísticas
func[0] = func[0]/sp                           #Erro Absoluto Médio (m3/s)
func[7] = 1.0 - func[1]/func[7]                #Coeficiente de Nash-Sutcliffe (adim.)
func[1] = (func[1]/sp)**0.5                    #Raiz do Erro Quadrático Médio (m3/s)
func[2] = func[2]/func[3]                      #Erro Positivo Médio (m3/s)
func[3] = func[3]*100.0/sp                     #Frequência de Erros Positivos (%)
func[4] = func[4]/func[5]                      #Erro Negativo Médio (m3/s)
func[5] = func[5]*100.0/sp                     #Frequência de Erros Negativos (%)
func[6] = func[6]/(desvp[0]*desvp[1]*sp)       #Correlação linear de Pearson (adim.)
func[8] = 0.30 * abs(media[2])/media[0] + 0.70 * desvp[2]/desvp[0]   #Coeficiente sem nome (adim.)


#Gravando resultados em arquivo
arq = open(arqOUT, 'w')
arq.write('# %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %9.6f %9.6f %9.6f\n' % tuple(func))

for i in range(len(qmod)):
    arq.write('%s %7.2f\n' % (datas[i].strftime('%Y %m %d %H'), qmod[i]))

print('Simulacao HORARIA concluida!')
#1230123012301230123012301230123012301230123012301230123012301230123012301230123012301230123012301230123012301230123012301230123
