import sys
sys.path.append('/Users/arlan/github/modelos/')
import propagacao
import matplotlib.pyplot as plt

# Teste - exercicio de propagacao hidrologica do Chow (1988) p. 262
# Dados
t = [0, 3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 60]
Qmon = [0, 100, 300, 200, 100, 0, 0, 0, 0, 0, 0, 0] # m3/s (Qmon = chuva, neste caso)
Qobs = [0, 0, 10, 70, 165, 180, 142, 79, 38, 13, 3, 0] # m3/s
dt = 6   # horas
k = 3.14 # horas
n = 5.31
# Executa funcoes de propagacao
Qprop_analitica = propagacao.rsvs_lineares_nash_sol_analitica(Qmon, k, n, dt=dt)
Qprop_diferencial = propagacao.rsvs_lineares_nash_sol_diferencial(Qmon, k, n, dt=dt)
# Plota
plt.plot(t, Qmon, label='Chuva (m3/s)', color='black')
plt.plot(t, Qobs, label='Vazão observada (m3/s)', color='black')
plt.plot(t, Qprop_analitica, label='Vazão simulada - Nash solução analítica (m3/s)')
plt.plot(t, Qprop_diferencial, label='Vazão simulada - Nash solução diferencial (m3/s)')
plt.legend()
plt.show()
