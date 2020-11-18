import sys
sys.path.append('/Users/arlan/github/modelos/')
import propagacao
import matplotlib.pyplot as plt

# Teste muskingum
# Exercicio de propagacao hidrologica
# Ven te Chow, pg. 262
t = [i for i in range(1,21)]
Qmon = [93,137,208,320,442,546,630,678,691,675,634,571,477,390,329,247,184,134,108,90]
Qresp = [85,91,114,159,233,324,420,509,578,623,642,635,603,546,479,413,341,274,215,170]
k = 2.3
x = 0.15
Qmeu = propagacao.muskingum(Qmon, k, x, qini=85, dt=1)
# Plotar
plt.plot(t, Qmon, label='Inflow (m3/s)', color='blue')
plt.plot(t, Qresp, label='Outflow - exercicio (m3/s)', color='black')
plt.plot(t, Qmeu, label='Outflow - obtido (m3/s)', color='red', linestyle=':')
plt.legend()
plt.show()
