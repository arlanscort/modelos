from matplotlib import pyplot as plt
from datetime import datetime

Q, peso = {}, {}
arq = open('ecv_simul_zdr_pc.txt', 'r')
for l in arq:
    if l[0] == '#': continue
    l = l.split()
    dt = datetime(int(l[0]), int(l[1]), int(l[2]), int(l[3]), 0, 0)
    Q[dt] = [float(l[7]), None, None]
    peso[dt] = float(l[8])
arq.close()

arq = open('result_fortran.txt', 'r')
for l in arq:
    if l[0] == '#': continue
    l = l.split()
    dt = datetime(int(l[0]), int(l[1]), int(l[2]), int(l[3]), 0, 0)
    Q[dt][1] = float(l[4])
arq.close()

arq = open('result_python.txt', 'r')
for l in arq:
    if l[0] == '#': continue
    l = l.split()
    dt = datetime(int(l[0]), int(l[1]), int(l[2]), int(l[3]), 0, 0)
    Q[dt][2] = float(l[4])
arq.close()

x = sorted(Q)
obs = [Q[i][0] for i in x]
fortran = [Q[i][1] for i in x]
python = [Q[i][2] for i in x]

plt.plot(x, obs, 'k-', label="Observado")
plt.plot(x, fortran, 'b-', label="Fortran")
plt.plot(x, python, 'g-', label="Python")
plt.grid()
plt.legend()
plt.show()