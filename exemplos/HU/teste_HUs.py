import sys
sys.path.append('/Users/arlan/github/modelos/')
import HU
import matplotlib.pyplot as plt

# Exemplo do artigo de Jeng e Coon (2003)
t = [i for i in range(0,114,3)]
ERH = [0] * len(t)
ERH[0] = 0
ERH[1] = 51.7915
ERH[2] = 99.9585
ERH[3] = 235.8793
DRH1_ex = [0, 2.1566,  8.8020, 23.7430, 35.5430, 36.4216, 34.7098, 31.9036, 28.6874, 25.4162, 22.2785, 19.3709, 16.7364, 14.3867, 12.3151, 10.5049, 8.9341, 7.5787, 6.4146, 5.4188, 4.5696, 3.8476, 3.2352, 2.7168, 2.2790, 1.9098, 1.5988, 1.3374, 1.1178, 0.9336, 0.7792, 0.6499, 0.5418, 0.4514, 0.3759, 0.3129, 0.2603, 0.2164]
DRH2_ex = [0, 3.4846, 11.5534, 30.5041, 37.5314, 39.3604, 38.0968, 35.1316, 31.3636, 27.3562, 23.4469, 19.8238, 16.5780, 13.7397, 11.3022,  9.2382, 7.5099, 6.0761, 4.8956, 3.9300, 3.1445, 2.5087, 1.9961, 1.5844, 1.2549, 0.9920, 0.7827, 0.6165, 0.4848, 0.3807, 0.2986, 0.2339, 0.1830, 0.1430, 0.1116, 0.0870, 0.0678, 0.0528]
# Metodo 1 - HUI de Nash
DRH1_meu = HU.HUI_nash(ERH, 1.49, 15.1, dt=3)
# Metodo 2 - HUI de Jeng e Coon
DRH2_meu = HU.HUI_3rsv_p_dist(ERH, 10.535, 0.1728, 0.7867, dt=3)
# Visualizar
# plt.plot(t, DRH1_ex, label='Nash ex', color='black')
# plt.plot(t, DRH1_meu, label='Nash meu')
plt.plot(t, DRH2_ex, label='3RSV ex', color='black')
plt.plot(t, DRH2_meu, label='3RSV meu')
plt.legend()
plt.grid()
plt.show()
