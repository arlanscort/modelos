import matplotlib.pyplot as plt

def plotar_hidro(idx, PME, ETP, Qobs, Qmon=None, Qsims=None):

    fig, (ax1, ax2) = plt.subplots(2, 1,sharex='all',gridspec_kw={'height_ratios': [1, 3]})

    ax1.bar( idx, PME, label='PME', color='blue')
    ax1.plot(idx, ETP, label='etp', color='red')
    ax1.invert_yaxis()
    ax1.set_ylabel('Altura de precipitação (mm)', fontsize=8)
    ax1.legend(loc='upper right', fontsize=8)


    ax2.plot(idx, Qobs, label='Qobs', color='black')
    ax2.set_ylabel('Vazão (m3/s)')
    if Qmon is not None:
        ax2.plot(idx, Qmon, label='Qmon', color='black', linestyle='--')

    if Qsims is not None:
        for chave in Qsims.keys():
            ax2.plot(idx, Qsims[chave], label=chave)
            ax2.legend(loc='upper right', fontsize=8)

    return fig
