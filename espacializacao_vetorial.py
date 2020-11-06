'''
--------------------------------------------------------------------------------
ALGORITMOS DE ESPACIALIZACAO
--------------------------------------------------------------------------------
Implementacao - Arlan Scortegagna, maio/2019
Ultima atualizacao - Arlan Scortegagna, nov/2020
Revisao - PENDENTE
--------------------------------------------------------------------------------
Descricao:
    Este script contem os algoritmos que permitem o calculo de Estimativas
    Quantititativas de Precipitacao (QPEs) para os pontos de uma grade, a partir
    dos dados pontuais de medicoes pluviometricas.
--------------------------------------------------------------------------------
Metodos:
    IDW - Inverse Distance Weighting
    ...
--------------------------------------------------------------------------------
Entradas:
    grade   - GeoDataFrame contendo a geometria da grade de pontos
    path_pd - path da pasta contendo os arquivos .pd (todos nela serao usados)
    EPSG    - codigo int do SRC projetadas (ex: 31983 - SIRGAS 2000 / UTM 23S)
--------------------------------------------------------------------------------
Saidas:
    DF_grade - DataFrame contendo a precipitacao interpolada nos pontos da grade
    PME      - Serie de dados de Precipitacao Media Espacial
--------------------------------------------------------------------------------
Observacoes:
    Os arquivos .pd devem estar no formato preconizado pelo padrao de dados
    hidrologicos do Simepar.
    As coordenadas nos arquivos .pd DEVEM ser geograficas (lat/long em graus
    decimais e em DATUM WGS84).
--------------------------------------------------------------------------------
'''

import geopandas as gpd
import pandas as pd
import glob
import numpy as np

def tem(valor):#Função para verificar se a estação tem dados na data
    if np.isnan(valor):
        return 0
    return 1

def idw(grade, path_pd, EPSG):

    # 1 - Converter os pontos da geometria para o SRC projetado
    grade = grade.to_crs('EPSG:{}'.format(EPSG))

    # 2 - Inicializar os DataFrames
    DF_postos = pd.DataFrame()
    DF_dists = pd.DataFrame()

    # 3 - Abrir os arquivos .pd individualmente
    arquivos_pd = glob.glob(path_pd+'*.pd')

    for arquivo_pd in arquivos_pd:
        # 4 - Capturar as coordenadas do posto no arquivo .pd
        nome_posto = arquivo_pd.split('/')[-1].split('.pd')[0]
        lat_long = pd.read_csv(arquivo_pd, nrows=1, header=None).iloc[0]
        lat  = lat_long[0]
        long = lat_long[1]

        # 5 - Converter a coordenada de WGS84 (padrao) para o SRC projetado
        df_temp = pd.DataFrame({'posto':[nome_posto],
                                'Lat':[lat],
                                'Long':[long]})
        geometry = gpd.points_from_xy(x=df_temp.Long, y=df_temp.Lat)
        gdf_temp = gpd.GeoDataFrame(df_temp, geometry=geometry)
        gdf_temp.set_crs('EPSG:4326', inplace=True)
        gdf_temp = gdf_temp.to_crs('EPSG:{}'.format(EPSG))
        posto_x = gdf_temp.loc[0].geometry.x
        posto_y = gdf_temp.loc[0].geometry.y

        # 6 - Calcular a distancia (km) entre o posto e todos os pontos da grade
        ponto_x=grade.geometry.x
        ponto_y=grade.geometry.y
        dist = ((posto_x-ponto_x)**2+(posto_y-ponto_y)**2)**(1/2)
        DF_dists[nome_posto]=dist/1000

        # 7 - Coletar a serie de dados do posto e conceatenar em DF_postos
        sr_posto = pd.read_csv(arquivo_pd, skiprows=1, parse_dates=True, index_col='data')['p']
        sr_posto = sr_posto.rename(nome_posto)
        DF_postos = DF_postos.join(sr_posto, how='outer')

    # 8 - Calcular a precipitacao interpolada para cada ponto de grade
    L = len(grade)
    no_postos = len(DF_postos.columns)

    DF_grade = pd.DataFrame(columns=grade.index)
    DF_grade.index.name = 'data'
    for pi in grade.index:
        print(f'processando grade {pi}')
        D=DF_dists.loc[pi].values#vetor de distancia
        W=1/(D**2)# vetor de pesos
        aux=pd.DataFrame({'WD':(DF_postos*W).sum(axis=1),'soma':(DF_postos.apply(np.vectorize(tem))*W).sum(axis=1)})
        aux['res']=aux['WD']/aux['soma']
        DF_grade[pi]=aux['res'].round(decimals=2)

    # 9 - Calcular a PME na grade
    PME = DF_grade.mean(axis=1, skipna=True)

    return DF_grade, PME
