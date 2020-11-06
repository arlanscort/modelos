'''
Arquivo de exemplo para executar a espacializacao da precipitacao
'''

# 1 - Importar as bibliotecas necessarias, principalmente geopandas
import geopandas as gpd

# 2 - Adicionar o path do repositorio "modelos"
import sys
sys.path.append('/Users/arlan/github/modelos/')

# 3 - Importar a funcao de interpolacao desejada do modulo de espacializacao
from espacializacao import idw

# 4 - Coletar as entradas
grade = gpd.read_file('grade.gpkg')
path_pd = 'dados/'
EPSG = 31983

# 5 - Executar IDWx
DF_grade, PME = idw(grade, path_pd, EPSG)

# 6 - Exportar a PME no padrao preconizado nas diretrizes da Hidrologia - Simepar
PME.rename('p').round(2).to_csv('pme_teste.pd', index_label='data')
# Obs: falta salvar o lat,long no arquivo .pd, o ideal seria o centroide da grade
