import pandas as pd
from wordcloud import WordCloud
import matplotlib.pyplot as plt
import numpy

data={'FUSE': 19917,
'ADAS': 942,
'BoundaryPlasmaModels': 488,
'CHEASE': 258,
'CoordinateConventions': 392,
'EPEDNN': 564,
'FiniteElementHermite': 503,
'Fortran90Namelists': 523,
'FuseUtils': 45,
'FusionMaterials': 529,
'FXP': 149,
'IMAS': 14528,
'IMASDD': 4037,
'MXHEquilibrium': 2626,
'MeshTools': 370,
'MillerExtendedHarmonic': 1281,
'NEO': 1042,
'NNeutronics': 238,
'QED': 587,
'SimulationParameters': 1831,
'TEQUILA': 3035,
'TGLFNN': 1271,
'TJLF': 12809,
'VacuumFields': 2003,
'XSteamTP': 4218,
'ThermalSystemModels': 8197}

avg = sum([v for v in data.values()]) / len(data)

data['TJLF'] -= 5000
data['IMAS'] -= 5000
data['FUSE'] -= 5000
del data['Fortran90Namelists']
for key in data:
    data[key] = numpy.log(data[key]+avg / 10)

# Convert the data into a suitable format
df = pd.DataFrame(list(data.items()), columns=['Package', 'Lines of Code'])

# Generate the word cloud
wordcloud = WordCloud(width=600, height=800, background_color='white').generate_from_frequencies(data)

# Display the word cloud
plt.figure(figsize=(6, 8))
plt.imshow(wordcloud, interpolation='bilinear')
plt.axis('off')
plt.show()