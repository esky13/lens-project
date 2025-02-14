import corner
import numpy as np
import matplotlib.pyplot as plt

folder_name = 'test'

figure = plt.figure(figsize=(8,8))

samples_data = {}
names = []
#O1 : x (arcsec)
#O1 : y (arcsec)
#O1 : emass
#O1 : theta (deg)
#O1 : rc (arcsec)
#O1 : sigma (km/s)
with open('../lens-data/RXJ2248/a_model_test/bayes.dat', 'r') as file:
    lines = file.readlines()
    chain_num = 0
    for line in lines:
        if line.startswith('#'):
            name = line.strip().lstrip('#')
            samples_data[name] = []
            names.append(name)
        else:
            sample_values = line.split()
            for i, value in enumerate(sample_values):
                samples_data[names[i]].append(float(value))
    corner.corner(np.column_stack((samples_data['O1 : x (arcsec)'],samples_data['O1 : y (arcsec)'],\
                                   samples_data['O1 : emass'],samples_data['O1 : theta (deg)'],\
                                   samples_data['O1 : rc (arcsec)'],samples_data['O1 : sigma (km/s)']))\
                                    ,labels=['x center (arcsec)','y center (arcsec)','ellipticity','theta (deg)','r cut (arcsec)','sigma(km/s)'],color='grey',fig=figure)
plt.savefig(f'/home/lyumx/work/lens-project/{folder_name}.png')
plt.show()
