import corner
import numpy as np
import matplotlib.pyplot as plt

true_r_cut = 38.39/4.9213  # we need to calculate D_L for conversion from kpc to arcsecond
true_sigma = 282.28

folder_name = '20250218'

figure = plt.figure(figsize=(8,8))
# fig_sigma, axes_sigma = plt.subplots(3,1,constrained_layout=True)
# fig_r, axes_r = plt.subplots(2,1,constrained_layout=True)

def smooth(samples,window_size = 100, step_size = 50):
    sample_num = samples.shape[-1]
    sm_samples = []
    for start in range(0,sample_num-window_size+1,step_size):
        end = start + window_size
        window_samples = samples[:,start:end]
        sm_sample = np.mean(window_samples,axis=1)
        sm_samples.append(sm_sample)
    return np.array(sm_samples)

samples_data = {}
samples_data_chain = {}
names = []
with open('../lens-data/RXJ2248/origin_data/bayes.dat', 'r') as file:
    lines = file.readlines()
    chain_num = 0
    for line in lines:
        if line.startswith('#'):
            name = line.strip().lstrip('#')
            samples_data[name] = []
            samples_data_chain[name] = []
            for i in range(10):
                samples_data_chain[name].append([])
            names.append(name)
        else:
            sample_values = line.split()
            for i, value in enumerate(sample_values):
                samples_data[names[i]].append(float(value))
                samples_data_chain[names[i]][chain_num].append(float(value))
            chain_num = (chain_num+1)%10
    corner.corner(np.column_stack((samples_data['Pot0 rcut (arcsec)'],samples_data['Pot0 sigma (km/s)'])),labels=['r_cut(arcsec)','sigma(km/s)'],truths=[true_r_cut,true_sigma],color='grey',fig=figure)
    # sm_sigma_chain = smooth(np.array(samples_data_chain['Pot0 sigma (km/s)']))
    # axes_sigma[2].plot(sm_sigma_chain)

# samples_data = {}
# samples_data_chain = {}
# names = []
# with open(f'../lens-data/RXJ2248/without_prior/{folder_name}/bayes.dat', 'r') as file:
#     lines = file.readlines()
#     chain_num = 0
#     for line in lines:
#         if line.startswith('#'):
#             name = line.strip().lstrip('#')
#             samples_data[name] = []
#             samples_data_chain[name] = []
#             for i in range(10):
#                 samples_data_chain[name].append([])
#             names.append(name)
#         else:
#             sample_values = line.split()
#             for i, value in enumerate(sample_values):
#                 samples_data[names[i]].append(float(value))
#                 samples_data_chain[names[i]][chain_num].append(float(value))
#             chain_num = (chain_num+1)%10
#     corner.corner(np.column_stack((samples_data['Pot0 rcut (arcsec)'],samples_data['Pot0 sigma (km/s)'])),labels=['r_cut(arcsec)','sigma(km/s)'],truths=[true_r_cut,true_sigma],color='red',fig=figure)
#     # sm_sigma_chain = smooth(np.array(samples_data_chain['Pot0 sigma (km/s)']))
#     # axes_sigma[0].plot(sm_sigma_chain)

samples_data = {}
samples_data_chain = {}
names = []
with open(f'../lens-data/RXJ2248/with_prior/{folder_name}/bayes.dat', 'r') as file:
    lines = file.readlines()
    chain_num = 0
    for line in lines:
        if line.startswith('#'):
            name = line.strip().lstrip('#')
            samples_data[name] = []
            samples_data_chain[name] = []
            for i in range(10):
                samples_data_chain[name].append([])
            names.append(name)
        else:
            sample_values = line.split()
            for i, value in enumerate(sample_values):
                samples_data[names[i]].append(float(value))
                samples_data_chain[names[i]][chain_num].append(float(value))
            chain_num = (chain_num+1)%10
    corner.corner(np.column_stack((samples_data['Pot0 rcut (arcsec)'],samples_data['Pot0 sigma (km/s)'])),labels=['r_cut(arcsec)','sigma(km/s)'],truths=[true_r_cut,true_sigma],color='blue',fig=figure)
    # sm_sigma_chain = smooth(np.array(samples_data_chain['Pot0 sigma (km/s)']))
    # axes_sigma[1].plot(sm_sigma_chain)

i = np.argmin(samples_data['Chi2'])

sigma_mean = samples_data['Pot0 sigma (km/s)'][i]
rcut_mean = samples_data['Pot0 rcut (arcsec)'][i]
print(f'average sigma={sigma_mean}, average rcut={rcut_mean}')

plt.plot([],[],color='grey',label='origin fitting')
plt.plot([],[],color='red',label='my fitting without prior')
plt.plot([],[],color='blue',label='my fitting with prior')
plt.legend(loc='upper right', bbox_to_anchor=(1,1.5))
plt.savefig(f'/home/lyumx/work/lens-project/{folder_name}.png')
# plt.show()

plt.figure()
sm_chi2_chain = smooth(np.array(samples_data_chain['Chi2']))
plt.plot(sm_chi2_chain[:,0])
# plt.show()

plt.figure()
plt.plot(np.array(samples_data_chain['Chi2'][0]))
plt.show()