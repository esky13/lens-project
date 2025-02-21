from scipy.optimize import curve_fit
from pyLensLib import lenstool as lst
import numpy as np
import matplotlib.pyplot as plt

VISUALIZATION = False

cluster_name = 'RXJ2248'
# read galaxy and cluster information from this file.
filename = f'/home/lyumx/work/lens-data/{cluster_name}/with_prior/20250220_withscatter_lenstool/scatter.par'

co = lst.getCosmo(filename)
df = lst.getClMembers(filename)
fov = lst.getFoV(filename)
z_source = 1
z_lens = df.loc[0].z_lens

sigma0 = 282.284649
rcut0 = 38.386646*1e-3/co.angular_diameter_distance(z_lens).value*3600/np.pi*180
rcore0 = 0.05
m0 = 16.176001
slope = 3
vdslope = 3.740000
# fit the velocity dipersion-m relation and visualize
v = []
l = []
r = []
rcore = []
for i in range(len(df)):
    gal = df.loc[i]
    if len(gal.gal_id) < 3:
        continue
    v.append(gal.v_disp)
    r.append(gal.cut_radius)
    rcore.append(gal.core_radius)

# read magnitude from members.cat
with open(f'/home/lyumx/work/lens-data/{cluster_name}/members.cat',"r") as f:
    lines = f.readlines()
    for line in lines[1:]:
        numbers = line.strip().split()
        l.append(float(numbers[6]))
v = np.array(v)
r = np.array(r)
rcore = np.array(rcore)
l = np.array(l)
l_line = np.linspace(np.min(l),np.max(l),100)

if VISUALIZATION:
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    plt.scatter(l,v)
    plt.plot(l_line,sigma0*10**(0.4*(m0-l_line)/vdslope))
    plt.xlabel('magnitude')
    plt.ylabel('velocity dispersion (km/s)')
    plt.show()

    fig,ax = plt.subplots(1,1,figsize=(10,10))
    plt.scatter(l,r)
    plt.plot(l_line,rcut0*10**(0.4*(m0-l_line)*2/slope))
    plt.xlabel('magnitude')
    plt.ylabel('cut radius (arcsec)')
    plt.show()

# fit slope sigma0
def v_m_func(l,sigma0):
    m0 = 16.176001
    return sigma0*10**(0.4*(m0-l)/vdslope)

def rcut_m_func(l,rcut0):
    return rcut0*10**(0.4*(m0-l)*2/slope)

popt_v,_ = curve_fit(v_m_func,l,v)
popt_rcut,_ = curve_fit(rcut_m_func,l,r)
print(f'fit velocity sigma0={popt_v[0]}')
print(f'given velocity slope={vdslope}, sigma0={sigma0}')
print(f'fit cut radius rcut0={popt_rcut[0]}')
print(f'given cut radius slope={slope}, rcut0={rcut0}')
