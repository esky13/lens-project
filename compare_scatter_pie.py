from pyLensLib import lenstool as lst
import numpy as np
import matplotlib.pyplot as plt
from pyLensLib.sersic_numba import sersic
from pyLensLib.pointsrc import pointsrc
import string
from pyLensLib.deflector import deflector
from pyLensLib.piemd import piemd
import astropy.constants as const
from scipy.optimize import curve_fit
from tkinter import Tk, Entry, Button, Label, StringVar
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def my_create_deflector(number,fov,dpie,z_lens,z_source):
    theta1 = np.linspace(fov[0],fov[1],number)
    theta2 = np.linspace(fov[2],fov[3],number)
    theta1, theta2 = np.meshgrid(theta1,theta2)

    angx = np.zeros((number,number))
    angy = np.zeros((number,number))
    for gal_i in dpie:  
        ax, ay = gal_i.angle(theta1,theta2)
        angx += ax
        angy += ay
    kwargs_def = {'zl': z_lens, 'zs': z_source}
    df = deflector(co, angx=angx, angy=angy, **kwargs_def)
    df.setGrid(theta=np.linspace(fov[0],fov[1],number), compute_potential=False)
    return df

# three files are required:
# best.par
# RXJ2248_angx.fits, RXJ2248_angy.fits

cluster_name = 'RXJ2248'
filename = f'/home/lmx/work/lens-data/{cluster_name}/origin_data/best.par'

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
# generate scatter velocity and cut radius
dpie_scatter = []
for i in range(len(df)):
    gal = df.loc[i]
    e = 1-np.sqrt((1-gal.ellipticite)/(1+gal.ellipticite))
    q = 1-e
    if len(gal.gal_id) < 3:
        noise_scale = 0
        offset = 0
    else:
        noise_scale = 0/20
        offset = 0/10
    kwargs = {'zl': gal.z_lens,
            'zs': z_source,
            'sigma0': gal.v_disp*np.sqrt(3/2)+np.random.randn()*sigma0*noise_scale+sigma0*offset,
            'q': q,
            'theta_c': gal.core_radius,
            'theta_t': gal.cut_radius,
            'pa': gal.angle_pos*np.pi/180.0+np.pi/2,
            'x1': gal.x_centre,
            'x2': gal.y_centre}
    dpie_scatter.append(piemd(co, **kwargs))

# generate scatter df
number = 2000
theta1 = np.linspace(fov[0],fov[1],number)
theta2 = np.linspace(fov[2],fov[3],number)
theta1, theta2 = np.meshgrid(theta1,theta2)
df_scatter = my_create_deflector(number,fov,dpie_scatter,z_lens,z_source)
kappa_scatter = df_scatter.kappa(theta1,theta2)

# get ref df
df_ref = lst.create_deflector(f'/home/lmx/work/lens-data/{cluster_name}/origin_data/best.par',
                          filex=f'/home/lmx/work/lens-data/{cluster_name}/origin_data/{cluster_name}_angx.fits',
                          filey=f'/home/lmx/work/lens-data/{cluster_name}/origin_data/{cluster_name}_angy.fits')

print('relative kappa difference: max:',np.max(abs(kappa_scatter-df_ref.ka)/df_ref.ka),'min:',np.min(abs(kappa_scatter-df_ref.ka)/df_ref.ka))
print('deflection angle difference: max:',np.max(df_scatter.angx-df_ref.angx),'min:',np.min(df_scatter.angx-df_ref.angx))


fig,ax = plt.subplots(1,2,figsize=(30,15))

ax[0].imshow(kappa_scatter,origin='lower',vmax = 3,vmin=0,extent=fov)
# critical curve
tl = df_scatter.tancl()
cl = df_scatter.getCaustics(tl)
for t in tl:
    x,y = df_scatter.getCritPoints(t)
    ax[0].plot(x,y,'-',color='white')
for c in cl:
    x,y = df_scatter.getCausticPoints(c)
    ax[0].plot(x,y,'-',color='red')
ax[0].set_title('scattered galaxies generated by pyLensLib',fontsize=30)

ax[0].imshow(df_ref.ka,origin='lower',vmax = 3,vmin=0,extent=fov)
# critical curve
tl = df_ref.tancl()
cl = df_ref.getCaustics(tl)
for t in tl:
    x,y = df_ref.getCritPoints(t)
    ax[0].plot(x,y,'-',color='yellow')
for c in cl:
    x,y = df_ref.getCausticPoints(c)
    ax[0].plot(x,y,'-',color='orange')
ax[0].set_title('origin critical and cautics',fontsize=30)
plt.savefig('./compare_scattered_pie.png')
plt.show()