from pyLensLib import lenstool as lst
import numpy as np
import matplotlib.pyplot as plt
from pyLensLib.sersic_numba import sersic
from pyLensLib.pointsrc import pointsrc
import string

cluster_name = 'RXJ2248'

df = lst.create_deflector(f'/home/lmx/work/lens-data/{cluster_name}/origin_data/best.par',
                          filex=f'/home/lmx/work/lens-data/{cluster_name}/origin_data/{cluster_name}_angx.fits',
                          filey=f'/home/lmx/work/lens-data/{cluster_name}/origin_data/{cluster_name}_angy.fits')

fov = lst.getFoV(f'/home/lmx/work/lens-data/{cluster_name}/origin_data/best.par')

tl = df.tancl()
cl = df.getCaustics(tl)
ggsl_cs = df.ggslCrossSection(clt=tl,minsize=0.5,maxsize=3.0)

fig,ax = plt.subplots(1,1,figsize=(10,10))
ax.imshow(df.ka,vmax=2.0,vmin=0.0,origin='lower',extent=fov)
for t in tl:
    x,y = df.getCritPoints(t)
    ax.plot(x,y,'-',color='white')

for c in cl:
    x,y = df.getCausticPoints(c)
    ax.plot(x,y,'-',color='red')

kwargs = {
    'n': 1.0,
    'q': 0.8,
    'pa': np.pi/7.0,
    're': 1.0,
    'flux': 1.0,
    'zs': df.zs
}

# generate random point sources
ps_imgx_lst = []
ps_imgy_lst = []

ref_RA, ref_DEC = lst.getRef_RA_DEC(f'/home/lmx/work/lens-data/{cluster_name}/origin_data/best.par')  
with open(f'../lens-data/{cluster_name}/generate_origin_sources/obs_arcs.dat',"r") as image_file:
    current_group = 1
    xi = []
    yi = []
    for line in image_file:
        if line.startswith('#'):
            continue
        data = line.strip().split()
        group = ''.join(filter(str.isdigit, data[0]))
        if group != current_group:
            current_group = group
            ps_imgx_lst.append(xi)
            ps_imgy_lst.append(yi)
            xi = []
            yi = []
        else:
            xi.append((float(data[1])-ref_RA)*3600)
            yi.append((float(data[2])-ref_DEC)*3600)

cmap = plt.get_cmap('plasma')
for i, (xi,yi) in enumerate(zip(ps_imgx_lst,ps_imgy_lst)):
    ax.plot(xi,yi,'o',color=cmap(i/len(ps_imgx_lst)),markersize=3)
    for j, (x,y) in enumerate(zip(xi,yi)):
        plt.text(x, y, f'{i}-{j}', fontsize=8, ha='right')

ps_x_lst = []
ps_y_lst = []
with open(f'../lens-data/{cluster_name}/generate_origin_sources/source.dat',"r") as source_file:
    current_group = 1
    xi = []
    yi = []
    for line in source_file:
        if line.startswith('#'):
            continue
        data = line.strip().split()
        group = ''.join(filter(str.isdigit, data[0]))
        if group != current_group:
            current_group = group
            ps_x_lst.append(np.mean(xi))
            ps_y_lst.append(np.mean(yi))
            xi = []
            yi = []
        else:
            xi.append(float(data[1]))
            yi.append(float(data[2]))

cmap = plt.get_cmap('plasma')
for i, (xi,yi) in enumerate(zip(ps_x_lst,ps_y_lst)):
    ax.plot(xi,yi,'o',color=cmap(i/len(ps_x_lst)),markersize=3)
    plt.text(xi, yi, f'{i}', fontsize=8, ha='right', color='b')
plt.savefig('/home/lmx/work/lens-project/lensing-image.png')
 