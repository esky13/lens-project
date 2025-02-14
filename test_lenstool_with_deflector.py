from pyLensLib import lenstool as lst
import numpy as np
import matplotlib.pyplot as plt
from pyLensLib.sersic_numba import sersic
from pyLensLib.pointsrc import pointsrc
import string

cluster_name = 'RXJ2248'

df = lst.create_deflector(f'/home/lyumx/work/lens-data/{cluster_name}/origin_data/best.par',
                          filex=f'/home/lyumx/work/lens-data/{cluster_name}/origin_data/{cluster_name}_angx.fits',
                          filey=f'/home/lyumx/work/lens-data/{cluster_name}/origin_data/{cluster_name}_angy.fits')
df2 = lst.create_deflector(f'/home/lyumx/work/lens-data/{cluster_name}/with_prior/20241204/best.par',
                          filex=f'/home/lyumx/work/lens-data/{cluster_name}/with_prior/20241204/{cluster_name}_angx.fits',
                          filey=f'/home/lyumx/work/lens-data/{cluster_name}/with_prior/20241204/{cluster_name}_angy.fits')
fov = lst.getFoV(f'/home/lyumx/work/lens-data/{cluster_name}/origin_data/best.par')
print(fov)
# df = lst.create_deflector('/home/lyumx/work/lens-data/M0416_B22.par',
#                           filex='/home/lyumx/work/lens-data/M0416_B22_angx.fits',
#                           filey='/home/lyumx/work/lens-data/M0416_B22_angy.fits')

# fov = lst.getFoV('/home/lyumx/work/lens-data/M0416_B22.par')
# print(fov)
# calculate galaxy lensing cross section
fig,ax = plt.subplots(2,1,figsize=(10,10))

tl = df.tancl()
cl = df.getCaustics(tl)
ax[0].imshow(df.ka,vmax=2.0,vmin=0.0,origin='lower',extent=fov)
for t in tl:
    x,y = df.getCritPoints(t)
    ax[0].plot(x,y,'-',color='white')
for c in cl:
    x,y = df.getCausticPoints(c)
    ax[0].plot(x,y,'-',color='red')

tl = df2.tancl()
cl = df2.getCaustics(tl)
ax[1].imshow(df.ka,vmax=2.0,vmin=0.0,origin='lower',extent=fov)
for t in tl:
    x,y = df2.getCritPoints(t)
    ax[1].plot(x,y,'-',color='white')
for c in cl:
    x,y = df2.getCausticPoints(c)
    ax[1].plot(x,y,'-',color='red')
# plt.show()
plt.savefig('/home/lyumx/work/lens-project/lensing-critical.png')

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
ps_imgmu_lst = []
ps_sourcex_lst = []
ps_sourcey_lst = []

x_lst = [-26.8,17.01,-37.4,20.79,-22.036,-2.781,-24.6,9.3934,-4.3,3.051]
y_lst = [17.6,-15.64,30.4,33.73,-26.387,32.628,27.33,-36.3891,-1.7,-0.228]
# x_lst = [-26.8,17.01,-37.4,20.79,-22.036,3.8867246,8.7515315,-15.86,3.051,-4.3,-2.781,-24.6,9.3934]
# y_lst = [17.6,-15.64,30.4,33.73,-26.387,-5.0425294,-6.8338928,11.38,-0.228,-1.7,32.628,27.33,-36.3891]
for x,y in zip(x_lst,y_lst):
    ps1 = pointsrc(size=200, Npix=2000, gl=df, ys1=x, ys2=y, **kwargs)
    xi,yi, mui = ps1.find_images()
    print('origin:',xi,yi)
    ps2 = pointsrc(size=200, Npix=2000, gl=df2, ys1=x, ys2=y, **kwargs)
    xi,yi, mui = ps2.find_images()
    print('sampled:',xi,yi)
    ps_imgx_lst.append(xi)
    ps_imgy_lst.append(yi)
    ps_imgmu_lst.append(mui)
    ps_sourcex_lst.append(x)
    ps_sourcey_lst.append(y)

# cmap = plt.get_cmap('Pastel1')
# for i, (xi,yi,mui,x,y) in enumerate(zip(ps_imgx_lst,ps_imgy_lst,ps_imgmu_lst,ps_sourcex_lst,ps_sourcey_lst)):
#     ax.plot(xi,yi,'o',color=cmap(i),markersize=2)
#     ax.plot(x,y,'*',color=cmap(i),markersize=10)

# # plt.show()
# plt.savefig('/home/lyumx/work/lens-project/lensing-image.png')

# # save image position to file
# # x y a b theta z mu

# ref_RA, ref_DEC = lst.getRef_RA_DEC(f'/home/lyumx/work/lens-data/{cluster_name}/origin_data/best.par')
# with open(f'/home/lyumx/work/lens-data/{cluster_name}/my_ggsl_images.dat','w') as file:
# # ref_RA, ref_DEC = lst.getRef_RA_DEC('/home/lyumx/work/lens-data/fit_galaxies/ground_true/bestopt.par')
# # with open('/home/lyumx/work/lens-data/fit_galaxies/without_velocity_prior/ggsl_images.dat','w') as file:
#     file.write(f"#REFERENCE 3 {ref_RA} {ref_DEC}\n")
#     for i, (xi,yi) in enumerate(zip(ps_imgx_lst,ps_imgy_lst)):
#         for j, (xii,yii) in enumerate(zip(xi,yi)):
#             file.write(f"{i+1}.{j+1} {xii:.7f} {yii:.7f} {0.481567004332} {0.481567004332} {90.0} {kwargs['zs']:.4f} {25.0}\n")    
# with open(f'/home/lyumx/work/lens-data/{cluster_name}/my_ggsl_sources.dat','w') as file:
#     file.write(f"#REFERENCE 3 {ref_RA} {ref_DEC}\n")
#     for i, (xi,yi) in enumerate(zip(ps_sourcex_lst,ps_sourcey_lst)):
#         file.write(f"{i+1} {xi:.7f} {yi:.7f} {0.481567004332} {0.481567004332} {90.0} {kwargs['zs']:.4f} {25.0}\n")    