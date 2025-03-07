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

VISUALIZE_SCATTER = 0

###############################
# generate scattered lens mass
# select sources based on the mass calculated by PyLensLib
###############################

# update the caustic and critics according to the new z source
def update_zs():
    try:
        z = float(entry.get())
        ax.cla()
        ax.imshow(kappa_scatter,origin='lower',vmax = 3,vmin=0,extent=fov)
        df_scatter.change_redshift(z)
        # critical curve
        tl = df_scatter.tancl()
        cl = df_scatter.getCaustics(tl)
        for t in tl:
            x,y = df_scatter.getCritPoints(t)
            ax.plot(x,y,'-',color='white')
        for c in cl:
            x,y = df_scatter.getCausticPoints(c)
            ax.plot(x,y,'-',color='red')
        canvas.draw()
        kwargs['zs'] = z
    except ValueError:
        print('Please enter a valid number.')

# recorde the source postion        
def on_click(event):
    if event.inaxes == ax:
        x,y = event.xdata, event.ydata
        ps1 = pointsrc(size=200, Npix=2000, gl=df_scatter, ys1=x, ys2=y, **kwargs)
        xi,yi, mui = ps1.find_images()
        if len(mui)>1:
            print('images at:',xi,yi)
            ps_imgx_lst.append(xi)
            ps_imgy_lst.append(yi)
            ps_imgmu_lst.append(mui)
            ps_sourcex_lst.append(x)
            ps_sourcey_lst.append(y)
            ps_z_lst.append(kwargs['zs'])

# create deflector
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

def on_close():
    root.quit()
    root.destroy()

# zoom the figure 
def zoom(event):
    global zoom_factor
    mouse_x, mouse_y = event.xdata, event.ydata 

    if mouse_x is None or mouse_y is None:
        return
    if event.button == 'up': 
        zoom_factor += zoom_step
    else: 
        zoom_factor -= zoom_step
    a = max(0.1, zoom_factor)
    b = min(10,zoom_factor)
    zoom_factor = min(max(0.1, zoom_factor),10)
    update_zoomed_view(mouse_x, mouse_y)

def update_zoomed_view(mouse_x, mouse_y):
    global zoom_factor

    center_x = mouse_x
    center_y = mouse_y

    x_min, x_max, y_min, y_max = ax.get_xlim()[0], ax.get_xlim()[1], ax.get_ylim()[0], ax.get_ylim()[1]

    new_width = (x_max - x_min) / zoom_factor
    new_height = (y_max - y_min) / zoom_factor

    new_x_min = center_x - new_width / 2
    new_x_max = center_x + new_width / 2
    new_y_min = center_y - new_height / 2
    new_y_max = center_y + new_height / 2

    ax.set_xlim(new_x_min, new_x_max)
    ax.set_ylim(new_y_min, new_y_max)

    canvas.draw()

# two files are needed:
# best.par: provide information of cosmology and cluster
# members.cat: provide maginitude of the galaxies
# three files are generated:
# scatter.par: record the scatter galaxies
# scatter_ggsl_images/sources.dat: the position of images/sources, images are used for sampling

cluster_name = 'RXJ2248'

# read galaxy and cluster information from this file.
filename = f'/home/lyumx/work/lens-data/{cluster_name}/origin_data/best.par'

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
        noise_scale = 1/30
        offset = 0/10
    kwargs = {'zl': gal.z_lens,
            'zs': z_source,
            'sigma0': gal.v_disp*np.sqrt(3/2)+np.random.randn()*sigma0*noise_scale+sigma0*offset,
            'q': q,
            'theta_c': gal.core_radius,
            'theta_t': gal.cut_radius+np.random.randn()*rcut0*noise_scale+rcut0*offset,
            'pa': gal.angle_pos*np.pi/180+np.pi/2,
            'x1': gal.x_centre,
            'x2': gal.y_centre}
    dpie_scatter.append(piemd(co, **kwargs))

# generate scatter df
number = 2000    # this parameter =2000 to have the same resolution with deflection angle generated by lenstool
theta1 = np.linspace(fov[0],fov[1],number)
theta2 = np.linspace(fov[2],fov[3],number)
theta1, theta2 = np.meshgrid(theta1,theta2)
df_scatter = my_create_deflector(number,fov,dpie_scatter,z_lens,z_source)
kappa_scatter = df_scatter.kappa(theta1,theta2)

# save scatter galaxies to file scatter.par
with open(f'/home/lyumx/work/lens-data/{cluster_name}/scatter.par','w') as f_scatter:
    with open(filename,'r') as f_origin:
        lines = f_origin.readlines()
        case_id = 0
        for line in lines:
            if 'potentiel' in line or 'fini' in line:
                case_id = 1
            match case_id:
                case 1:
                    pass
                case 0:
                    if 'source' in line:
                        f_scatter.write('\tsource 1 scatter_ggsl_sources.dat\n')
                    else:
                        f_scatter.write(line)
            if case_id == 1 and 'end' in line:
                case_id = 0
    for i in range(len(df)):
        gal = dpie_scatter[i]
        f_scatter.write(f'potentiel {df.loc[i].gal_id}\n')
        f_scatter.write(f'\tprofil 81\n')
        f_scatter.write(f'\tx_centre {gal.x1}\n')
        f_scatter.write(f'\ty_centre {gal.x2}\n')
        f_scatter.write(f'\tellipticite {df.loc[i].ellipticite}\n')
        f_scatter.write(f'\tangle_pos {df.loc[i].angle_pos}\n')
        f_scatter.write(f'\tcore_radius {gal.theta_c}\n')
        f_scatter.write(f'\tcut_radius {gal.theta_t}\n')
        f_scatter.write(f'\tv_disp {gal.sigma0*np.sqrt(2/3)}\n')
        f_scatter.write(f'\tz_lens {gal.zl}\n')
        f_scatter.write(f'end\n')
    f_scatter.write('fini')

if VISUALIZE_SCATTER:
    # fit the velocity dipersion-m relation and visualize
    v = []
    l = []
    r = []
    rcore = []
    for i in range(len(dpie_scatter)):
        gal = df.loc[i]
        if len(gal.gal_id) < 3:
            continue
        v.append(dpie_scatter[i].sigma0*np.sqrt(2/3))
        r.append(dpie_scatter[i].theta_t)
        rcore.append(dpie_scatter[i].theta_c)
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
    def v_m_func(l,vdslope,sigma0):
        m0 = 16.176001
        return sigma0*10**(0.4*(m0-l)/vdslope)

    def rcut_m_func(l,slope,rcut0):
        return rcut0*10**(0.4*(m0-l)*2/slope)

    popt_v,_ = curve_fit(v_m_func,l,v)
    popt_rcut,_ = curve_fit(rcut_m_func,l,r)
    print(f'fit velocity slope={popt_v[0]}, sigma0={popt_v[1]}')
    print(f'given velocity slope={vdslope}, sigma0={sigma0}')
    print(f'fit cut radius slope={popt_rcut[0]}, rcut0={popt_rcut[1]}')
    print(f'given cut radius slope={slope}, rcut0={rcut0}')

# generate image using scatter lens
ps_imgx_lst = []
ps_imgy_lst = []
ps_imgmu_lst = []
ps_sourcex_lst = []
ps_sourcey_lst = []
ps_z_lst = []
kwargs = {
    'n': 1.0,
    'q': 0.8,
    'pa': np.pi/7.0,
    're': 1.0,
    'flux': 1.0,
    'zs': z_source
}

# visualize and select sources
zoom_factor = 1.0
zoom_step = 0.1
root = Tk()
root.protocol("WM_DELETE_WINDOW", on_close)

label = Label(root, text="Enter zs:")
label.pack()
entry = Entry(root)
entry.insert(0,str(1.0))
entry.pack()
button = Button(root,text="confirm", command = update_zs)
button.pack()

fig,ax = plt.subplots(figsize=(10,10))
canvas = FigureCanvasTkAgg(fig, master=root)
update_zs()
canvas.get_tk_widget().pack()
canvas.mpl_connect('button_press_event',on_click)
canvas.mpl_connect('scroll_event', zoom)
root.mainloop()

# save image position to file
# x y a b theta z mu
ref_RA, ref_DEC = lst.getRef_RA_DEC(f'/home/lyumx/work/lens-data/{cluster_name}/origin_data/best.par')
with open(f'/home/lyumx/work/lens-data/{cluster_name}/scatter_ggsl_images.dat','w') as file:
    file.write(f"#REFERENCE 3 {ref_RA} {ref_DEC}\n")
    for i, (xi,yi,zs) in enumerate(zip(ps_imgx_lst,ps_imgy_lst,ps_z_lst)):
        for j, (xii,yii) in enumerate(zip(xi,yi)):
            file.write(f"{i+1}.{j+1} {xii:.7f} {yii:.7f} {0.481567004332} {0.481567004332} {90.0} {zs:.4f} {25.0}\n")    
with open(f'/home/lyumx/work/lens-data/{cluster_name}/scatter_ggsl_sources.dat','w') as file:
    file.write(f"#REFERENCE 3 {ref_RA} {ref_DEC}\n")
    for i, (xi,yi,zs) in enumerate(zip(ps_sourcex_lst,ps_sourcey_lst,ps_z_lst)):
        file.write(f"{i+1} {xi:.7f} {yi:.7f} {0.481567004332} {0.481567004332} {90.0} {zs:.4f} {25.0}\n")    
