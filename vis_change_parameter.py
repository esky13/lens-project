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
from tkinter import Tk, Entry, Button, Label, StringVar, Listbox, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from utils import import_sources, visualize_sources, import_images, calculate_images, visualize_images

def plot_source_image():
    visualize_images(ax, ps_imgx_lst, ps_imgy_lst, ps_shown_lst)
    visualize_sources(ax,ps_sourcex_lst,ps_sourcey_lst,ps_shown_lst)

def plot_lensing():
    global ps_imgx_lst
    global ps_imgy_lst
    ax.clear()
    ax.cla()
    ax.imshow(kappa_scatter,origin='lower',vmax = 3,vmin=0,extent=fov)
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

    ps_imgx_lst,ps_imgy_lst = calculate_images(df_scatter,ps_sourcex_lst,ps_sourcey_lst,ps_sourcez_lst,kwargs)
    plot_source_image()

# update the caustic and critics according to the new z source
def update_zs():
    global z_source
    try:
        z_source = float(entry_zs.get())
        df_scatter.change_redshift(z_source)
        plot_lensing()
    except ValueError:
        print('Please enter a valid redshift.')

def update_sigma0():
    global sigma0
    global kappa_scatter
    global df_scatter
    try:
        prev_sigma0 = sigma0
        sigma0 = float(entry_v.get())
        dpie_scatter = []
        for i in range(len(df)):
            gal = df.loc[i]
            e = 1-np.sqrt((1-gal.ellipticite)/(1+gal.ellipticite))
            q = 1-e
            if len(gal.gal_id) < 3:
                v_disp = gal.v_disp
            else:
                v_disp = gal.v_disp*sigma0/prev_sigma0
            kwargs = {'zl': gal.z_lens,
                    'zs': z_source,
                    'sigma0': v_disp*np.sqrt(3/2),
                    'q': q,
                    'theta_c': gal.core_radius,
                    'theta_t': gal.cut_radius,
                    'pa': gal.angle_pos*np.pi/180+np.pi/2,
                    'x1': gal.x_centre,
                    'x2': gal.y_centre}
            dpie_scatter.append(piemd(co, **kwargs))
        df_scatter = my_create_deflector(grid_number,fov,dpie_scatter,z_lens,z_source)
        kappa_scatter = df_scatter.kappa(theta1,theta2)
        plot_lensing()
    except ValueError:
        print('Please enter a valid dispersion velocity')

def update_rcut0():
    global rcut0
    global kappa_scatter
    global df_scatter
    try:
        prev_rcut0 = rcut0
        rcut0 = float(entry_v.get())
        dpie_scatter = []
        for i in range(len(df)):
            gal = df.loc[i]
            e = 1-np.sqrt((1-gal.ellipticite)/(1+gal.ellipticite))
            q = 1-e
            if len(gal.gal_id) < 3:
                cut_radius = gal.cut_radius
            else:
                cut_radius = gal.cut_radius*rcut0/prev_rcut0
            kwargs = {'zl': gal.z_lens,
                    'zs': z_source,
                    'sigma0': gal.v_disp*np.sqrt(3/2),
                    'q': q,
                    'theta_c': gal.core_radius,
                    'theta_t': cut_radius,
                    'pa': gal.angle_pos*np.pi/180+np.pi/2,
                    'x1': gal.x_centre,
                    'x2': gal.y_centre}
            dpie_scatter.append(piemd(co, **kwargs))
        df_scatter = my_create_deflector(grid_number,fov,dpie_scatter,z_lens,z_source)
        kappa_scatter = df_scatter.kappa(theta1,theta2)
        plot_lensing()
    except ValueError:
        print('Please enter a valid cut radius')

# recorde the source postion        
def on_click(event):
    global ps_imgx_lst
    global ps_imgy_lst
    if event.inaxes == ax:
        x,y = event.xdata, event.ydata
        ps_sourcex_lst.append(x)
        ps_sourcey_lst.append(y)
        ps_sourcez_lst.append(z_source)
        ps_imgx_lst,ps_imgy_lst = calculate_images(df_scatter,ps_sourcex_lst,ps_sourcey_lst,ps_sourcez_lst,kwargs)
        plot_source_image()

# create deflector
def my_create_deflector(grid_number,fov,dpie,z_lens,z_source):
    theta1 = np.linspace(fov[0],fov[1],grid_number)
    theta2 = np.linspace(fov[2],fov[3],grid_number)
    theta1, theta2 = np.meshgrid(theta1,theta2)

    angx = np.zeros((grid_number,grid_number))
    angy = np.zeros((grid_number,grid_number))
    for gal_i in dpie:  
        ax, ay = gal_i.angle(theta1,theta2)
        angx += ax
        angy += ay
    kwargs_def = {'zl': z_lens, 'zs': z_source}
    df = deflector(co, angx=angx, angy=angy, **kwargs_def)
    df.setGrid(theta=np.linspace(fov[0],fov[1],grid_number), compute_potential=False)
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
    zoom_factor = min(max(0.1, zoom_factor), 10)
    update_zoomed_view(mouse_x, mouse_y)

def update_zoomed_view(mouse_x, mouse_y):
    global zoom_factor

    center_x = mouse_x
    center_y = mouse_y

    x_min, x_max, y_min, y_max = fov[0], fov[1], fov[2], fov[3]

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
    kwargs = {'zl': gal.z_lens,
            'zs': z_source,
            'sigma0': gal.v_disp*np.sqrt(3/2),
            'q': q,
            'theta_c': gal.core_radius,
            'theta_t': gal.cut_radius,
            'pa': gal.angle_pos*np.pi/180+np.pi/2,
            'x1': gal.x_centre,
            'x2': gal.y_centre}
    dpie_scatter.append(piemd(co, **kwargs))

# generate scatter df
grid_number = 500    # this parameter =2000 to have the same resolution with deflection angle generated by lenstool
theta1 = np.linspace(fov[0],fov[1],grid_number)
theta2 = np.linspace(fov[2],fov[3],grid_number)
theta1, theta2 = np.meshgrid(theta1,theta2)
df_scatter = my_create_deflector(grid_number,fov,dpie_scatter,z_lens,z_source)
kappa_scatter = df_scatter.kappa(theta1,theta2)

# generate image using scatter lens
ps_imgx_lst = []
ps_imgy_lst = []
ps_imgmu_lst = []
ps_sourcex_lst = []
ps_sourcey_lst = []
ps_sourcez_lst = []
ps_shown_lst = []
kwargs = {
    'n': 1.0,
    'q': 0.8,
    'pa': np.pi/7.0,
    're': 1.0,
    'flux': 1.0,
    'zs': z_source
}

ref_RA, ref_DEC = lst.getRef_RA_DEC(f'/home/lyumx/work/lens-data/{cluster_name}/origin_data/best.par')  
ps_sourcex_lst, ps_sourcey_lst, ps_sourcez_lst = import_sources(f'../lens-data/{cluster_name}/generate_origin_sources/source.dat')
ps_shown_lst = np.ones((len(ps_sourcex_lst)))
# ps_imgx_lst_ref, ps_imgy_lst_ref = import_images(f'../lens-data/{cluster_name}/generate_origin_sources/obs_arcs.dat', ref_RA, ref_DEC)
ps_imgx_lst,ps_imgy_lst = calculate_images(df_scatter,ps_sourcex_lst,ps_sourcey_lst,ps_sourcez_lst,kwargs)

# visualize
zoom_factor = 1.0
zoom_step = 0.1

root = Tk()
root.protocol("WM_DELETE_WINDOW", on_close)

# Add input fields and buttons for zs, v, and r
# zs
label_zs = Label(root, text="Enter zs:")
label_zs.pack()
entry_zs = Entry(root)
entry_zs.insert(0, str(z_source))
entry_zs.pack()
button_zs = Button(root, text="Confirm zs", command=update_zs)
button_zs.pack()

# v
label_v = Label(root, text="Enter sigma0 (km/s):")
label_v.pack()
entry_v = Entry(root)
entry_v.insert(0, str(sigma0))
entry_v.pack()
button_v = Button(root, text="Confirm sigma0", command=update_sigma0)
button_v.pack()

# r
label_r = Label(root, text="Enter rcut0 (arcsecond):")
label_r.pack()
entry_r = Entry(root)
entry_r.insert(0, str(rcut0))
entry_r.pack()
button_r = Button(root, text="Confirm rcut0", command=update_rcut0)
button_r.pack()

fig,ax = plt.subplots(figsize=(10,10))
canvas = FigureCanvasTkAgg(fig, master=root)

ps_imgx_lst,ps_imgy_lst = calculate_images(df_scatter,ps_sourcex_lst,ps_sourcey_lst,ps_sourcez_lst,kwargs)
update_zs()
# visualize_images(ax, ps_imgx_lst_ref, ps_imgy_lst_ref, c='r')

canvas.get_tk_widget().pack()
canvas.mpl_connect('button_press_event',on_click)
canvas.mpl_connect('scroll_event', zoom)

import tkinter as tk

# Function to update the state of a selected item
def toggle_state():
    global ps_shown_lst
    selected = listbox.curselection()  # Get selected item index
    if not selected:
        print("No Source Selection", "Please select a source from the list.")
        return
    index = selected[0]  # Get the first selected index
    ps_shown_lst[index] = 1 if ps_shown_lst[index] == 0 else 0  # Toggle state
    update_listbox()  # Refresh the listbox

# Function to set all states to 1 or 0
def set_all_states(state):
    global ps_shown_lst
    ps_shown_lst = np.ones_like(ps_shown_lst)*state
    update_listbox()  # Refresh the listbox

# Function to refresh the listbox
def update_listbox():
    global ps_shown_lst
    listbox.delete(0, tk.END)  # Clear the listbox
    for i, (x,y,z,shown) in enumerate(zip(ps_sourcex_lst,ps_sourcey_lst,ps_sourcez_lst,ps_shown_lst)):
        listbox.insert(tk.END, f"{'S' if shown == 1 else ' '} source {i}-z: {z:.2f}-({x:.3f},{y:.3f})")
    plot_lensing()

# Create a Listbox to display items
listbox = Listbox(root)
listbox.pack()

# Populate the Listbox with initial data
update_listbox()

# Create buttons for toggling and setting states
toggle_button = Button(root, text="Toggle Selected Source", command=toggle_state)
toggle_button.pack()
enable_all_button = Button(root, text="Visible All", command=lambda: set_all_states(1))
enable_all_button.pack()
disable_all_button = Button(root, text="Invisible All", command=lambda: set_all_states(0))
disable_all_button.pack()

# Start the Tkinter main loop
root.mainloop()
