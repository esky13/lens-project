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

# update the caustic and critics according to the new z source
def update_zs():
    try:
        z = float(entry.get())
        ax.cla()
        ax.imshow(df.ka,origin='lower',vmax = 3,vmin=0,extent=fov)
        df.change_redshift(z)
        # critical curve
        tl = df.tancl()
        cl = df.getCaustics(tl)
        for t in tl:
            x,y = df.getCritPoints(t)
            ax.plot(x,y,'-',color='white')
        for c in cl:
            x,y = df.getCausticPoints(c)
            ax.plot(x,y,'-',color='red')
        canvas.draw()
        kwargs['zs'] = z
    except ValueError:
        print('Please enter a valid number.')

# recorde the source postion        
def on_click(event):
    if event.inaxes == ax:
        x,y = event.xdata, event.ydata
        ps1 = pointsrc(size=200, Npix=2000, gl=df, ys1=x, ys2=y, **kwargs)
        xi,yi, mui = ps1.find_images()
        if len(mui)>1:
            print('images at:',xi,yi)
            ps_imgx_lst.append(xi)
            ps_imgy_lst.append(yi)
            ps_imgmu_lst.append(mui)
            ps_sourcex_lst.append(x)
            ps_sourcey_lst.append(y)
            ps_z_lst.append(kwargs['zs'])

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


cluster_name = 'RXJ2248'


df = lst.create_deflector(f'/home/lyumx/work/lens-data/{cluster_name}/generate_images/scatter.par',
                          filex=f'/home/lyumx/work/lens-data/{cluster_name}/generate_images/{cluster_name}_angx.fits',
                          filey=f'/home/lyumx/work/lens-data/{cluster_name}/generate_images/{cluster_name}_angy.fits',
                          zl=0.348)
fov = lst.getFoV(f'/home/lyumx/work/lens-data/{cluster_name}/generate_images/scatter.par')

z_source = 1.5
z_lens = df.zl

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
entry.insert(0,str(z_source))
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
