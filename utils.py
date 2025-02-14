import numpy as np
import matplotlib.pyplot as plt
from pyLensLib.pointsrc import pointsrc
from pyLensLib import lenstool as lst
import math

def import_sources(source_filename):
    ps_x_lst = []
    ps_y_lst = []
    ps_z_lst = []
    with open(source_filename,"r") as source_file:
        Is_starting = True
        xi = []
        yi = []
        z = []
        for line in source_file:
            if line.startswith('#'):
                continue
            data = line.strip().split()
            group = ''.join(filter(str.isdigit, data[0]))
            if Is_starting:
                current_group = group
                Is_starting = False
            if group != current_group:
                current_group = group
                ps_x_lst.append(np.mean(xi))
                ps_y_lst.append(np.mean(yi))
                ps_z_lst.append(np.mean(z))
                xi = []
                yi = []
                z = []
                xi.append(float(data[1]))
                yi.append(float(data[2]))
                z.append(float(data[6]))
            else:
                xi.append(float(data[1]))
                yi.append(float(data[2]))
                z.append(float(data[6]))
    return ps_x_lst, ps_y_lst, ps_z_lst

def import_images(image_filename, ref_RA, ref_DEC):
    with open(image_filename,"r") as image_file:
        Is_starting = True
        xi = []
        yi = []
        ps_imgx_lst = []
        ps_imgy_lst = []
        for line in image_file:
            if line.startswith('#'):
                continue
            data = line.strip().split()
            group = ''.join(filter(str.isdigit, data[0]))
            if Is_starting:
                current_group = group
                Is_starting = False
            if group != current_group:
                current_group = group
                ps_imgx_lst.append(xi)
                ps_imgy_lst.append(yi)
                xi = []
                yi = []
                trans_x, trans_y = convertXY(ref_RA, ref_DEC, float(data[1]), float(data[2]))
                xi.append(trans_x)
                yi.append(trans_y)
            else:
                trans_x, trans_y = convertXY(ref_RA, ref_DEC, float(data[1]), float(data[2]))
                xi.append(trans_x)
                yi.append(trans_y)
    return ps_imgx_lst, ps_imgy_lst

def visualize_sources(ax, ps_x_lst, ps_y_lst, ps_shown_lst):
    cmap = plt.get_cmap('plasma')
    for i, (xi,yi,shown) in enumerate(zip(ps_x_lst,ps_y_lst,ps_shown_lst)):
        if not shown:
            continue
        ax.plot(xi, yi, 'o', color=cmap(i/len(ps_x_lst)), markersize=3)
        ax.text(xi, yi, f'{i}', fontsize=8, ha='right', color='b')

def visualize_images(ax, ps_imgx_lst, ps_imgy_lst, ps_shown_lst, c='k'):
    cmap = plt.get_cmap('plasma')
    for i, (xi,yi,shown) in enumerate(zip(ps_imgx_lst,ps_imgy_lst,ps_shown_lst)):
        if not shown:
            continue
        ax.plot(xi,yi,'o',color=cmap(i/len(ps_imgx_lst)),markersize=3)
        for j, (x,y) in enumerate(zip(xi,yi)):
            plt.text(x, y, f'{i}-{j}', fontsize=8, ha='right', color=c)

def calculate_images(df_scatter,ps_x_lst,ps_y_lst,ps_z_lst,kwargs):
    ps_imgx_lst = []
    ps_imgy_lst = []
    for (x,y,z) in zip(ps_x_lst, ps_y_lst, ps_z_lst):
        kwargs['zs']=z
        ps1 = pointsrc(size=200, Npix=2000, gl=df_scatter, ys1=x, ys2=y, **kwargs)
        xi,yi, mui = ps1.find_images()
        ps_imgx_lst.append(xi)
        ps_imgy_lst.append(yi)
    return ps_imgx_lst, ps_imgy_lst

def convertXY(ref_RA, ref_DEC, x, y):
    x -= ref_RA
    if x > 180: x -= 360
    elif x < -180: x += 360
    x *= -3600*np.cos(math.radians(ref_DEC))
    y -= ref_DEC
    y *= 3600
    return x, y
