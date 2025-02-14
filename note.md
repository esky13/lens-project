## Steps
- Preliminary steps:
  - Install pyLensLib
  - Install lenstool
  - use lenstool to create maps of the deflection angle components
- Step 1:
  - generate a lens model with lens tool
  - generate random sources in the background of the cluster and find their multiple images (lenstool or pyLensLib)
  - fit the multiple images
    - using lenstool without stellar kinematics priors
    - using lenstool with stellar kinematics priors
- Step 2:

## Notes
### Install lenstool
download compile version from https://projets.lam.fr/projects/lenstool/wiki. Add path to environment:
```
gedit ~/.bashrc
export LENSTOOL_DIR="~/lenstool-7.1-linux64"
export PATH="${LENSTOOL_DIR}:${PATH}"
```
Manual can be found in source code, in the directory `doc`.

### create deflection angle maps
open command line:
```
lenstool filename
```
enter `:wq` after modifying the file. 

### generate random sources
#### with pyLensLib
save the position of images in ```images.cat```, the format is $\{i\ x\ y\ a\ b\ \theta\ z_s\ \mu\}$. Identifier $i$ is $i.j$ representing $j^{th}$ image of $i^{th}$ source.  
#### with lenstool
position of source can be defined in ```source.dat``` following the same format of image, like
```
1 -1.0193879 -2.2174569 0.183664 0.111620 54.34383 1.1424 27.58 
```
then add source parameter in runmode:
```
source    1 source.dat
```
Source can also be added randomly by ```source``` identifiers:
```
source
    random -1
    n_source 100
    elip_max 0.8
    dist_z 1
    z_source 8
    2z_source 9
    end
```
But in this way, many sources have only one image.
### fit multiple images
edit .par file for data to fit like:
```
runmode
    reference     3  0.0 0.0
    dpl       1 2000 1.0 fit_ggls_angx.fits fit_ggls_angy.fits
    poten     1 2000 1.0 fit_ggls_pot.fits 
    image     1 ggls_images.cat
    inverse   3 0.5   # options for fitting
    end
image
    multfile  1  ggls_images.cat
    mult_wcs  1
    sigposArcsec 0.1
    form     0
    end
grid
    number    100
    polar     0
    nlens     3   # number of potential to fit
    nlens_opt 3   # same as above
    end
potential 1
    ...
limit 1
    ...
    v_disp  1 200. 1000.  0.10000   # mark the parameter want to fit with one, and provide limits
    end
...
fini
```
#### with stellar kinematics priors
use potfile to define $\sigma_{ref}$ and $r_{cut}$ for fitting, galaxies defined in file ```cgal.cat```.
```
potfile
	filein	3 cgal.cat
	type 81
	zlens  0.3
	core   0.027
	cutkpc  1 1.0 50.0
	sigma   3 150 50.
	slope 	0 4
	vdsplot 0 4
	mag0 16.7
	end
```

#### result
| parameters | cluster $\sigma_0$ (km/s)  |$\sigma_0^{ref}$(km/s) | $r_{cut}^{ref}$ (kpc) |
| ---------- | --- | --------- | ---------|
| true value | 1411.34 | 157.47 | 2.59 |
| with kinematics prior | 1411.62 | 156.74 (7.54)| 2.61 (0.33)|
| without kinematics prior | 1411.80 | 156.72 (7.55)| 2.61 (0.32)|



- Step 2:

    - use galaxy cluster obtained from high-resolution numerical hydro simulation
    - use the outputs of Subfind to model the scaling relations  for the lenstool code
    - repeat: generating sources & fit model (only with stellar kinematic priors)
  
## simple simulation
use lenstool (+velocity dispersion or not)
Subfind, pyLensLib
## fit other mass distribution

difference between simulation and obs.
how to model matters.
