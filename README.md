# Alignement

## description :

Alignement functions for Imagers, SPICE or EIS with imagers or synthetic rasters of imagers. 
Co-align image using the Pearson's coefficient or the residus method. The alignement can be performed with both images: 

- in helioprojective coordinates (Alignement)
- in carrington coordinates
- by slicing over the large FOV image pixels, complitely neglecting the headers

The code also allows to create synthetic rasters of SPICE or EIS image, using a list of imager paths.

## installation
either use pip install + git in your virtual environment :

```shell
pip install git+https://github.com/adolliou/Alignement
```
to create a virtual environement do in your local python folder :
```shell
python -m venv env
source env/bin/activate # write "deactivate" in shell to go out of your virtual environement. 
```
or locally while in the local folder:

```shell
cd path/to/git/repo
pip install .
```


## Usage

See the test folder for multiple examples.
We show here a typical example to align SPICE data with a synthetic raster created from FSI 304 files. 

### Creation of a SPICE synthetic raster 
First of all, we need to create a synthetic raster of the SPICE raster using the reference SPICE raster and a list of the imager fits files. 
```python
from SPICE_alignment.synras.map_builder import SPICEComposedMapBuilder
from glob import glob
import astropy.units as u


path_spice = "path/to/spice/l2.fits"
path_to_imager_list = glob("path/to/fsi304/*l2.fits")
window_spice = "Ly-gamma-CIII group (Merged)" # int or str: window in the HDUList used.depends on the alignmenet you want to do. 
window_imager = -1 # same for imagers in imager_list
threshold_time = u.Quantity(30, "s") # maximum threshold time you want
output_L3_fits = "path/to/output/synthetic_raster_folder"

C = SPICEComposedMapBuilder(path_to_spectro=path_spice, list_imager_paths=path_to_imager_list,
                               window_imager=window_imager, window_spectro=window_spice,
                               threshold_time=threshold_time)
C.process(path_output=output_L3_fits)
```
### Align SPICE raster with the created synthetic raster

Create a SPICE pseudo raster by spectrally summing over the chosen HDUList window. Then perform a cross-correlation algorithm by shifting the headers values.
It returns a cross-correlation matrix that can be used to determine the optimal shift to apply to the header values. The header values that can be shifted are CRVAL1, CRVAL2, CROTA, CDELT1 (not recommended as of now) and CDELT2 (not recommanded as of now).

```python
import numpy as np
from SPICE_alignment.plot.plot import PlotFunctions
from SPICE_alignment.utils.Util import SpiceUtil
from SPICE_alignment.hdrshift.alignement_spice import AlignmentSpice


path_to_synthetic_raster_fits = "path/to/input/synthetic_raster.fits"
path_spice_input = "path/to/spice/l2.fits"
window_spice = "Ly-gamma-CIII group (Merged)" # int or str: window in the HDUList used.depends on the alignmenet you want to do. 
window_sr = -1 # same for imagers in imager_list
path_save_figure= "path/to/output/figures/folder"

lag_crval1 = np.arange(-30, -15, 4) # lag crvals in the headers, in arcsec
lag_crval2 = np.arange(30, 51, 4)  # in arcsec
lag_crota = np.array([0]) # in degrees
lag_cdelt1 = np.array([0]) # in arcsec
lag_cdelt2 = np.array([0]) # in arcsec
parallelism = True

A = AlignmentSpice(large_fov_known_pointing=path_to_synthetic_raster_fits, small_fov_to_correct=path_spice_input,
                        lag_crval1=lag_crval1, lag_crval2=lag_crval2, lag_crota=lag_crota, use_tqdm=True,
                        lag_cdelta1=lag_cdelt1, lag_cdelta2=lag_cdelt2, parallelism=parallelism,
                        large_fov_window=-1, small_fov_window=window_sr,
                        path_save_figure=path_save_figure,)

corr = A.align_using_helioprojective(method='correlation')
PlotFunctions.plot_correlation(corr, lag_crval1, lag_crval2,)
SpiceUtil.write_corrected_fits(path_spice_l2_input=path_spice_input, 
                               path_spice_l2_output="path/where/to/save/corrected/fits", lag_crval1=lag_crval1, 
                               lag_crval2=lag_crval2, corr=corr, window_spice=window_spice)

PlotFunctions.plot_spice_co_alignment(imager_window=-1, large_fov_fits_path=path_to_synthetic_raster_fits,
                                           corr=corr, raster_window= window_spice, levels_percentile=[80, 90],
                                           results_folder=None, spice_raster_path=path_spice_input, show=True,
                                           lag_crval1=lag_crval1, lag_crval2=lag_crval2)


```
Example of a results for co-alignment between SPICE and FSI 304, from plot_spice_co_alignement:
![Example of a results for co-alignment between SPICE and FSI 304, from plot_spice_co_alignement](co_alignment_SPICE_FSI.png)


## credits

- carrington transform: [F. Auchère](https://github.com/frederic-auchere)
- SPICE utils: [G. Pelouze](https://github.com/gpelouze)
- matric transform: [F. Auchère](https://github.com/frederic-auchere)

## Contact

Author: Antoine Dolliou (antoine.dolliou@universite-paris-saclay.fr)

## Acknowledgment

Please acknowledge Dolliou et al, 2024 (in prep)