# SPICE_Alignement

## Description :

This package provide tools to co-align an image with another one (called "reference" image).  You can provide an imager as input, using the Alignment class. You can also provide a SPICE raster, using the specific SpiceAlignment class. 
In addition, the pakcage provides tools to create a synthethic raster corresponding to a given SPICE raster, using a time sequence of imager. The obtained synthetic raster can then be used as the reference image to co-align the SPICE raster.

It is advised to use an imager with a Full Sun field of view as the reference image, where Limb fitting has been previously applied. Example of them include 
the L2 FITS files of the FSI 174 and 304 imagers. 
The co-alignment itself is performed using a cross-correlation tehcnique, through the Pearson's coefficient, or the residus method. The alignement can be done in the following coordinates: 

- Helioprojective coordinates.
- Carrington coordinates. In that case, you have to provide the information to build a pixel grid in a carrington frame.
- by slicing over the image pixels, complitely neglecting the headers informations.

**Warning** As of now, the code works for day to day cases, but has not been thouroughly tested. Please verify the results with the plot_co_alignment method.
Report any bug you encounter with Github or by e-mail to the author. 

## Installation
To create a virtual environement, write the following command in the shell :
```shell
python -m venv env
source env/bin/activate # write "deactivate" in shell to go out of your virtual environement. 
```

You can use the pip install + git command, while you are in virtual environment :

```shell
pip install git+https://github.com/adolliou/SPICE_Alignement
```
or clone the SPICE_alignment repository locally. Then, while you are in your virtual environment:

```shell
cd path/to/SPICE_alignment
pip install .
```



## Usage

We show here a typical example to align SPICE data with a synthetic raster created from FSI 304 files. 

### Creation of a SPICE synthetic raster 
First of all, we need to create a synthetic raster for the SPICE raster using a list of FSI 304 FITS files.
```python
from SPICE_alignment.synras.map_builder import SPICEComposedMapBuilder
from glob import glob
import astropy.units as u


path_spice = "path/to/spice/l2.fits"
path_to_imager_list = glob("path/to/fsi304/*l2.fits")
window_spice = "Ly-gamma-CIII group (Merged)" # The window of the HDULIST for the SPICE FITS file. 
window_imager = -1 # The widow of the HDULIST of the imagers FITS files
threshold_time = u.Quantity(30, "s") # maximum threshold between the SPICE acquisition time, and the closest FSI 304 image. If the code can't any FSI below the threshold, it returns an error 
output_L3_fits = "path/to/output/synthetic_raster_folder"

C = SPICEComposedMapBuilder(path_to_spectro=path_spice, list_imager_paths=path_to_imager_list,
                               window_imager=window_imager, window_spectro=window_spice,
                               threshold_time=threshold_time)
C.process(path_output=output_L3_fits)
```
### Align SPICE raster with the synthetic raster

The code first creates a SPICE pseudo raster by spectrally summing over the chosen HDUList window (here, the C III window). It then cross-correlates the SPICE image with the synthetic raster, while shifting the header values of the SPICE image.
It returns a cross-correlation matrix that can be used to determine the optimal shift to apply to the header values.
The header values that can be shifted are CRVAL1, CRVAL2, CROTA, CDELT1  and CDELT2.

```python
import numpy as np
from SPICE_alignment.hdrshift.alignement_spice import AlignmentSpice
from SPICE_alignment.plot.plot import PlotFunctions
from SPICE_alignment.utils.Util import AlignSpiceUtil


path_to_synthetic_raster_fits = "path/to/input/synthetic_raster.fits"
path_spice_input = "path/to/spice/l2.fits"
window_spice_to_align =  "Ly-gamma-CIII group (Merged)" # the window used for the co-alignment, here the one which includes the C III line.
windows_spice = ["Mg IX 706 - Peak", # The windows where the pointing will be corrected. It is adviced to correct the shift in all of the spectral windows. 
            "Ne VIII 770 - Peak",
            "S V 786 / O IV 787 - Peak",
            "Ly-gamma-CIII group (Merged)",
            "LyB- FeX group (Merged)",
            "O VI 1032 - Peak"] 
window_sr = -1 # the HDULIST index for the synthetic raster. 
path_save_figure= "path/to/output/figures/folder"

param_alignement = {
    "lag_crval1": np.arange(-30, -15, 4), # lag crvals in the headers, in arcsec
    "lag_crval2": np.arange(30, 51, 4),  # in arcsec
    "lag_crota": np.array([0]), # in degrees
    "lag_cdelt1": np.array([0]), # in arcsec
    "lag_cdelt2": np.array([0]), # in arcsec
}

parallelism = True

A = AlignmentSpice(large_fov_known_pointing=path_to_synthetic_raster_fits, small_fov_to_correct=path_spice_input,
                         use_tqdm=True,
                   parallelism=parallelism, counts_cpu_max=10,
                        large_fov_window=-1, small_fov_window=window_sr,
                        path_save_figure=path_save_figure,
                   **param_alignement)

corr = A.align_using_helioprojective(method='correlation')
PlotFunctions.plot_correlation(corr, **param_alignement)
AlignSpiceUtil.write_corrected_fits(path_spice_l2_input=path_spice_input, 
                               path_spice_l2_output="path/where/to/save/corrected/fits", corr=corr,
                                    window_spice_list=windows_spice, 
                                    **param_alignement)

PlotFunctions.plot_co_alignment(large_fov_window=-1, large_fov_path=path_to_synthetic_raster_fits,
                                           corr=corr, small_fov_window= window_spice_to_align, 
                                levels_percentile=[80, 90], 
                                           results_folder=None, small_fov_path=path_spice_input, show=True,
                                           **param_alignement)


```
Example of a results for co-alignment between a SPICE C III image and a FSI 304 synthetic raster, obtained with plot_co_alignment :
![Example of a results for co-alignment between SPICE and FSI 304, from plot_spice_co_alignement](co_alignment_SPICE_FSI.png)


## Credits

- Carrington transform: [F. Auchère](https://github.com/frederic-auchere)
- SPICE utils: [G. Pelouze](https://github.com/gpelouze)
- matrix transform: [F. Auchère](https://github.com/frederic-auchere)

## Contact

Author: Antoine Dolliou (antoine.dolliou@universite-paris-saclay.fr)

## Acknowledgment

If you use this code in your publication, please cite Dolliou et al, 2024 (in prep). 