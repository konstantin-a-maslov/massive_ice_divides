# Ice divide reconstuction algorithms (MASSIVE)

This repository contains a set of algorithms implemented/tested/employed within [the MASSIVE project](https://www.mn.uio.no/geo/english/research/projects/massive/index.html).


## Installation

We recommend using the [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) Python distributions. 
After installing one of them, one can use the `conda` package manager to install the required libraries in a new environment called `massive-dem` and activate it by running

```
conda env create -f env.yml
conda activate massive-dem
```


## Kienholz et al., 2013

This script implements a routine for ice divide delineation proposed by [Kienholz et al., 2013](https://doi.org/10.3189/2013JoG12J138). 
Our implementation is open source as it is based upon [WhiteboxTools](https://github.com/jblindsay/whitebox-tools) and [pysheds](https://github.com/mdbartos/pysheds). 
Note, however, that because of the different software utilised for hydrological calculations and all calculations handled in the raster space, the outputs may be different from the original algorithm. 
Make sure that the input rasters are defined on the same grid before running. 
Also, for some erroneous data and bad parameter choices, the script can yield artifacts and provide wrong ice divides. 
Always inspect the outputs visually as their correctness is not guaranteed. 

The manual for the script is provided below, use accordingly. 
Please check the paper for a more detailed description of the parameters. 

```
usage: kienholz.py [-h] [--fill_holes] [--dem_conditioned] [--dem_conditioning_method {fill,breach}]
                   [--dem_smoothing DEM_SMOOTHING] [--pixel_size PIXEL_SIZE]
                   [--buffer1_size BUFFER1_SIZE] [--buffer2_size BUFFER2_SIZE]
                   [--gutter_size GUTTER_SIZE] [--l_gutter L_GUTTER] [--routing {d8,mfd,dinf}]
                   [--sortby {accumulation,elevation}] [--a A] [--b B] [--c C]
                   [--percentile PERCENTILE] [--skip_missing_artifacts] [--a_sliver_1 A_SLIVER_1]
                   [--a_sliver_2 A_SLIVER_2] [--c_sliver C_SLIVER]
                   input_outlines_path input_dem_path output_path

positional arguments:
  input_outlines_path   Path to input outlines .tif file
  input_dem_path        Path to input DEM .tif file
  output_path           Path to output .shp file

options:
  -h, --help            show this help message and exit
  --fill_holes          Fill holes in the outlines
  --dem_conditioned     Skip DEM conditioning
  --dem_conditioning_method {fill,breach}
                        DEM conditioning method
  --dem_smoothing DEM_SMOOTHING
                        DEM smoothing parameter
  --pixel_size PIXEL_SIZE
                        Pixel size in CRS units
  --buffer1_size BUFFER1_SIZE
                        Outer buffer size in CRS units
  --buffer2_size BUFFER2_SIZE
                        Inner buffer size in CRS units
  --gutter_size GUTTER_SIZE
                        Gutter size in pixels
  --l_gutter L_GUTTER   Gutter depth
  --routing {d8,mfd,dinf}
                        Routing method for hydrological calculations
  --sortby {accumulation,elevation}
                        Sorting criteria for pour points
  --a A                 Radii calculation parameter a
  --b B                 Radii calculation parameter b
  --c C                 Maximum radius to merge pour points
  --percentile PERCENTILE
                        Percentile to filter out pour point with low accumulation/elevation
  --skip_missing_artifacts
                        Skip pixels that were missed because of poor calibration/bad data/software
                        errors
  --a_sliver_1 A_SLIVER_1
                        Area threshold for filtering out slivers
  --a_sliver_2 A_SLIVER_2
                        Area threshold for filtering out prolongated slivers
  --c_sliver C_SLIVER   Compactness threshold for filtering out slivers
```

You can also access the manual page of the script by running

```
(massive-dem) python kienholz.py -h
```

To cite the original algorithm, please use the following bib entry.

```
@article{kienholz2013,
    title = {A new semi-automatic approach for dividing glacier complexes into individual glaciers},
    author = {Kienholz, Christian and Hock, Regine and Arendt, Anthony},
    journal = {Journal of Glaciology},
    volume = {59},
    year = {2013},
    month = {10},
    pages = {925-937},
    doi = {10.3189/2013JoG12J138},
}
```


## RGI ice divide copying

**In development.**


## License

This software is licensed under the [GNU General Public License v2](LICENSE).


<br/>

> If you notice any inaccuracies, mistakes or errors, feel free to submit a pull request or kindly email the authors.
