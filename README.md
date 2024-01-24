# oibeam
This repository contains the description on how to compute the size of the so-called "beam" in Optical Long Baseline Interferometry. It is the size of the central lobe of the impulse response (aka  the dirty beam or the point spread function) and sets the size of a resolution element of the system.

This computation is described in  detail in the follwoing document:
[![there](https://github.com/JMMC-OpenDev/oibeam/blob/gh-pages/beam.svg)](https://github.com/JMMC-OpenDev/oibeam/blob/gh-pages/beam.pdf)
### Notebook
Example of this comptation are shown in the following Julia [notebooks](https://github.com/JMMC-OpenDev/oibeam/blob/notebooks/BeamExample.ipynb).
As Nbviewer cannot show plotly plots properly, it can be directly seen [here](https://jovian.com/ferreols/beamexample). There you can zoom in to compare the image of the beam and the estimated ellipse.
 
### BeamEstimation package:
This repository contains also a small Julia tool shipped in the package `BeamEstimation` to compute this beam. Please look at the notebook how to use this tool. 
To install it:
```julia
using Pkg
Pkg.add(url="https://github.com/JMMC-OpenDev/oibeam.git")
```
