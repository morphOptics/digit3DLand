# R/digit3DLand: R tool to digitize 3D landmarks on triangular mesh
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

[Nicolas Navarro](http://nnavarro.free.fr) and [Rémi Laffont](http://biogeosciences.u-bourgogne.fr)

[MorphOptics platform](http://biogeosciences.u-bourgogne.fr/fr/services-analytiques/160-imagerie-et-morphometrie-morphoptics.html)

R/digit3DLand is an [R](http://www.r-project.org) package to digitize 3D landmarks 
on triangular meshes. 

This tool digitizes 3D landmarks using a multiresolution approach. First, the landmark is selected roughly on a low resolution mesh 
and then a zoom with the local area around this point is open in full resolution. The landmakrk is taken as the closest vertex on the 
full resolution mesh. A template can be used and is adjusted using at least four landmarks. 

#### Installation

Install R/digit3DLand from its [GitHub repository](https://github.com/nnavarro/digit3DLand).

##### Install prerequisites
Install [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
```

Then install R/digit3DLand 

```r
	require(devtools)
	install_github("morphOptics/digit3DLand", local=FALSE)
```
