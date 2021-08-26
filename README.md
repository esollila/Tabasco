# Tabasco MATLAB toolbox
Tabasco toolbox provides an easy to use function, `tabasco.m`, to compute the TABASCO (TApered or BAnded Shrinkage COvariance) estimator. We also include a demo example that reproduces Figure 3 (upper left corner plot) in our paper.

### Compatibility

The code is tested on MATLAB R2020a, but should work on other versions of Matlab with no or little changes as well as all platforms (Windows 64-bit, Linux 64-bit, or Mac 64-bit).

### Dependency

This toolbox does not depend on any other toolbox. However, the demo example `tabasco_demo.m` uses RegularizedSCM toolbox for computing the Ledoit-Wolf estimator (LWE) proposed by Ledoit and Wolf (2004). 

Hence, in order to run `tabasco_demo.m` please install the RegularizedSCM toolbox from http://users.spa.aalto.fi/esollila/regscm/

### Installation for MATLAB version >= 2014b

Simply double click the `Tabasco.mltbx`  and the MATLAB® installs the toolbox. If it does not work follow the instructions below for installation for MATLAB version < 2014b.

### Installation for MATLAB version < 2014b

Download the files to a local folder, called _Tabasco_ (say). 
Add the folder to the MATLAB search path as follows.
-  Start MATLAB and go to the _Tabasco_ folder, and execute the lines:

> `addpath(pwd) %<-- Add the toolbox to the MATLAB path`

> `save path  %<-- Save the path %`

### How to cite

If you use this toolbox, please cite the the publication:

[1] Esa Ollila and Arnaud Breloy, "Regularized Tapered or Banded Sample Covariance Matrix", ArXiv preprint, submitted for publication, Sept. 2021.

Now you are good to go!

### Getting started

- Simply type  `help tabasco`

- To have a quick acces to the provided demo example in the toolbox simply type `demo Tabasco`

### Authors
* Esa Ollila, esa DOT ollila AT aalto.fi
* Arnaud Breloy, abreloy AT parisnanterre DOT fr
