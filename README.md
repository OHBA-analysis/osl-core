### OSL2

OHBA software library

### Dependencies

- Matlab R2014b and newer fully supported
- Matlab R2012b-R2014a have basic testing for processing, known graphics incompatibilities

### Setup

For turnkey operation, download a `.zip` distribution release of OSL.

To set up from GitHub, perform the following

- Clone this repository, using the branch workshop
- Clone `https://github.com/OHBA-analysis/ohba-external`
- Clone `https://github.com/OHBA-analysis/HMM-MAR`
- Clone `https://github.com/OHBA-analysis/GLEAN`
- Clone `https://github.com/OHBA-analysis/MEG-ROI-nets`
- Download `spm12` [here](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- Download OSL standard masks etc. from the OSL wiki

Your directory structure should then look like

	- some_directory
		- osl2
		- ohba-external
		- HMM-MAR
		- GLEAN
		- MEG-ROI-nets
		- spm12
		- layouts
		- std_masks
		- parcellations

Add the `osl2` folder to your path e.g.

	addpath('some_directory/osl2')

Run

	osl2_startup

### Troubleshooting

##### FSL

Some OSL functions require FSL to be installed. First, make sure you have installed FSL - see [https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation). 

When you start up OSL for the first time, a file `fsl_location.txt` will be created in your OSL directory (i.e. the path returned by the function `osldir()`). This file contains a list of locations for components of FSL, specifying where the binaries, libraries, and Matlab utilities are located. By default, OSL will look in some standard locations for FSL and use those if FSL is present. If FSL is in an unexpected location, an error will be raised, and you will need to specify the location of FSL manually. Edit `fsl_location.txt` to specify the directories on your system. 

You can also use this facility to specify which version of FSL you want to use, if you have multiple versions of FSL on your system.

##### FieldTrip mex errors

To recompile Fieldtrip, start up OSL, and then run

	osl_recompile_fieldtrip

