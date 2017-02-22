### OSL2

OHBA software library

### Dependencies

- Matlab R2014b and newer fully supported
- Matlab R2012b-R2014a have basic testing for processing, known graphics incompatibilities

### Setup

For turnkey operation, download a `.zip` distribution release of OSL.

To set up from GitHub, perform the following

- Clone this repository
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
		- spm12
		- layouts
		- std_masks
		- parcellations

Add the `osl2` folder to your path e.g.

	addpath('some_directory/osl2')

Run

	osl2_startup

### Troubleshooting

##### FieldTrip mex errors

To recompile Fieldtrip, start up OSL, and then run

	osl_recompile_fieldtrip

