### OSL2

OHBA software library

### Dependencies

- Matlab R2014b and newer fully supported
- Matlab R2012b-R2014a have basic testing for processing, known graphics incompatibilities

### Setup

For turnkey operation, download a `.zip` distribution release of OSL.

To set up from GitHub, perform the following

- Clone this repository
- Clone `https://github.com/OHBA-analysis/osl-external`
- Download `spm12` [here](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- Download OSL standard masks etc. from the OSL wiki

Your directory structure should then look like

	- some_directory
		- osl2
		- osl-external
		- spm12
		- layouts
		- std_masks

Add the `osl2` folder to your path e.g.

	addpath('some_directory/osl2')

Run

	osl2_startup
	
