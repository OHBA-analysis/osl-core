### OSL2

OHBA software library

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

	osl_startup
	
