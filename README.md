# OSL - OHBA SOFTWARE LIBRARY

To set up OSL, follow the directions below:

### Prerequisites

In order to use OSL, you must have FSL already installed - see [https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation). 

When you start up OSL for the first time, a file `osl.conf` will be created in your OSL directory. Among other settings, this file contains a list of locations for components of FSL, specifying where the binaries, libraries, and Matlab utilities are located. By default, OSL will look in some standard locations for FSL and use those if FSL is present. If FSL is in an unexpected location, an error will be raised, and you will need to specify the location of FSL manually. Edit `osl.conf` to specify the directories on your system. 

You can also use this facility to specify which version of FSL you want to use, if you have multiple versions of FSL on your system.

There are some OSL functions that rely on Workbench. As with FSL, after you have run OSL for the first time,  you can edit `osl.conf` to specify the location of Workbench on your system.

##### MATLAB version

Matlab R2014b and newer are fully supported. Matlab R2012b-R2014a have basic testing for processing, but there are known graphics incompatibilities. We have not tested functionality for versions of Matlab below R2012b.

##### SPM12 version

OSL includes its own version of SPM12 because we have made some changes to functions within SPM. Some functions in OSL may not work correctly if you use an unsupported version of SPM. You can test different versions of SPM by changing the `SPMDIR` specified in `osl.conf`

### Getting started

You should have the following directory structure

	- osl
		- osl-core
		- ohba-external
		- HMM-MAR
		- GLEAN
		- MEG-ROI-nets
		- spm12
		- layouts
		- std_masks
		- parcellations

To start using OSL, add the `osl-core` folder to your path, and run `osl_startup`.

We refer to the 'osl' folder above as `OSLDIR`, which you can identify in Matlab by typing `osldir` after initializing OSL.  

To check if there are any compatibility issues, add the `osl-core` folder to your path and run `osl_check_installation`

#### Example data

To run the tutorial examples, download the example data tar file, and extract it. Place the resulting folder inside your OSL directory, so that you have `osl/example_data`

#### Shutting down OSL

Sometimes you may wish to use OSL for part of your pipeline, but then switch back to something else e.g. FieldTrip. When you run `osl_startup`, your Matlab path is automatically backed up. If you run `osl_shutdown`, your path will be restored. 

#### Configuration options

When OSL is first started, a configuration file `osl.conf` will be written. This file defines the following variables

- `FSLDIR` - the root directory for FSL
- `FSLBIN` - the `bin` directory for FSL. On some platforms, this is not simply `FSLDIR/bin`
- `FSLLIB` - the `lib` directory for FSL
- `WORKBENCH` - the root directory for Workbench (i.e. the folder containing either `bin_macosx64` or `bin_linux64`)
- `FREESURFER` - root directory for Freesurfer (currently not used)
- `SPMDIR` - root directory for SPM, use this to specify a version of SPM12 (otherwise, the included `OSLDIR/spm12` will be used)
- `PATH_BACKUP` - automatically set, this is a snapshot of the Matlab path when `osl_startup` is called, which will be restored if `osl_shutdown` is run

#### In-place upgrading

It is possible to update the code by downloading the latest version from GitHub. This can be accomplished by opening a terminal, going into the `osl-core` folder, and running `upgrade.sh`. This will overwrite the code folders with the latest versions from GitHub. 

**Note that your local copies of `osl-core` `ohba-external` `HMM-MAR` `GLEAN` and `MEG-ROI-nets` will be deleted, so any new files you have placed in these folders will be lost**. 

An OSL release consists of a matched set of the files in the code repositories, and the extra content such as `std_masks`, `spm12` etc. These folders will not be affected by the upgrade process. It is possible that there may be incompatibilities between the latest code, and your old copies of these supplementary files. If you encounter unexpected behaviour after upgrading your release, you may instead need to download the complete package again. 

### Troubleshooting

##### FSL paths

Some OSL functions require FSL to be installed. First, make sure you have installed FSL - see [https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation). 

When you start up OSL for the first time, a file `osl.conf` will be created in your OSL directory (i.e. the path returned by the function `osldir()`). Among other settings, this file contains a list of locations for components of FSL, specifying where the binaries, libraries, and Matlab utilities are located. By default, OSL will look in some standard locations for FSL and use those if FSL is present. If FSL is in an unexpected location, an error will be raised, and you will need to specify the location of FSL manually. Edit `osl.conf` to specify the directories on your system. 

You can also use this facility to specify which version of FSL you want to use, if you have multiple versions of FSL on your system.

##### FSLView doesn't start on Mac OS 10.12

Update XQuartz to the latest version

##### FSLVIEW appears offscreen on Mac OS

It's possible for FSLView to appear off the screen, if an external monitor was previously used and is now not present (or vice versa). If the title bar is not visible, it may not be possible to reposition the window, rendering FSLView unusable. We have included an app to fix this. 

1. Open the Mac OS 'System Preferences', click 'Security and Privacy', and click 'Accessibility' in the list on the left side of the window
2. Unlock the preferences by clicking the lock icon in the bottom left of the window (may not be necessary depending on your machine)
3. Press the '+' button to add a new app to the list of permitted apps. Navigate to your `osl2` directory, then go into `util` and select `fix_fslview.app`
4. Open FSL, either by typing `fslview` in the terminal if your system paths are set up correctly, or in Matlab, type

		p = parcellation(8)
		p.fslview

5. Once FSLView is opened (i.e. the icon has appeared in the dock, even if the FSLView window is located offscreen), double click `fix_fslview.app`
6. The FSLView window should now be repositioned on your screen

##### FieldTrip mex errors

To recompile Fieldtrip, start up OSL, and then run

	osl_recompile_fieldtrip

If the problem is in `ft_getopt` it is also fine to just delete the MEX file i.e.

	osl/spm12/external/fieldtrip/utilities/ft_getopt.mexmaci64

On recent versions of Matlab, there is essentially no performance advantage from using the compiled version of this function.

##### SPM crashes or hangs

If you experience unexpected behaviour in SPM e.g. it hangs when you try to load an MEEG object, you may need to recompile the MEX files. To do this, in a terminal go into your SPM directory and then

	cd src
	make distclean
	make && make install
	make external-distclean
	make external && make external-install

This will also recompile FieldTrip
