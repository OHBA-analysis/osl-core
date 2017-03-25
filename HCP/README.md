## OSL HCP tools

OSL is designed to facilitate working with HCP data, either in the same analysis framework but different voxel coordinates, or by conversion into OSL standard brain space. This folder contains functions to convert FieldTrip sourcemodel information into the OSL standard NIFTI representation, and to map data onto OSL standard space. 

### Working with HCP data in OSL

The HCP provides lead fields but not structure, which means that the HCP data can _only_ be beamformed onto standard masks from FieldTrip. Although these masks have the same grid resolution as some of the standard masks, they have different sizes and also different offsets. 

Starting with the raw HCP data, it is first beamformed onto one of the FieldTrip masks - these are in the `std_masks` folder with filenames `ft_*.nii.gz`. From this point, there are several options

1. Apply standard parcellation pipelines, but using parcellations defined on the FieldTrip mask
2. Convert the data onto one of the OSL standard masks, and then use existing parcellations

For the first option, you need a parcellation defined on a FieldTrip mask. There are some parcellations already defined on this grid. However, if you only have a parcellation on an OSL standard mask, it can be resampled onto the FieldTrip grid

	osl_resample_nii_matlab('original_parcellation.nii.gz','ft_8mm.nii.gz','new_parcellation','interptype','nearest')

This will

