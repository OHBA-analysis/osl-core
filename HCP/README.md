## OSL HCP tools

OSL is designed to facilitate working with HCP data, either in the same analysis framework but different voxel coordinates, or by conversion into OSL standard brain space. This folder contains functions to convert FieldTrip sourcemodel information into the OSL standard NIFTI representation, and to map data onto OSL standard space. 

### Conversion of sourcemodel

See `HCP_sourcemodel_to_nii`. The workflow for creating OSL standard masks based on sourcemodel structs is provided in `HCP_create_ft_masks`. This file

1. Converts a sourcemodel struct to a binary mask
2. Uses interpolation to convert the OSL standard brain images to the new grid

