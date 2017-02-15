### Parcellation object

A parcellation maps voxels to one or more parcels. They can be represented in a variety of equivalent (or semi-equivalent) ways. Different analyses may expect different input formats. This class provides a simple interface for common operations involving parcellations, including:

- Loading in parcels in a range of different input formats, and associating them with an input object
- Visualizing the parcellation
- Representing the parcellation in a standard format for OSL

Example usage

```matlab
	p = parcellation('parcellation_with_PCC_reduced8mm.nii') % Load the parcellation
	p.fslview % Plot the parcellation using fslview
	p.to_matrix(p.weight_mask) % Return the parcellation matrix used by osl_apply_parcellation()
	p.parcelflag % Return the parcelflag variable required by ROInets
```

#### Constructing a parcellation

Use `parcellation(input_mask)` where input mask

- Is the name of a NIFTI file
- Is the name of a Matlab file containing a variable 'mask'
- Is a matrix

The input formats supported are

- XYZ x Parcels (4D matrix)
- Vox x Parcels (2D matrix)
- Vox x 1 (1D matrix, value is parcel)
- XYZ x 1 (3D matrix, value is parcel)

Internally, all of these are converted to XYZ x Parcels. 

Whether the parcellation is weighted or binary is determined automatically. A parcellation is _not_ weighted if either

- All of the values in the input mask are `1` or `0`
- All of the values in the input mask are integers _and_ the input mask is either `Vox x 1` or `XYZ x 1` (if the values are not integers, it is assumed that this is instead a parcellation containing a single weighted parcel)

A parcel is overlapping if the weight mask has more than one parcel assigned to a voxel. These properties are computed automatically, so you can use the parcellation object to inspect these properties of a parcellation contained within a file.

As part of the construction, the `parcellation` object will try and guess which standard mask corresponds to the input mask, based on the size of the input. This generally seems to work surprisingly well. The template mask, the voxel coordinates, and the mask filename, are all stored in the parcellation as well. You can quickly load a standard mask just by specifying the resolution e.g.

    p = parcellation(8)

Will load in an 8mm mask, where `p.template_mask` has the greyscale image, and `p.weight_mask` is a binary mask. 

#### Common tasks

##### Reshaping

To convert XYZ to a matrix (e.g. for 8mm, `23×27×23` to `3559x1`) use

	mat = p.to_matrix(vol)

For the reverse

    vol = p.to_vol(mat)

In the latter case, the matrix could be `voxels x frame` or `parcels x frame` (or `frame x parcels` - if the second dimension matches the number of parcels, the matrix will be transposed). 

##### Binarizing

To remove weights, use

	weight_mask = p.remove_weights()

To remove overlap, use

	weight_mask = p.remove_overlap()

To do both, use

	weight_mask = p.binarize()

Note that because `remove_overlap` assigns voxels to parcels based on their weights, this needs to be done prior to removing the weights. Thus `binarize()` is equivalent to `remove_overlap()` followed by `remove_weights()`, in that order. 