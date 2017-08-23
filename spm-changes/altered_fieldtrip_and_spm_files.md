## List of OSL SPM changes

WARNING - `spm-changes` should NOT be added to the path, `osl_startup` will move this files into their required locations.

#### Beamforming addons

SPM12 doesn't come with the `spm-beamforming-toolbox`. This will be copied into SPM. Note that the version in `osl-core` has OSL-specific changes and is thus different to the public distribution of `spm-beamforming-toolbox`.

#### Masking of triggers

For binary masking of triggers to work for elekta data
i.e. corresponds to the `S.fixoxfordneuromag` flag passed into `spm_eeg_convert_4osl` from `osl_convert_script`

	osl/spm_eeg_convert_4osl.m
	spm12/external/fieldtrip/fileio/ft_read_event_4osl.m
	spm12/external/fieldtrip/fileio/private/read_trigger_4osl.m

#### Default option in `spm_eeg_inv_mesh_ui.m`

**Problem** If `sMRI` is empty then the program `spm_eeg_inv_mesh_ui.m` prompts the user for input

**Solution** `spm_eeg_inv_mesh_ui.m` uses template in this case without waiting for user input

#### Incorrect coil mapping in `ft_headmodel_localspheres`

*Problem*: Fieldtrip's implementation of local spheres does not correctly map coils to channels when tra matrix is modified. This was previously solved by reloading the raw tra matrix prior to running the forward model, but this was not robust to splitting the coregistration and forward model stages, since the sensor information is saved into D.inv during coregistration, and was an somewhat cumbersome hack.

**Solution:** save new field coilchan in D.sensors('MEG') in `osl_convert_script.m` (which will in turn be a field in grad in `ft_headmodel_localspheres`). Modified `ft_headmodel_localspheres.m` to use this field instead of `grad.tra` if present.

	osl2.0/spm12/external/fieldtrip/forward/ft_headmodel_localspheres.m

####  TOG correction in `undobalancing.m`

**Problem:** Fieldtrip undoes TOG corrections on CTF data during source recon and data conversion

**Solution:** replace undobalancing.m with a file which takes no action, but issues an informative warning. 

	spm12/external/fieldtrip/forward/private/undobalancing.m
	spm12/external/fieldtrip/plotting/private/undobalancing.m
	spm12/external/fieldtrip/fileio/private/undobalancing.m
	spm12/external/fieldtrip/utilities/private/undobalancing.m
	spm12/external/fieldtrip/private/undobalancing.m

#### Montage bug

_Temporarily disabled in case this has now been fixed in the latest SPM12_

**Problem:** `spm_eeg_montage` has a bug that occurs when the case of strings in D.chantype and D.sensors are inconsistent.

**Solution:** copy patched `spm_eeg_montage` (provided by Vlad) to the SPM12 directory via spm-changes

	osl2.0/spm12/spm_eeg_montage.m

#### gca, gcf fieldtrip bug

_Temporarily disabled to see if these are fixed now (probably already changed upstream)_

**Problem:** topoplot_common, ft_select_range, ft_singleplotER, ft_singleplotTFR still thiink gca and gcf will return a double

**Solution:** copy local copy to the SPM12/fieldtrip directory

