## List of OSL SPM changes

#### Beamforming addons

`bf_*` in this directory are placed in `spm12/toolbox/spm-beamforming-toolbox`

### Fixes

#### Masking of triggers

For binary masking of triggers to work for elekta data
i.e. corresponds to the `S.fixoxfordneuromag` flag passed into `spm_eeg_convert_4osl` from `osl_convert_script`

	osl2.0/osl/spm_eeg_convert_4osl.m
	osl2.0/spm12/external/fieldtrip/fileio/ft_read_event_4osl.m
	osl2.0/spm12/external/fieldtrip/fileio/private/read_trigger_4osl.m


#### Mex error in `ft_multiplotER`

**Problem:**  Calls to `ft_multiplotER` were producing an error in calling the mex file `spm12/external/fieldtrip/utilities/private/mxSerialize`

**Solution:** inserted all `fieldtrip/utilities/private/*.mexmaci6`4 files from spm12b (beta) 

	osl2.0/spm12/external/fieldtrip/utilities/private/*.mexmaci64

**Update:** This fix did not work me (Adam) and should therefore not be in place as a general *solution.* Instead I would recommend reinstalling SPM12 and recompiling the fieldtrip mex files using `ft_compile_mex(true)`. This eventually worked for me without having to copy any SPM12b files. The original *solution *has been removed.

#### Incorrect coil mapping

*Problem*: Fieldtrip's implementation of local spheres does not correctly map coils to channels when tra matrix is modified. This was previously solved by reloading the raw tra matrix prior to running the forward model, but this was not robust to splitting the coregistration and forward model stages, since the sensor information is saved into D.inv during coregistration, and was an somewhat cumbersome hack.

**Solution:** save new field coilchan in D.sensors('MEG') in `osl_convert_script.m` (which will in turn be a field in grad in `ft_headmodel_localspheres`). Modified `ft_headmodel_localspheres.m` to use this field instead of `grad.tra` if present.

	osl2.0/spm12/external/fieldtrip/forward/ft_headmodel_localspheres.m

#### Badsamples changes

**Problem:** `badsamples.m` has been updated on the SPM12 SVN, this change may not have been updated.

**Solution:** copy local copy to the SPM12 directory

	osl2.0/spm12/@meeg/badsamples.m

####  TOG correction

**Problem:** Fieldtrip undoes TOG corrections on CTF data during source recon and data conversion

**Solution:** replace undobalancing.m with a file which takes no action, but issues an informative warning. 

	osl2.0/spm12/external/fieldtrip/forward/private/undobalancing.m
	osl2.0/spm12/external/fieldtrip/plotting/private/undobalancing.m
	osl2.0/spm12/external/fieldtrip/fileio/private/undobalancing.m
	osl2.0/spm12/external/fieldtrip/utilities/private/undobalancing.m
	osl2.0/spm12/external/fieldtrip/private/undobalancing.m

#### Path function problem

**Problem:** path.m does not return "this" when the path is set. 

**Solution:** copy local correct copy to the SPM12 directory via spm-changes

	osl2.0/spm12/@meeg/path.m

#### Montage bug

**Problem:** `spm_eeg_montage` has a bug that occurs when the case of strings in D.chantype and D.sensors are inconsistent.

**Solution:** copy patched `spm_eeg_montage` (provided by Vlad) to the SPM12 directory via spm-changes

	osl2.0/spm12/spm_eeg_montage.m

#### Subsref bug

**Problem:** @meeg/subsref.m has been updated on the SPM12 SVN to fix a bug when indexing badsamples when online montages are used.

**Solution:** copy local copy to the SPM12 directory

	osl2.0/spm12/@meeg/subsref.mm

#### gca, gcf fieldtrip bug

**Problem:** topoplot_common, ft_select_range and ft_singleplotER still thiink gca and gcf will return a double

**Solution:** copy local copy to the SPM12/fieldtrip directory

