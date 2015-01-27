% MEGSIM
%
% Files
%   AB_get_source_timecourses                 - Reconstruct source time courses for vector or scalar OSL beamformer
%   assign_inputs                             - parses the inputs to parent function OSL_simulate_MEG_data
%   check_filename_is_clear                   - checks for and removes existing objects
%   component_topoplot                        - plots a field over channels on a head shape. 
%   create_new_meeg_object                    - creates new meeg object by inheritance
%   do_source_recon_using_osl                 - performs source reconstruction using OSL tools
%   find_MEG_chantypes                        - finds MEG modalities in spm object robustly
%   find_nearest_coordinate                   - finds closest points on mesh to specified locations
%   get_oil_dipole_magnitudes                 - extract dipole magnitudes from source recon
%   get_orientations_from_mesh                - Find orthogonal directions to a surface mesh
%   inverse_affine_transform_points           - computes inverse of affine transformation
%   inverse_transform_dipole_orientation      - applies affine transformation in reverse direction
%   invert_affine_transformation              - inverts an affine transformation matrix
%   isposdef                                  - Test for positive definite matrix.
%   LCMV_beamformer                           - simple LCMV beamformer applied to demeaned, filtered data
%   make_debugging_plots                      - generate plots for debugging purposes
%   megchannels                               - returns indices of meeg channels
%   mwpinv                                    - Pseudoinverse with reduced rank
%   randnorm                                  - Sample from multivariate normal.
%   setup_source_recon_params                 - sets oil structure for source recon pipeline
%   setup_source_recon_params_for_blank_recon - sets oil structure for setting up simulation
%   simulate_MEG_signal                       - calculates a simulated signal in MEG sensors
%   transform_dipole_orientation              - applies affine transformation to vector
