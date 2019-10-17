%% Example of using Workbench to plot surface data
%
% This example shows how to render surface data from volume data using OSL and 
% workbench. First, we need some volume data to plot. For this example, we will
% calculate the alpha power using a parcellation. First, load an SPM object
% and a parcellation

D = spm_eeg_load(fullfile(osldir,'example_data','roinets_example','subject_1.mat'));
p = parcellation(fullfile(osldir,'parcellations','fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'));

%%
% Then, use ROInets to compute parcel timecourses from the data, including
% leakage correction. Don't forget to select the correct montage first
D = D.montage('switch',2);
D = ROInets.get_node_tcs(D,p.parcelflag(true),'PCA');
D = ROInets.remove_source_leakage(D,'symmetric');

%% 
% For each parcel, compute the power spectrum 
V = D(:,:,:).';
for j = 1:size(V,2)
    [spec_power(:,j),f] = pwelch(detrend(V(:,j)),10*D.fsample,[],[],D.fsample);
end

%%
% Now compute the power in the alpha band for each parcel
P_alpha = trapz(f(f>=8&f<=12),spec_power(f>=8&f<=12,:));

%%
% This parcel-wise data can be expanded into a volume, and saved to a .nii file
fname = p.savenii(p.to_vol(P_alpha),fullfile(osldir,'practical','alpha_power'));

%%
% Then we can display the .nii file using osl_render4D
osl_render4D(fname);

%%
% This command will write .gii and CIFTI files from the .nii file, and will open
% them using Workbench. The .gii and CIFTI files are saved in the same directory
% as the input file. You can optionally specify a directory to put them in

% osl_render4D('alpha_power.nii.gz','savedir','test_4d') 