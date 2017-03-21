D = spm_eeg_load(fullfile(osldir,'example_data','roinets_example','subject_1.mat'));
D = D.montage('switch',2);
p = parcellation(fullfile(osldir,'parcellations','fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'));
D = ROInets.get_node_tcs(D,p.parcelflag(true),'PCA');
D = ROInets.remove_source_leakage(D,'symmetric');
V = D(:,:,:).';
figure
plot(V)
for j = 1:size(V,2)
	[spec_power(:,j),f] = pwelch(detrend(V(:,j)),10*D.fsample,[],[],D.fsample);
end
figure
loglog(f,spec_power)
P_alpha = trapz(f(f>=8&f<=12),spec_power(f>=8&f<=12,:))

p.savenii(p.to_vol(P_alpha),'alpha_power');
osl_render4D('alpha_power.nii.gz','savedir','test_4d')

% Or in a single step in the current directory
osl_render4D(p.savenii(p.to_vol(P_alpha),'alpha_power'))