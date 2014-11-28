function V = osl_source_variance(D)
% Computes the variance of source reconstructed data
% V = osl_source_variance(D)
dat = D.montage('switch');

tbad = get_bad_sections(D,'logical');
V = var(dat(:,find(~tbad),:),[],2); %#ok - logical indexing didn't work...
V = D.montage('getmontage').tra * V;

end