function V = osl_source_variance(D)
% Computes the variance of source reconstructed data
% V = osl_source_variance(D)
dat = D.montage('switch');

tbad = osl_bad_sections(D,'logical');
C = cov(dat(:,find(~tbad),:)'); %#ok - logical indexing didn't work...
V = diag(D.montage('getmontage').tra * C * D.montage('getmontage').tra');

end