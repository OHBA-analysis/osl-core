function V = osl_source_variance(D)
% Computes the variance of source reconstructed data
% V = osl_source_variance(D)
currentMontage = montage(D,'getindex');

if currentMontage == 0
    error('No virtual montage applied!')
else  
    D = D.montage('switch');
    C = osl_cov(D);
    V = diag(D.montage('getmontage',currentMontage).tra * C * D.montage('getmontage',currentMontage).tra');
end

end