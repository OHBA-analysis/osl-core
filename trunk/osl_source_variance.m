function V = osl_source_variance(D)
% Computes the temporal variance of data in MEEG object
% V = osl_source_variance(D)
currentMontage = montage(D,'getindex');

D = D.montage('switch');
C = osl_cov(D);

if currentMontage ~= 0
    tra = D.montage('getmontage',currentMontage).tra;
else
    tra = eye(size(C));
end

V = diag(tra * C * tra');

end

