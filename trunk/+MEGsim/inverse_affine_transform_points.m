function new = inverse_affine_transform_points(M, old)
%INVERSE_AFFINE_TRANSFORM_POINTS computes inverse of affine transformation

% Giles Colclough 28 October 2013

Minv = MEGsim.invert_affine_transformation(M);

new  = spm_eeg_inv_transform_points(Minv, old);