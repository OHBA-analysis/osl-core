function beta = beta_update(prior,lambda,Y,X)
%performs the VB update equations for the coefficients of a General Linear
%Model
beta.Sigma = inv(inv(prior.Sigma0) + lambda.b*lambda.c*X'*X);
beta.mu = beta.Sigma*(inv(prior.Sigma0)*prior.mu0 + lambda.b*lambda.c*X'*Y);
end