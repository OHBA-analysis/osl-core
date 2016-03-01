function lambda = lambda_update(prior,beta,Y,X)

lambda.b = (prior.b^-1+0.5*(Y'*Y - 2*Y' * X * beta.mu + beta.mu'*X'*X*beta.mu+ ...
        trace(beta.Sigma*X'*X)))^-1;
lambda.c = 0.5*(size(Y,1)*size(Y,2)) + prior.c;

end