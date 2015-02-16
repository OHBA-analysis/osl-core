function [D] = wishart_kl (B_q,B_p,alpha_q,alpha_p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [D] = wishart_kl (B_q,B_p,alpha_q,alpha_p)
%
%   computes the divergence 
%                /
%      D(q||p) = | q(x)*log(q(x)/p(x)) dx
%               /
%   between two k-dimensional Wishart propability density for C  given
%   scale matrix B and shape parameter alpha
%
%            alpha
%         |B|              alpha-(k+1)/2
%   p(x)= ------------- |C|              exp (-tr(BC))
%         Gamma (alpha)
%              k
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4,
  error('Incorrect number of input arguements');
end;

if size(B_q)~=size(B_p),
  error('Distributions must have equal dimensions');
end;

K=size(B_p,1);

Lq = -logdet(B_q,'chol') + K * log(2);
Lp = -logdet(B_p,'chol') + K * log(2);

lZq = log(2) * (alpha_q*K/2)  - Lq * (-alpha_q/2) + K*(K-1)/4 * log(pi); 
lZp = log(2) * (alpha_p*K/2)  - Lp * (-alpha_p/2) + K*(K-1)/4 * log(pi); 

Lq = Lq + K * log(2);
Lp = Lp + K * log(2);

for k=1:K
    lZq = lZq + gammaln(alpha_q/2+0.5-0.5*k);
    lZp = lZp + gammaln(alpha_p/2+0.5-0.5*k);
    Lq = Lq + psi(alpha_q/2+0.5-0.5*k);
    Lp = Lp + psi(alpha_p/2+0.5-0.5*k);
end

D = (alpha_q/2-0.5-0.5*K)*Lq - (alpha_p/2-0.5-0.5*K)*Lp ...
    - alpha_q * K / 2 + alpha_q * trace(B_p*inv(B_q)) / 2 + lZp - lZq;

return;

