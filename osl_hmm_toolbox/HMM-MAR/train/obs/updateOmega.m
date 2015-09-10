function hmm = updateOmega(hmm,Gamma,Gammasum,residuals,orders,Tres,XX,XXGXX,XW)

K = length(hmm.state); ndim = size(XW,3);
Sind = hmm.train.Sind==1; S = hmm.train.S==1;
if ~hmm.train.zeromean, Sind = [true(1,size(hmm.state(1).W.Mu_W,2)); Sind]; end
regressed = sum(S,1)>0;
if isempty(orders), 
    order = 0; 
else
    order = orders(end);
end

switch hmm.train.covtype,
    case 'uniquediag'
        hmm.Omega.Gam_rate(regressed) = hmm.prior.Omega.Gam_rate(regressed);
        for k=1:K
            XWk = permute(XW(k,:,:),[2 3 1]);
            e = (residuals(:,regressed) - XWk(:,regressed)).^2;
            swx2 = zeros(Tres,ndim);
            if order>0 || hmm.train.zeromean==0
                for n=1:ndim
                    if ~regressed(n), continue; end
                    tmp = XX(:,Sind(:,n)) * permute(hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)),[2 3 1]);
                    swx2(:,n) = sum(tmp .* XX(:,Sind(:,n)),2);
                end;
            end
            hmm.Omega.Gam_rate(regressed) = hmm.Omega.Gam_rate(regressed) + ...
                0.5 * sum( repmat(Gamma(:,k),1,sum(regressed)) .* (e + swx2(:,regressed) ) );
        end;
        hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres / 2;
        
    case 'diag'
        for k=1:K
            XWk = permute(XW(k,:,:),[2 3 1]);
            e = (residuals(:,regressed) - XWk(:,regressed)).^2;
            swx2 = zeros(Tres,ndim);
            if order>0 || hmm.train.zeromean==0
                for n=1:ndim
                    if ~regressed(n), continue; end
                    tmp = XX(:,Sind(:,n)) * permute(hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)),[2 3 1]);
                    swx2(:,n) = sum(tmp .* XX(:,Sind(:,n)),2);
                end;
            end
            hmm.state(k).Omega.Gam_rate(regressed) = hmm.state(k).prior.Omega.Gam_rate(regressed) + ...
                0.5* sum( repmat(Gamma(:,k),1,sum(regressed)) .* (e + swx2(:,regressed) ) );
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + Gammasum(k) / 2;
        end
        
    case 'uniquefull'
        hmm.Omega.Gam_rate(regressed,regressed) = hmm.prior.Omega.Gam_rate(regressed,regressed);
        for k=1:K
            XWk = permute(XW(k,:,:),[2 3 1]);
            e = (residuals(:,regressed) - XWk(:,regressed));
            e = (e' .* repmat(Gamma(:,k)',sum(regressed),1)) * e;
            swx2 =  zeros(ndim,ndim);
            if order>0 || hmm.train.zeromean==0
                for n1=find(regressed)
                    for n2=find(regressed)
                        if n2<n1, continue, end;
                        index1 = (0:length(orders)*ndim+(~hmm.train.zeromean)-1) * ndim + n1;
                        index2 = (0:length(orders)*ndim+(~hmm.train.zeromean)-1) * ndim + n2;
                        index1 = index1(Sind(:,n1)); index2 = index2(Sind(:,n2));                       
                        swx2(n1,n2) = sum(sum(hmm.state(k).W.S_W(index1,index2) .* XXGXX(Sind(:,n1),Sind(:,n2),k)));
                        swx2(n2,n1) = swx2(n1,n2);
                    end
                end
            end
            hmm.Omega.Gam_rate(regressed,regressed) = hmm.Omega.Gam_rate(regressed,regressed) + (e + swx2(regressed,regressed));
        end
        hmm.Omega.Gam_irate(regressed,regressed) = inv(hmm.Omega.Gam_rate(regressed,regressed));
        hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres;
        
    case 'full'
        for k=1:K
            XWk = permute(XW(k,:,:),[2 3 1]);
            e = (residuals(:,regressed) - XWk(:,regressed));
            e = (e' .* repmat(Gamma(:,k)',sum(regressed),1)) * e;
            swx2 =  zeros(ndim,ndim);
            if order>0 || hmm.train.zeromean==0
                for n1=find(regressed)
                    for n2=find(regressed)
                        if n2<n1, continue, end;
                        index1 = (0:length(orders)*ndim+(~hmm.train.zeromean)-1) * ndim + n1;
                        index2 = (0:length(orders)*ndim+(~hmm.train.zeromean)-1) * ndim + n2;
                        index1 = index1(Sind(:,n1)); index2 = index2(Sind(:,n2));
                        swx2(n1,n2) = sum(sum(hmm.state(k).W.S_W(index1,index2) .* XXGXX(Sind(:,n1),Sind(:,n2),k)));
                        swx2(n2,n1) = swx2(n1,n2);
                    end
                end
            end
            hmm.state(k).Omega.Gam_rate(regressed,regressed) = hmm.state(k).prior.Omega.Gam_rate(regressed,regressed) + ...
                (e + swx2(regressed,regressed));
            hmm.state(k).Omega.Gam_irate(regressed,regressed) = inv(hmm.state(k).Omega.Gam_rate(regressed,regressed));
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + Gammasum(k);
        end
end;

end