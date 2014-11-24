function [ hmm, block, frenbest ] = run_multistart_hmm( S )

data=S.data;
NK=S.NK;
num_starts=S.num_starts;
frenbest=inf;

if ~isfield(S,'force_zero_means'),
    S.force_zero_means=1;
end;

failed=0;
for st=1:num_starts,
    
    hmm=struct('K',NK);
    hmm=hmminit(data,hmm,'diag'); 
    hmm.obsmodel='Gauss';
    hmm.train.obsupdate=ones(1,hmm.K);    % update observation models ?        
    hmm.train.cyc=40;
    hmm.train.init=1;         % Yes, we've already done initialisation
    % hmm.prior.Dir2d_alpha=ones(size(hmm.prior.Dir2d_alpha));
    % hmm.prior.Dir2d_alpha=hmm.prior.Dir2d_alpha+1000000*diag(diag(hmm.prior.Dir2d_alpha));
    % hmm.prior.Dir2d_alpha=hmm.prior.Dir2d_alpha/100;

    if(S.force_zero_means)
        for ndx=1:length(hmm.state),
            hmm.state(ndx).Mu=0*hmm.state(ndx).Norm_Mu;  
            hmm.state(ndx).constraint=1;
        end;
    end;
    
    try,
        hmm=hmmtrain(data,size(data,1),hmm);
            
        if(st==1 | hmm.train.FrEn(end)<frenbest),
            hmmbest=hmm;
            frenbest=hmm.train.FrEn(end);
        end;

        frens(st)=hmm.train.FrEn(end);
    catch,
        frens(st)=inf;
        failed=failed+1;
    end;

    
end;

if isinf(frenbest),    
   error('All HMMs have failed')   
end;

if(failed)
    warning([num2str(failed) ' out of ' num2str(num_starts) ' have failed']);
end;

[aval aind]=min(frens);

hmm=hmmbest;
clear hmmbest;

hmm.train.cyc=40;
T=length(data);
[block]=hmmdecode(data,T,hmm);

stop=1;
while ~stop,

    NK,
    [y x]=hist(block(1).q_star,1:NK);

    y2=y/sum(y);
    y2

    [aval aind]=min(y2);

    if(0),
    %if(aval<0.1),            

        NK=NK-1;       
        
        if(NK==1)
            stop=1;
        else,
            
            hmm=struct('K',NK);                    
            hmm=hmminit(data,hmm,'diag'); 
            hmm.obsmodel='Gauss';
            hmm.train.obsupdate=ones(1,hmm.K);    % update observation models ?        
            hmm.train.cyc=40;
            hmm.train.init=1;         % Yes, we've already done initialisation

            if(S.force_zero_means)
                for ndx=1:length(hmm.state),
                    hmm.state(ndx).Mu=0*hmm.state(ndx).Norm_Mu;  
                    hmm.state(ndx).constraint=1;
                end;
            end;
            
            hmm=hmmtrain(data,size(data,1),hmm);
        
            T=length(data);
            [block]=hmmdecode(data,T,hmm);


        end;

    else
        stop=1;
    end;
end;


end

