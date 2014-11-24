function x=setup_beamformer_designmatrix(S)

% creates a GLM design matrix
%
% Need S to contain:
% - Xsummary: a parsimonious description of the design matrix 
% contains values Xsummary{reg,cond}, where reg is regressor no. and cond is condition no.  
%
% - trialtypes: vector num_trials long indicating which condition each trial in the MEG data belongs to 
% OR
% - Xsummary is the path to the text file containing the design matrix. That text file
% is num_trials x num_regressors, where the num_trials (and trial order) assumes that the D.badtrials trials have been removed. 
% outputs the design matrix, x, which is num_trials x num_regressors
% - trial_rejects (optional) is a text file containing a list of trial indices ( indexed via the trial order in
% the loaded in design matrix ) to indicate any further trials ( above and beyond the D.badtrials trials ) 
% that you do not want to include in the analysis (e.g. for behavioural reasons), and which will get set to 0 in the design matrix.
%
% MWW 2011

if(isfield(S,'design_matrix')),
    x=S.design_matrix;
    
else,
    
    if(ischar(S.Xsummary)),    

        x=load(S.Xsummary);

        if(isfield(S,'trial_rejects')),
            trial_rejects=load(S.trial_rejects);

            x(trial_rejects,:)=0;   
        end;

    else,

        if(~isfield(S,'trialtypes'))
            if(isstr(S.res))
                tmp=load(S.res);

                res=tmp.res;
                if(~isfield(res,'trialtypes'))
                    res=res{1};
                end;
            else
                res=S.res;

            end;
            S.trialtypes=res.trialtypes;
        end;

        % Set up design matrix 
        x=[];
        for d=1:size(S.Xsummary,2), % indexes regressor

            regtmp=[];          

            for t=1:length(S.Xsummary{d}), % indexes trial type        
                regtmp(find(S.trialtypes==t))=S.Xsummary{d}(t);
            end;

            x(:,d)=regtmp;
        end;

    end;

end;