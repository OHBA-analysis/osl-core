function osl_hmm_plotstatepath(hmm,options)
% Plots the HMM statepath
% osl_hmm_plotstatepath(hmm)
% osl_hmm_plotstatepath(hmm,options)
% AB 2013

if nargin < 2
  options=[];
end;

try 
    options = ft_checkopt(options,'mode','char',{'separate','overlay'});
catch
    options.mode='separate';
    disp(['options.mode not set correctly. Defaulting to options.mode=''' options.mode '''']);
end; % stack all states on separate lines on y-axis or overlay
mode=options.mode; 

try t0=options.t0; catch t0 = 0; end; % start time point in secs
try bc=options.bc; catch bc = []; end; % baseline correct epoched states using this time range
try win=options.win; catch win = []; end; % smoothing 
try epoched_statepath_sub=options.epoched_statepath_sub; catch epoched_statepath_sub = []; end; 

if ~isfield(hmm,'fsample') || isempty(hmm.fsample)
    fs = 1;
else
    fs = hmm.fsample;
end

cla(gca); hold(gca,'on')

k = max(hmm.statepath);

col = colormap(gca,'lines');

leg=[];
for s = 1:k
  
  if isempty(epoched_statepath_sub)
    tmp=double(logical(hmm.statepath==s));
  else
    count=0;
    tmp=zeros(1,size(epoched_statepath_sub{1},2));
    for sub = 1:length(epoched_statepath_sub)
        stateinds=double(logical(epoched_statepath_sub{sub}==s));
        tmp=tmp+sum(stateinds,3);
        count=count+size(stateinds,3);
    end;
    tmp=tmp/count;
  end;
  
  if ~isempty(win)
    FO = conv(tmp,rectwin(fs*win),'same')./ (fs*win);
  else
    FO=tmp;
  end;
  
  ts=linspace(t0,t0+length(FO)/fs,length(FO));

  if ~isempty(bc)
    FO=FO-mean(FO(intersect(ts>bc(1),ts < bc(2))))
  end;
  
  if strcmp(mode,'separate'),      
    plot(ts,0.8*FO + s - 0.4,'color',col(s,:));
    xlabel('Time (s)'); ylabel('State #')
  else
    plot(ts,FO,'color',col(s,:));  
    xlabel('Time (s)'); ylabel('FO')
    leg=[leg; 's' num2str(s)];
    
  end;  
  
  %plot((1:length(hmm.statepath))/fs, 0.8*double(hmm.statepath==s)+ s - 0.4,'color',col(s,:));
end

legend(leg);

if strcmp(mode,'separate'),
    ylabel('State #')
 
    set(gca,'ylim',[0 hmm.K]+0.5)
else
    ylabel('FO')
end;

if ~isfield(hmm,'fsample') || isempty(hmm.fsample)
  xlabel('Samples');    
else
  xlabel('Time (s)');
end
 


end
