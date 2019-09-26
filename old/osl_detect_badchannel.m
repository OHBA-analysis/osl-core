function [list_bch]=osl_detect_badchannel(S)

modality=S.modality;
rel_thres=S.rel_thres;
Dx=S.D;

ica_nonlinearity='tanh';

Fsnew=S.fs_to_use;
Fs=fsample(Dx);


fname=fnamedat(Dx);
fname=fname(1:end-4);

S = [];
S.D = [path(Dx) filesep fname];
S.fsample_new = Fsnew;
Dd = spm_eeg_downsample(S);

triggers=events(Dx);
if isempty(triggers)
    timescale_dec=[1:size(Dd,2)];
else
    t_in=triggers(1).time;
    t_fin=triggers(end).time;
    timescale_dec=[ceil(Fsnew*t_in):fix(Fsnew*t_fin)];
end


list_chan=find(strcmp(chantype(Dd),modality));
data=Dd(list_chan,timescale_dec);

[nchan leng]=size(data);

% for i=1:nchan
%     data(i,:)=(data(i,:)-mean(data(i,:)))/std(data(i,:));
% end

data=1000*data/max(data(:));

timeline=1/Fsnew*[timescale_dec(1):timescale_dec(1)+leng-1];

flag=0;

cont=0;

list_gch=[1:nchan];

while flag==0
    
    clear A;
    clear W;
    clear IC;
    clear power;
    clear order;
    
    cont=cont+1;
    
    sigs=data(list_gch,:);
    disp(['Unworking channel identification: ICA run ' num2str(cont)]);
    randn('state',1);
    [IC, A, W] = fastica(sigs,'approach','defl','numOfIC',nchan,'g',ica_nonlinearity,'verbose','on','displayMode','off','maxNumIterations',75);
    
    [nchan,Nc]=size(A);
    
    [Nc,leng]=size(IC);
    
    if Nc > 0
        % order ICs by var(A_i) without demeaning - WHY are we doing this?
        clear power;
        for i=1:Nc
            somma=0;
            for k=1:nchan
                somma=somma+A(k,i)^2;
            end
            power(i)=somma/nchan;
        end
        
        [power,order]=sort(power);
        
        power=fliplr(power);
        
        order=fliplr(order);
        
        A=A(:,order);
        
        W=W(order,:);
        
        IC=IC(order,:);
        
    end
    
    % ICA gives A nchans x nIC 
    %   IC: nIC x ntpts
    badchannels=[];
    list_ic=[];
    
    for i=1:Nc % loop over components
        [vt,order]=sort(abs(A(:,i))); % sort by magnitude of each chan in this component
        vt=flipud(vt);
        order=flipud(order);
               
        if vt(1)/vt(2) > rel_thres
            disp('found one!!!');
            disp(vt(1)/vt(2));
            bc=order(1);
            badchannels = [badchannels bc];
            list_ic=[list_ic i];
        end
    end
    
    
    [badchannels,order]=sort(badchannels);
    list_ic=list_ic(order);
    
    base=timeline;
    
    for i=1:length(badchannels)
        disp(['channel ' num2str(list_gch(badchannels(i)))]);
        sig=IC(list_ic(i),:);
        der=diff(sig);
        der_abs=abs(der);
        thres=10*std(der_abs);
        [peaks]=trig_dev(der_abs,thres);
        time_art=base(peaks);
        for k=length(time_art):-1:2
            if time_art(k)-time_art(k-1) < 0.5
                time_art(k)=[];
                peaks(k)=[];
            end
        end
        time_art=0.01*round(100*time_art);
        if ~isempty(time_art)
            disp(['time instants : ' num2str(time_art) ' sec']);
        else
            disp('no peaks detected');
        end
        
        
    end
    
    badchannels=unique(badchannels);
    if isempty(badchannels)
        flag=1;
    else
        list_gch(badchannels)=[];
    end
    
end

list_bch=[1:size(data,1)];
list_bch(list_gch)=[];
list_bch=list_chan(list_bch);

delete(Dd);

disp(['NWC processing completed!']);