function Dnew = osl_merge_meeg(filelist,newfile)
% Merge continuous MEEG objects in filelist to a new single object newfile
% Adam Baker 2015

D = cell(size(filelist));
for i = 1:length(filelist)
    D{i} = spm_eeg_load(filelist{i});
end
   
num_samples = cellfun(@nsamples,D);

end_samples = cumsum(num_samples);
start_samples = [1 end_samples(1:end-1)+1];

% Copy data:
Dnew = clone(D{1},newfile,[D{1}.nchannels,sum(num_samples),1]);
for i = 1:length(filelist)
    Dnew(:,start_samples(i):end_samples(i),:) = D{i}(:,:,:);
end

% Ensure bad channels are consistent:
Dnew = badchannels(Dnew,1:Dnew.nchannels,0);
Dnew = badchannels(Dnew,unique(cell2mat(cellfun(@badchannels,D,'uniformoutput',0))),1);

% Correct events:
ev = cellfun(@events,D,'uniformoutput',0);
for i = 1:length(filelist)
    for j = 1:numel(ev{i})
        ev{i}(j).time = ev{i}(j).time + (start_samples(i) - 1)/Dnew.fsample;
    end
end
ev = cat(1,ev{:});
Dnew = events(Dnew,1,ev);

Dnew.save;

end