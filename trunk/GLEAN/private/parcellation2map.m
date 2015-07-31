function pmap = parcellation2map(map,parcellation_fname,mask_fname)
% converts a parcels x components parcellation map to a full spatial map 
% using parcellation_fname as a mask

parcellation = readnii(parcellation_fname,mask_fname);

if size(parcellation,2) > 1
    % Create non-overlapping parcellation
    [~,assignments] = max(parcellation,[],2);
    assignments(all(parcellation==0,2)) = 0;
else
    assignments = logical(parcellation);
end

pmap = zeros(size(assignments,1),size(map,2));
for p = 1:size(map,1)
    for comp = 1:size(map,2)
        pmap(assignments == p,comp) = map(p,comp);
    end
end
    
end