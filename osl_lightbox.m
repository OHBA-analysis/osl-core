function osl_lightbox(in,spat_res,thresh,comp2use)

if nargin<4
    comp2use=0;
end
if nargin<3
    thresh=0;
end
if nargin<2
    spat_res=2;
end
if numel(thresh)<2
    thresh=[thresh -thresh];
end

OSLDIR = getenv('OSLDIR');

%% Read in Standard Brain
stdbrain=nii.load([OSLDIR '/std_masks/MNI152_T1_' num2str(spat_res) 'mm_brain.nii.gz']);
stdbrain=uint8(stdbrain/max(stdbrain(:))*255);
x_range=[min(find(sum(sum(stdbrain>0,2),3)>0)) max(find(sum(sum(stdbrain>0,2),3)>0))];
y_range=[min(find(sum(sum(stdbrain>0,1),3)>0)) max(find(sum(sum(stdbrain>0,1),3)>0))];
z_range=[min(find(sum(sum(stdbrain>0,1),2)>0)) max(find(sum(sum(stdbrain>0,1),2)>0))];

plot_sel=[1 2 7 8 13 14 19 20 25 26];

x_ind=round(x_range(1):diff(x_range)/10:x_range(2));
y_ind=round(y_range(1):diff(y_range)/10:y_range(2));
z_ind=round(z_range(1):diff(z_range)/10:z_range(2));

%% Read in images
if strcmp(class(in),'char')
    try, vol_in=nii.load(in); catch, error(['Unable to read nifti file ' in]); end
else
    if numel(size(in))<3
        try, vol_in = matrix2vols(in,stdbrain); catch, error('Unable to convert input matrix into volume');end
    else
        vol_in = in;
    end
end

fsldir=getenv('FSLDIR');
[cdata_p,map_p] = imread([fsldir '/etc/luts/ramp.gif'],'gif');
[cdata_n,map_n] = imread([fsldir '/etc/luts/ramp2.gif'],'gif');
cbar_p = permute(ind2rgb( cdata_p, map_p ),[2 1 3]);
cbar_p=uint8(255*cbar_p(end:-1:1,:,:));
cbar_p=cbar_p(end:-1:1,:,:);
cbar_n = permute(ind2rgb( cdata_n, map_n ),[2 1 3]);
cbar_n=uint8(255*cbar_n(end:-1:1,:,:));

if comp2use~=0
    vol_in=vol_in(:,:,:,comp2use);
end

%% Generate plots
for I=1:size(vol_in,4)
    vol=vol_in(:,:,:,I);
    volp=vol;volp(volp<thresh(1))=0;volp=volp*255/max(abs(vol(:)));
    voln=-vol;voln(voln<-thresh(2))=0;voln=voln*255/max(abs(vol(:)));
       
    ind=0;
    figure('Position',[100 100 1000 700]);
    for plot_ind=plot_sel
        ind=ind+1;
        subplot(5,6,plot_ind);
        brain_slice=flipud(permute(stdbrain(x_ind(ind),:,:),[2 3 1])');
        p_slice=flipud(permute(volp(x_ind(ind),:,:),[2 3 1])');
        n_slice=flipud(permute(voln(x_ind(ind),:,:),[2 3 1])');
        
        slice_image=zeros([size(brain_slice) 3]);slice_image(:,:,1)=brain_slice;slice_image(:,:,2)=brain_slice;slice_image(:,:,3)=brain_slice;
        
        for x=1:size(slice_image,2)
            for y=1:size(slice_image,1);
                if p_slice(y,x)>0
                    cind=ceil(size(cbar_p,1)*p_slice(y,x)/255);
                    slice_image(y,x,:)=cbar_p(cind,1,:);
                end
                if n_slice(y,x)>0
                    cind=ceil(size(cbar_n,1)*n_slice(y,x)/255);
                    slice_image(y,x,:)=cbar_n(cind,1,:);
                end
            end
        end
        slice_image=uint8(slice_image);
        image(slice_image);
        axis equal;
        axis tight;
    end
    
    ind=0;
    for plot_ind=plot_sel+2
        ind=ind+1;
        subplot(5,6,plot_ind);
        brain_slice=flipud(permute(stdbrain(:,y_ind(ind),:),[1 3 2])');
        p_slice=flipud(permute(volp(:,y_ind(ind),:),[1 3 2])');
        n_slice=flipud(permute(voln(:,y_ind(ind),:),[1 3 2])');
        
        slice_image=zeros([size(brain_slice) 3]);slice_image(:,:,1)=brain_slice;slice_image(:,:,2)=brain_slice;slice_image(:,:,3)=brain_slice;
        
        for x=1:size(slice_image,2)
            for y=1:size(slice_image,1);
                if p_slice(y,x)>0
                    cind=ceil(size(cbar_p,1)*p_slice(y,x)/255);
                    slice_image(y,x,:)=cbar_p(cind,1,:);
                end
                if n_slice(y,x)>0
                    cind=ceil(size(cbar_n,1)*n_slice(y,x)/255);
                    slice_image(y,x,:)=cbar_n(cind,1,:);
                end
            end
        end
        
        slice_image=uint8(slice_image);
        image(slice_image);
        axis equal;
        axis tight;
    end
    
    ind=0;
    for plot_ind=plot_sel+4
        ind=ind+1;
        subplot(5,6,plot_ind);
        brain_slice=flipud(permute(stdbrain(:,:,z_ind(ind)),[1 2 3])');
        p_slice=flipud(permute(volp(:,:,z_ind(ind)),[1 2 3])');
        n_slice=flipud(permute(voln(:,:,z_ind(ind)),[1 2 3])');
        
        slice_image=zeros([size(brain_slice) 3]);slice_image(:,:,1)=brain_slice;slice_image(:,:,2)=brain_slice;slice_image(:,:,3)=brain_slice;
        
        for x=1:size(slice_image,2)
            for y=1:size(slice_image,1);
                if p_slice(y,x)>0
                    cind=ceil(size(cbar_p,1)*p_slice(y,x)/255);
                    slice_image(y,x,:)=cbar_p(cind,1,:);
                end
                if n_slice(y,x)>0
                    cind=ceil(size(cbar_n,1)*n_slice(y,x)/255);
                    slice_image(y,x,:)=cbar_n(cind,1,:);
                end
            end
        end
       
        slice_image=uint8(slice_image);
        image(slice_image);
        axis equal;
        axis tight;
    end
end
end