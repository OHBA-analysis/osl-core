function [movie_fname vols2surf_results] = osl_make_surf_movie(S)

% osl_make_surf_movie(S)
%
% S=[];
% S.vol='/Users/woolrich/homedir/matlab/osl_testdata_dir/faces_subject1_data/spm8_meg1.mat_wideband_beamform.oat/session1_wholebrain_first_level_dir/tstat3_2mm.nii.gz'
% S.minmax=[3, 7]; S.time_range=[0 0.2]; S.times=stats.times;
%
% Mark Woolrich 2013

S.view=0;
S.type='freesurfer';
S.time_indices=intersect(find(S.times>=S.time_range(1)), find(S.times<=S.time_range(2)));
times_included=S.times(S.time_indices);

if ~isfield(S,'vols2surf_results'),
    vols2surf_results = osl_render_vols_to_surf(S);
else
    vols2surf_results = S.vols2surf_results;
end;

movie_fname=[vols2surf_results.workingdir '/surface_frames/movie'];

runcmd(['rm -f ' movie_fname '.avi']);

avi=avifile(movie_fname);
no_frame=length(times_included);
Movie = cell(no_frame,1);
 
% make timer:

ii=no_frame;
IML=imread(vols2surf_results.left_tiff_names{ii});
IMR=imread(vols2surf_results.right_tiff_names{ii});

x_timer_buffer=100;
y_timer_buffer=8;

htimer=figure;
set(htimer,'position', [100,100, 800, 100]);
axes('position',[0 0 1 1]);
timer=zeros(round(size([IML IML],1)/50),size([IML IML],2),3);
imagesc(timer);

set(gca, 'visible', 'off') 
set(gcf, 'color', 'k'); 
timerf=getframe(htimer);

IML=imread(vols2surf_results.left_tiff_names{ii});
IMR=imread(vols2surf_results.right_tiff_names{ii});
%IM=double(cat(1,[IML IMR], uint16(timerf.cdata)));
%IM=double(timerf.cdata)/double(intmax('uint8'));
%IM=double(cat(1,[IML IMR], uint16(timerf.cdata)));
IM=double(timerf.cdata)/double(intmax('uint8'));
IM2=cat(1,double([IML IMR])/double(intmax('uint16')),IM);
figure;image(IM2);

for ii=1:no_frame;

    htimer=figure;
    set(htimer,'position', [100,100, 800, 100]);
    axes('position',[0 0 1 1]);
    timer=zeros(round(size([IML IML],1)/50),size([IML IML],2),3);
    timer(y_timer_buffer:end-y_timer_buffer,x_timer_buffer:x_timer_buffer+round(ii*(size([IML IML],2)-2*x_timer_buffer)/no_frame),:)=1;
    imagesc(timer);

    set(gca, 'visible', 'off') 
    set(gcf, 'color', 'k'); 
    timerf=getframe(htimer);

    IML=imread(vols2surf_results.left_tiff_names{ii});
    IMR=imread(vols2surf_results.right_tiff_names{ii});
    %IM=double(cat(1,[IML IMR], uint16(timerf.cdata)));
    %IM=double(timerf.cdata)/double(intmax('uint8'));
    %IM=double(cat(1,[IML IMR], uint16(timerf.cdata)));
    IM=double(timerf.cdata)/double(intmax('uint8'));
    IM2=cat(1,double([IML IMR])/double(intmax('uint16')),IM);
    figure;image(IM2);

   Movie=im2frame(IM2);
 %  Movie=im2frame(IM);
   avi=addframe(avi,Movie);
    ca
end 
avi=close(avi);

movie_fname=[movie_fname '.avi'];
['!open ' movie_fname ' &']
