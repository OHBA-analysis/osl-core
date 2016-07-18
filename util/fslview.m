% fslview.m
%
% Matlab wrapper for calling FSLVIEW to view .nii files.
%
% Syntax: fslview(fnames, thresholds, colour_maps,anatomical)
%
% e.g. fslview('/home/myfiles/tstats');
% e.g. fslview({'tstats.nii.gz'; 'copes.nii.gz'}, [2.3 5; 0.2 0.5], [1 2]);
% e.g. fslview({'tstats.nii.gz'; 'copes.nii.gz'}, [2.3 5; 0.2 0.5], {'Red-Yellow';'Blue-Lightblue'});
% e.g. fslview({'tstats.nii.gz'; 'copes.nii.gz'}, [2.3 5; 0.2 0.5], 'mni_brain');%
%
% thresholds can be left empty, i.e. [], to allow adaptive setting of min max
%
% Version 1.0
% 16.01.13

function fslview(fnames,thresholds,colour_maps,anatomical)

global OSLDIR

available_maps = {'Red-Yellow';
                  'Blue-Lightblue'; 
                  'Red'; 
                  'Green';
                  'Blue'; 
                  'Yellow'; 
                  'Pink'; 
                  'Hot'; 
                  'Cool'; 
                  'Copper'};


if nargin<4
    anatomical='mni_brain';
end;

if nargin<3
   colour_maps = {'Red-Yellow'};
end

if nargin<2
   thresholds = [2.3 5];
end

if nargin==0 || isempty(fnames)
    [filename, pathname]=uigetfile({'*.nii; *.nii.gz', 'Nifti Files (.nii , .nii.gz)';'*.*', 'All Files'},'Select a nifti file');
    fnames=fullfile(pathname,filename);
end

warning off

if ischar(fnames)
  fnames = {fnames};
end

% Convert index list into colour map names
if isnumeric(colour_maps)
    tmp=colour_maps;colour_maps=[];
    for f=1:numel(tmp)
        colour_maps{f}=available_maps{tmp(f)};
    end
    clear tmp;
end

% Pad out missing colour maps to match number of volumes
if numel(colour_maps)<numel(fnames)
    for f=numel(colour_maps):numel(fnames)
        colour_maps{f}=colour_maps{end};
    end
end

% Pad out missing threshold values to match number of volumes
if ~isempty(thresholds)
    if size(thresholds,1)<numel(fnames)
        for f=size(thresholds,1):numel(fnames)
            thresholds(f,:)=thresholds(end,:);
        end
    end
end;

fnames_formatted = [];
for f=1:numel(fnames)
  [~,~,scales,~,~] = read_avw(fnames{f});
  gs(f) = scales(1);
  colmap=[' -l "' colour_maps{f} '" '];
  if ~isempty(thresholds)
    vals = [' -b ' num2str(thresholds(f,1)) ',' num2str(thresholds(f,2)) ' '];
  else
    vals = '';
  end;
  
  fnames_formatted = [fnames_formatted fnames{f} colmap vals];
end

if numel(unique(gs)) ~= 1
  error('spatial maps must be of equal sizes')
end

switch anatomical
    case 'mni_brain'
        anatomical_fname=[OSLDIR '/std_masks/MNI152_T1_' num2str(gs(1)) 'mm_brain'];
    otherwise
        anatomical_fname=anatomical;    
end;

runcmd(['fslview ' fnames_formatted ' ' anatomical_fname  '&']);
 
warning on
end

% HELP TEXT FROM FSLVIEW
%
% fslview [-m 3d|ortho|lightbox] <baseimage> [-l lutname] [-b low,hi]
% 	[ <overlay> [-l lutname] [-b low,hi] ] ...
% fslview -m ortho,lightbox filtered_func_data thresh_zstat1 -t 0.5 thresh_zstat2 -l "Cool" -t 0.5
% 
% Optional arguments (You may optionally specify one or more of):
% 	-V,--verbose	switch on diagnostic messages
% 	-h,--help	display this message
% 	-m,--mode	Initial viewer mode. Comma separated list of: 3d; single, ortho; lightbox
% 
% 
% Per-image options
% 
% Usage: 
% image [-l GreyScale] [-t 0.1] [-b 2.3,6]
% 	-l,--lut	Lookup table name. As per GUI, one of: Greyscale;
% 			"Red-Yellow"; "Blue-Lightblue"; Red; Green;
% 			Blue; Yellow; Pink; Hot; Cool; Copper, etc.
% 	-b,--bricon	Initial bricon range, e.g., 2.3,6
% 	-t,--trans	Initial transparency, e.g., 0.2
