function component_topoplot(D,comp,modality)
%COMPONENT_TOPOPLOT plots a field over channels on a head shape. 
% subfunction taken from identify_artefactual_components_manual.m in OSL

%	Copyright 2014 OHBA
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 213 $
%	$LastChangedDate: 2014-07-24 12:39:09 +0100 (Thu, 24 Jul 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 28-Oct-2013 

cfg=[];
data=[];

OSLDIR = getenv('OSLDIR');

chan_inds=setdiff(find(any([strcmp(D.chantype,'EEG');strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)),D.badchannels);
map_inds(find(any([strcmp(D.chantype,'EEG');strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1))) = 1:numel(find(any([strcmp(D.chantype,'EEG');strcmp(D.chantype,'MEGMAG');strcmp(D.chantype,'MEGPLANAR');strcmp(D.chantype,'MEGGRAD')],1)));


comp_full(map_inds(chan_inds),:)=comp;
comp_full(D.badchannels) = 0;
comp = comp_full;

if (strcmp(modality,'MEGMAG')),
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = [OSLDIR '/layouts/neuromag306mag.lay'];
    chanind = strmatch(cfg.channel, D.chantype);
    
    comp2view=zeros(102,1);
    for ind=1:length(chanind),
        indp=(ind-1)*3+3;
        comp2view(ind,:)=comp(indp);
    end;
    data.dimord='chan_comp';
    data.topo=comp2view;
    data.topolabel = D.chanlabels(chanind);
    data.time = {1};
    
elseif (strcmp(modality,'MEGPLANAR')),
    cfg.channel     = {'MEGMAG'};
    cfg.layout      = [OSLDIR '/layouts/neuromag306mag.lay'];
    chanind = strmatch(cfg.channel, D.chantype);
    comp2view=zeros(102,1);
    for ind=1:length(chanind),
        indp=(ind-1)*3+1;
        comp2view(ind,:)=comp(indp)+comp(indp+1);
    end;
    data.dimord='chan_comp';
    data.topo=comp2view;
    data.topolabel = D.chanlabels(chanind);
    data.time = {1};
    warning('The planar grads are being overlayed on the magnetometer locations and using the magnetometer labels');

elseif (strcmp(modality,'MEGGRAD')),
    cfg.channel     = {'MEG'};
    cfg.layout      = [OSLDIR '/layouts/CTF275.lay'];
    chanind = strmatch(cfg.channel, D.chantype);
    comp2view=comp;
    data.dimord='chan_comp';
    data.topo=comp2view;
    data.topolabel = D.chanlabels(chanind);
    data.time = {1};
    
elseif (strcmp(modality,'EEG')),
    cfg.channel     = {'EEG'};
    A=sensors(D,'EEG');
    pnt=A.chanpos;
    elec = [];
    elec.pnt = pnt(chan_inds,:);
    elec.label = A.label(chan_inds);
    cfg.elec  = elec;
    %cfg.layout = ft_prepare_layout(cfg);
    cfg.channel={'all'};
    
    comp2view=comp;
    data.dimord='chan_comp';
    data.topo=comp2view;
    data.topolabel = A.label(chan_inds);
    data.label = A.label(chan_inds);
    data.time = {1};
else
    
    error('Unsupported modality');
end

cfg.component = 1;
cfg.interactive = 'no';
cfg.comment     = '';
ft_topoplotIC(cfg,data);

end