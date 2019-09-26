function osl_braingraph(G, Glims, NodeValues, NodeLims, mnipos, labels, thresh, spherecols, scaleEdges)
% Pretty graph plotting
% osl_braingraph(G, Glims, NodeValues, NodeLims, mnipos, labels, thresh, spherecols, edgeScales)
% --------------------------------------------------------------
% G           - a matrix describing the graph connections (or [] for no connections)
% Glims       - threshold for G colourmap (or [] for default size)
% NodeValues  - scalar value for each node size
% NodeLims    - size limits for NodeValues
% mnipos      - mni coordinates of nodes
% labels      - text labels for nodes (or [] for no labels)
% thresh      - percentile threshold above which graph connections are displayed
% spherecols  - [r g b] colours for each node
% edgeScales  - limits for linewidths when drawing edges - leave empty to
%               have constant thickness, or set a range such as [2 15].
% --------------------------------------------------------------
% Adam Baker 2012

OSLDIR = getenv('OSLDIR');

cla, hold on

% prettiness = 0.1; % not very pretty
prettiness = 0.5; % quite pretty
% prettiness = 1; % very pretty

brain_shading = repmat(0.5,1,3); brain_alpha = 0.1;

mesh = export(gifti([OSLDIR '/spm12/canonical/cortex_5124.surf.gii']), 'spm');
mesh = struct('faces',mesh.face,'vertices',mesh.vert);
mesh = reducepatch(mesh,prettiness);

trisurf(mesh.faces,mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3),'facecolor',brain_shading,'edgecolor','none','FaceAlpha',brain_alpha)

if ~isempty(G)
  G(logical(eye(size(G)))) = 0;

  if islogical(G)
      G = single(G);
      G(G==0) = nan;
  else
      if exist('thresh','var') && ~isempty(thresh)
          Gt = G(~tril(ones(size(G))));
          G(abs(G)<percentile(abs(Gt),thresh))=nan;
      end
  end
  
  if isempty(Glims)
    Glims = [-max(abs(G(:))) max(abs(G(:)))];
  end
  
  if Glims(1) >= Glims(2)
    Glims(2) = Glims(1)+0.1;
  end
  
end

if isempty(NodeLims)
  NodeLims = [0 max(abs(NodeValues(:)))];
end

NodeValues = (NodeValues - NodeLims(1)) ./ (NodeLims(2) - NodeLims(1));

NodeValues = NodeValues*8;

cmap      = colormap(bluewhitered(256));
% allow conversion to one-colour
if Glims(1) >= 0,
    cmap = cmap(129:end,:);
    isSingleColour = true;
elseif Glims(2) <= 0,
    cmap = cmap(1:128,:);
    isSingleColour = true;
else
    
    isSingleColour = false;
end
if ~exist('scaleEdges', 'var') || isempty(scaleEdges),
    scaleEdges = [5 5];
end
weightMap = linspace(scaleEdges(1), scaleEdges(2), size(cmap,1));

if ~isempty(G)
  [i,j] = find(~isnan(G));
  for p=length(i):-1:1,
    colorInd(p) = closest(G(i(p),j(p)),linspace(Glims(1),Glims(2),size(cmap,1)));
  end
  
  for p=1:length(i),
    edgecolour  = cmap(colorInd(p),:);
    
    if exist('scaleEdges', 'var') && ~isempty(scaleEdges),
        if isSingleColour,
            weightInd  = closest(colorInd(p), linspace(min(colorInd), max(colorInd), size(cmap, 1)));
        else
            weightInd = 2*closest(abs(colorInd(p) - round(size(cmap,1)./2)), linspace(0, round(max(colorInd)./2), round(size(cmap, 1)./2)));
        end
        edgeWeight = weightMap(weightInd);
    else
        edgeWeight = 3;
    end


    line(mnipos([i(p) j(p)],1),mnipos([i(p) j(p)],2),mnipos([i(p) j(p)],3),'color',edgecolour,'linewidth',edgeWeight)
  end
  
  set(gca,'clim',Glims);
  colormap(bluewhitered(256));
  hc = colorbar;
  FONTSIZE = 14;
  if isSingleColour, 
      YTicks = [Glims(1) Glims(2)];
  else
      YTicks = [Glims(1) 0 Glims(2)];
  end%if
  set(hc, ...
      'FontName', 'Helvetica', ...
      'FontSize', FONTSIZE, ...
      'Box', 'on', ...
      'TickDir', 'in', ...
      'XColor', [0.3 0.3 0.3], ...
      'YColor', [0.3 0.3 0.3], ...
      'LineWidth', 2, ...
      'YTick', YTicks, ...
      'TickLength', [0 0]);
  
% these bits don't work in Matlab 2015
%       'XMinorTick', 'off', ...
%       'YMinorTick', 'off', ...
%       'YGrid', 'off', ...
%       'XGrid', 'off', ...

YTL = get(hc,'yticklabel');
% set(hc,'yticklabel',[repmat(' ',size(YTL,1),1), YTL]);
end

% Draw coloured spheres for nodes
if ~exist('spherecols','var')
  spherecols = [];
end

[X,Y,Z] = sphere(10);

for n = 1:length(NodeValues)
  x = abs(NodeValues(n))*X; y = abs(NodeValues(n))*Y; z = abs(NodeValues(n))*Z;
  if isempty(spherecols)
    if NodeValues(n) < 0, spherecolour = [0 1 1]; else spherecolour = [1 0 1]; end
  else spherecolour = spherecols(n,:);
  end
  surf(x+mnipos(n,1),y+mnipos(n,2),z+mnipos(n,3),'edgecolor','none','facecolor',spherecolour,'facealpha',0.2);
end

if ~isempty(labels)
  text(1.1*mnipos(:,1),1.1*mnipos(:,2),1.1*mnipos(:,3), labels, 'FontSize',12,'HorizontalAlignment','center')
end

set(gcf,'renderer','openGL')

handles = get(gca,'children');
is_text  = strcmp(get(handles,'type'),'text');
arrayfun(@uistack,handles(is_text))

set(gca,'xColor','w')
set(gca,'yColor','w')
set(gca,'zColor','w')
set(gcf,'color','w')
axis image
axis vis3d
hold off
rotate3d('on')
view(-40,80);
end

function i = closest(a,k)
%CLOSEST finds index of vector a closest to k
assert(isscalar(k) | isscalar(a));

[~,i] = min(abs(a-k));
end%closest

function newmap = bluewhitered(m)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.


if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

% RcolorBrewer RdBu
bottom    = [5 113 176] / 255;
botmiddle = [146 197 222] / 255;
middle    = [247 247 247] / 255;
topmiddle = [244 165 130] / 255;
top       = [202   0  32] / 255;

% Find middle
lims = get(gca, 'CLim');

% Find ratio of negative to positive
if (lims(1) < 0) && (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % And put 'em together
    newmap = [newmap1; newmap];
    
elseif lims(1) >= 0
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
else
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
end
end%bluewhitered