function M = rhino_icp(data1,data2,Ninits,ax)
% Run Iterative Closest Point with multiple initialisations
% M = rhino_icp(data1,data2,Ninits,do_plots
%
% REQUIRED INPUTS:
%
% data1      - [3 x N] locations of the set of static points 
%              (e.g. MRI scalp surface)
% data2      - [3 x N] locations of the set of moving points 
%              (e.g. Polhemus headshape points)
%
% OPTIONAL INPUTS:
%
% Ninits     - Number of random initialisations to perform (default 10)
% ax         - Axes handle or [] to plot iterations
%
% OUTPUTS:
%
% M          - [4 x 4] rigid transformation matrix mapping data2 to data1
%
% Based on ICP.m by Martin Kjer and Jakob Wilm, Technical University of 
% Denmark
%
% Adam Baker 2014

if nargin < 3
    Ninits = 10;
end

if nargin < 4 || isempty(ax)
    do_plots = 0;
else
    do_plots = 1;
    if ~ishandle(ax)
        ax = axes;
    end
end


            
err_old = Inf;
err = zeros(1,Ninits);

Mr = eye(4);
data2r  = data2;

for init = 1:Ninits

    [TR,TT,e] = icp(data1,data2r,50,'Minimize','plane','Matching','kDtree','ReturnAll',true);
    
    [~,emin] = min(e);
    e  = e(1:emin);
    TR = TR(:,:,1:emin);
    TT = TT(:,:,1:emin);
    
    Mi = [TR(:,:,end) TT(:,:,end);0 0 0 1];

    err(init) = e(end);
 
    if err(init) < err_old
        err_old = err(init);
        M     = Mi*Mr;
    end
        
    if do_plots
        hs = plot(ax,[]);
        xlim = get(ax,'xlim'); ylim = get(ax,'ylim'); zlim = get(ax,'zlim');
        for i = 1:length(e)
            delete(hs);
            Mi = [TR(:,:,i) TT(:,:,i);0 0 0 1];
            Mi = Mi*Mr;
            dataM = spm_eeg_inv_transform_points(Mi,data2')';
            hold(ax,'on')
            hs = plot3(ax,dataM(1,:),dataM(2,:),dataM(3,:),'.b');
            set(ax,'xlim',xlim); set(ax,'ylim',ylim); set(ax,'zlim',zlim);
            title(ax,num2str(e(i)))
            drawnow
        end
        delete(hs)
    end
    

    
    % Give the registration a kick...
    % 15 degrees rotation about each axis, 10 mm translation
    a = (rand-0.5)*pi/6; b = (rand-0.5)*pi/6; c = (rand-0.5)*pi/6;
    Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
    Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
    Rz = [cos(c) -sin(c) 0; sin(c) cos(c) 0; 0 0 1];
    T  = 10*[(rand-0.5) (rand-0.5) (rand-0.5)];
    Mr = [Rx*Ry*Rz T'; 0 0 0 1];
    
    data2r  = Mr*[data2; ones(1,size(data2,2))];
    data2r  = data2r(1:3,:);
   
end


end


