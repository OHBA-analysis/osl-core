function osl_mlg_plot
    % This file is OSL2's 1000th commit!

    if isempty(get(0,'children'))
        p = parcellation(8);
        p.plot
    end

    figure(gcf)
    
    fh = gcf;
    ax = gca;

    c = get(fh,'Color');
    v = get(ax,'View');
    u = get(ax,'Units');
    vis = get(ax,'Visible');
    set(ax,'Units','normalized')
    pos = get(ax,'Position');
    set(ax,'Visible','off')
    set(ax,'View',[  200    6.0000])

    front = [1     1     1     1     1     1     1     1     1     1     1     1     1   1     1     1     1     1     1     1     1     1     1     1     1     1;...
    1     1     0     1     0     1     1     1     1     1     1     1     0   0     1     0     1     0     1     1     1     1     1     1     1     1;...
    0     1     1     0     1     0     1     1     1     1     1     1     0   0     1     1     0     1     0     1     1     1     1     1     1     0;...
    0     0     1     1     0     1     0     1     1     1     1     1     0   0     0     1     1     0     1     0     1     1     1     1     1     0;...
    0     0     0     1     1     1     1     1     1     1     1     0     0   0     0     0     1     1     1     1     1     1     1     1     0     0;...
    0     0     0     0     1     1     1     1     1     1     0     0     0   0     0     0     0     1     1     1     1     1     1     0     0     0];
    front = repmat(+~front,1,1,3);

    front_a = [1     1     1     1     1     1     1     1     1     1     1     1     1   1     1     1     1     1     1     1     1     1     1     1     1     1;...
    1     1     1     1    1     1     1     1     1     1     1     1     0   0     1     1     1     1     1     1     1     1     1     1     1     1;...
    0     1     1     1     1     1     1     1     1     1     1     1     0   0     1     1     1     1     1     1     1     1     1     1     1     0;...
    0     0     1     1     1     1     1     1     1     1     1     1     0   0     0     1     1     1     1     1     1     1     1     1     1     0;...
    0     0     0     1     1     1     1     1     1     1     1     0     0   0     0     0     1     1     1     1     1     1     1     1     0     0;...
    0     0     0     0     1     1     1     1     1     1     0     0     0   0     0     0     0     1     1     1     1     1     1     0     0     0];

    set(ax,'YLim',[-105 80])
    s1 = surface(-[-50,50;-50,50],80*ones(2),[20,20;0,0]);
    set(s1,'CData',front,'alphadata',0*front_a,'facecolor','texture','edgealpha',0,'facealpha','texturemap','alphadatamapping','none');
    s2 = surface([50 60;50 60],[80 0;80 0],[20,20;14,14],'FaceAlpha',1,'alphadatamapping','none');
    s3 = surface(-[50 60;50 60],[80 0;80 0],[20,20;14,14],'FaceAlpha',1,'alphadatamapping','none');
    set([s2 s3],'CData',zeros(2,2,3));

    s = [s1 s2 s3];

    cl = onCleanup(@() restore(fh,ax,pos,c,v,u,vis,s));

    offsets = linspace(6*pi,0,100);
    position(s,offsets(1));

    alphafade = linspace(0,1,50);
    for j = 1:length(alphafade)
        set(s(1),'alphadata',alphafade(j)*front_a);
        set(s(2:3),'FaceAlpha',alphafade(j));
        drawnow
    end

    pause(0.8)

    for j = 1:length(offsets)
        position(s,offsets(j));
        drawnow
    end

    pause(0.8)

    while 1 && ishandle(fh) % if figure is not deleted
        j=j+1;
        if mod(j,4)==0
            set(fh,'Color',rand(3,1));
        end
        set(ax,'View', [200+5*sin(j/20) 5*cos(2*j/25)])
        set(ax,'Position',pos + [0 0 0.2 0.2]*sin(j/10) + [0.1 0.1 0 0]*sin(pi+j/10))
        drawnow
    end

    function restore(fh,ax,pos,c,v,u,vis,s)
        if ~ishandle(fh)
            return
        end
        set(ax,'Position',pos);
        set(fh,'Color',c);
        set(ax,'View',v);
        set(ax,'Units',u);
        set(ax,'Visible',vis);
        delete(s)

    function position(s,z)
        set(s(1),'ZData',[z+20 z+20;z z])
        set(s(2:3),'ZData',[z+20 z+20;z+14 z+14])

