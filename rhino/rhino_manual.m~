function rhino_manual(D)

MainFig = figure('Name',                  'RHINO'   ,...
                 'NumberTitle',           'off'     ,...
                 'Menubar',               'none'    ,...
                 'DockControls',          'off'     ,...
                 'KeyPressFcn',           @key_press);          
             
uitools.toolbar = uitoolbar;
rotateicon = load(fullfile(matlabroot,'toolbox','matlab','icons','rotate.mat'));
uitools.rotate = uipushtool(uitools.toolbar,'ClickedCallback',@togglerot,'CData',rotateicon.cdata,'TooltipString','Toggle rotation');

rhino_display(D,MainFig)
%set(MainFig,'busyaction','queue','KeyPressFcn',@key_press);
%togglerot


    function key_press(~,evnt)

        p = D.inv{1}.mesh.tess_rhino.vert';
        
        inc_rot = 0.05;
        inc_trans = 1;
        
        switch evnt.Key
            case 'q'
                p = rhino_rigidtransform(p,inc_rot,0,0,0,0,0);
            case 'e'
                p = rhino_rigidtransform(p,-inc_rot,0,0,0,0,0);
            case 's'
                p = rhino_rigidtransform(p,0,inc_rot,0,0,0,0);
            case 'w'
                p = rhino_rigidtransform(p,0,-inc_rot,0,0,0,0);
            case 'd'
                p = rhino_rigidtransform(p,0,0,inc_rot,0,0,0);
            case 'a'
                p = rhino_rigidtransform(p,0,0,-inc_rot,0,0,0);
                
            case 'uparrow'
                p = rhino_rigidtransform(p,0,0,0,inc_trans,0,0);
            case 'downarrow'
                p = rhino_rigidtransform(p,0,0,0,-inc_trans,0,0);
            case 'leftarrow'
                p = rhino_rigidtransform(p,0,0,0,0,inc_trans,0);
            case 'rightarrow'
                p = rhino_rigidtransform(p,0,0,0,0,-inc_trans,0);
            case '['
                p = rhino_rigidtransform(p,0,0,0,0,0,inc_trans);
            case ']'
                p = rhino_rigidtransform(p,0,0,0,0,0,-inc_trans);
        end
        
        D.inv{1}.mesh.tess_rhino.vert = p';
        cla, rhino_display(D,MainFig)
        rotate3d off

    end

    function togglerot(~,evnt)
        rotate3d
    end


end