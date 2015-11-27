function rhino_manual(D)
% Interactive tool for manual tweaking of RHINO coregistration.
%
% rhino_manual(D) - displays the coregistration for MEEG object D in a 
%                   new figure
% 
% Adam Baker 2015

try 
    D.inv{1}.mesh.tess_rhino;
catch
    error('RHINO has not been run!')
end


MainFig = figure('NumberTitle',           'off'     ,...
                 'Menubar',               'none'    ,...
                 'DockControls',          'off'     ,...
                 'CloseRequestFcn',       @closefig ,...
                 'KeyPressFcn',           @key_press);


% Display error
err = sum(rhino_cperror(D.inv{1}.mesh.tess_rhino.vert,D.inv{1}.datareg(1).fid_eeg.pnt));
set(MainFig,'name',['RHINO - error = ' num2str(err)])


uitools.toolbar = uitoolbar;
rotateicon = load(fullfile(matlabroot,'toolbox','matlab','icons','rotate.mat'));
uitools.rotate = uitoggletool(uitools.toolbar,'ClickedCallback',@togglerot,'CData',rotateicon.cdata,'TooltipString','Toggle rotation','state','on');

uicontrol('Style','text',...
          'String',sprintf('Translate using cursor keys and N/M keys. Rotate using Q/E,W/S,A/D keys'),...
          'Units','Normalized',...
          'Position',[0 0 1 0.05],...
          'Parent',MainFig);

rhino_display(D,MainFig)
%set(MainFig,'busyaction','queue','KeyPressFcn',@key_press);
%togglerot

inv_org = D.inv{1};
is_modified = 0;

    function key_press(~,evnt)
        is_modified = 1;
        
        p = D.inv{1}.mesh.tess_rhino.vert';
        
        inc_rot = 0.03;
        inc_trans = 1;
        
        switch evnt.Key
            case 'w'
                [p,Mr] = rhino_rigidtransform(p,inc_rot,0,0,0,0,0);
            case 's'
                [p,Mr] = rhino_rigidtransform(p,-inc_rot,0,0,0,0,0);
            case 'q'
                [p,Mr] = rhino_rigidtransform(p,0,inc_rot,0,0,0,0);
            case 'e'
                [p,Mr] = rhino_rigidtransform(p,0,-inc_rot,0,0,0,0);
            case 'd'
                [p,Mr] = rhino_rigidtransform(p,0,0,inc_rot,0,0,0);
            case 'a'
                [p,Mr] = rhino_rigidtransform(p,0,0,-inc_rot,0,0,0);
                
            case 'rightarrow'
                [p,Mr] = rhino_rigidtransform(p,0,0,0,inc_trans,0,0);
            case 'leftarrow'
                [p,Mr] = rhino_rigidtransform(p,0,0,0,-inc_trans,0,0);
            case 'n'
                [p,Mr] = rhino_rigidtransform(p,0,0,0,0,inc_trans,0);
            case 'm'
                [p,Mr] = rhino_rigidtransform(p,0,0,0,0,-inc_trans,0);
            case 'downarrow'
                [p,Mr] = rhino_rigidtransform(p,0,0,0,0,0,inc_trans);
            case 'uparrow'
                [p,Mr] = rhino_rigidtransform(p,0,0,0,0,0,-inc_trans);
        end
        
  
        % APPLY NEW REGISTRATION TO SPM OBJECT
        D.inv{1}.mesh.tess_rhino.vert      = p';
        D.inv{1}.datareg(1).fid_mri.fid.pnt   = spm_eeg_inv_transform_points(Mr, D.inv{1}.datareg(1).fid_mri.fid.pnt);
        D.inv{1}.datareg(1).fid_mri.pnt       = spm_eeg_inv_transform_points(Mr, D.inv{1}.datareg(1).fid_mri.pnt);
        D.inv{1}.datareg(1).fromMNI           = Mr * D.inv{1}.datareg(1).fromMNI;
        D.inv{1}.datareg(1).toMNI             = inv(D.inv{1}.datareg(1).fromMNI);
        %D.save;
 
        % DISPLAY NEW COREGISTRATION
        cla, rhino_display(D,MainFig)
        rotate3d off
        
        % Display error
        err = sum(rhino_cperror(D.inv{1}.mesh.tess_rhino.vert,D.inv{1}.datareg(1).fid_eeg.pnt));
        set(MainFig,'name',['RHINO - error = ' num2str(err)])
        
    end

    function togglerot(~,evnt)
        if strcmp(get(uitools.rotate,'state'),'on')
            rotate3d ON
        else
            rotate3d off
        end
    end

    function closefig(src,evnt)
        if is_modified
            selection = questdlg('Keep new position?','Exit?');
            switch selection,
                case 'Yes'
                    % Do nothing
                case 'No'
                    D.inv{1} = inv_org;
                case 'Cancel'
                    return
            end
        end
        D.save;
        delete(MainFig)
    end


end