function h = coreg(D)

    if strcmp(D.inv{1}.comment,'rhino')
        h = figure('name','RHINO coregistration','tag','rhino_coreg');
        rhino_display(D,h);
    end

