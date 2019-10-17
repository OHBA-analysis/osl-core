function openfif(filename)
    % Open nii files in Matlab using osleyes when double clicking
    [~,fname] = fileparts(filename);
    fprintf(2,'%s.fif is a binary file, cannot open by double clicking\n',fname);
    