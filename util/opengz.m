function o = opengz(filename)
    % Use osleyes to open .gz files
    % Just open all gz files because user is unlikely to be 
    % double-clicking tar.gz files (for example) and because it 
    % doesn't seem possible to provide open functions for
    % 'nii.gz'
    o = osleyes(filename);