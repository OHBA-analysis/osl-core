function sens = undobalancing(sens)
    % In order to retain the third-order gradiometers, this file disables UndoBalancing entirely
    warning([mfilename ':Disabled'],'UndoBalancing is disabled. Your analysis will proceed using the third-order gradiometer-corrected signals.\n');
