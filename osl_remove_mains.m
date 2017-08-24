function D2 = osl_remove_mains(D,ac_freq)
    % Apply a notch filter at given AC frequency and first harmonic if present
    % As a filtering step, this makes a new D object and saves it to disk

    if nargin < 2 || isempty(ac_freq) 
        ac_freq = 50;
    end
    
    D2 = osl_filter(D,-1*(ac_freq+[-2 2]));

    if D2.fsample > 2*(ac_freq+2) % If second harmonic is present
        D2 = osl_filter(D,-1*(2*ac_freq+[-2 2]));
    end
    
end
