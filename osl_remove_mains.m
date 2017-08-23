function D2 = osl_remove_mains(D,ac_freq)
    % Apply a notch filter at given AC frequency and first harmonic if present
    % As a filtering step, this makes a new D object and saves it to disk

    if nargin < 2 || isempty(ac_freq) 
        ac_freq = 50;
    end
    
    S3              = [];
    S3.D            = D;
    S3.type          = 'butterworth';
    S3.freq         = ac_freq+[-2 2];
    S3.band         = 'stop';
    S3.dir          = 'twopass';
    S3.order        = 5;
    D2 = spm_eeg_filter(S3);

    if D2.fsample > 2*(ac_freq+2) % If second harmonic is present
        S3              = [];
        S3.D            = D2;
        S3.type          = 'butterworth';
        S3.freq         = 2*ac_freq+[-2 2];
        S3.band         = 'stop';
        S3.dir          = 'twopass';
        S3.order        = 5;
        D2 = spm_eeg_filter(S3);
    end
    
end
