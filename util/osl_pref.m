function varargout = osl_pref(key,val)
    if nargin == 1
        operation = 'get';
    else
        operation = 'set';
        assert(ischar(val),'Value to write to preference file must be a string')
    end


    fname = fullfile(osldir,'preferences.txt');

    if ~osl_util.isfile(fname)
        if strcmp(operation,'get')
            varargout{1} = [];
            return
        else
            fclose(fopen(fname,'w'));
        end
    end

    % Get the conf struct
    f = fopen(fname,'r');
    conf = struct();
    while ~feof(f)
        l = fgetl(f); 
        if ~ischar(l)
            break
        end
        [n,v] = strtok(l,'=');
        conf.(n) = v(2:end);
    end
    fclose(f);
    keys = fields(conf);

    if strcmp(operation,'get') && ~ismember(key,keys)
        varargout{1} = [];
    elseif strcmp(operation,'get')
        varargout{1} = conf.(key);
    else
        conf.(key) = val;

        f = fopen(fname,'w');
        f_list = sort(fields(conf));
        for j = 1:length(f_list)
            fprintf(f,'%s=%s\n',f_list{j},conf.(f_list{j}));
        end
        fclose(f);

    end


