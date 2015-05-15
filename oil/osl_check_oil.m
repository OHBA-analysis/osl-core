function oil = osl_check_oil(oilin)

% oil = osl_check_oil(oil)
%
% Checks an OIL (OHBA's ICA pipeLine) struct has all the appropriate 
% settings, which can then be passed to osl_run_oil to do an OIL 
% analysis. Throws an error if any required inputs are missing, fills other 
% settings with default values.
%
% OIL can analysis task-positive and resting state data. However, the user
% must specify the type of analysis using:
%
% oil.paradigm = 'task'; OR oil.paradigm = 'rest';
%
% OIL has six stages with the following required inputs
%
% For oil.source_recon:
%
% oil.source_recon.conditions; list of conditions to include in the
% analysis, e.g. oil.source_recon.conditions={'Motorbike','Neutral face'}; 
% oil.source_recon.D_continuous; list of continuous time SPM MEEG object
% file paths to run the analysis on (list order should correspond to
% oil.source_recon.mri and oil.source_recon.D_epoched fields if provided) 
% AND/OR
% oil.source_recon.D_epoched; list of epoched SPM MEEG object file paths to
% run the analysis on (list order should correspond to oil.source_recon.mri
% and oil.source_recon.D_epoched fields if provided) 
%
% e.g. oil.source_recon.D_epoched{1}='subject1'; oil.source_recon.D_epoched{2}='subject2'
% oil.source_recon.time_range; Time range (from to) within trial, in secs,
% need to all be the same duration and one for each condition.   
% oil.source_recon.freq_range; Frequency range (from to) in Hz
%
% For oil.enveloping:
%
% NONE
%
% For oil.concat_subs:
%
% NONE
%
% For oil.ica:
%
% NONE
%
% For oil.ica_first_level:
%
% oil.ica_first_level.design_matrix_summary & oil.ica.contrast are required
% to run first level statistics on ICA outputs. See osl_example_oil.m for
% example inputs.
%
% For oil.ica_group_level:
%
% NONE
%
% Optional inputs:
%
% See inside this function (e.g. use "type osl_check_oil") to see the other
% optional settings
%
% HL 250313
%
% Version 1.1
%
% henry.luckhoo@trinity.ox.ac.uk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% overall settings:
try oil.to_do    = oilin.to_do;    oilin = rmfield(oilin,'to_do');     catch, oil.to_do=[1 1 1 1 0 0];                                                                            end;
try oil.do_plots = oilin.do_plots; oilin = rmfield(oilin,'do_plots');  catch, oil.do_plots=0;                                                                                     end;
try oil.paradigm = oilin.paradigm; oilin = rmfield(oilin,'paradigm');  catch, error('User must specify the paradigm type: oil.paradigm, It can be either ''rest'' or ''task''.'); end;
% allow for a name
if isfield(oilin, 'fname'), oil.fname = oilin.fname; oilin = rmfield(oilin, 'fname'); end

%% Source recon from oat
oil.source_recon = oilin.source_recon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Enveloping Settings 

oil.enveloping=[]; if ~isfield(oilin,'enveloping'); oilin.enveloping = struct; end
try oil.enveloping.window_length            = oilin.enveloping.window_length;            oilin.enveloping = rmfield(oilin.enveloping,'window_length');            catch, oil.enveloping.window_length             = 1;                                                                                           end;    
try oil.enveloping.ss                       = oilin.enveloping.ss;                       oilin.enveloping = rmfield(oilin.enveloping,'ss');        catch, oil.enveloping.ss         = 4;                                                                                           end;
try oil.enveloping.gridstep                 = oilin.enveloping.gridstep;                 oilin.enveloping = rmfield(oilin.enveloping,'gridstep');                 catch, oil.enveloping.gridstep                  = 8;                                                                                           end;
try oil.enveloping.timewindow               = oilin.enveloping.timewindow;               oilin.enveloping = rmfield(oilin.enveloping,'timewindow');               catch, oil.enveloping.timewindow                = 'all';                                                                                       end;
try oil.enveloping.name                     = oilin.enveloping.name;                     oilin.enveloping = rmfield(oilin.enveloping,'name');                     catch, oil.enveloping.name                      = [num2str(oil.enveloping.window_length) 's_win_' num2str(oil.enveloping.gridstep) 'mm_grid']; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Concatenation Settings 

oil.concat_subs=[]; if ~isfield(oilin,'concat_subs'); oilin.concat_subs = struct; end
try oil.concat_subs.sessions_to_do     = oilin.concat_subs.sessions_to_do;      oilin.concat_subs = rmfield(oilin.concat_subs,'sessions_to_do');      catch, oil.concat_subs.sessions_to_do = oilin.source_recon.sessions_to_do;  end;
try oil.concat_subs.normalise_subjects = oilin.concat_subs.normalise_subjects;  oilin.concat_subs = rmfield(oilin.concat_subs,'normalise_subjects');  catch, oil.concat_subs.normalise_subjects  = 0;                           end;
try oil.concat_subs.demean_vox         = oilin.concat_subs.demean_vox;          oilin.concat_subs = rmfield(oilin.concat_subs,'demean_vox');          catch, oil.concat_subs.demean_vox          = 1;                           end;
try oil.concat_subs.name               = oilin.concat_subs.name;                oilin.concat_subs = rmfield(oilin.concat_subs,'name');                catch, oil.concat_subs.name                = 'all_subs';                  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ICA Settings 

oil.ica=[]; if ~isfield(oilin,'ica'); oilin.ica = struct; end
try oil.ica.temp_or_spat       = oilin.ica.temp_or_spat;        oilin.ica = rmfield(oilin.ica,'temp_or_spat');        catch, oil.ica.temp_or_spat  = 'temporal';                                                                                          end;
try oil.ica.use_gm_mask        = oilin.ica.use_gm_mask;         oilin.ica = rmfield(oilin.ica,'use_gm_mask');         catch, oil.ica.use_gm_mask   = 0;                                                                                                   end;
try oil.ica.num_ics            = oilin.ica.num_ics;             oilin.ica = rmfield(oilin.ica,'num_ics');             catch, oil.ica.num_ics       = 25;                                                                                                  end;
try oil.ica.icasso_its         = oilin.ica.icasso_its;          oilin.ica = rmfield(oilin.ica,'icasso_its');          catch, oil.ica.icasso_its    = 0;                                                                                                   end;
try oil.ica.last_eig           = oilin.ica.last_eig;            oilin.ica = rmfield(oilin.ica,'last_eig');            catch, oil.ica.last_eig      = oil.ica.num_ics;                                                                                     end;
try oil.ica.other_oils         = oilin.ica.other_oils;          oilin.ica = rmfield(oilin.ica,'other_oils');          catch,                                                                                                                              end;
try oil.ica.nuisance_oil_names = oilin.ica.nuisance_oil_names;  oilin.ica = rmfield(oilin.ica,'nuisance_oil_names');  catch,                                                                                                                              end;
try oil.ica.normalise_vox      = oilin.ica.normalise_vox;       oilin.ica = rmfield(oilin.ica,'normalise_vox');       catch, oil.ica.normalise_vox = 0;                                                                                                   end;
try oil.ica.initGuess          = oilin.ica.initGuess;           oilin.ica = rmfield(oilin.ica,'initGuess');           catch, oil.ica.initGuess = 0;                                                                                                       end;
try oil.ica.name               = oilin.ica.name;                oilin.ica = rmfield(oilin.ica,'name');                catch, oil.ica.name          = [oil.ica.temp_or_spat '_ica_' num2str(oil.ica.num_ics) 'comps_' num2str(oil.ica.icasso_its) '_its']; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  First Level ICA Stats Settings 

oil.ica_first_level=[]; if ~isfield(oilin,'ica_first_level'); oilin.ica_first_level = struct; end
% required settings:
if strcmp(oil.paradigm,'task')
try oil.ica_first_level.design_matrix_summary = oilin.ica_first_level.design_matrix_summary; oilin.ica_first_level = rmfield(oilin.ica_first_level,'design_matrix_summary'); catch, error('design_matrix_summary not specified');                     end;
try oil.ica_first_level.contrast              = oilin.ica_first_level.contrast;              oilin.ica_first_level = rmfield(oilin.ica_first_level,'contrast');              catch, error('contrast not specified');                                  end;
end
% optional settings
try oil.ica_first_level.name                  = oilin.ica_first_level.name;                  oilin.ica_first_level = rmfield(oilin.ica_first_level,'name');                  catch, oil.ica_first_level.name           = 'ica_first_level';           end;
try oil.ica_first_level.time_range            = oilin.ica_first_level.time_range;            oilin.ica_first_level = rmfield(oilin.ica_first_level,'time_range');            catch, oil.ica_first_level.time_range     = oil.source_recon.time_range; end;
try oil.ica_first_level.use_robust_glm        = oilin.ica_first_level.use_robust_glm;        oilin.ica_first_level = rmfield(oilin.ica_first_level,'use_robust_glm');        catch, oil.ica_first_level.use_robust_glm = 0;                           end; % do robust GLM (uses bisquare via Matlab robustfit fn)
try oil.ica_first_level.cope_type             = oilin.ica_first_level.cope_type;             oilin.ica_first_level = rmfield(oilin.ica_first_level,'cope_type');             catch, oil.ica_first_level.cope_type      = 'cope';                      end; % cope type to input to group GLM (from first level), set to 'coape', 'cope', or 'acope'
try oil.ica_first_level.spatial_basis_set     = oilin.ica_first_level.spatial_basis_set;     oilin.ica_first_level = rmfield(oilin.ica_first_level,'spatial_basis_set');     catch,                                                                   end; % Spatial Basis Set from alternative analysis e.g. weights-normalised MEG or fMRI.
% try oil.ica_first_level.use_nonorm_data       =
% oilin.ica_first_level.use_nonorm_data;       oilin.ica_first_level =
% rmfield(oilin.ica_first_level,'use_nonorm_data');       catch,
% oil.ica_first_level.use_nonorm_data = 1;                          end; %
% To get single subject masks if needed for parcellation - this field no
% longer has an effect. Normalised data are used by default. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Group Level ICA Stats Settings 

oil.ica_group_level=[]; if ~isfield(oilin,'ica_group_level'); oilin.ica_group_level = struct; end
try oil.ica_group_level.group_design_matrix               = oilin.ica_group_level.group_design_matrix;               oilin.ica_group_level = rmfield(oilin.ica_group_level,'group_design_matrix');               catch, oil.ica_group_level.group_design_matrix               = ones(length(oil.concat_subs.sessions_to_do),1); end;
try oil.ica_group_level.group_contrast                    = oilin.ica_group_level.group_contrast;                    oilin.ica_group_level = rmfield(oilin.ica_group_level,'group_contrast');                    catch, oil.ica_group_level.group_contrast{1}                 = [1];                                            end;
try oil.ica_group_level.use_randomise                     = oilin.ica_group_level.use_randomise;                     oilin.ica_group_level = rmfield(oilin.ica_group_level,'use_randomise');                     catch, oil.ica_group_level.use_randomise                     = 0;                                              end;
try oil.ica_group_level.Npermutations                     = oilin.ica_group_level.Npermutations;                     oilin.ica_group_level = rmfield(oilin.ica_group_level,'Npermutations');                     catch, oil.ica_group_level.Npermutations                     = 5000;                                           end;
try oil.ica_group_level.name                              = oilin.ica_group_level.name;                              oilin.ica_group_level = rmfield(oilin.ica_group_level,'name');                              catch, oil.ica_group_level.name                              = 'ica_group_level';                              end;
try oil.ica_group_level.demean_designmatrix               = oilin.ica_group_level.demean_designmatrix;               oilin.ica_group_level = rmfield(oilin.ica_group_level,'demean_designmatrix');               catch, oil.ica_group_level.demean_designmatrix               = 0;                                              end;
try oil.ica_group_level.group_varcope_spatial_smooth_fwhm = oilin.ica_group_level.group_varcope_spatial_smooth_fwhm; oilin.ica_group_level = rmfield(oilin.ica_group_level,'group_varcope_spatial_smooth_fwhm'); catch, oil.ica_group_level.group_varcope_spatial_smooth_fwhm = 0;                                              end;
try oil.ica_group_level.cluster_threshold                 = oilin.ica_group_level.cluster_threshold;                 oilin.ica_group_level = rmfield(oilin.ica_group_level,'cluster_threshold');                 catch, oil.ica_group_level.cluster_threshold                 = 2.3;                                            end;
try oil.ica_group_level.comps2use                         = oilin.ica_group_level.comps2use;                         oilin.ica_group_level = rmfield(oilin.ica_group_level,'comps2use');                         catch, oil.ica_group_level.comps2use                         = 1:oil.ica.num_ics;                              end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% copy any results
try oil.results                     = oilin.results;                     oilin                 = rmfield(oilin,'results');                                          catch, end;
% Source Recon Results
try oil.source_recon.results_fnames = oilin.source_recon.results_fnames; oilin.source_recon    = rmfield(oilin.source_recon,'results_fnames');                      catch, end;
% Enveloping Results
try oil.enveloping.results          = oilin.enveloping.results;          oilin.enveloping      = rmfield(oilin.enveloping,'results');                               catch, end;
% Concatenation Results
try oil.concat_subs.results         = oilin.concat_subs.results;         oilin.concat_subs     = rmfield(oilin.concat_subs,'results');                              catch, end;
% ICA Results
try oil.ica.results                 = oilin.ica.results;                 oilin.ica             = rmfield(oilin.ica,'results');                                      catch, end;
% First-level Results
try oil.ica_first_level.results     = oilin.ica_first_level.results;     oilin.ica_first_level = rmfield(oilin.ica_first_level,'results');                          catch, end;
% Group-level Results
try oil.ica_group_level.results     = oilin.ica_group_level.results;     oilin.ica_group_level = rmfield(oilin.ica_group_level,'results');                          catch, end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check people haven't set any weird fields

weirdfields = fieldnames(oilin.enveloping);
if ~isempty(weirdfields)
    displayedWeirdFields = [];
    for iprint = 1:numel(weirdfields)
        displayedWeirdFields = [displayedWeirdFields weirdfields{iprint} ' '];
    end
    error([mfilename ':UnrecognisedFields'], ...
          ['The following oil.enveloping settings were not recognized ', ...
           'by osl_check_oil: \n ' displayedWeirdFields '\n']);
end % if ~isempty(weirdfields)

weirdfields = fieldnames(oilin.concat_subs);
if ~isempty(weirdfields)
    displayedWeirdFields = [];
    for iprint = 1:numel(weirdfields)
        displayedWeirdFields = [displayedWeirdFields weirdfields{iprint} ' '];
    end
    error([mfilename ':UnrecognisedFields'], ...
          ['The following oil.concat_subs settings were not recognized ', ...
           'by osl_check_oil: \n ' displayedWeirdFields '\n']);
end % if ~isempty(weirdfields)

weirdfields = fieldnames(oilin.ica);
if ~isempty(weirdfields)
    displayedWeirdFields = [];
    for iprint = 1:numel(weirdfields)
        displayedWeirdFields = [displayedWeirdFields weirdfields{iprint} ' '];
    end
    error([mfilename ':UnrecognisedFields'], ...
          ['The following oil.ica settings were not recognized ', ...
           'by osl_check_oil: \n ' displayedWeirdFields '\n']);
end % if ~isempty(weirdfields)

weirdfields = fieldnames(oilin.ica_first_level);
if ~isempty(weirdfields)
    displayedWeirdFields = [];
    for iprint = 1:numel(weirdfields)
        displayedWeirdFields = [displayedWeirdFields weirdfields{iprint} ' '];
    end
    error([mfilename ':UnrecognisedFields'], ...
          ['The following oil.ica_first_level settings were not recognized ', ...
           'by osl_check_oil: \n ' displayedWeirdFields '\n']);
end % if ~isempty(weirdfields)

weirdfields = fieldnames(oilin.ica_group_level);
if ~isempty(weirdfields)
    displayedWeirdFields = [];
    for iprint = 1:numel(weirdfields)
        displayedWeirdFields = [displayedWeirdFields weirdfields{iprint} ' '];
    end
    error([mfilename ':UnrecognisedFields'], ...
          ['The following oil.ica_group_level settings were not recognized ', ...
           'by osl_check_oil: \n ' displayedWeirdFields '\n']);
end % if ~isempty(weirdfields)

oilin = rmfield(oilin,'source_recon');
oilin = rmfield(oilin,'enveloping');
oilin = rmfield(oilin,'concat_subs');
oilin = rmfield(oilin,'ica');
oilin = rmfield(oilin,'ica_first_level');
oilin = rmfield(oilin,'ica_group_level');

try 
    oilin = rmfield(oilin,'fname'); 
catch
    % do nothing
end;

weirdfields = fieldnames(oilin);
if ~isempty(weirdfields)
    displayedWeirdFields = [];
    for iprint = 1:numel(weirdfields)
        displayedWeirdFields = [displayedWeirdFields weirdfields{iprint} ' '];
    end
    error([mfilename ':UnrecognisedFields'], ...
          ['The following oil settings were not recognized ', ...
           'by osl_check_oil: \n ' displayedWeirdFields '\n']);
end % if ~isempty(weirdfields)
