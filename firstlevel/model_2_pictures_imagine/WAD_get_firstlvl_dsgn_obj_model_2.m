%% WAD_get_firstlvl_dsgn_obj_model_2
%
% This script contains a function that defines a number of fields of 
% the CANlab style first level DSGN structure array for the WAD study 
% on neck-related pictures of movements (taken from the PFActS-C) that are 
% presented and vary in the degree to which they are perceived as fearful.
% This function is used in WAD_first_m2_s1_spm_fit_firstlvl_models.m to
% run first level analysis using CANlab tools
%
% IMPORTANT NOTE: function and script are study-specific!
%
% See canlab_glm_single_subject('dsgninfo')
% OR Github\CanlabCore\CanlabCore\GLM_Batch_tools\canlab_glm_dsgninfo.txt 
% for details on first level analysis using CANlab tools, 
% including options and defaults
%
% This script is adapted by @lukasvo76 from the scripts
% 1) get_firslvl_dsgn_obj.m and get_single_trial_dsgn_obj.m by @bogpetre on
% Google Drive\CANlab\CANLAB Lab Member Documents\GLM_batch_tools\
% bogdan_paingen\classic_glm_contrasts
% 2) MPA2_set_design_model1_blanca.m by @martaceko on
% Google Drive\CANlab\CANLAB Lab Member Documents\GLM_batch_tools\
% Marta_MPA2\MPA2code_1stlevel\Code
% contact @lukasvo76 if you need those original scripts
% 
% DEPENDENCIES ON YOUR MATLAB PATH
% a) SPM12
% b) CANlab tools cloned from Github (see canlab.github.io)
% 
% INPUTS 
% none - you need to adapt the function script below to your study
%
% OUTPUT
% CANlab style first level DSGN structure array in your Matlab workspace, 
% to be used by WAD_first_m2_s1_spm_fit_firstlvl_models.m
%
%__________________________________________________________________________
%
% authors: Lukas Van Oudenhove & Iris Coppieters
%
% date:   April, 2021
%
%__________________________________________________________________________
% @(#)% WAD_get_firstlvl_dsgn_obj_model_2.m         v1.0       
% last modified: 2021/04/16
%
%
%% function code
function DSGN = WAD_get_firstlvl_dsgn_obj_model_2()
    % INPUT
    % required fields
    DSGN.metadata = "WAD first level analysis model 2, i.e. modeling 3 conditions for high, moderate, neutral fear, and 3 conditions for imagine periods following high, moderate, and neutral fear pictures, all as 3 second events"; % field for annotation with study info, or whatever you like
    DSGN.modeldir = 'C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI\firstlevel\model_2_pictures_imagine'; % directory where you want to write first level results
    DSGN.subjects = {}; % sets up empty cell array field for subjects in structure array DSGN
        fnames = dir('C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI\derivatives\fmriprep\sub*'); % get subject names from root directory with subject data 
        fnames=fnames([fnames.isdir]'); % duplicate names of subjects in fmriprep dir because of .html files with same subject names as folders
            for i = 1:size(fnames,1)
                this_f = fnames(i);
                DSGN.subjects = [DSGN.subjects, [this_f.folder, '\', this_f.name]]; % cell array of subject directories (absolute paths)
            end
    DSGN.funcnames = {'ses-1\func\run-A\s6*.nii',...
        'ses-1\func\run-B\s6*.nii'}; % cell array (one cell per session) of paths to functional files, relative to absolute path specific in DSGN.subjects
    % optional fields
    DSGN.concatenation = {}; % default: none; cell array of arrays of runs to concatenate; see documentation for when to concatenate, and how it works exactly
    DSGN.allowmissingfunc = true; % default: false; true will prevent erroring out when functional file is missing for at least one run is missing for at least one subject
    DSGN.customrunintercepts = {}; % default: none; will only work if DSGN.concatenation is specified; cell array of vectors specifying custom intercepts
    
    % PARAMETERS
    DSGN.tr = 1; % repetition time (TR) in seconds
    DSGN.hpf = 128; % high pass filter in seconds; SPM default is 128, CANlab default is 180 since the brain response to pain stimuli last long and variance may be lost at shorter lengths, use scn_spm_design_check output, and review the SPM.mat in spm for diagnostics; 
    % STUDY-SPECIFIC: in this study, with rather short non-pain events and
    % relatively short ITI, we stick with the SPM default
    DSGN.fmri_t = 56; % microtime resolution - t=number of slices since we did slice timing; spm (and CANlab) default 16 can be kept for multiband w/o slice timing
    DSGN.fmri_t0 = 28; % microtime onset - reference slice used in slice timing correction; spm (and CANlab) default 1 can be kept for multiband w/o slice timing
    
    % MODELING
    % required fields: cell array (one cell per session (i.e. run in our case)) of cell arrays (one cell per condition) of MAT-file names; if only one session is specified, it will be applied to all sessions
    c=0;
    c=c+1;DSGN.conditions{c}={'high_fear' 'moderate_fear' 'neutral_fear' 'imagine_high' 'imagine_moderate' 'imagine_neutral'};
    c=c+1;DSGN.conditions{c}={'high_fear' 'moderate_fear' 'neutral_fear' 'imagine_high' 'imagine_moderate' 'imagine_neutral'};
    % optional fields
%     DSGN.pmods = {{}}; % cell array (one cell per session) of cell arrays (one cell per condition) of cell arrays (one cell per modulator) of MAT-file names
%     DSGN.convolution; default hrf.derivs = [0 0]; structure specifying the convolution to use for conditions different fields required depending on convolution type; 
%     DSGN.ar1 = false; % autoregressive AR(1) to model serial correlations; SPM default is true, CANlab default is false, Tor recommends turning autocorrelation off, because this algorithm pools across the whole brain, and does not perform well in some situations; if you are performing a group analysis, the autocorrelation problem is not as concerning
    DSGN.notimemod = true; % default: false; if true, turn off time modulation of conditions, i.e. when you do not expect linear trends over time
%     DSGN.singletrials = {{}}; % a cell array (1 cell per session) of cell arrays (1 cell per condition) of (corresponding to DSGN.conditions) of true/false values indicating whether to convert specified condition to set of single trial conditions
%     DSGN.singletrialsall = false; % default: false; if true, set DSGN.singletrials to true for all conditions
    DSGN.modelingfilesdir = 'model_2_pictures_imagine'; % name of subfolder which will be created within directory containing functional files where .mat files containing fields of DSGN structure will be saved
%     DSGN.allowemptycond = false; % default:false; if true, allow empty conditions
%     DSGN.allowmissingcondfiles = false; % default:false; if true, throw warning instead of error when no file(s) are found corresponding to a MAT-file name/wildcard
    DSGN.multireg = 'noise_regs'; % specify name for matfile with noise parameters you want to save
    
    % CONTRASTS
    % required fields
    % cell array (one cell per contrast) of contrast definitions
    c=0;
    % picture only contrasts
    c=c+1;DSGN.contrasts{c} = {{'high_fear'}}; % CON_0001
    DSGN.contrastnames{c} = 'high_fear'; % not needed strictly, because this will be automatically generated for standard contrasts like this (high fear versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'moderate_fear'}}; % CON_0002
    DSGN.contrastnames{c} = 'moderate_fear'; % not needed strictly, because this will be automatically generated for standard contrasts like this (moderate fear versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'neutral_fear'}}; % CON_0003
    DSGN.contrastnames{c} = 'neutral_fear'; % not needed strictly, because this will be automatically generated for standard contrasts like this (neutral versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'high_fear'} {'neutral_fear'}}; % CON_0004
    DSGN.contrastnames{c} = 'high_vs_neutral_fear'; % not needed strictly, because this will be automatically generated for standard contrasts like this (high fear versus neutral)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'moderate_fear'} {'neutral_fear'}}; % CON_0005
    DSGN.contrastnames{c} = 'moderate_vs_neutral_fear'; % not needed strictly, because this will be automatically generated for standard contrasts like this (moderate fear versus neutral)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'high_fear'} {'moderate_fear'}}; % CON_0006
    DSGN.contrastnames{c} = 'high_vs_moderate_fear'; % not needed strictly, because this will be automatically generated for standard contrasts like this (high fear versus moderate fear)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'high_fear'} {'moderate_fear'} {'neutral_fear'}}; % CON_0007
    DSGN.contrastnames{c} = 'fear (high + moderate) versus neutral'; % needed because this is not a standard contrast
    DSGN.contrastweights{c} = [1 1 -2];
    % imagine only contrasts
    c=c+1;DSGN.contrasts{c} = {{'imagine_high'}}; % CON_0008
    DSGN.contrastnames{c} = 'imagine_high'; % not needed strictly, because this will be automatically generated for standard contrasts like this (imagine high versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'imagine_moderate'}}; % CON_0009
    DSGN.contrastnames{c} = 'imagine_moderate'; % not needed strictly, because this will be automatically generated for standard contrasts like this (imagine moderate versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'imagine_neutral'}}; % CON_0010
    DSGN.contrastnames{c} = 'imagine_neutral'; % not needed strictly, because this will be automatically generated for standard contrasts like this (imagine neutral versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'imagine_high'} {'imagine_neutral'}}; % CON_0011
    DSGN.contrastnames{c} = 'high_vs_neutral_imagine'; % not needed strictly, because this will be automatically generated for standard contrasts like this (high fear versus neutral)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'imagine_moderate'} {'imagine_neutral'}}; % CON_0012
    DSGN.contrastnames{c} = 'moderate_vs_neutral_imagine'; % not needed strictly, because this will be automatically generated for standard contrasts like this (imagine high versus imagine neutral)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'imagine_high'} {'imagine_moderate'}}; % CON_0013
    DSGN.contrastnames{c} = 'high_vs_moderate_imagine'; % not needed strictly, because this will be automatically generated for standard contrasts like this (imagine high versus imagine moderate)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'imagine_high'} {'imagine_moderate'} {'imagine_neutral'}}; % CON_0014
    DSGN.contrastnames{c} = 'imagine (high + moderate) versus neutral'; % needed because this is not a standard contrast
    DSGN.contrastweights{c} = [1 1 -2];
    % picture + imagine contrasts
    c=c+1;DSGN.contrasts{c} = {{'high_fear'} {'imagine_high'}}; % CON_0015
    DSGN.contrastnames{c} = 'fear_imagine_high'; % not needed strictly, because this will be automatically generated for standard contrasts like this (high fear versus baseline)
    DSGN.contrastweights{c} = [1 1];
    c=c+1;DSGN.contrasts{c} = {{'moderate_fear'} {'imagine_moderate'}}; % CON_0016
    DSGN.contrastnames{c} = 'fear_imagine_moderate'; % not needed strictly, because this will be automatically generated for standard contrasts like this (moderate fear versus baseline)
    DSGN.contrastweights{c} = [1 1];
    c=c+1;DSGN.contrasts{c} = {{'neutral_fear'} {'imagine_neutral'}}; % CON_0017
    DSGN.contrastnames{c} = 'fear_imagine_neutral'; % not needed strictly, because this will be automatically generated for standard contrasts like this (neutral versus baseline)
    DSGN.contrastweights{c} = [1 1];
    c=c+1;DSGN.contrasts{c} = {{'high_fear'} {'imagine_high'} {'neutral_fear'} {'imagine_neutral'}}; % CON_0018
    DSGN.contrastnames{c} = 'high_vs_neutral_fear_imagine'; % not needed strictly, because this will be automatically generated for standard contrasts like this (high fear versus neutral)
    DSGN.contrastweights{c} = [1 1 -1 -1];
    c=c+1;DSGN.contrasts{c} = {{'moderate_fear'} {'imagine_moderate'} {'neutral_fear'} {'imagine_neutral'}}; % CON_0019
    DSGN.contrastnames{c} = 'moderate_vs_neutral_fear_imagine'; % not needed strictly, because this will be automatically generated for standard contrasts like this (moderate fear versus neutral)
    DSGN.contrastweights{c} = [1 1 -1 -1];
    c=c+1;DSGN.contrasts{c} = {{'high_fear'} {'imagine_high'} {'moderate_fear'} {'imagine_moderate'}}; % CON_0020
    DSGN.contrastnames{c} = 'high_vs_moderate_fear_imagine'; % not needed strictly, because this will be automatically generated for standard contrasts like this (high fear versus moderate fear)
    DSGN.contrastweights{c} = [1 1 -1 -1];
    c=c+1;DSGN.contrasts{c} = {{'high_fear'} {'imagine_high'} {'moderate_fear'} {'imagine_moderate'} {'neutral_fear'} {'imagine_neutral'}}; % CON_0021
    DSGN.contrastnames{c} = 'fear_imagine (high + moderate) versus neutral'; % needed because this is not a standard contrast
    DSGN.contrastweights{c} = [1 1 1 1 -2 -2];
end