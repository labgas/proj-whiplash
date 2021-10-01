%% WAD_prep_s3b_correct_onset_files_model_2
%
% This is a quick and (hopefully not) dirty script to overwrite the
% onsets.mat files for model 2 of Iris' study, with separate event
% categories for imagine events following high, moderate, and neutral fear
% pictures
%
%__________________________________________________________________________
%
% authors: Lukas Van Oudenhove & Iris Coppieters
% date:   April, 2021
%
%__________________________________________________________________________
% @(#)% WAD_prep_s3b_correct_onset_files_model_2.m         v1.0        
% last modified: 2021/04/15


%% DEFINE AND CREATE DIRECTORIES, IDENTIFY SUBJECTS, AND SET OPTIONS
%--------------------------------------------------------------------------

% DIRECTORY STRUCTURE SETUP
% existing directories
rootdir = 'C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI';
rawdir = fullfile(rootdir,'rawdata');
derivdir = fullfile(rootdir,'derivatives\fmriprep');

models = {'model_1_pictures';'model_2_pictures_imagine'};

% create subject list from derivdir, and check whether it matches subjdirs
% in rawdir, otherwise throw error
subjs = dir(fullfile(derivdir,'sub-*'));
idx = [subjs.isdir]';
subjs = {subjs(idx).name}';

subjs_raw = dir(fullfile(rawdir,'sub-*'));
idx_raw = [subjs_raw.isdir]';
subjs_raw = {subjs_raw(idx_raw).name}';

if ~isequal(subjs,subjs_raw)
    error('subject directories in derivatives and rawdata do not match; check and correct before proceeding');
end

clear idx idx_raw subjs_raw

% define runs
runs = {'run-A';'run-B'};


%% LOOP OVER SUBJECTS
%--------------------------------------------------------------------------

for sub=1:size(subjs,1)
    
    % DEFINE SUBJECT LEVEL DIRS
    subjrawdir=fullfile(rawdir,subjs{sub},'\ses-1\func');
    subjderivdir=fullfile(derivdir,subjs{sub},'\ses-1\func');
    
    % DEFINE SUBJECT LEVEL FILENAMES
    fmriprep_noisefiles = dir(fullfile(subjderivdir,'*pic*desc-confounds*.json'));
    fmriprep_noisefiles = {fmriprep_noisefiles(:).name}';
    runnames = char(fmriprep_noisefiles);
    runnames = runnames(:,1:24); % this is study-specific, depends on length of subjectname and taskname - could be adapted based on regexp, but easy enough to adapt per study

    % CALCULATE AND EXTRACT CONFOUND REGRESSORS AND ONSETS, AND WRITE TO
    % FILE
    for run=1:2
        
        % define subdir for this run
        rundir = fullfile(subjderivdir,runs{run});

        % EVENTS FILES
        % read events.tsv files with onsets, durations, and trial type
        eventsfiles = dir(fullfile(subjrawdir,strcat(runnames(run,:),'*events.tsv')));
        eventsfiles = {eventsfiles(:).name}';
        
        O=readtable(fullfile(subjrawdir,eventsfiles{2}),'FileType', 'text', 'Delimiter', 'tab');
        O.trial_type = categorical(O.trial_type);
        
        % save events file as .mat file
        filename_events = strcat(rundir,'\onsets_',runnames(run,:),'_',models{2});
        save(filename_events,'O');
        clear O filename_events
        
    end % for loop runs
    
end % for loop subjects 