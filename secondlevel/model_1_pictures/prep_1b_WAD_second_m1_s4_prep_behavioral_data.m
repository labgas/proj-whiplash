%% About this script
%
% Optional: Run these load and attach behavioral data from files (e.g., from Excel)            
%
% This script is an example script only.  You should modify it to fit your
% needs, which will depend on which types of behavioral/non-imaging data
% you have and what variables you want to store and analyze. The basics are
% desribed here:
%
% - Store behavioral data tables in any ad hoc format in DAT.BEHAVIOR.
% This can be a useful reference if you want to change/add custom analyses
% relating brain to behavior. You can create custom scripts that pull data 
% from .BEHAVIOR and use it in analyses.
%
% - Store a between-person grouping variable (e.g., patient vs. control,
% etc.) in DAT.BETWEENPERSON.group. This should be coded with values of 1
% and -1. Also add fields (see below) for names and colors associated with
% each group, and a description of what the 1/-1 codes mean.
% Some analyses consider this variable and run between-group contrasts
% and/or control for them in analyses of the entire sample (e.g.,
% "signature response" analyses).  
% SVM analyses can also be run that use the .group variable. See:
% prep_3d_run_SVM_betweenperson_contrasts and 
% c2b_SVM_betweenperson_contrasts  
%
% - If you have no binary group variable,  it is OK to leave the .group
% field empty. 
%
% - If you have continuous variable(s) instead of a binary group variable,
% you can enter a continuous variable in .group (for now!) -- this script
% uses that continuous variable:  (it may cause problems with other scripts
% that assume binary .group data, and may be changed in future versions):
% prep_3a_run_second_level_regression_and_save
%
% - Instead of a single .group variable to be tested with all
% conditions/contrasts, you can also enter different variables for each
% condition and/or contrast.  This is useful if you want to correlate each
% contrast with a different behavioral variable (maybe, e.g., for each contrast,
% reaction time differences for the same contrast). 
% If so, enter DAT.BETWEENPERSON.conditions and/or
% DAT.BETWEENPERSON.contrasts.  These should be cell arrays with one cell
% per condition or contrast.  Each cell contains a matlab "table" object
% with the data and possibly subject IDs (see below).
% 
% - You can run this script as part of a workflow (prep_1...
% prep_2...prep_3 etc)
% You can also run the script AFTER you've prepped all the imaging data, just
% to add behavioral data to the existing DAT structure.  If so, make sure
% you RELOAD the existing DAT structure with b_reload_saved_matfiles.m
% before you run this script.  Otherwise, if you create a new DAT
% structure, important information saved during the data prep (prep_2...,
% prep_3...) process will be missing, and you will need to re-run the whole
% prep sequence.
%
% @LUKASVO76 NOTES
% - this is a custom script based on the two prep_1b CANlab example
%   scripts


%% CUSTOM CODE: READ behavioral data from files. These data will be put into standard structure compatible with analyses
% phenotype.xlsx BIDS file
rootdir = 'C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI';
behavioral_data_filename = 'phenotype_subjects.xlsx';
phenodir = fullfile(rootdir,'phenotype');
behavioral_fname_path = fullfile(phenodir, behavioral_data_filename);

if ~exist(behavioral_fname_path, 'file'), fprintf(1, 'CANNOT FIND FILE: %s\n',behavioral_fname_path); end

behavioral_data_table = readtable(behavioral_fname_path,'TreatAsEmpty','n/a'); % read .xlsx file into Matlab table format
behavioral_data_table = sortrows(behavioral_data_table,'subject_id'); % sort table on subject_id variable
% behavioral_data_table((behavioral_data_table.included == 0),:) = []; % delete rows of excluded subjects
behavioral_data_table.id = [1:height(behavioral_data_table)]'; % create consecutively numbered subject idx

% calculate contrasts between ratings, and zscore them
% NOTE: respect the order of the raw rating variables in the table - first NA, then symptoms)
behavioral_data_table.pain_high_neu = zscore((behavioral_data_table.Pain_highfear - behavioral_data_table.Pain_control),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.pain_mod_neu = zscore((behavioral_data_table.Pain_moderatefear - behavioral_data_table.Pain_control),0,'omitnan');
behavioral_data_table.pain_high_mod = zscore((behavioral_data_table.Pain_highfear - behavioral_data_table.Pain_moderatefear),0,'omitnan'); 
behavioral_data_table.pain_fear_neu = zscore((((behavioral_data_table.Pain_highfear + behavioral_data_table.Pain_moderatefear)./2) - behavioral_data_table.Pain_control),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.worry_high_neu = zscore((behavioral_data_table.Worry_highfear - behavioral_data_table.Worry_control),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.worry_mod_neu = zscore((behavioral_data_table.Worry_moderatefear - behavioral_data_table.Worry_control),0,'omitnan');
behavioral_data_table.worry_high_mod = zscore((behavioral_data_table.Worry_highfear - behavioral_data_table.Worry_moderatefear),0,'omitnan'); 
behavioral_data_table.worry_fear_neu = zscore((((behavioral_data_table.Worry_highfear + behavioral_data_table.Worry_moderatefear)./2) - behavioral_data_table.Worry_control),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.anxious_high_neu = zscore((behavioral_data_table.anxious_highfear - behavioral_data_table.anxious_control),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.anxious_mod_neu = zscore((behavioral_data_table.anxious_moderatefear - behavioral_data_table.anxious_control),0,'omitnan');
behavioral_data_table.anxious_high_mod = zscore((behavioral_data_table.anxious_highfear - behavioral_data_table.anxious_moderatefear),0,'omitnan'); 
behavioral_data_table.anxious_fear_neu = zscore((((behavioral_data_table.anxious_highfear + behavioral_data_table.anxious_moderatefear)./2) - behavioral_data_table.anxious_control),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.avoidance_high_neu = zscore((behavioral_data_table.avoidance_highfear - behavioral_data_table.avoidance_control),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.avoidance_mod_neu = zscore((behavioral_data_table.avoidance_moderatefear - behavioral_data_table.avoidance_control),0,'omitnan');
behavioral_data_table.avoidance_high_mod = zscore((behavioral_data_table.avoidance_highfear - behavioral_data_table.avoidance_moderatefear),0,'omitnan'); 
behavioral_data_table.avoidance_fear_neu = zscore((((behavioral_data_table.avoidance_highfear + behavioral_data_table.avoidance_moderatefear)./2) - behavioral_data_table.avoidance_control),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1

% zscore raw ratings too
% NOTE: order does not matter here as variables already exist in table
% behavioral_data_table.NA_neg = zscore(behavioral_data_table.NA_neg,0,'omitnan');
% behavioral_data_table.NA_neu = zscore(behavioral_data_table.NA_neu,0,'omitnan');
% behavioral_data_table.NA_pos = zscore(behavioral_data_table.NA_pos,0,'omitnan');
% behavioral_data_table.symptoms_neg = zscore(behavioral_data_table.symptoms_neg,0,'omitnan');
% behavioral_data_table.symptoms_neu = zscore(behavioral_data_table.symptoms_neu,0,'omitnan');
% behavioral_data_table.symptoms_pos = zscore(behavioral_data_table.symptoms_pos,0,'omitnan');

% Add to DAT for record, and flexible use later
DAT.BEHAVIOR.behavioral_data_table = behavioral_data_table;


%% INITIALIZE GROUP VARIABLE

% Initialize empty variables
DAT.BETWEENPERSON = [];

% Single group variable, optional, for convenience
% These fields are mandatory, but they can be empty
% If group variable and between-person variables vary by
% condition/contrast, leave these empty
% -------------------------------------------------------------------------
DAT.BETWEENPERSON.group = [];
DAT.BETWEENPERSON.groupnames = {};
DAT.BETWEENPERSON.groupcolors = {};


%% INITIALIZE CONDITION/CONTRAST-SPECIFIC BETWEEN-PERSON DATA TABLES

% Cell array with table of group (between-person) variables for each
% condition, and for each contrast.
% If variables are entered:
% 1. they will be controlled for in analyses of the
% overall condition/contrast effects.
% 2. Analyses relating conditions/contrasts to these variables will be
% performed.
% If no variables are entered (empty elements), only
% within-person/whole-group effects will be analyzed.

% Between-person variables influencing each condition
% Table of [n images in condition x q variables]
% names in table can be any valid name.

DAT.BETWEENPERSON.conditions = cell(1, length(DAT.conditions));
[DAT.BETWEENPERSON.conditions{:}] = deal(table());  % empty tables

% Between-person variables influencing each condition
% Table of [n images in contrast x q variables]
% names in table can be any valid name.

DAT.BETWEENPERSON.contrasts = cell(1, length(DAT.contrastnames));
[DAT.BETWEENPERSON.contrasts{:}] = deal(table());  % empty tables


%% CUSTOM CODE: TRANSFORM INTO between_design_table
%
% Create a table for each condition/contrast with the between-person design, or leave empty.
%
% a variable called 'id' contains subject identfiers.  Other variables will
% be used as regressors.  Variables with only two levels should be effects
% coded, with [1 -1] values.

id = DAT.BEHAVIOR.behavioral_data_table.id;
group = DAT.BEHAVIOR.behavioral_data_table.group;
covs = DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames(contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,'pain','IgnoreCase',true) | contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,'worry','IgnoreCase',true) | contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,'anxious','IgnoreCase',true) | contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,'avoidance','IgnoreCase',true)); % we want to be able to use the calculated (contrasts between) ratings as covariates for the respective brain conditions/contrasts flexibly

for i = 1:length(DAT.conditions)

DAT.BETWEENPERSON.conditions{i}.group = group;
DAT.BETWEENPERSON.conditions{i}.pain_rating = DAT.BEHAVIOR.behavioral_data_table.(covs{i}); % we include the symptom ratings for each condition here, to be able to include these as covariate in analyses on conditions!
DAT.BETWEENPERSON.conditions{i}.worry_rating = DAT.BEHAVIOR.behavioral_data_table.(covs{3+i});
DAT.BETWEENPERSON.conditions{i}.anxious_rating = DAT.BEHAVIOR.behavioral_data_table.(covs{6+i});
DAT.BETWEENPERSON.conditions{i}.avoidance_rating = DAT.BEHAVIOR.behavioral_data_table.(covs{9+i});

end

for j = 1:length(DAT.contrasts)

DAT.BETWEENPERSON.contrasts{j}.group = group;
DAT.BETWEENPERSON.contrasts{j}.delta_pain_rating = DAT.BEHAVIOR.behavioral_data_table.(covs{(length(DAT.conditions)*4)+j}); % we include the calculated contrasts between conditions on symptom ratings here, they come after the NA and symptom ratings for conditions, and the NA ratings for contrasts, to include those as covariate in our analysis on contrasts
DAT.BETWEENPERSON.contrasts{j}.delta_worry_rating = DAT.BEHAVIOR.behavioral_data_table.(covs{(length(DAT.conditions)*4)+3+j});
DAT.BETWEENPERSON.contrasts{j}.delta_anxious_rating = DAT.BEHAVIOR.behavioral_data_table.(covs{(length(DAT.conditions)*4)+6+j});
DAT.BETWEENPERSON.contrasts{j}.delta_avoidance_rating = DAT.BEHAVIOR.behavioral_data_table.(covs{(length(DAT.conditions)*4)+9+j});

end

DAT.BETWEENPERSON.group = group; % @lukasvo76: ALWAYS define this group variable if you have different groups of subjects, even when you have already defined them in the conditions{i}.group or contrasts{j}.group fields above!
DAT.BETWEENPERSON.groupnames = {'patient' 'control'}; % @lukasvo76: IMPORTANT to put the name of the group coded 1 before the name of the group coded -1 for correct labeling in later scripts comparing groups!
DAT.BETWEENPERSON.groupcolors = {[.7 .3 .5] [.3 .5 .7]};

%% Check DAT, print warnings, save DAT structure

if ~isfield(DAT, 'conditions') 
    printhdr('Incomplete DAT structure');
    disp('The DAT field is incomplete. Run prep_1_set_conditions_contrasts_colors before running prep_1b...')
end
    
if isfield(DAT, 'SIG_conditions') && isfield(DAT, 'gray_white_csf')
    % Looks complete, we already have data, no warnings 
else
    printhdr('DAT structure ready for data prep');
    disp('DAT field does not have info from prep_2, prep_3, or prep_4 sequences');
    disp('prep_2/3/4 scripts should be run before generating results.');
end

printhdr('Save results');

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, 'DAT');

