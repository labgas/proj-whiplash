% This script runs single-level mediation models for Iris' proj-WAD, 
% with the following variables
% X: group (patients versus controls)
% M: different options
% 1. signature responses
% a. NPS response (including subregions) for each contrast
% b. FAPS response for each contrast
% 2. voxel-wise brain mediator
% 3. Principal Directions of Mediation (PDM) multivariate brain mediator
% Y: somatic symptom ratings for each contrast
%
% It requires all study-specific and most core prep_ scripts for proj-WAD 
% and CANlab_help examples (LaBGAS fork, https://github.com/labgas/CANlab_help_examples) 
% to be run first, or load your saved .mat files if you ran these
% scripts earlier!
%
% Options for masking and scaling are currently still hard-coded, but 
% can easily be set as options in a2_emosymp_m1_s2_set_default_options, 
% or as part of this script at a later stage if desired
%
% This script is adapted from the following CANlab scripts:
%
% https://www.dropbox.com/s/uvfhdstfveve393/s8_mediation_model_FAPS.m?dl=0 
% by @pkragel
% https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/Multivariate_Mediation_ExampleScript.m
% by @tor and @martin
% 
% see also 
% https://canlab.github.io/_pages/mediation_example_script_1/mediation_example_script_1.html
% https://canlab.github.io/_pages/mediation_brain_sample_report/mediation_brain_results_report.html
%
% @lukasvo76 @Dartmouth, March-April 2021
% 
% WORK IN PROGRESS
% FIGURES, TITLES, AND STRUCTURE OF HTML OUTPUT CAN BE IMPROVED
% LOOPING OVER SUBPATTERNS/REGIONS CAN BE IMPLEMENTED


%% DEFINE PATHS AND VARIABLES TO BE USED IN ALL ANALYSES
% -------------------------------------------------------------------------
addpath(genpath('C:\Users\lukas\Documents\GitHub\MediationToolbox')); % add CANlab mediation toolbox to your path
addpath(genpath('C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI\scripts_windows')); % add WAD scripts folder to your path - still on Dropbox rather than Github for now

a_WAD_second_m1_s1_set_up_paths_always_run_first

underscores = '_____________________________________________________________________';

load(fullfile(resultsdir,'image_names_and_setup.mat'));
load(fullfile(resultsdir,'data_objects.mat'));

mediationdir = fullfile(resultsdir,'mediation_analysis'); % a_ script needs to be run first (always!) so resultsdir is known
brain_mediationdir = fullfile(mediationdir,'brain_voxelwise');
pdm_mediationdir = fullfile(mediationdir,'brain_pdm');
if ~isfolder(mediationdir)
    mkdir(mediationdir);
end
if ~isfolder(brain_mediationdir)
    mkdir(brain_mediationdir);
end
if ~isfolder(pdm_mediationdir)
    mkdir(pdm_mediationdir);
end
cd(mediationdir);

% note: you do not need to remove subjects with NaNs for behavioral data
% here as the mediation function called below contains a nanremove command
idx = logical(DAT.BEHAVIOR.behavioral_data_table.included);
X = DAT.BETWEENPERSON.group; % define predictor
Ydat_pain = [DAT.BEHAVIOR.behavioral_data_table.pain_high_neu, DAT.BEHAVIOR.behavioral_data_table.pain_mod_neu, DAT.BEHAVIOR.behavioral_data_table.pain_high_mod, DAT.BEHAVIOR.behavioral_data_table.pain_fear_neu]; % define matrix of outcomes (subjects * contrasts)
Ydat_worry = [DAT.BEHAVIOR.behavioral_data_table.worry_high_neu, DAT.BEHAVIOR.behavioral_data_table.worry_mod_neu, DAT.BEHAVIOR.behavioral_data_table.worry_high_mod, DAT.BEHAVIOR.behavioral_data_table.worry_fear_neu];
Ydat_anxious = [DAT.BEHAVIOR.behavioral_data_table.anxious_high_neu, DAT.BEHAVIOR.behavioral_data_table.anxious_mod_neu, DAT.BEHAVIOR.behavioral_data_table.anxious_high_mod, DAT.BEHAVIOR.behavioral_data_table.anxious_fear_neu];
Ydat_avoidance = [DAT.BEHAVIOR.behavioral_data_table.avoidance_high_neu, DAT.BEHAVIOR.behavioral_data_table.avoidance_mod_neu, DAT.BEHAVIOR.behavioral_data_table.avoidance_high_mod, DAT.BEHAVIOR.behavioral_data_table.avoidance_fear_neu];
Ydat = {Ydat_pain, Ydat_avoidance, Ydat_worry, Ydat_anxious};
Ynames = {'pain','avoidance','worry','anxious'};
contrastnames = DAT.SIG_contrasts.raw.dotproduct.conditionnames; % get names of contrasts


%% NPS & FAPS
% -------------------------------------------------------------------------

%%
% <html><h3>PAIN</h3></html>

for i = 1:size(contrastnames,2) % loop over contrasts
    
    % define outcome for contrast i
    Y{i} = Ydat_pain(:,i);
    
    % define mediators
    M{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPS(:,i));
%     Mpos{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSpos(:,i));
%     Mneg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSneg(:,i));
    MFAPS{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.FAPS(:,i));
    MVIFS{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.VIFS(:,i));
    Mreddan{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.reddan(:,i));
    
    % run mediation models for each of the mediators
    [paths, toplevelstats, ~] = mediation(X,Y{i},M{i},'names',{'group',contrastnames{i},'NPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPS{i} = paths;
    toplevelstats_NPS{i} = toplevelstats;
    clear paths toplevelstats;
    
%     [paths, toplevelstats, ~] = mediation(X,Y{i},Mpos{i},'names',{'group',contrastnames{i},'NPSpos response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_NPSpos{i} = paths;
%     toplevelstats_NPSpos{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},Mneg{i},'names',{'group',contrastnames{i},'NPSneg response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_NPSneg{i} = paths;
%     toplevelstats_NPSneg{i} = toplevelstats;
%     clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MFAPS{i},'names',{'group',contrastnames{i},'FAPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_FAPS{i} = paths;
    toplevelstats_FAPS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MVIFS{i},'names',{'group',contrastnames{i},'VIFS response'},'boottop','plots');
    drawnow;snapnow;
    paths_VIFS{i} = paths;
    toplevelstats_VIFS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},Mreddan{i},'names',{'group',contrastnames{i},'Reddan threat response'},'boottop','plots');
    drawnow;snapnow;
    paths_reddan{i} = paths;
    toplevelstats_reddan{i} = toplevelstats;
    clear paths toplevelstats;
    
end


%%
% <html><h3>AVOIDANCE</h3></html>

for i = 1:size(contrastnames,2) % loop over contrasts
    
    % define outcome for contrast i
    Y{i} = Ydat_avoidance(:,i);
    
    % define mediators
    M{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPS(:,i));
%     Mpos{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSpos(:,i));
%     Mneg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSneg(:,i));
    MFAPS{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.FAPS(:,i));
    MVIFS{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.VIFS(:,i));
    Mreddan{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.reddan(:,i));
    
    % run mediation models for each of the mediators
    [paths, toplevelstats, ~] = mediation(X,Y{i},M{i},'names',{'group',contrastnames{i},'NPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPS{i} = paths;
    toplevelstats_NPS{i} = toplevelstats;
    clear paths toplevelstats;
    
%     [paths, toplevelstats, ~] = mediation(X,Y{i},Mpos{i},'names',{'group',contrastnames{i},'NPSpos response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_NPSpos{i} = paths;
%     toplevelstats_NPSpos{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},Mneg{i},'names',{'group',contrastnames{i},'NPSneg response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_NPSneg{i} = paths;
%     toplevelstats_NPSneg{i} = toplevelstats;
%     clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MFAPS{i},'names',{'group',contrastnames{i},'FAPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_FAPS{i} = paths;
    toplevelstats_FAPS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MVIFS{i},'names',{'group',contrastnames{i},'VIFS response'},'boottop','plots');
    drawnow;snapnow;
    paths_VIFS{i} = paths;
    toplevelstats_VIFS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},Mreddan{i},'names',{'group',contrastnames{i},'Reddan threat response'},'boottop','plots');
    drawnow;snapnow;
    paths_reddan{i} = paths;
    toplevelstats_reddan{i} = toplevelstats;
    clear paths toplevelstats;
    
end


%%
% <html><h3>WORRY</h3></html>

for i = 1:size(contrastnames,2) % loop over contrasts
    
    % define outcome for contrast i
    Y{i} = Ydat_worry(:,i);
    
    % define mediators
    M{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPS(:,i));
%     Mpos{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSpos(:,i));
%     Mneg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSneg(:,i));
    MFAPS{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.FAPS(:,i));
    MVIFS{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.VIFS(:,i));
    Mreddan{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.reddan(:,i));
    
    % run mediation models for each of the mediators
    [paths, toplevelstats, ~] = mediation(X,Y{i},M{i},'names',{'group',contrastnames{i},'NPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPS{i} = paths;
    toplevelstats_NPS{i} = toplevelstats;
    clear paths toplevelstats;
    
%     [paths, toplevelstats, ~] = mediation(X,Y{i},Mpos{i},'names',{'group',contrastnames{i},'NPSpos response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_NPSpos{i} = paths;
%     toplevelstats_NPSpos{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},Mneg{i},'names',{'group',contrastnames{i},'NPSneg response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_NPSneg{i} = paths;
%     toplevelstats_NPSneg{i} = toplevelstats;
%     clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MFAPS{i},'names',{'group',contrastnames{i},'FAPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_FAPS{i} = paths;
    toplevelstats_FAPS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MVIFS{i},'names',{'group',contrastnames{i},'VIFS response'},'boottop','plots');
    drawnow;snapnow;
    paths_VIFS{i} = paths;
    toplevelstats_VIFS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},Mreddan{i},'names',{'group',contrastnames{i},'Reddan threat response'},'boottop','plots');
    drawnow;snapnow;
    paths_reddan{i} = paths;
    toplevelstats_reddan{i} = toplevelstats;
    clear paths toplevelstats;
    
end


%%
% <html><h3>ANXIOUS</h3></html>

for i = 1:size(contrastnames,2) % loop over contrasts
    
    % define outcome for contrast i
    Y{i} = Ydat_anxious(:,i);
    
    % define mediators
    M{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPS(:,i));
%     Mpos{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSpos(:,i));
%     Mneg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSneg(:,i));
    MFAPS{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.FAPS(:,i));
    MVIFS{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.VIFS(:,i));
    Mreddan{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.reddan(:,i));
    
    % run mediation models for each of the mediators
    [paths, toplevelstats, ~] = mediation(X,Y{i},M{i},'names',{'group',contrastnames{i},'NPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPS{i} = paths;
    toplevelstats_NPS{i} = toplevelstats;
    clear paths toplevelstats;
    
%     [paths, toplevelstats, ~] = mediation(X,Y{i},Mpos{i},'names',{'group',contrastnames{i},'NPSpos response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_NPSpos{i} = paths;
%     toplevelstats_NPSpos{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},Mneg{i},'names',{'group',contrastnames{i},'NPSneg response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_NPSneg{i} = paths;
%     toplevelstats_NPSneg{i} = toplevelstats;
%     clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MFAPS{i},'names',{'group',contrastnames{i},'FAPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_FAPS{i} = paths;
    toplevelstats_FAPS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MVIFS{i},'names',{'group',contrastnames{i},'VIFS response'},'boottop','plots');
    drawnow;snapnow;
    paths_VIFS{i} = paths;
    toplevelstats_VIFS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},Mreddan{i},'names',{'group',contrastnames{i},'Reddan threat response'},'boottop','plots');
    drawnow;snapnow;
    paths_reddan{i} = paths;
    toplevelstats_reddan{i} = toplevelstats;
    clear paths toplevelstats;
    
end


%% SIMPLE REGRESSIONS OF SIGNATURE RESPONSES ON RATINGS
% ---------------------------------------------------------------------------

%%
% <html><h3>FAPS</h3></html>
for i = 1:size(contrastnames,2)
    
    printhdr(contrastnames{i});
    
    for j = 1:size(Ydat,2)
        
        printstr(Ynames{j});
        
        X_reg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.FAPS(:,i));
        Y{i,j} = Ydat{j}(:,i);
        [b_rob,stat_rob] = robustfit(X_reg{i},Y{i,j});
        betas_rob{i,j} = b_rob;
        stats_rob{i,j} = stat_rob;
        
        [b_ols,bint_ols,r_ols,rint_ols,stat_ols] = regress(Y{i,j},[ones(size(X_reg{i},1),1) X_reg{i}]);
        betas_ols{i,j} = b_ols;
        betaints_ols{i,j} = bint_ols;
        resids_ols{i,j} = r_ols;
        residints_ols{i,j} = rint_ols;
        stats_ols{i,j} = stat_ols;
        
        figure;
        title(contrastnames{i});
        scatter(X_reg{i},Y{i,j},'filled'); grid on; hold on
        plot(X_reg{i},b_ols(1)+b_ols(2)*X_reg{i},'r','LineWidth',2);
        plot(X_reg{i},b_rob(1)+b_rob(2)*X_reg{i},'g','LineWidth',2)
        legend('Data','Ordinary Least Squares','Robust Regression');
        xlabel(['FAPS response']);
        ylabel([Ynames{j} ' rating']);
        hold off;
        drawnow,snapnow;
        
        fprintf('\nROBUST BETA COEFFICIENTS\n%s\n\t\tbeta\t\tse\t\tt-value\t\tp-value\n',underscores);
        fprintf('intercept\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\n',betas_rob{i,j}(1,1),stats_rob{i,j}.se(1,1),stats_rob{i,j}.t(1,1),stats_rob{i,j}.p(1,1));
        fprintf('FAPS\t\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\n',betas_rob{i,j}(2,1),stats_rob{i,j}.se(2,1),stats_rob{i,j}.t(2,1),stats_rob{i,j}.p(2,1));
        fprintf('%s\n\n',underscores);
        
        clear b_rob stat_rob b_ols bint_ols r_ols rint_ols stat_ols;
        
    end
    
end 

%%
% <html><h3>VIFS</h3></html>
for i = 1:size(contrastnames,2)
    
    printhdr(contrastnames{i});
    
    for j = 1:size(Ydat,2)
        
        X_reg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.VIFS(:,i));
        Y{i,j} = Ydat{j}(:,i);
        [b_rob,stat_rob] = robustfit(X_reg{i},Y{i,j});
        betas_rob{i,j} = b_rob;
        stats_rob{i,j} = stat_rob;
        
        [b_ols,bint_ols,r_ols,rint_ols,stat_ols] = regress(Y{i,j},[ones(size(X_reg{i},1),1) X_reg{i}]);
        betas_ols{i,j} = b_ols;
        betaints_ols{i,j} = bint_ols;
        resids_ols{i,j} = r_ols;
        residints_ols{i,j} = rint_ols;
        stats_ols{i,j} = stat_ols;
        
        figure;
        scatter(X_reg{i},Y{i,j},'filled'); grid on; hold on
        plot(X_reg{i},b_ols(1)+b_ols(2)*X_reg{i},'r','LineWidth',2);
        plot(X_reg{i},b_rob(1)+b_rob(2)*X_reg{i},'g','LineWidth',2)
        legend('Data','Ordinary Least Squares','Robust Regression');
        xlabel(['VIFS response']);
        ylabel([Ynames{j} ' rating']);
        hold off;
        drawnow,snapnow;
        
        fprintf('\nROBUST BETA COEFFICIENTS\n%s\n\t\tbeta\t\tse\t\tt-value\t\tp-value\n',underscores);
        fprintf('intercept\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\n',betas_rob{i,j}(1,1),stats_rob{i,j}.se(1,1),stats_rob{i,j}.t(1,1),stats_rob{i,j}.p(1,1));
        fprintf('VIFS\t\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\n',betas_rob{i,j}(2,1),stats_rob{i,j}.se(2,1),stats_rob{i,j}.t(2,1),stats_rob{i,j}.p(2,1));
        fprintf('%s\n\n',underscores);
        
        clear b_rob stat_rob b_ols bint_ols r_ols rint_ols stat_ols;
        
    end
    
end 

%%
% <html><h3>REDDAN</h3></html>
for i = 1:size(contrastnames,2)
    
    printhdr(contrastnames{i});
    
    for j = 1:size(Ydat,2)
        
        X_reg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.reddan(:,i));
        Y{i,j} = Ydat{j}(:,i);
        [b_rob,stat_rob] = robustfit(X_reg{i},Y{i,j});
        betas_rob{i,j} = b_rob;
        stats_rob{i,j} = stat_rob;
        
        [b_ols,bint_ols,r_ols,rint_ols,stat_ols] = regress(Y{i,j},[ones(size(X_reg{i},1),1) X_reg{i}]);
        betas_ols{i,j} = b_ols;
        betaints_ols{i,j} = bint_ols;
        resids_ols{i,j} = r_ols;
        residints_ols{i,j} = rint_ols;
        stats_ols{i,j} = stat_ols;
        
        figure;
        scatter(X_reg{i},Y{i,j},'filled'); grid on; hold on
        plot(X_reg{i},b_ols(1)+b_ols(2)*X_reg{i},'r','LineWidth',2);
        plot(X_reg{i},b_rob(1)+b_rob(2)*X_reg{i},'g','LineWidth',2)
        legend('Data','Ordinary Least Squares','Robust Regression');
        xlabel(['Reddan threat response']);
        ylabel([Ynames{j} ' rating']);
        hold off;
        drawnow,snapnow;
        
        fprintf('\nROBUST BETA COEFFICIENTS\n%s\n\t\tbeta\t\tse\t\tt-value\t\tp-value\n',underscores);
        fprintf('intercept\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\n',betas_rob{i,j}(1,1),stats_rob{i,j}.se(1,1),stats_rob{i,j}.t(1,1),stats_rob{i,j}.p(1,1));
        fprintf('Reddan\t\t%5.4f\t\t%5.4f\t\t%5.4f\t\t%5.4f\n',betas_rob{i,j}(2,1),stats_rob{i,j}.se(2,1),stats_rob{i,j}.t(2,1),stats_rob{i,j}.p(2,1));
        fprintf('%s\n\n',underscores);
        
        clear b_rob stat_rob b_ols bint_ols r_ols rint_ols stat_ols;
        
    end
    
end 


% %% SELECTED NPS SUBREGIONS
% % -------------------------------------------------------------------------
% 
% for i = 1:size(contrastnames,2) 
%     
%     Y{i} = Ydat(:,i);
%     
%     MrIns{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx,2);
%     MrdpIns{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx,6);
%     MrS2_Op{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx,7);
%     MdACC{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx,8);
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MrIns{i},'names',{'group',contrastnames{i},'rIns response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_rIns{i} = paths;
%     toplevelstats_rIns{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MrdpIns{i},'names',{'group',contrastnames{i},'rdpIns response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_rdpIns{i} = paths;
%     toplevelstats_rdpIns{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MrS2_Op{i},'names',{'group',contrastnames{i},'rS2_Op response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_rS2_Op{i} = paths;
%     toplevelstats_rS2_Op{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MdACC{i},'names',{'group',contrastnames{i},'dACC response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_dACC{i} = paths;
%     toplevelstats_dACC{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MrIns{i},'M',[MrdpIns{i},MrS2_Op{i},MdACC{i}],'names',{'group',contrastnames{i},'rIns response','rdpIns response','rS2_Op response','dACC response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_mult_med{i} = paths;
%     toplevelstats_mult_med{i} = toplevelstats;
%     clear paths toplevelstats;
%     
% end


% %% VOXEL-BASED MEDIATION
% %--------------------------------------------------------------------------
% 
% for i = 1:size(contrastnames,2)
%     
%     brain_mediationdir_contrast = fullfile(brain_mediationdir, contrastnames{i});
%     if ~isfolder(brain_mediationdir_contrast)
%         mkdir(brain_mediationdir_contrast);
%     end
%     cd(brain_mediationdir_contrast);
%     
%     Y{i} = Ydat_pain(:,i);
%     
%     M_brain{i} = DAT.imgs{i}(:);
%     
%     mediation_brain(X,Y{i},char(M_brain{i}),'mask',maskname_glm,'names',{'group',contrastnames{i},'brain'});
%     
%     cd(brain_mediationdir);
%     
% end
% 
% % now run the following script from each results directory
% 
% publish_mediation_report;


%% MULTIVARIATE MEDIATION
%--------------------------------------------------------------------------
% based on
% https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/Multivariate_Mediation_ExampleScript.m
% and personal communication with Martin Lindquist and Xiaochun Han
% kudos to Bogdan Petre for help with source reconstruction code

X_c = num2cell(X(idx)); % convert predictor to cell array and exclude subject without behavioral data

for i = 1:size(contrastnames,2)
    
    pdm_mediationdir_contrast = fullfile(pdm_mediationdir, contrastnames{i});
        if ~isfolder(pdm_mediationdir_contrast)
            mkdir(pdm_mediationdir_contrast);
        end
    cd(pdm_mediationdir_contrast);
    
    for j = 1:size(Ydat,2)
        pdm_mediationdir_contrast_behav = fullfile(pdm_mediationdir_contrast, Ynames{j});
            if ~isfolder(pdm_mediationdir_contrast_behav)
                mkdir(pdm_mediationdir_contrast_behav);
            end
        cd(pdm_mediationdir_contrast_behav);
    
        if ~isfile(strcat('PDMresults_',contrastnames{i},'_',Ynames{j},'.mat'))

            Y{i,j} = Ydat{j}(idx,i);
            Y_c{i,j} = num2cell(Y{i,j});

            M_brain{i} = DAT.imgs{i}(idx);
            names = M_brain{i};

            mask = which('gray_matter_mask.nii');

            for k=1:size(X_c,1)
                dat = fmri_data(names{k},mask); 
                m{k} = dat.dat; 
            end

            save(strcat('data_objects_',contrastnames{i},'_',Ynames{j},'.mat'),'dat','-v7.3');

            pdm = multivariateMediation(X_c,Y_c{i,j},m,'B',20,'svd','plots'); % run up to here first to decide on number of pdms to retain based on plots
            pdmfull = pdm;
            pdm = multivariateMediation(pdm,'nPDM',2);
            pdm = multivariateMediation(pdm,'noPDMestimation','bootPDM',1:2,'bootjPDM','Bsamp',5000,'save2file',strcat('PDMresults_',contrastnames{i},'_',Ynames{j},'.mat'));
            save(strcat('PDMresults_',contrastnames{i},'_',Ynames{j},'.mat'),'pdmfull','-append');

        else

            load(strcat('PDMresults_',contrastnames{i},'_',Ynames{j},'.mat'));
            pdm = out;
            clear out;
            load(strcat('data_objects_',contrastnames{i},'_',Ynames{j},'.mat'));

        end % if loop
        
        % source reconstruction
        % "source reconstruction is thus no more than a covariance map
        % which shows how much each voxel covaries with model predictions
        % across images" - Bogdan Petre
        for n = 1:size(pdm.Wfull,2)
            
            X_source{i,j} = fmri_data(DAT.imgs{i}(idx),which('gray_matter_mask.nii'));
            Xz_source{i,j} = rescale(X_source{i,j},'centervoxels');
            M_source{i,j,n} = pdm.Wfull{1,n};
            P_source{i,j,n} = M_source{i,j,n}'*Xz_source{i,j}.dat;
            source{i,j,n} = (P_source{i,j,n}*Xz_source{i,j}.dat');
            source{i,j,n} = source{i,j,n}'./(size(P_source{i,j,n},2)-1);
            source_obj{i,j,n} = dat;
            source_obj{i,j,n}.dat = source{i,j,n};
            source_obj{i,j,n}.image_names = 'source reconstruction';
            source_obj{i,j,n}.fullpath = '';
            source_obj{i,j,n}.history = {['source reconstructed from con images and pdm_',num2str(n)]};
            source_obj_j{n} = source_obj{i,j,n};
            
        end % source recon loop
        
        save(strcat('PDM_source_recon_',contrastnames{i},'_',Ynames{j},'.mat'),'source_obj_j');
    
    end % for loop outcomes
        
    cd(pdm_mediationdir);
    
end % for loop contrasts

%% now run the following script from each results directory

publish_multivariate_mediation_report;
