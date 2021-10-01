% This script calculates the signature response for the FAPS (fear in the
% anticipation of pain signature, Kragel, Van Oudenhove, et al, in
% preparation) to Iris' WAD dataset, but can easily be
% adapted/functionalized for any dataset in the future
% it is adapted from prep_4_apply_signatures_and_save and should be
% incorporated into that script once the FAPS pattern becomes publicly
% available on CANlab Github


%% LOAD FAPS PATTERN INTO FMRI_DATA OBJECT, AND APPLY TO CONDITIONS AND CONTRASTS IN DAT
%---------------------------------------------------------------------------------------

FAPS_weights = fmri_data('C:\Users\lukas\Dropbox (Personal)\V_visceral_common\projects\anticipation_visceral\scripts\FAPS.nii');

% CONDITIONS
%-----------
k = size(DAT.conditions,2);

for i = 1:k
        
    DAT.SIG_conditions.raw.dotproduct.FAPS(:,i) = apply_mask(DATA_OBJ{i},FAPS_weights,'pattern_expression','ignore_missing');
    DAT.SIG_conditions.raw.cosine_sim.FAPS(:,i) = apply_mask(DATA_OBJ{i},FAPS_weights,'pattern_expression','ignore_missing','cosine_similarity');
    
    DAT.SIG_conditions.scaled.dotproduct.FAPS(:,i) = apply_mask(DATA_OBJsc{i},FAPS_weights,'pattern_expression','ignore_missing');
    DAT.SIG_conditions.scaled.cosine_sim.FAPS(:,i) = apply_mask(DATA_OBJsc{i},FAPS_weights,'pattern_expression','ignore_missing','cosine_similarity');

end

DAT.SIG_conditions.raw.dotproduct.FAPS = array2table(DAT.SIG_conditions.raw.dotproduct.FAPS,'VariableNames',DAT.SIG_conditions.raw.dotproduct.conditionnames);
DAT.SIG_conditions.raw.cosine_sim.FAPS = array2table(DAT.SIG_conditions.raw.cosine_sim.FAPS,'VariableNames',DAT.SIG_conditions.raw.cosine_sim.conditionnames);

DAT.SIG_conditions.scaled.dotproduct.FAPS = array2table(DAT.SIG_conditions.scaled.dotproduct.FAPS,'VariableNames',DAT.SIG_conditions.scaled.dotproduct.conditionnames);
DAT.SIG_conditions.scaled.cosine_sim.FAPS = array2table(DAT.SIG_conditions.scaled.cosine_sim.FAPS,'VariableNames',DAT.SIG_conditions.scaled.cosine_sim.conditionnames);

% CONTRASTS
%----------

kc = size(DAT.contrasts,1);

for j = 1:kc
    
    DAT.SIG_contrasts.raw.dotproduct.FAPS(:,j) = apply_mask(DATA_OBJ_CON{j},FAPS_weights,'pattern_expression','ignore_missing');
    DAT.SIG_contrasts.raw.cosine_sim.FAPS(:,j) = apply_mask(DATA_OBJ_CON{j},FAPS_weights,'pattern_expression','ignore_missing','cosine_similarity');
    
    DAT.SIG_contrasts.scaled.dotproduct.FAPS(:,j) = apply_mask(DATA_OBJ_CONsc{j},FAPS_weights,'pattern_expression','ignore_missing');
    DAT.SIG_contrasts.scaled.cosine_sim.FAPS(:,j) = apply_mask(DATA_OBJ_CONsc{j},FAPS_weights,'pattern_expression','ignore_missing','cosine_similarity');
    
end

DAT.SIG_contrasts.raw.dotproduct.FAPS = array2table(DAT.SIG_contrasts.raw.dotproduct.FAPS,'VariableNames',DAT.SIG_contrasts.raw.dotproduct.conditionnames);
DAT.SIG_contrasts.raw.cosine_sim.FAPS = array2table(DAT.SIG_contrasts.raw.cosine_sim.FAPS,'VariableNames',DAT.SIG_contrasts.raw.cosine_sim.conditionnames);

DAT.SIG_contrasts.scaled.dotproduct.FAPS = array2table(DAT.SIG_contrasts.scaled.dotproduct.FAPS,'VariableNames',DAT.SIG_contrasts.scaled.dotproduct.conditionnames);
DAT.SIG_contrasts.scaled.cosine_sim.FAPS = array2table(DAT.SIG_contrasts.scaled.cosine_sim.FAPS,'VariableNames',DAT.SIG_contrasts.scaled.cosine_sim.conditionnames);


%% LOAD VIFS PATTERN INTO FMRI_DATA OBJECT, AND APPLY TO CONDITIONS AND CONTRASTS IN DAT
%---------------------------------------------------------------------------------------

VIFS_weights = fmri_data('C:\Users\lukas\Dropbox (Personal)\V_visceral_common\projects\anticipation_visceral\scripts\VIFS.nii');

% CONDITIONS
%-----------
k = size(DAT.conditions,2);

for i = 1:k
        
    DAT.SIG_conditions.raw.dotproduct.VIFS(:,i) = apply_mask(DATA_OBJ{i},VIFS_weights,'pattern_expression','ignore_missing');
    DAT.SIG_conditions.raw.cosine_sim.VIFS(:,i) = apply_mask(DATA_OBJ{i},VIFS_weights,'pattern_expression','ignore_missing','cosine_similarity');
    
    DAT.SIG_conditions.scaled.dotproduct.VIFS(:,i) = apply_mask(DATA_OBJsc{i},VIFS_weights,'pattern_expression','ignore_missing');
    DAT.SIG_conditions.scaled.cosine_sim.VIFS(:,i) = apply_mask(DATA_OBJsc{i},VIFS_weights,'pattern_expression','ignore_missing','cosine_similarity');

end

DAT.SIG_conditions.raw.dotproduct.VIFS = array2table(DAT.SIG_conditions.raw.dotproduct.VIFS,'VariableNames',DAT.SIG_conditions.raw.dotproduct.conditionnames);
DAT.SIG_conditions.raw.cosine_sim.VIFS = array2table(DAT.SIG_conditions.raw.cosine_sim.VIFS,'VariableNames',DAT.SIG_conditions.raw.cosine_sim.conditionnames);

DAT.SIG_conditions.scaled.dotproduct.VIFS = array2table(DAT.SIG_conditions.scaled.dotproduct.VIFS,'VariableNames',DAT.SIG_conditions.scaled.dotproduct.conditionnames);
DAT.SIG_conditions.scaled.cosine_sim.VIFS = array2table(DAT.SIG_conditions.scaled.cosine_sim.VIFS,'VariableNames',DAT.SIG_conditions.scaled.cosine_sim.conditionnames);

% CONTRASTS
%----------

kc = size(DAT.contrasts,1);

for j = 1:kc
    
    DAT.SIG_contrasts.raw.dotproduct.VIFS(:,j) = apply_mask(DATA_OBJ_CON{j},VIFS_weights,'pattern_expression','ignore_missing');
    DAT.SIG_contrasts.raw.cosine_sim.VIFS(:,j) = apply_mask(DATA_OBJ_CON{j},VIFS_weights,'pattern_expression','ignore_missing','cosine_similarity');
    
    DAT.SIG_contrasts.scaled.dotproduct.VIFS(:,j) = apply_mask(DATA_OBJ_CONsc{j},VIFS_weights,'pattern_expression','ignore_missing');
    DAT.SIG_contrasts.scaled.cosine_sim.VIFS(:,j) = apply_mask(DATA_OBJ_CONsc{j},VIFS_weights,'pattern_expression','ignore_missing','cosine_similarity');
    
end

DAT.SIG_contrasts.raw.dotproduct.VIFS = array2table(DAT.SIG_contrasts.raw.dotproduct.VIFS,'VariableNames',DAT.SIG_contrasts.raw.dotproduct.conditionnames);
DAT.SIG_contrasts.raw.cosine_sim.VIFS = array2table(DAT.SIG_contrasts.raw.cosine_sim.VIFS,'VariableNames',DAT.SIG_contrasts.raw.cosine_sim.conditionnames);

DAT.SIG_contrasts.scaled.dotproduct.VIFS = array2table(DAT.SIG_contrasts.scaled.dotproduct.VIFS,'VariableNames',DAT.SIG_contrasts.scaled.dotproduct.conditionnames);
DAT.SIG_contrasts.scaled.cosine_sim.VIFS = array2table(DAT.SIG_contrasts.scaled.cosine_sim.VIFS,'VariableNames',DAT.SIG_contrasts.scaled.cosine_sim.conditionnames);


%% LOAD REDDAN THREAT PATTERN INTO FMRI_DATA OBJECT, AND APPLY TO CONDITIONS AND CONTRASTS IN DAT
%---------------------------------------------------------------------------------------

reddan_weights = fmri_data(which('IE_ImEx_Acq_Threat_SVM_nothresh.nii.gz'));

% CONDITIONS
%-----------
k = size(DAT.conditions,2);

for i = 1:k
        
    DAT.SIG_conditions.raw.dotproduct.reddan(:,i) = apply_mask(DATA_OBJ{i},reddan_weights,'pattern_expression','ignore_missing');
    DAT.SIG_conditions.raw.cosine_sim.reddan(:,i) = apply_mask(DATA_OBJ{i},reddan_weights,'pattern_expression','ignore_missing','cosine_similarity');
    
    DAT.SIG_conditions.scaled.dotproduct.reddan(:,i) = apply_mask(DATA_OBJsc{i},reddan_weights,'pattern_expression','ignore_missing');
    DAT.SIG_conditions.scaled.cosine_sim.reddan(:,i) = apply_mask(DATA_OBJsc{i},reddan_weights,'pattern_expression','ignore_missing','cosine_similarity');

end

DAT.SIG_conditions.raw.dotproduct.reddan = array2table(DAT.SIG_conditions.raw.dotproduct.reddan,'VariableNames',DAT.SIG_conditions.raw.dotproduct.conditionnames);
DAT.SIG_conditions.raw.cosine_sim.reddan = array2table(DAT.SIG_conditions.raw.cosine_sim.reddan,'VariableNames',DAT.SIG_conditions.raw.cosine_sim.conditionnames);

DAT.SIG_conditions.scaled.dotproduct.reddan = array2table(DAT.SIG_conditions.scaled.dotproduct.reddan,'VariableNames',DAT.SIG_conditions.scaled.dotproduct.conditionnames);
DAT.SIG_conditions.scaled.cosine_sim.reddan = array2table(DAT.SIG_conditions.scaled.cosine_sim.reddan,'VariableNames',DAT.SIG_conditions.scaled.cosine_sim.conditionnames);

% CONTRASTS
%----------

kc = size(DAT.contrasts,1);

for j = 1:kc
    
    DAT.SIG_contrasts.raw.dotproduct.reddan(:,j) = apply_mask(DATA_OBJ_CON{j},reddan_weights,'pattern_expression','ignore_missing');
    DAT.SIG_contrasts.raw.cosine_sim.reddan(:,j) = apply_mask(DATA_OBJ_CON{j},reddan_weights,'pattern_expression','ignore_missing','cosine_similarity');
    
    DAT.SIG_contrasts.scaled.dotproduct.reddan(:,j) = apply_mask(DATA_OBJ_CONsc{j},reddan_weights,'pattern_expression','ignore_missing');
    DAT.SIG_contrasts.scaled.cosine_sim.reddan(:,j) = apply_mask(DATA_OBJ_CONsc{j},reddan_weights,'pattern_expression','ignore_missing','cosine_similarity');
    
end

DAT.SIG_contrasts.raw.dotproduct.reddan = array2table(DAT.SIG_contrasts.raw.dotproduct.reddan,'VariableNames',DAT.SIG_contrasts.raw.dotproduct.conditionnames);
DAT.SIG_contrasts.raw.cosine_sim.reddan = array2table(DAT.SIG_contrasts.raw.cosine_sim.reddan,'VariableNames',DAT.SIG_contrasts.raw.cosine_sim.conditionnames);

DAT.SIG_contrasts.scaled.dotproduct.reddan = array2table(DAT.SIG_contrasts.scaled.dotproduct.reddan,'VariableNames',DAT.SIG_contrasts.scaled.dotproduct.conditionnames);
DAT.SIG_contrasts.scaled.cosine_sim.reddan = array2table(DAT.SIG_contrasts.scaled.cosine_sim.reddan,'VariableNames',DAT.SIG_contrasts.scaled.cosine_sim.conditionnames);


%% SAVE
% -------------------------------------------------------------------------

printhdr('Save results');

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, 'DAT', '-append');