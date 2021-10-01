%% DEFINE PATHS AND VARIABLES TO BE USED IN ALL ANALYSES
% -------------------------------------------------------------------------
addpath(genpath('C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI\scripts_windows')); % add WAD scripts dir to your path
addpath(genpath('C:\Users\lukas\Documents\GitHub\RainCloudPlots')); % add RainCloudPlots Github repo to your path
addpath(genpath('C:\Users\lukas\Documents\GitHub\Robust_Statistical_Toolbox')); % add Robust Statistical Toolbox Github repo to your path
addpath(genpath('C:\Users\lukas\Documents\MATLAB\cbrewer')); % add colorbrewer to your path for more color options

a_WAD_second_m1_s1_set_up_paths_always_run_first

figspubdir = fullfile(resultsdir,'figures_publication');
if ~isfolder(figspubdir)
    mkdir(figspubdir);
end

load(fullfile(resultsdir,'image_names_and_setup.mat'));

try
    % get nice colours from colorbrewer
    % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
    [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
catch
    % if you don't have colorbrewer, accept these far more boring colours
    cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
end

cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);

fig_position = [200 200 600 400]; % coordinates for figures


%% BEHAVIORAL DATA
%--------------------------------------------------------------------------
% CONDITIONS
idx = logical(DAT.BEHAVIOR.behavioral_data_table.included);
behdat = DAT.BEHAVIOR.behavioral_data_table(idx,:);
behdat = sortrows(behdat,'group','ascend');
behdat.group(behdat.group == -1) = 2;
group = [behdat.group; behdat.group; behdat.group];
condition = [ones(height(behdat),1); 2.*ones(height(behdat),1); 3.*ones(height(behdat),1)];
pain = [behdat.Pain_highfear; behdat.Pain_moderatefear; behdat.Pain_control];
worry = [behdat.Worry_highfear; behdat.Worry_moderatefear; behdat.Worry_control];
anxiety = [behdat.anxious_highfear; behdat.anxious_moderatefear; behdat.anxious_control];
avoidance = [behdat.avoidance_highfear; behdat.avoidance_moderatefear; behdat.avoidance_control];
outcomes = {pain,worry,anxiety,avoidance};
outcome_names = {'pain','worry','anxiety','avoidance'};

for o = 1:size(outcomes,2)
    
    D{o} = [outcomes{o},condition,group];
    for i = 1:3
        for j = 1:2
        data{i,j} = D{o}(D{o}(:, 2) == i & D{o}(:, 3) ==j);
        end
    end
    
    data_all{o} = data;
    
end

for o = 1:size(outcomes,2)
    
    f  = figure('Position', fig_position,'WindowState','maximized');
    h   = rm_raincloud(data_all{o}, cl, 0, 'rash');
    ax{o} = gca;
    ax{o}.FontSize = 14;
    ax{o}.FontName = 'Cambria';
    ax{o}.XAxis.LineWidth = 1;
    ax{o}.YAxis.LineWidth = 1;
    xlim([0 10]);
    xlabel({(strcat(outcome_names{o},' rating')),''},'FontSize',24,'FontWeight','bold');
    ylabel({'','condition'},'FontSize',24,'FontWeight','bold');
    yticklabels({'\fontsize{20} \bf control','\fontsize{20} \bf moderate fear','\fontsize{20} \bf high fear'});
    legend([h.l(1,1) h.l(1,2)],{'whiplash patients','healthy controls'},'Location','northeast','FontSize',24,'FontWeight','bold','Box','off');

    for i = 1:3
        for j = 1:2
            h.s{i, j}.SizeData = 150;
        end
    end
    
    % save
    print(f,fullfile(figspubdir,strcat('behav_',outcome_names{o},'.png')),'-dpng','-r600');
    
    f_all{o}= f;
    h_all{o} = h;
    
    clear f h;
    
end

% CONTRASTS - CLASSIC RAINCLOUD PLOT

for k = 1:max(unique(behdat.group))
    
    pain_contrast{k} = [behdat.pain_high_neu(behdat.group==k),behdat.pain_mod_neu(behdat.group==k)];
    worry_contrast{k} = [behdat.worry_high_neu(behdat.group==k),behdat.worry_mod_neu(behdat.group==k)];
    anxiety_contrast{k} = [behdat.anxious_high_neu(behdat.group==k),behdat.anxious_mod_neu(behdat.group==k)];
    avoidance_contrast{k} = [behdat.avoidance_high_neu(behdat.group==k),behdat.avoidance_mod_neu(behdat.group==k)];
    
end

outcomes_contrast = {pain_contrast,worry_contrast,anxiety_contrast,avoidance_contrast};
contrast_names = {'high fear versus neutral','moderate fear versus neutral'};

for o = 1:size(outcomes,2)
    
    for m = 1:size(pain_contrast,2)
        
        f2 = figure('Position', fig_position,'WindowState','maximized');
        h1 = raincloud_plot(outcomes_contrast{o}{1}(:,m), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
             'box_col_match', 1, 'line_width', 3);
        h2 = raincloud_plot(outcomes_contrast{o}{2}(:,m), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1, 'line_width', 3);
        h1{2}.SizeData = 50;
        h2{2}.SizeData = 50;
        h1{1}.EdgeColor = 'none';
        h2{1}.EdgeColor = 'none';
        h1{2}.SizeData = 50;
        h2{2}.SizeData = 50;
        ax3{o,m} = gca;
        ax3{o,m}.FontSize = 14;
        ax3{o,m}.FontName = 'Cambria';
        ax3{o,m}.YAxisLocation = 'origin';
        ax3{o,m}.YTick = [];
        ax3{o,m}.LineWidth = 0.25;
        xlabel({'',[outcome_names{o},' rating']},'FontSize',24,'FontWeight','bold');
        legend([h1{1} h2{1}], {'whiplash patients', 'healthy controls'},'Location','best','FontSize',24,'FontWeight','bold','Box','off');
        title(contrast_names{m},'FontSize',28,'FontWeight','bold');
        ylim([(h2{3}.Position(2)+0.10.*h2{3}.Position(2)) (max([h1{1}.YData h2{1}.YData])+0.05.*max([h1{1}.YData h2{1}.YData]))]);
        box off
    
    
        % save
        print(f2,fullfile(figspubdir,strcat('behav_contrasts_',outcome_names{o},'_',contrast_names{m},'.png')),'-dpng','-r600');

        f2_all{o,m}= f2;
        h1_all{o,m} = h1;
        h2_all{o,m} = h2;

        clear f2 h1 h2;
        
    end
    
end


%% SIGNATURE RESPONSES
%--------------------------------------------------------------------------
sigdat = DAT.SIG_contrasts.raw.dotproduct;
behdat_full = DAT.BEHAVIOR.behavioral_data_table;
behdat_full.group(behdat_full.group == -1) = 2;

for n = 1:max(unique(behdat_full.group))
    
    NPS{n} = [sigdat.NPS.high_fear_v_neutral(behdat_full.group==n),sigdat.NPS.moderate_fear_v_neutral(behdat_full.group==n)];
    NPSpos{n} = [sigdat.NPSpos.high_fear_v_neutral(behdat_full.group==n),sigdat.NPSpos.moderate_fear_v_neutral(behdat_full.group==n)];
    NPSneg{n} = [sigdat.NPSneg.high_fear_v_neutral(behdat_full.group==n),sigdat.NPSneg.moderate_fear_v_neutral(behdat_full.group==n)];
    PINES{n} = [sigdat.PINES.high_fear_v_neutral(behdat_full.group==n),sigdat.PINES.moderate_fear_v_neutral(behdat_full.group==n)];
    FAPS{n} = [sigdat.FAPS.high_fear_v_neutral(behdat_full.group==n),sigdat.FAPS.moderate_fear_v_neutral(behdat_full.group==n)];
    VIFS{n} = [sigdat.VIFS.high_fear_v_neutral(behdat_full.group==n),sigdat.VIFS.moderate_fear_v_neutral(behdat_full.group==n)];
    reddan{n} = [sigdat.reddan.high_fear_v_neutral(behdat_full.group==n),sigdat.reddan.moderate_fear_v_neutral(behdat_full.group==n)];
    
end

signatures = {NPS, NPSpos, NPSneg, PINES, FAPS, VIFS, reddan};
signature_names = {'NPS', 'NPS positive', 'NPS negative', 'PINES', 'FAPS', 'VIFS', 'Reddan threat'};
contrast_names = {'high fear versus neutral','moderate fear versus neutral'};

for s = 1:size(signatures,2)
    
    for p = 1:size(NPS,2)
        
        f3 = figure('Position', fig_position,'WindowState','maximized');
        h3 = raincloud_plot(signatures{s}{1}(:,p), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
             'box_col_match', 1, 'line_width', 3);
        h4 = raincloud_plot(signatures{s}{2}(:,p), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1, 'line_width', 3);
        h3{2}.SizeData = 50;
        h4{2}.SizeData = 50;
        h3{1}.EdgeColor = 'none';
        h4{1}.EdgeColor = 'none';
        h3{2}.SizeData = 50;
        h4{2}.SizeData = 50;
        ax4{s,p} = gca;
        ax4{s,p}.FontSize = 14;
        ax4{s,p}.FontName = 'Cambria';
        ax4{s,p}.YAxisLocation = 'origin';
        ax4{s,p}.YTick = [];
        ax4{s,p}.LineWidth = 0.25;
        xlabel({'',[signature_names{s},' response']},'FontSize',24,'FontWeight','bold');
        legend([h3{1} h4{1}], {'whiplash patients', 'healthy controls'},'Location','northeast','FontSize',24,'FontWeight','bold','Box','off');
        title(contrast_names{p},'FontSize',28,'FontWeight','bold');
        ylim([(h4{3}.Position(2)+0.10.*h4{3}.Position(2)) (max([h3{1}.YData h4{1}.YData])+0.05.*max([h3{1}.YData h4{1}.YData]))]);
        box off
    
    
        % save
        print(f3,fullfile(figspubdir,strcat(signature_names{s},'_',contrast_names{p},'.png')),'-dpng','-r600');

        f3_all{s,p}= f3;
        h3_all{s,p} = h3;
        h4_all{s,p} = h4;

        clear f3 h3 h4;
        
    end
    
end