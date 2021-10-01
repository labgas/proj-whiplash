%% WAD_second_m1_copy_con_imgs_CANlab_folder_struct
%
% This simple script copies con1-3 images (high_fear, moderate_fear, neutral_fear) 
% from rootdir\firstlevel\model_x_yyy\sub-zz\ to
% rootdir\secondlevel\model_x_yyy\con_images\sub-zz
% to prepare second level analysis with CANlab batch scripts
%
%__________________________________________________________________________
%
% authors: Lukas Van Oudenhove & Iris Coppieters
% date:   April, 2021
%__________________________________________________________________________
% @(#)% WAD_second_m1_copy_con_imgs_CANlab_folder_struct     v1.0        
% last modified: 2021/04/13


%% DEFINE MODEL AND DIRECTORIES
%--------------------------------------------------------------------------
conimgs2copy=[1:3]; % numbers of con images to copy
model='model_1_pictures';
rootdir='C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI';
firstleveldir=fullfile(rootdir,'firstlevel\',model);
secondleveldir=fullfile(rootdir,'secondlevel\',model);
conimgsdir=fullfile(secondleveldir,'\data');
if ~exist(conimgsdir,'dir')
    mkdir(conimgsdir);
end


%% GET SUBJECT DIRECTORY NAMES FROM FIRSTLEVELDIR
%--------------------------------------------------------------------------
subjs = dir([firstleveldir,'\sub-*']);
subjs = {subjs(:).name}';


%% WRITE SUBJECT DIRECTORIES IN CONIMGSDIR
%--------------------------------------------------------------------------
cd(conimgsdir);
sm=@(x)spm_mkdir(x); % defines spm_mkdir as an anonymous function sm
cellfun(sm,subjs); % applies function sm to all cells of subjs


%% LOOP OVER SUBJECTS TO COPY THE SPECIFIED CON IMAGES OVER
%--------------------------------------------------------------------------
for i=1:size(subjs,1)
    subjfirstleveldir=fullfile(firstleveldir,subjs{i});
    subjconimgsdir=fullfile(conimgsdir,subjs{i});
    conimgs=ls(fullfile(subjfirstleveldir,'con*.nii'));
    conimgs=conimgs(conimgs2copy,:);
        for j=1:size(conimgs,1)
            copyfile(fullfile(subjfirstleveldir,conimgs(j,:)),fullfile(subjconimgsdir,conimgs(j,:)));
        end
end