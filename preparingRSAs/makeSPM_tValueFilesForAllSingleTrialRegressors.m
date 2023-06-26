%% create t-images for trialwise betas

clear all
clear matlabbatch;

%% initialize spm batch
% spm needs to be in path
% addpath 'PATH_TO_SPM' 
spm('defaults','fmri');
spm_jobman('initcfg');

%% initialize paths and sjs
subjects =  [1:4]; % subject numbers
data_dir = '.\data\mriData\'; % path to subject data

%% run in spm batch
for sub = subjects
    
    %try
    clear matlabbatch;

    VPFolder = dir(fullfile(data_dir, strcat('sj',sprintf('%03d', sub))));
    
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(data_dir, VPFolder.name, '\trialwiseGLM\SPM.mat')};
    
    for tr = 1:420 %number of trials
        matlabbatch{1}.spm.stats.con.consess{tr}.tcon.name = ['trial', num2str(tr)];
        matlabbatch{1}.spm.stats.con.consess{tr}.tcon.weights = [zeros(1,tr-1),1];
        matlabbatch{1}.spm.stats.con.consess{tr}.tcon.sessrep = 'none';
    end
    matlabbatch{1}.spm.stats.con.delete = 1;
    
    spm_jobman('run', matlabbatch)
end


