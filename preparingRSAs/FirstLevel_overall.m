%%  Fist level examplew without MCF files anf or overall condition comparison (not single trial)

clear all; close all; clc;

batch_counter = 1;
clear batch

su =  [1:99] ;

for sub =1:length(su)
    clear matlabbatch
    
    subj = su(sub);
    folderName=sprintf('%02.0f',subj);
    
    Func1 = {};
    Functionals1 ={};
    Func2 = {};
    Functionals2 ={};
    mkdir(['M:\YOURPATH\DATA\Subj_',folderName,'\First_Level_normalized_allDays_unsmoothed\']);
    
    step = 1;
    
    matlabbatch{step}.spm.stats.fmri_spec.dir = {['M:\YOURPATH\DATA\Subj_',folderName,'\First_Level_normalized_allDays_unsmoothed\']};
    matlabbatch{step}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{step}.spm.stats.fmri_spec.timing.RT = 2;
    matlabbatch{step}.spm.stats.fmri_spec.timing.fmri_t = 62;
    matlabbatch{step}.spm.stats.fmri_spec.timing.fmri_t0 = 15;
    
    
    a = dir(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\Bold\Block_1\','wauf*']);
    b = dir(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\Bold\Block_2\','wauf*']);
    
    FuncRetD1={};
    FunctionalsRetD1B1 ={};
    FunctionalsRetD1B2 ={};
    
    for i = 1:length(a)
        FunctionalsRetD1B1{i,1} = [ a(i).folder,'\',a(i).name];
    end
    
    for i = 1:length(b)
        FunctionalsRetD1B2{i,1} = [ b(i).folder,'\',b(i).name];
    end
    
    FuncRetD1 =[(FunctionalsRetD1B1);(FunctionalsRetD1B2)];
    
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).scans = cellstr(FuncRetD1);

    
    b = dir(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\','movement_nuisance*.txt']);
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).multi_reg(1) = {[b(1).folder,'\',b(1).name]};
   
    
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_reac_rememb1.mat']);
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_no_reac_rememb1.mat']);
    
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_reac_rememb2.mat']);
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_no_reac_rememb2.mat']);
    
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_reac_forgottenb1.mat'])
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_no_reac_forgottenb1.mat']);
    
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_reac_forgottenb2.mat'])
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_no_reac_forgottenb2.mat']);
    
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_CR_b1.mat']);
    load(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\first_Retrieval\onsets_CR_b2.mat']);
    
    first_block_length =  length(FunctionalsRetD1B1)*2;
    
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(1).name = 'Reac_rememD1';
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(1).onset = [onsets_reac_rememb1,onsets_reac_rememb2+first_block_length];
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(1).duration =2;
    %
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(2).name = 'Reac_forgottenD1';
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(2).onset = [onsets_reac_forgottenb1,onsets_reac_forgottenb2+first_block_length];
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(2).duration =2;
    
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(3).name = 'No_Reac_RememD1';
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(3).onset =  [onsets_no_reac_rememb1,onsets_no_reac_rememb2+first_block_length];
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(3).duration =2;
    %
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(4).name = 'No_Reac_ForgottenD1';
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(4).onset =  [onsets_no_reac_forgottenb1,onsets_no_reac_forgottenb2+first_block_length];
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(4).duration =2;
    
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(5).name = 'CRD1';
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(5).onset =  [onsets_CR_b1,onsets_CR_b2+first_block_length];
    matlabbatch{step}.spm.stats.fmri_spec.sess(1).cond(5).duration =2;
    
    %%%%%
    matlabbatch{step}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{step}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{step}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{step}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{step}.spm.stats.fmri_spec.mthresh = 0.6;
    matlabbatch{step}.spm.stats.fmri_spec.cvi = 'AR(1)';
    step = step+1;
    
    %%Model Estiamtion
    matlabbatch{step}.spm.stats.fmri_est.spmmat(1) = cellstr(['I:\STRESSNET\MRI\MRI_DATA\StressNet_MRI_converted_nii\Subj_',folderName,'\First_Level_normalized_allDays_unsmoothed\SPM.mat']);
    matlabbatch{step}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{step}.spm.stats.fmri_est.method.Classical = 1;
    
    
    batch{batch_counter} = matlabbatch;
    
    batch_counter = batch_counter+1;
end

poolobj = gcp('nocreate');
delete(poolobj);
parfor i = 1:length(batch)
    spm_jobman('run', batch{i})
    disp(['Finished jo for subj.', i]);
end
%%



