%% COMPUTE FIRST LEVEL WITH EACH TRIAL AS ONE REGRESSOR

clear all; close all; clc;

su = [1:99] ;

batch_counter = 1; %counter especially for parallelization
clear batch

for sub =1:length(su)
    
    subj = su(sub);
    folderName=sprintf('%02.0f',subj);
    trl = 1;
    
    clear matlabbatch
    
    mkdir(['M:\YOURPATH\DATA\Subj_',folderName,'\firstlevel_native']) %make firstlevel dir define firstlevel
    matlabbatch{1}.spm.stats.fmri_spec.dir = {['M:\YOURPATH\DATA\Subj_',folderName,'\firstlevel_native']}; %define firstlevel output dir
    
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 62;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 15;
    
    %for each encoding run, a Multiple Conditions file (MCF) and the
    %corresponding movement regressors are put in
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(spm_select('FPList',['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\Encoding\Bold\Block_1'],'^rarf2*'));
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\Encoding\MCF_EncodeB1_singleTrial_subj_',folderName,'.mat']};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    x = dir(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\Encoding\Bold\Block_1\','rp_f*']); 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {[x(1).folder,'\',x(1).name]};%movement regressors
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = cellstr(spm_select('FPList',['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\Encoding\Bold\Block_2'],'^rarf2*'));
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\Encoding\MCF_EncodeB2_singleTrial_subj_',folderName,'.mat']};
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
    x = dir(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\Encoding\Bold\Block_2\','rp_f*']);%movement regressors
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {[x(1).folder,'\',x(1).name]};
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(3).scans = cellstr(spm_select('FPList',['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\Encoding\Bold\Block_3'],'^rarf2*'));
    matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi = {['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\Encoding\MCF_EncodeB3_singleTrial_subj_',folderName,'.mat']};
    matlabbatch{1}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
    x = dir(['M:\YOURPATH\DATA\Subj_',folderName,'\Day1\Encoding\Bold\Block_3\','rp_f*']);
    matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi_reg = {[x(1).folder,'\',x(1).name]};%movement regressors
    matlabbatch{1}.spm.stats.fmri_spec.sess(3).hpf = 128;
    
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.6;
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; 
    
    batch{batch_counter} = matlabbatch;
        batch_counter = batch_counter+1;
end


 poolobj = gcp('nocreate');
delete(poolobj);
parfor i = 1:length(batch)
    spm_jobman('run', batch{i})
    disp(['Finished jo for subj.', i]);
end