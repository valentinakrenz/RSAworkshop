%% back transformation of ROIs for multivariate analyses of fMRI data in native space
% requires SPM12, marsbar toolbox parfor toolbox (otherwise change to for loops)

% script by Hendrik Heinbockel, adapted by Valentina Krenz

%%
clear all;
su = [1:99] ;


%%  1   Segement structural T1 to get the inverse deformation field

% batch_counter =1; %if you don't want to parallelize
clear batch;

for sub =1:length(su)
    clear matlabbatch;
    
    subj = su(sub);
    folderName=sprintf('%02.0f',subj);
    
    if ismember(str2num(folderName), excluded_su)
        continue;
    end
    
    delete(['M:\YOURPATH\DATA\Subj_',folderName,'\Day2\struct\*.mat']);
    
    struct =  {spm_select('FPList',['M:\YOURPATH\DATA\Subj_',folderName,'\Day2\struct\'],'^rs20*')};
    
    matlabbatch{1}.spm.spatial.preproc.channel.vols = struct;
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
    %gray matter
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'C:\spm12\tpm\TPM.nii,1'};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    % white matter
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'C:\spm12\tpm\TPM.nii,2'};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    % CSF
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'C:\spm12\tpm\TPM.nii,3'};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'C:\spm12\tpm\TPM.nii,4'};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'C:\spm12\tpm\TPM.nii,5'};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'C:\spm12\tpm\TPM.nii,6'};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
    matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
    matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
        NaN NaN NaN];
    
    spm_jobman('run', matlabbatch); % comment out if you want to parallelize
    
% %    if you want to parallelized
%      batch{batch_counter} = matlabbatch;
        
%     batch_counter =  batch_counter +1;
%     parfor i = 1:length(batch)
%     spm_jobman('run', batch{i})
%     disp(['Finished specify for job no.', i]);
%     end

end


%% MNI space ROI gets backwards-normalized to native space

clear all; close all; clc;
always = 'M:\YOURPATH\DATA\CONVERTEDROIS';

% e.g. two ROis from the Harvard oxford atlas (which are now in MNI
%  space) and are supposed to be warped to native space

ROIS={
    'yourROI1.nii', 'yourROI2.nii' };

batch_counter = 1;
clear batch;

for sub =1:length(su) %loop over subjects
    batch_counter = 1;
    clear batch
    
    for R = 1:length(ROIS)  %loop over all (both) ROis
        clear matlabbatch
        
        subj = su(sub);
        folderName=sprintf('%02.0f',subj);
        
        %copy the MNI roi in a subject specific folder
        copyfile('W:\PATHTOATLASROI...',[always,'\Subj_',folderName,'\Day2\struct\']);
        
        %Path to the inverse deformation field
        x = dir(['M:\YOURPATH\DATA\Subj_',folderName,'\Day2\struct\','iy_rs*']);
        
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr([x(1).folder,'\',x(1).name]);
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr([always,'\Subj_',folderName,'\Day2\struct\',ROIS{R}]);
        
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-100 -150 -100
            100 100 100];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 2;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'subjnorm';
        
        spm_jobman('run', matlabbatch);
%         batch{batch_counter} = matlabbatch;
        
%         batch_counter =  batch_counter +1;
%     end
%     parfor i = 1:length(batch)
%         spm_jobman('run', batch{i})
%         disp(['Finished specify for job no.', i]);
%     end
    end 
end

%% coregister ROIs to individual beta image

for sub =1:length(su)
    batch_counter = 1;
    clear batch
    
    for R = 1:length(ROIS)
        clear matlabbatch
        
        subj = su(sub);
        folderName=sprintf('%02.0f',subj);
        
        
        if ismember(str2num(folderName), excluded_su)
            continue;
        end
        
        m = dir([always,'\Subj_',folderName,'\Overall_First_Level_native\','beta*']);  %I coregistered to a beta image, because i later used the native space masks on beta images in RSA and MVPA
        beta = {};
        beta = {strcat(m(3).folder,'\',m(3).name)};
        
        
        matlabbatch{1}.spm.spatial.coreg.write.ref = beta;
        matlabbatch{1}.spm.spatial.coreg.write.source = cellstr([always,'\Subj_',folderName,'\Day2\struct\subjnorm',ROIS{R}]);
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 7;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
        
        spm_jobman('run', matlabbatch);
%         batch{batch_counter} = matlabbatch;
%         
%         batch_counter =  batch_counter +1;
    end
%     parfor i = 1:length(batch)
%         spm_jobman('run', batch{i})
%         disp(['Finished specify for job no.', i]);
%     end
    
end


%% Now marsbar is uses  grey matter mask which comes also from
%%segmenting to cut out every voxel from the mask, that is not really
%%covering grey matter (but white or csf)

marsbar %open marsbar UI

for sub =1:length(su)
    clear matlabbatch;
    
    for R = 1:length(ROIS)
        clear matlabbatch
        
        subj = su(sub);
        folderName=sprintf('%02.0f',subj);
        
        
        % For selecting images, later
        img_flt = mars_veropts('get_img_ext');
        d = [];
        imgname = ['M:\YOURPATH\DATA\Subj_',folderName,'\Day2\struct\rsubjnorm',ROIS{R}];
        [p, f e] = fileparts(imgname);
        func = '';
        d = f; l = f;
        d = [d ' - binarized'];
        l = [l '_bin'];
        o = maroi_image(struct('vol', spm_vol(imgname), 'binarize',1));
        o = maroi_matrix(o);
        o = descrip(o,d);
        o = label(o,l);
        
        saveroi(o, ['M:\YOURPATH\DATA\Subj_',folderName,'\Day2\struct\rm',ROIS{R}(1:end-4),'.mat'])
        save_as_image(o, ['M:\YOURPATH\DATA\Subj_',folderName,'\Day2\struct\rm',ROIS{R}])
        
        
        roilist = {o,om};
        for i = 1:length(roilist)
            eval(sprintf('r%d = roilist{%d};', i, i));
        end
        eval(['o=' ' r1 & r2 ' ';']);
        save_as_image(o, ['M:\YOURPATH\DATA\Subj_',folderName,'\Day2\struct\rmask',ROIS{R}])
        
    end
end


%% coregistrating this mask

for sub =1:length(su)
    batch_counter = 1;
    clear batch
    
    for R = 1:length(ROIS)
        clear matlabbatch
        
        subj = su(sub);
        folderName=sprintf('%02.0f',subj);
        
        
        if ismember(str2num(folderName), excluded_su)
            continue;
        end
        
        m = dir([always,'\Subj_',folderName,'\Overall_First_Level_native\','beta*']);
        beta = {};
        beta = {strcat(m(3).folder,'\',m(3).name)};
        
        
        matlabbatch{1}.spm.spatial.coreg.write.ref = beta;
        matlabbatch{1}.spm.spatial.coreg.write.source = cellstr([always,'\Subj_',folderName,'\Day2\struct\', ROIS{R}]);
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
        
        spm_jobman('run', matlabbatch);
%         batch{batch_counter} = matlabbatch;
%         
%         batch_counter =  batch_counter +1;
    end
%     parfor i = 1:length(batch)
%         spm_jobman('run', batch{i})
%         disp(['Finished specify for job no.', i]);
%     end
    
end
