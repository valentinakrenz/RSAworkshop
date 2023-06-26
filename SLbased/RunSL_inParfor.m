%% RUNS SEARCHLIGHT BASED RSA
% set paths, variable names and indices in compute_SL
% function is in the same directory as script %otherwise:
% add directory of function that you want to run in your (par)for-loop

clear all
delete(gcp('nocreate'))

cd(pwd) %set directory to current script
addpath('..\NiftiTools\') % go to parent folder and add NiftiTools

%% choose your SL radius
radius = 2; %in voxel
fisher_transform = 1; %1 = yes, 0=0
sub = 1:4;

%% run SL RSA
% change to for, if you can't use parfor loop

parfor s =1:length(nSjs)
    sj = sub(s);
    compTime= RSA_SL_EES(sj,radius, fisher_transform) 
end
delete(gcp('nocreate'))