%% RUNS SEARCHLIGHT BASED RSA
% set paths, variable names and indices in compute_SL

% %% prepare parallelization
% % adapt to your avaible number of workes
% try
%     delete(workers);
% catch
% end
% workers=parpool(3);%change to number of available kernels

%% choose your SL radius
radius = 2; %in voxel
fisher_transform = 1;

%% run SL RSA
% % change to for, if you can't use parfor loop
% 
% parfor sj=1:4
%     compTime= RSA_SL_ERS(sj,radius) 
% end
% delete(workers);

for sj=1:1
    compTime= RSA_SL_EES(sj,radius, fisher_transform) 
end
% delete(workers);