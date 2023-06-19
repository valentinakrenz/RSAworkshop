% ROI-BASED REPRESENTATIONAL SIMILARIRY ANALAYSIS
% CORRELATES ACTIVIATION PATTERNS OF EACH ITEM AT ENCODING WITH
% ACTIVATION PATTERNS WITH ITEMS DURING MEMORY TESTING IN A SPECIFIC ROI

% Valentina Krenz 2023

%% PREPARE AND READ IN DATA PATHS
clear all 
addpath NiftiTools\

% output name and directory
fileName = 'ERS';
outputDir = '.\results';
betaDir= 'data\mriData\'; % directory with subject folders
ROIdir= 'ROIs\'; % directory with ROI files
SubFolder='\trialwiseGLM\'; % name of subject sub-folder were beta files are stored
nrBetaMaps=420; % number of beta images
SJfolders=dir([betaDir 'sj*']); % get each subject folder
RoiNames=dir([ROIdir  '*.nii']); % take ROI names from ROI path

% Check if the directory exists
if ~isfolder(outputDir)
    % Create the directory if it doesn't exist
    mkdir(outputDir);
end
filePath = fullfile(outputDir, [fileName '.xlsx']); % Construct the full file path with the desired extension

nSjs=numel(SJfolders); % number of subject files

%% RSM indices per stimulus type
% beta images are sorted by condition

% 1:30 encoding1 negative
% 31:60 encoding1 neutral

% 61:90 encoding2 negative
% 91:120 encoding2 neutral

% 121:150 encoding3 negative
% 151:180 encoding3 neutral

% 181:210 old negative
% 211:240 old neutral

% 241:270 perc negative
% 271:300 perc neutral

% 301:330 sem negative
% 331:360 sem neutral

% 361:390 unrelared negative
% 391:420 unrelated neutral

% initialize struct for different kinds of contrasts, i.e. within-subject conditions
Contrasts=struct(); 
% contrast names should be 'RSAtype_emotion' as they are later split into 
% 2 variables by parts = strsplit(AllContrastNames{ctr}, '_') for the
% output file
Contrasts.ERS_negative = NaN(nrBetaMaps);% Encoding-Retrieval-Similarity negative
Contrasts.ERS_neutral = NaN(nrBetaMaps);% Encoding-Retrieval-Similarity neutral
AllContrastNames = fields(Contrasts);% save contrast names

% create array to map stimulus indices to different contrasts
ImageInds=NaN(420);
% create 180x420 matrix, where each row of the matrix contains the numbers 
% 1 to 30, twice, and this pattern is repeated 420 times along the columns
ImageInds(1:180,:)=repmat(repmat([1:30 1:30]',3,1),1,420); 

%% DEFINE WHICH IMAGES ARE COMPARED

for EncIm=1:30
    % get encoding activity of one item over all encoding runs
    ImInds_negative = [EncIm EncIm+60 EncIm+120];
    ImInds_neutral = [EncIm+30 EncIm+90 EncIm+150];
    
%   % image indices if you want only encoding run 1
%     ImInds_negative = EncIm; % run 1 neg
%     ImInds_neutral = EncIm+30; % run 1 neut
    
%   % image indices if you want only encoding run 2
%     ImInds_negative=EncIm+60; % run 2 neg
%     ImInds_neutral=EncIm+90; % run 2 neut

%   % image indices if you want only encoding run 3
%     ImInds_negative=EncIm+120; % run 3 neg
%     ImInds_neutral=EncIm+150; % run 3 neut

    Contrasts.ERS_negative(ImInds_negative,EncIm+180) = 1; % Encoding-Retrieval-Similarity negative
    Contrasts.ERS_neutral(ImInds_neutral,EncIm+210) = 1; % Encoding-Retrieval-Similarity neutral
    
end

%% PREPARE ROI FILES AND REPRESENTATIONAL-SIMILARITY-MATRICES (RSMs)

nRois=numel(RoiNames); % number of ROIs
ROIdata=struct(); % initialize struct to store ROI data information
ROInames=cell(nRois,1); % get ROI names

for roi=1:nRois %loop through each ROI
    ROI=load_nii([ROIdir RoiNames(roi).name]); % read in nifti with corresponding ROI
    ROIdata.([RoiNames(roi).name(1:end-4) '_Inds']) = find(ROI.img); % extract the voxel indices where the ROI is present 
    ROIdata.([RoiNames(roi).name(1:end-4) '_RSMs']) = zeros(numel(SJfolders),nrBetaMaps,nrBetaMaps,'single'); % initialize RSM for each subject
    ROInames{roi} = RoiNames(roi).name(1:end-4); % name according to file name
end

%% COMPUTE REPRESENTATIONAL SIMILARITY MATRICES
% read in the beta images for each subject, organize them into a 4D array, 
% and compute a representational similarity matrix (RSM) for each ROI by 
% correlating the voxel-wise activity patterns across images 
% Store RSMs in a structure (ROI.data)

for sj=1:nSjs   
    sj % prints current subject to command window
    tic % stopwatch timer
    betaPath=[betaDir,SJfolders(sj).name,SubFolder]; % constructs path to the current sj's beta images
    betaFiles=dir([betaPath, '*.nii']); % read in all beta files from the current sj's folder
    for b=1:numel(betaFiles) % iterate over all beta images of the current sj
        betaDat=load_nii([betaPath,betaFiles(b).name]); % load current beta image into workspace
        if b==1
            BetaMaps=zeros([numel(betaFiles),size(betaDat.img)],'single'); % if its the first beta image, initialize a 4D array to hold all beta images of the current sj with dimensions: number of beta images x size of 2d beta image
        end
        BetaMaps(b,:,:,:)= betaDat.img; % add current beta image to the previously initialized 4D array at the appropriate index
    end     
    for r=1:nRois
        RoiBetas=BetaMaps(:,ROIdata.([RoiNames(r).name(1:end-4) '_Inds'])); %extract beta values of the voxels that fall within the current ROI for all beta images
        % remove voxels within the current ROI that have NaN values
        RoiBetas(:,isnan(mean(RoiBetas,1)))=[]; % computing the mean beta value across images for each voxel, checking if result is NaN and removing the columns (voxels) where result is NaN
        % COMPUTE RSM FOR THE CURRENT ROI AND SJ  
        ROIdata.([RoiNames(r).name(1:end-4) '_RSMs'])(sj,:,:)=corr(squeeze(RoiBetas)'); %correlates beta values across images for each voxel within the ROI
    end
    toc % stop stopwatch timer and display elapsed time
end

%% COMPUTE AND EXPORT FISHER Z-TRANSFORMED PEARSON'S R 
% compute mean similarity per sj, ROI, contrast, item number and saved the
% fisher z-transformed value in a table 
% apply some restructuring of the table before exporting to excel

% Initialize empty arrays to store the data
data = [];
sjs = [];
rois = [];
contrasts = [];
RSAtypes = [];
emotions = [];
stims = [];
items = [];

allSjs = {SJfolders(:).name}'; % collect all subject names

for roi=1:nRois
    for ctr=1:numel(AllContrastNames)
        
        % Extract RSA type and emotion from contrast name
        parts = strsplit(AllContrastNames{ctr}, '_'); % Split the string on the underscore
        rsa_type = parts{1}; % The RSA type is the first part
        emotion = parts{2}; % The emotion is the second part
        
        for stim=1:30
            % get RSM data from current ROI
            roiRSMs=ROIdata.([ROInames{roi} '_RSMs']);
            % compute mean and Fisher z-transform r-values for the
            % specific sub condition
            currData = atanh(mean(roiRSMs(:,Contrasts.(AllContrastNames{ctr})==1&ImageInds==stim),2));
            
            item = strcat(emotion, '_', sprintf('%02d', stim)); % create item string
            
            % Append the data and variables to the respective arrays
            data = [data; currData];
            sjs = [sjs; allSjs];
            rois = [rois; repmat(ROInames(roi), numel(currData), 1)];
            emotions = [emotions; repmat({emotion}, numel(currData), 1)];
            RSAtypes = [RSAtypes; repmat({rsa_type}, numel(currData), 1)];
            items = [items; repmat({item}, numel(currData), 1)];
            contrasts = [contrasts; repmat(AllContrastNames(ctr), numel(currData), 1)];
            
        end
    end
end

% Create the table and name headers
OPinTable = table(sjs, rois, contrasts, RSAtypes, emotions, items, data, ...
    'VariableNames', {'sj', 'ROI', 'contrast', 'RSAtype', 'emotion', 'item', 'corr'});

% Write the table to a file
writetable(OPinTable, filePath);