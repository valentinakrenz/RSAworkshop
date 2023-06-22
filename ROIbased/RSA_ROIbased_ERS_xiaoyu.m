% ROI-BASED REPRESENTATIONAL SIMILARIRY ANALAYSIS
% CORRELATES ACTIVIATION PATTERNS OF EACH ITEM AT ENCODING WITH
% ACTIVATION PATTERNS WITH ITEMS DURING MEMORY TESTING IN A SPECIFIC ROI

% Valentina Krenz 2023

%% PREPARE AND READ IN DATA PATHS
% Goals: to get the corr for each subject, each condition, each image in
% each ROI

clear all 
addpath NiftiTools\

% output name and directory
fileName = 'ESmethod2';
outputDir = '.\results';
betaDir= 'data\mriData\'; % directory with subject folders
ROIdir= 'ROIs\'; % directory with ROI files
SubFolder='\trialwiseGLM\'; % name of subject sub-folder were beta files are stored
nrBetaMaps=420; % number of beta images
SJfolders=dir([betaDir 'sj*']); % get each subject folder
RoiNames=dir([ROIdir  '*.nii']); % take ROI names from ROI path

% % Check if the directory exists
% if ~isfolder(outputDir)
%     % Create the directory if it doesn't exist
%     mkdir(outputDir);
% end
filePath = fullfile(outputDir, [fileName '.xlsx']); % Construct the full file path with the desired extension

nSjs=numel(SJfolders); % number of subject files

%% RSM indices per stimulus type
% beta images are sorted by condition

% 1:30 encoding1 negative 11xx
% 31:60 encoding1 neutral 12xx

% 61:90 encoding2 negative 21xx
% 91:120 encoding2 neutral 22xx

% 121:150 encoding3 negative 31xx
% 151:180 encoding3 neutral 32xx

% 181:210 old negative 41xx
% 211:240 old neutral 42xx

% 241:270 perc negative
% 271:300 perc neutral

% 301:330 sem negative
% 331:360 sem neutral

% 361:390 unrelared negative
% 391:420 unrelated neutral

% imgType=NaN(nrBetaMaps,1);
% runSeq=NaN(nrBetaMaps,1);
% emotionType=NaN(nrBetaMaps,1);
% imageID=NaN(nrBetaMaps,1);

%% use four digits to represent each image
% first digit: runs 1-3
% second digit: emotion type
% third and fourth digit: image index 01-30
% for sequence 1-4  4xxx=recogntion
Sequence(1:240,1)=[ones(1,60)*1000 ones(1,60)*2000 ones(1,60)*3000 ones(1,60)*4000];
% for emotion type x1xx=negative x2xx=neutral
emotionType(1:240,1)=repmat([ones(1,30)*100,ones(1,30)*200],1,4);
% for image xxID
imageID(1:240,1)=repmat([1:30 1:30],1,4);
imgType=Sequence+emotionType+imageID; % denote each image

% e.g., the first negative image
% mod(imgType,1000)==101;
% (mod(imgType,1000)==101 & imgType<4000); %only during the encoding phase

% e.g., the first run
% fix(imgType/1000)==1;
% e.g., negative images in the first run
% fix(imgType/1000)==1 & fix(mod(imgType,1000)/100)==1;


%%

% initialize struct for different kinds of contrasts, i.e. within-subject conditions
% Contrasts=struct(); 
% contrast names should be 'RSAtype_emotion' as they are later split into 
% 2 variables by parts = strsplit(AllContrastNames{ctr}, '_') for the
% output file
% ! why we need Contrast? - for each condition
% Contrasts.ERS_negative = NaN(nrBetaMaps);% Encoding-Retrieval-Similarity negative
% Contrasts.ERS_neutral = NaN(nrBetaMaps);% Encoding-Retrieval-Similarity neutral
AllContrastNames = {'ERS_negative';'ERS_neutral'};% save contrast names

% create array to map stimulus indices to different contrasts
% ! why we need ImageInds? - to pick up the corr for a specific image
% for each image
% ImageInds=NaN(420);
% % create 180x420 matrix, where each row of the matrix contains the numbers 
% % 1 to 30, twice, and this pattern is repeated 420 times along the columns
% ImageInds(1:180,:)=repmat(repmat([1:30 1:30]',3,1),1,420); 

% ImageInds(1:240,:)=repmat(repmat([1:30 1:30]',4,1),1,420); % 
% not important. 180 rows are essential

%% DEFINE WHICH IMAGES ARE COMPARED

% for EncIm=1:30
%     % get encoding activity of one item over all encoding runs
%     ImInds_negative = [EncIm EncIm+60 EncIm+120];
%     ImInds_neutral = [EncIm+30 EncIm+90 EncIm+150];
%     
% %   % image indices if you want only encoding run 1
% %     ImInds_negative = EncIm; % run 1 neg
% %     ImInds_neutral = EncIm+30; % run 1 neut
%     
% %   % image indices if you want only encoding run 2
% %     ImInds_negative=EncIm+60; % run 2 neg
% %     ImInds_neutral=EncIm+90; % run 2 neut
% 
% %   % image indices if you want only encoding run 3
% %     ImInds_negative=EncIm+120; % run 3 neg
% %     ImInds_neutral=EncIm+150; % run 3 neut
% 
%     Contrasts.ERS_negative(ImInds_negative,EncIm+180) = 1; % Encoding-Retrieval-Similarity negative %原来的negative image和再认任务的image相同的时候，设置�?
%     Contrasts.ERS_neutral(ImInds_neutral,EncIm+210) = 1; % Encoding-Retrieval-Similarity neutral
%     
% end

%% PREPARE ROI FILES AND REPRESENTATIONAL-SIMILARITY-MATRICES (RSMs)

nRois=numel(RoiNames); % number of ROIs
ROIdata=struct(); % initialize struct to store ROI data information
ROInames=cell(nRois,1); % get ROI names

for roi=1:nRois %loop through each ROI
    ROI=load_nii([ROIdir RoiNames(roi).name]); % read in nifti with corresponding ROI
    ROIdata.([RoiNames(roi).name(1:end-4) '_Inds']) = find(ROI.img); % extract the voxel indices where the ROI is present 
    % participants*420trials*420trials
    ROIdata.([RoiNames(roi).name(1:end-4) '_RSMs']) = zeros(numel(SJfolders),nrBetaMaps,nrBetaMaps,'single'); % initialize RSM for each subject
    ROInames{roi} = RoiNames(roi).name(1:end-4); % name according to file name
end

%% COMPUTE REPRESENTATIONAL SIMILARITY MATRICES
% read in the beta images for each subject, organize them into a 4D array, 
% and compute a representational similarity matrix (RSM) for each ROI by 
% correlating the voxel-wise activity patterns across images 
% Store RSMs in a structure (ROI.data)

for sj=1:nSjs   
    sj; % prints current subject to command window
    tic % stopwatch timer
    % load each participant' beta files:
    betaPath=[betaDir,SJfolders(sj).name,SubFolder]; % constructs path to the current sj's beta images
    betaFiles=dir([betaPath, '*.nii']); % read in all beta files from the current sj's folder
    % aggregate all the trials' beta maps
    for b=1:numel(betaFiles) % iterate over all beta images of the current sj
        betaDat=load_nii([betaPath,betaFiles(b).name]); % load current beta image into workspace
        if b==1
            % 420*79*95*79
            BetaMaps=zeros([numel(betaFiles),size(betaDat.img)],'single'); % if its the first beta image, initialize a 4D array to hold all beta images of the current sj with dimensions: number of beta images x size of 2d beta image
        end
        BetaMaps(b,:,:,:)= betaDat.img; % add current beta image to the previously initialized 4D array at the appropriate index
    end     
    for r=1:nRois
        %420*766 trial*ROIvoxel
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
%             currData = atanh(mean(roiRSMs(:,Contrasts.(AllContrastNames{ctr})==1&ImageInds==stim),2));
            
            switch emotion
                case 'negative'
                    type2=100; % the second digit represents the emotion; x1xx = negative images
                case 'neutral'
                    type2=200;% x2xx = neutral images
            end
            % 
            stimuli_x_type=type2+stim; % to form the number for the specific image
            stimuli_y_type=4000+stimuli_x_type;% the number for recogntion image
            
            % corr_x_index:(mod(imgType,1000)==stimuli_x_type & imgType<4000 )
            % the negative/neutral images during the encoding phase
            
            % corr_y_index:imgType==stimuli_y_type
            % the specific image during the retrieval phase
            
            corr_roi=roiRSMs(:,(mod(imgType,1000)==stimuli_x_type & imgType<4000 ),imgType==stimuli_y_type);
            currData = atanh(mean(corr_roi,2));
            
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