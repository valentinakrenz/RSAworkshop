%% COMPARES MODEL-REPRESENTATIONAL-SIMILARITY-MATRIX (RSM) IN A NEURAL ROI  WITH A MODEL RSM

%BetaIndices per condition
% 1:30 encoding1 negative
% 31:60 encoding1 neutral
% 61:90 encoding2 negative
% 91:120 encoding2 neutral
% 121:150 encoding3 negative
% 151:180 encoding3 neutral
% 181:210 old negative at recognition
% 211:240 old neutral at recognition
% 241:270 perceptually related negative lures
% 271:300 perceptually related neutral lures
% 301:330 semantically related negative lures
% 331:360 semantically related neutral lures
% 361:390 unrelared negative lures
% 391:420 unrelated neutral lures

%% DEFINE YOUR PATHs
clear all
addpath NiftiTools\
dataFolder='.\data\mriData\'; %fmri data folder name
SubFolder='\trialwiseGLM\'; %firstlevel folder name
load('.\data\runDat.mat')  %load runDat struct
SJfolders=dir([dataFolder 'sj*']); %adapt to naming of your sj data folders
nSjs=numel(SJfolders); % how many sjs do you want to analyze
ROIfolder='.\ROIs\'; %path to your ROI files
RoiImNames=dir([ROIfolder '*.nii']); %read in all ROI files

%% DEFINE CONDITION INDICES
NegInds=[1:30,[1:30]+60, [1:30]+120,[1:30]+180];
NeutInds=[31:60,[31:60]+60, [31:60]+120,[31:60]+180];
OldInds=1:60;
PercInds=61:120;
SemInds=121:180;
NewInds=181:240;

%% DEFINE MODEL RSMs
% NEGATIVE ITEMS
% model 1: old items are distinct to all lures
OldAreSim_RestDiff_Neg=NaN(240);
OldAreSim_RestDiff_Neg(intersect(OldInds,NegInds),intersect(OldInds,NegInds))=.8; % old items are more similar to each other than all other items
OldAreSim_RestDiff_Neg(intersect(OldInds,NegInds),intersect(SemInds,NegInds))=0;
OldAreSim_RestDiff_Neg(intersect(OldInds,NegInds),intersect(PercInds,NegInds))=0;
OldAreSim_RestDiff_Neg(intersect(OldInds,NegInds),intersect(NewInds,NegInds))=0;
OldAreSim_RestDiff_Neg(intersect(SemInds,NegInds),intersect(OldInds,NegInds))=0;
OldAreSim_RestDiff_Neg(intersect(PercInds,NegInds),intersect(OldInds,NegInds))=0;
OldAreSim_RestDiff_Neg(intersect(NewInds,NegInds),intersect(OldInds,NegInds))=0;
% remove the diagonal!
OldAreSim_RestDiff_Neg(eye(size(OldAreSim_RestDiff_Neg))==1)=NaN;    
figure;imagesc(OldAreSim_RestDiff_Neg);

% model 2: old items and semantically related items are similar
OldAreSimToSem_Neg=NaN(240);
OldAreSimToSem_Neg(intersect(OldInds,NegInds),intersect(OldInds,NegInds))= .8; % old items are more similar to each other than all other items
OldAreSimToSem_Neg(intersect(OldInds,NegInds),intersect(SemInds,NegInds))= .5; % old items are also similar to semantically related items
OldAreSimToSem_Neg(intersect(SemInds,NegInds),intersect(OldInds,NegInds))= .5; % old items are also similar to semantically related items
OldAreSimToSem_Neg(intersect(OldInds,NegInds),intersect(PercInds,NegInds))= 0;
OldAreSimToSem_Neg(intersect(OldInds,NegInds),intersect(NewInds,NegInds))= 0;
OldAreSimToSem_Neg(intersect(PercInds,NegInds),intersect(OldInds,NegInds))= 0;
OldAreSimToSem_Neg(intersect(NewInds,NegInds),intersect(OldInds,NegInds))= 0;
% remove the diagonal
OldAreSimToSem_Neg(eye(size(OldAreSimToSem_Neg))==1)=NaN;
figure;imagesc(OldAreSimToSem_Neg);

% model 3: old items and perceptually related items are similar
OldAreSimToPerc_Neg=NaN(240);
OldAreSimToPerc_Neg(intersect(OldInds,NegInds),intersect(OldInds,NegInds))= .8; % old items are more similar to each other than all other items
OldAreSimToPerc_Neg(intersect(OldInds,NegInds),intersect(PercInds,NegInds))=.5; % old items are also similar to perceptually related items
OldAreSimToPerc_Neg(intersect(PercInds,NegInds),intersect(OldInds,NegInds))=.5; % old items are also similar to perceptually related items
OldAreSimToPerc_Neg(intersect(OldInds,NegInds),intersect(SemInds,NegInds))= 0;
OldAreSimToPerc_Neg(intersect(OldInds,NegInds),intersect(NewInds,NegInds))= 0;
OldAreSimToPerc_Neg(intersect(SemInds,NegInds),intersect(OldInds,NegInds))= 0;
OldAreSimToPerc_Neg(intersect(NewInds,NegInds),intersect(OldInds,NegInds))= 0;
% remove the diagonal
OldAreSimToPerc_Neg(eye(size(OldAreSimToPerc_Neg))==1)=NaN;
figure;imagesc(OldAreSimToPerc_Neg);

% NEUTRAL ITEMS
% model 1: old items are distinct to all lures
OldAreSim_RestDiff_Neut=NaN(240);
OldAreSim_RestDiff_Neut(intersect(OldInds,NeutInds),intersect(OldInds,NeutInds))=.8; % old items are more similar to each other than all other items
OldAreSim_RestDiff_Neut(intersect(OldInds,NegInds),intersect(SemInds,NegInds))=0;
OldAreSim_RestDiff_Neut(intersect(OldInds,NegInds),intersect(PercInds,NegInds))=0;
OldAreSim_RestDiff_Neut(intersect(OldInds,NegInds),intersect(NewInds,NegInds))=0;
OldAreSim_RestDiff_Neut(intersect(SemInds,NegInds),intersect(OldInds,NegInds))=0;
OldAreSim_RestDiff_Neut(intersect(PercInds,NegInds),intersect(OldInds,NegInds))=0;
OldAreSim_RestDiff_Neut(intersect(NewInds,NegInds),intersect(OldInds,NegInds))=0;
% remove the diagonal
OldAreSim_RestDiff_Neut(eye(size(OldAreSim_RestDiff_Neut))==1)=NaN;
figure;imagesc(OldAreSim_RestDiff_Neut);

% model 2: old items and semantically related items are similar
OldAreSimToSem_Neut=NaN(240);
OldAreSimToSem_Neut(intersect(OldInds,NeutInds),intersect(OldInds,NeutInds))=.8; % old items are more similar to each other than all other items
OldAreSimToSem_Neut(intersect(OldInds,NeutInds),intersect(SemInds,NeutInds))=.5; % old items are also similar to semantically related items
OldAreSimToSem_Neut(intersect(SemInds,NeutInds),intersect(OldInds,NeutInds))=.5; % old items are also similar to semantically related items
OldAreSimToSem_Neut(intersect(OldInds,NeutInds),intersect(PercInds,NeutInds))= 0;
OldAreSimToSem_Neut(intersect(OldInds,NeutInds),intersect(NewInds,NeutInds))= 0;
OldAreSimToSem_Neut(intersect(PercInds,NeutInds),intersect(OldInds,NeutInds))= 0;
OldAreSimToSem_Neut(intersect(NewInds,NeutInds),intersect(OldInds,NeutInds))= 0;
% remove the diagonal
OldAreSimToSem_Neut(eye(size(OldAreSimToSem_Neut))==1)=NaN;
figure;imagesc(OldAreSimToSem_Neut);

% model 3: old items and perceptually related items are similar
OldAreSimToPerc_Neut=NaN(240);%fill all items with NaNs
OldAreSimToPerc_Neut(intersect(OldInds,NeutInds),intersect(OldInds,NeutInds))= .8; % old items are more similar to each other than all other items
OldAreSimToPerc_Neut(intersect(OldInds,NeutInds),intersect(PercInds,NeutInds))=.5; % old items are also similar to perceptually related items
OldAreSimToPerc_Neut(intersect(PercInds,NeutInds),intersect(OldInds,NeutInds))=.5; % old items are also similar to perceptually related items
OldAreSimToPerc_Neut(intersect(OldInds,NeutInds),intersect(SemInds,NeutInds))= 0;
OldAreSimToPerc_Neut(intersect(OldInds,NeutInds),intersect(NewInds,NeutInds))= 0;
OldAreSimToPerc_Neut(intersect(SemInds,NeutInds),intersect(OldInds,NeutInds))= 0;
OldAreSimToPerc_Neut(intersect(NewInds,NeutInds),intersect(OldInds,NeutInds))= 0;
% remove the diagonal
OldAreSimToPerc_Neut(eye(size(OldAreSimToPerc_Neut))==1)=NaN;
figure;imagesc(OldAreSimToPerc_Neut);

%% PREPARE OUTPUT
% save all contrasts in a cell mat, and save the names

AllModelNames={'oldAreSimRestDiff_Neg','oldAreSimToSem_Neg','oldAreSimToPerc_Neg','oldAreSimRestDiff_Neu','oldAreSimToSem_Neu','oldAreSimToPerc_Neu'};
AllModels={OldAreSim_RestDiff_Neg,OldAreSim_RestDiff_Neut,OldAreSimToPerc_Neg,OldAreSimToPerc_Neut,OldAreSimToSem_Neg,OldAreSimToSem_Neut};

for i=1:numel(AllModelNames)
    subplot(2,4,i);
    imagesc(AllModels{i});
    title(AllModelNames{i});
end

OPcollumNames={};
count=1;
for r=1:numel(RoiImNames)
    for i=1:numel(AllModelNames)
        OPcollumNames{count}=[RoiImNames(r).name(1:end-4) '_' AllModelNames{i}];
        count=count+1;
    end
end

OPsimVals=zeros(nSjs,numel(OPcollumNames));
DataTable=array2table(OPsimVals);
DataTable.Properties.VariableNames=OPcollumNames;

%% RUN MODEL-COMPARISON RSA
for sj= 1:nSjs
    % compute runEffRSM
    SjNumber=str2num(SJfolders(sj).name(3:5)); %extract sj number from folder and convert it from string to number so it fits to our sj variables
    SjTrialTimingDat=runDat.(['vp' num2str(SjNumber)]); %acess trial timing data for the current sj from runDat file
    RunEffRSM=zeros(420); %initialize a RSM for the run-dependent similarity
    for tr1=1:240 %iterate through each pair of trial pairs (t1-t2)
        for tr2=1:240 
            trialRunCode=(SjTrialTimingDat(tr1,2)*10)+... %compute all possible run-combinations by taking the information in which trial each stimulus was presented from runDat
                SjTrialTimingDat(tr2,2);
            RunEffRSM(tr1+180,tr2+180)=trialRunCode; %consider only items from recognition task, i.e. only betas after 180 in RunEffRSM
        end
    end
    RunEffRSM(eye(420)==1)=0; % remove diagonal
    AllRetrRunCombos=unique(RunEffRSM); %extract the unique values of the run effect RSM representing different run codes
    AllRetrRunCombos=AllRetrRunCombos(AllRetrRunCombos>0); %filter out zero and negative values to keep only positive run effect codes
    RunEffRSM=RunEffRSM(181:end,181:end); %remove the first 180 rows and columns as the analysis concerns only recognition task
    
    sj
    tic
    
    % load in data for current sj
    betaPath=[dataFolder,SJfolders(sj).name,SubFolder]; %get beta path from this sj
    betaFiles=dir([betaPath, '*.nii']);%get beta images from this sj
    for b=181:420 %take only betas from memory test
        betaDat=load_nii([betaPath,betaFiles(b).name]); %load beta image for current trial
        if b==181
            BetaMaps=zeros([240,size(betaDat.img)],'single'); %initialize 4D matrix to store beta maps dor sj in first trial of memory test
        end
        BetaMaps(b-180,:,:,:)= betaDat.img; %fill 4D matrix with beta map for current trial
    end
    
    % define lower diagonal
    lowDiagSel=tril(ones(240),-1)==1; % take only the lower diagonal
    
    % run rsa for each roi
    for roi=1:numel(RoiImNames) %loop through rois
        roi
        % extract roi-patterns from each beta
        roiMap=load_nii([ROIfolder, RoiImNames(roi).name]);%load current roi
        roiInds=find(roiMap.img(:)==1); %load indices of current ROI
        ROI_Patts=BetaMaps(:,roiInds);%extract neural activation patterns in current ROI
        
        % correlate activation patterns between trials in specific ROI and
        % save in RSM with correlation of each pair of activation patterns
        RSM=corr(ROI_Patts','rows','complete'); %correlation is computed between complete rows, i.e. NaNs are discarted
        
        % remove retrieval run effect
        % compute mean correlation for a specific run-combination (e.g. all items from run1 and run2 correlated) and
        % subtract this run-related correlation from items in RSM that
        % belong to this run-combination, i.e. mean run1-run2 correlation
        % will be substraced from each correlation in RSM between items in
        % run 1 and run2
        for runcomp=1:numel(AllRetrRunCombos)% go through all run-combinations, i.e. run1-run1, run1-run2 etc.
            RSM(RunEffRSM==AllRetrRunCombos(runcomp))=RSM(RunEffRSM==AllRetrRunCombos(runcomp))-... 
                mean(RSM(RunEffRSM==AllRetrRunCombos(runcomp)));
        end
        
        % correlate neural RSM with model RSMs using only lower diagonal of
        % each RSM and spearman's rho and apply Fisher z-transformation
        for mdl=1:numel(AllModels)
            DataTable{sj, mdl+((roi-1)*numel(AllModels))}=atanh(corr(RSM(lowDiagSel),AllModels{mdl}(lowDiagSel),'rows','complete','type','spearman')); 
        end
    end
end
writetable(DataTable,'RetrRetrModelFits.xlsx');
