function compTime = RSA_SL_EES(sj,radius, fisher_transformed)

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
% Session constants=421:426;

%% DEFINE PATHS

addpath NiftiTools\
dataFolder='data\mriData\';%change accordingly
SubFolder='\trialwiseGLM\';
outputDir = 'results\RSA_SL';

%% READ IN PATHS AND FILES

% Check if the directory exists
if ~isfolder(outputDir)
    % Create the directory if it doesn't exist
    mkdir(outputDir);
end

SJfolders=dir([dataFolder 'sj*']);
betaPath=[dataFolder,SJfolders(sj).name,SubFolder]; % analyze sj that is defined by RunSL_inParfor
betaFiles=dir([betaPath, '*.nii']);
nii_template = load_nii('beta_0001.nii'); 

sj 
tic

for b=1:numel(betaFiles)
    betaDat=load_nii([betaPath,betaFiles(b).name]);
    if b==1
        BetaMaps=zeros([numel(betaFiles),size(betaDat.img)],'single'); % if b==1 BetaMaps is initialized
    end
    BetaMaps(b,:,:,:) = betaDat.img; %assign betaDat to the bth element of the BetaMaps variable
end

%% DEFINE COMPARISONS AND INDICES
emotions = {'Neg', 'Neu'};
Neg_SelDat = BetaMaps([1:30, 61:90, 121:150], :, :, :); % select only negative items at encoding
Neu_SelDat = BetaMaps([31:60, 91:120, 151:180], :, :, :); % select only neutral items at encoding
% define indices for each run
Enc1_Inds = 1:30; 
Enc2_Inds = 31:60;
Enc3_Inds = 61:90;
% save run indices
enc_runs = {Enc1_Inds, Enc2_Inds, Enc3_Inds};

%% RUN SEARCHLIGHT FOR EACH CONDITION AND SAVE NIFTI

for e = 1:length(emotions)
    emotion = emotions{e};
    
    if strcmp(emotion, 'Neg')
        SelDat = Neg_SelDat;
    elseif strcmp(emotion, 'Neu')
        SelDat = Neu_SelDat;
    end

    for i = 1:length(enc_runs)-1 % don't take the last encoding run to avoid duplicate
            % only place SLs within the brain
            SLmask = squeeze(std(SelDat)) ~= 0 & ~squeeze(isnan(std(SelDat)));
            [SLindsLinear, LinearInds, NrOfValidVoxelsPerSL] = GenerateAll_SL_Inds(radius, SLmask);
    
        for j = i+1:length(enc_runs)
            Inds_1 = enc_runs{i};
            Inds_2 = enc_runs{j};

            diffStim_PerSL = NaN(size(SLmask), 'single');
            sameStim_PerSL = NaN(size(SLmask), 'single');
            
            HypMat = eye(length(Inds_1)); % define diagonal

            for sl = 1:numel(LinearInds)
                %tic
                slInds = SLindsLinear(:, sl);
                slInds = slInds(~isnan(slInds));

                SL_Patts = SelDat(:, slInds);

                RSM = corr(SL_Patts(Inds_1, :)', SL_Patts(Inds_2, :)');
                sameStim_mean = mean(RSM(HypMat == 1)); % sim between same item (on-diagonal)
                diffStim_mean = mean(RSM(HypMat == 0)); % sim between different items(off-diagonal)
                sameStim_PerSL(LinearInds(sl)) = sameStim_mean;
                diffStim_PerSL(LinearInds(sl)) = diffStim_mean;
            end

            if fisher_transformed == 1
                
                sameStim_PerSL = atanh(sameStim_PerSL);
                diffStim_PerSL = atanh(diffStim_PerSL);
                value = 'FisherZ';
            
            elseif fisher_transformed == 0
                value = 'PearsonR';
            end 
            
            % Save each SL_PerSL_fishertransformed with the corresponding
            % emotion and run-combination in ffilename
            nii_template.img = sameStim_PerSL;
            save_nii(nii_template, [outputDir filesep 'EES_sameStim_PerSL_' emotion '_Enc' num2str(i) 'Enc' num2str(j) '_' value '_Rad' num2str(radius) '_SJ' sprintf('%03d', sj) '.nii']);

            nii_template.img = diffStim_PerSL;
            save_nii(nii_template, [outputDir filesep 'EES_diffStim_PerSL_' emotion '_Enc' num2str(i) 'Enc' num2str(j) '_' value '_Rad' num2str(radius) '_SJ' sprintf('%03d', sj) '.nii']);

        end
    end
end

compTime=toc