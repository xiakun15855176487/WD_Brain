%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_classification_analysis.m
%%%
%%% MATLAB script to perform SVM classify analysis using inhome matlab code
%%% Note: classify LWD and NWD
%%%       find the spectial brain patterns for these two groups
%%%
%%% Original: University of Science and Technology of China, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set path
output_path = './Results/SVM';

%% load data
%
load('LWD_data.mat')
load('NWD_data.mat')

%% construct train dataset
%
% get gray volume as train data
X_LWD = LWD_data.gv; 
X_NWD = NWD_data.gv;
X = [X_LWD;X_NWD];
X = zscore(X); % normalization

% set label for two groups
% LWD = 0, NWD = 1
Y = [zeros(size(X_LWD,1),1);ones(size(X_NWD,1),1)];

%% SVM analysis
%
% feature selected
features_idx = SVMRFE(X,Y);

% SVM train
X_Slected = X(:,features_idx);

[accu_mean,W_mean,cv,dist_total] = SVM_train(X_Slected,Y);

l_pad = Results.LWD.PAD_corrected;
n_pad = Results.NWD.PAD_corrected;
[R_l,P_l] = corr(l_pad,dist_total(1:78));
[R_n,P_n] = corr(n_pad,dist_total(79:end));

% permutation
accu_perm = zeros(1,10000);
W_perm = zeros(size(X_Slected,2),10000);
parfor i = 1:10000
    ordery = randperm(length(Y));
    Y_Perm = Y(ordery);

    [accu_perm(i),W_perm(:,i),~] = SVM_train(X_Slected,Y_Perm);
end

P_accu = length(find(accu_perm > accu_mean))/10000;


SVM_results_1113 = struct();
SVM_results_1113.features_idx = features_idx;
SVM_results_1113.W = W_mean;
SVM_results_1113.accu = accu_mean;
SVM_results_1113.accu_perm = accu_perm;
SVM_results_1113.P_accu = P_accu;
SVM_results_1113.CV = cv;
SVM_results_1113.dist_total = dist_total;

% save SVM results
save(fullfile(output_path,'SVM_results_1113.mat'),"SVM_results_1113")

%%  map weight back to 3D brain

% back to original feature size
W_full = zeros(size(X,2),1);
W_full(features_idx) = zscore(W_mean);

% load brain mask
gm_mask_file = fullfile('Data','Train_VBM','GM_mask_s.nii.gz');
gm_mask_nift = load_untouch_nii(gm_mask_file);
gm_mask_img = gm_mask_nift.img;
mask_idx = find(gm_mask_img);

W_3D = zeros(size(gm_mask_img));
W_3D(mask_idx) = W_full;

fname = 'SVM_weights.nii';
nifti_data = gm_mask_nift;
nifti_data.hdr.dime.datatype=16;
nifti_data.hdr.dime.bitpix=32;   
nifti_data.img = W_3D;
save_untouch_nii(nifti_data,fullfile(output_path,fname));


