%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_BA_prediction_analysis.m
%%%
%%% MATLAB script to how to do brain age prediction using inhome matlab code. 
%%% 
%%% Original: University of Science and Technology of China, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Narmalization 
load('Train_data.mat'); % Train_data.mat is a struct file including gray matter volume and demographic information.
[Normalized_Train_gv,PS] = mapstd(Train_data.gv');
Train_data.gv = Normalized_Train_gv';

% Apply to holdset, validation and patients data
load('Holdset_data.mat')
load('Validation_data_W2.mat')
load('Validation_data_W3.mat')
load('HC_data.mat')
load('CWD_data.mat')
load('HWD_data.mat')
load('LWD_data.mat')
load('NWD_data.mat')

% Holdset and Validation
Holdset_data.gv = mapstd('apply',Holdset_data.gv',PS)';
Validation_data_W2.gv = mapstd('apply',Validation_data_W2.gv',PS)';
Validation_data_W3.gv = mapstd('apply',Validation_data_W3.gv',PS)';

% WD
HC_data.gv = mapstd('apply',HC_data.gv',PS)';
CWD_data.gv = mapstd('apply',CWD_data.gv',PS)';
HWD_data.gv = mapstd('apply',HWD_data.gv',PS)';
LWD_data.gv = mapstd('apply',LWD_data.gv',PS)';
NWD_data.gv = mapstd('apply',NWD_data.gv',PS)';
%% Feature Dimension Reduction
%
[Coef_PCA,ScoreTrain,~,~,explained,PCAmu] = pca(Train_data.gv,'centered',true); 
Coef_PCA = Coef_PCA(:,1:500);

% Train
Train_data.gv = (Train_data.gv-PCAmu)*Coef_PCA;

% Holdset and Validation
Holdset_data.gv = (Holdset_data.gv-PCAmu)*Coef_PCA;
Validation_data_W2.gv = (Validation_data_W2.gv-PCAmu)*Coef_PCA;
Validation_data_W3.gv = (Validation_data_W3.gv-PCAmu)*Coef_PCA;

% WD
HC_data.gv = (HC_data.gv-PCAmu)*Coef_PCA;
CWD_data.gv = (CWD_data.gv-PCAmu)*Coef_PCA;
HWD_data.gv = (HWD_data.gv-PCAmu)*Coef_PCA;
LWD_data.gv = (LWD_data.gv-PCAmu)*Coef_PCA;
NWD_data.gv = (NWD_data.gv-PCAmu)*Coef_PCA;

%% Model Training
%
[ModelCoef_set,FitInfo] = lassoglm(Train_data.gv,Train_data.age,'normal','Link','identity','Alpha',0.5,'CV',5,'Lambda',[0.25:0.25:2.5],'MCReps',50);
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
Coef_Model = [FitInfo.Intercept(idxLambdaMinDeviance);ModelCoef_set(:,idxLambdaMinDeviance)];

%% Estimate brain age and PAD 
% Train
Train_data.brainage = glmval(Coef_Model,Train_data.gv,'identity');
Train_data.PAD = Train_data.brainage - Train_data.age;

% Holdset
Holdset_data.brainage = glmval(Coef_Model,Holdset_data.gv,'identity');
Holdset_data.PAD = Holdset_data.brainage - Holdset_data.age;

% Validation
Validation_data_W2.brainage = glmval(Coef_Model,Validation_data_W2.gv,'identity');
Validation_data_W2.PAD = Validation_data_W2.brainage - Validation_data_W2.age;

Validation_data_W3.brainage = glmval(Coef_Model,Validation_data_W3.gv,'identity');
Validation_data_W3.PAD = Validation_data_W2.brainage - Validation_data_W3.age;

% WD
HC_data.brainage = glmval(Coef_Model,HC_data.gv,'identity');
HC_data.PAD = HC_data.brainage - HC_data.age;

CWD_data.brainage = glmval(Coef_Model,CWD_data.gv,'identity');
CWD_data.PAD = CWD_data.brainage - CWD_data.age;

HWD_data.brainage = glmval(Coef_Model,HWD_data.gv,'identity');
HWD_data.PAD = HWD_data.brainage - HWD_data.age;

LWD_data.brainage = glmval(Coef_Model,LWD_data.gv,'identity');
LWD_data.PAD = LWD_data.brainage - LWD_data.age;

NWD_data.brainage = glmval(Coef_Model,NWD_data.gv,'identity');
NWD_data.PAD = NWD_data.brainage - NWD_data.age;


%% Estimate corrected brain age and PAD 
Coef_Correction = regress(Train_data.brainage,[Train_data.age,Train_data.sex,ones(length(Train_data.age),1)]);

% Train
Train_data.PAD_corrected = Train_data.brainage - sum(Coef_Correction'.*[Train_data.age,Train_data.sex,ones(size(Train_data.age,1),1)],2);
Train_data.brainage_corrected = Train_data.PAD_corrected + Train_data.age;

% Holdset
Holdset_data.PAD_corrected = Holdset_data.brainage - sum(Coef_Correction'.*[Holdset_data.age,Holdset_data.sex,ones(size(Holdset_data.age,1),1)],2);
Holdset_data.brainage_corrected = Holdset_data.PAD_corrected + Holdset_data.age;

% Validation
Validation_data_W2.PAD_corrected = Validation_data_W2.brainage - sum(Coef_Correction'.*[Validation_data_W2.age,Validation_data_W2.sex,ones(size(Validation_data_W2.age,1),1)],2);
Validation_data_W2.brainage_corrected = Validation_data_W2.PAD_corrected + Validation_data_W2.age;

Validation_data_W3.PAD_corrected = Validation_data_W3.brainage - sum(Coef_Correction'.*[Validation_data_W3.age,Validation_data_W3.sex,ones(size(Validation_data_W3.age,1),1)],2);
Validation_data_W3.brainage_corrected = Validation_data_W3.PAD_corrected + Validation_data_W3.age;

% WD
HC_data.PAD_corrected = HC_data.brainage - sum(Coef_Correction'.*[HC_data.age,HC_data.sex,ones(size(HC_data.age,1),1)],2);
HC_data.brainage_corrected = HC_data.PAD_corrected + HC_data.age;

CWD_data.PAD_corrected = CWD_data.brainage - sum(Coef_Correction'.*[CWD_data.age,CWD_data.sex,ones(size(CWD_data.age,1),1)],2);
CWD_data.brainage_corrected = CWD_data.PAD_corrected + CWD_data.age;

HWD_data.PAD_corrected = HWD_data.brainage - sum(Coef_Correction'.*[HWD_data.age,HWD_data.sex,ones(size(HWD_data.age,1),1)],2);
HWD_data.brainage_corrected = HWD_data.PAD_corrected + HWD_data.age;

LWD_data.PAD_corrected = LWD_data.brainage - sum(Coef_Correction'.*[LWD_data.age,LWD_data.sex,ones(size(LWD_data.age,1),1)],2);
LWD_data.brainage_corrected = LWD_data.PAD_corrected + LWD_data.age;

NWD_data.PAD_corrected = NWD_data.brainage - sum(Coef_Correction'.*[NWD_data.age,NWD_data.sex,ones(size(NWD_data.age,1),1)],2);
NWD_data.brainage_corrected = NWD_data.PAD_corrected + NWD_data.age;


%% Model Performance
Train_data.MAE = mae(Train_data.age,Train_data.brainage_corrected);
Train_data.R = corr(Train_data.age,Train_data.brainage_corrected);
Holdset_data.MAE = mae(Holdset_data.age,Holdset_data.brainage_corrected);
Holdset_data.R = corr(Holdset_data.age,Holdset_data.brainage_corrected);

%% Save Results and Parameters
Results = struct();

% Traind, Holdset and Validation
Results.Train = Train_data;
Results.Holdset = Holdset_data;
Results.Validation_W2 = Validation_data_W2;
Results.Validation_W3 = Validation_data_W3;

% WD
Results.HC = HC_data;
Results.CWD = CWD_data;
Results.HWD = HWD_data;
Results.LWD = LWD_data;
Results.NWD = NWD_data;

save('BAResults.mat','Results');