%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_visualization.m
%%%
%%% MATLAB script to visualize the results using inhome matlab code
%%%
%%% Original: University of Science and Technology of China, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load mask
% yeo 7 cortical network
yeo_file = fullfile('masks','Yeo_7Network_MNI152_2mm.nii.gz');
yeo_nifti = load_untouch_nii(yeo_file);
yeo_img = yeo_nifti.img;

% subcortex
subcort_file = fullfile('masks','subcort_atlas_2mm.nii.gz');
subcort_nifti = load_untouch_nii(subcort_file);
subcort_img = subcort_nifti.img;

% cerebellum
cerebellum_file = fullfile('masks','CRB_mask_2mm.nii.gz');
cerebellum_nifti = load_untouch_nii(cerebellum_file);
cerebellum_img = cerebellum_nifti.img;

%%  visualization of SVM weight
%
% load SVM weight
W_file = fullfile('Results','SVM','SVM_weights.nii');
W_nifti = load_untouch_nii(W_file);
W_img = W_nifti.img;

W_img_P = W_img; % positive 
W_img_P(W_img_P < 0) = 0;

W_img_N = W_img; % negative
W_img_N(W_img_N > 0) = 0;

vox_pec_P = vox_mpras_perc_net(yeo_img,subcort_img,cerebellum_img,W_img_P);
vox_pec_N = vox_mpras_perc_net(yeo_img,subcort_img,cerebellum_img,W_img_N);
vox_pec_All = vox_mpras_perc_net(yeo_img,subcort_img,cerebellum_img,W_img);
% Visaulization
labels = {'VIS','SMN','DAN','VAN','LIB','FNP','DMN','SUB','CRB'};
group1 = vox_pec_P';
group2 = vox_pec_N';
axes_limits = [zeros(1,length(vox_pec_N));zeros(1,length(vox_pec_N))+20];
radarPlot(group2,group1,labels,axes_limits)

%% visualization of spatial pattern of correlations between gray matter and PAD

% load data
NWD_file = fullfile('Results','GV_PAD_Corr','NWD_sigfdr.nii');
NWD_nifti = load_untouch_nii(NWD_file);
NWD_img = NWD_nifti.img;

LWD_file = fullfile('Results','GV_PAD_Corr','LWD_sigfdr.nii');
LWD_nifti = load_untouch_nii(LWD_file);
LWD_img = LWD_nifti.img;

vox_pec_NWD = vox_mpras_perc_net(yeo_img,subcort_img,cerebellum_img,NWD_img);
vox_pec_LWD = vox_mpras_perc_net(yeo_img,subcort_img,cerebellum_img,LWD_img);

% Visaulization
labels = {'VIS','SMN','DAN','VAN','LIB','FNP','DMN','SUB','CRB'};
group1 = vox_pec_NWD';
group2 = vox_pec_LWD';
axes_limits = [zeros(1,length(vox_pec_NWD));zeros(1,length(vox_pec_NWD))+40];
radarPlot(group2,group1,labels,axes_limits)

%% ksdensity plot -- Sample age
clear, clc;

age = NWD_data.age;
sex = NWD_data.sex;
male_id = find(sex ==1); female_id = find(sex ==2);

male_age   = age(male_id);
female_age = age(female_id);

mirror_density_plot_vertical(male_age, female_age);

%% scatter plot
load BAResults.mat

% Train
chronological_age = Results.Train.age;
predicted_age = Results.Train.brainage_corrected;
plot_brainage_scatter(chronological_age, predicted_age, 'Title', 'Training set');


% Hold-out
chronological_age = Results.Holdset.age;
predicted_age = Results.Holdset.brainage_corrected;
plot_brainage_scatter(chronological_age, predicted_age, 'Title', 'Hold-out set');

% Validation
chronological_age = Results.Validation_W2.brainage_corrected;
predicted_age = Results.Validation_W3.brainage_corrected;
plot_brainage_scatter(chronological_age, predicted_age, 'Title', ' ');

%% density_plot -- Brain predicted age difference

data_hc_pad = Results.HC.PAD_corrected;
data_wd_pad = Results.NWD.PAD_corrected;

custom_label_x = 'Brain predicted age difference (years)';
custom_label_y = 'Count density'; 
custom_color_wd = [0.36 0.64 0.34]; 
custom_color_hc = [0.6 0.6 0.6]; 

draw_density_plot_overlap(data_wd_pad, data_hc_pad, ...
    custom_label_x, custom_label_y, ...
    custom_color_wd, custom_color_hc);

%% density_plot -- SVM accuracy

load(fullfile('Results','SVM','SVM_results.mat'))
accu_perm = SVM_results.accu_perm;
accu = SVM_results.accu;

figure;
plot_density_with_p(accu_perm, accu);

%% gene weight plot

X_weights = subcort_spat_corr.NWD.rho;
custom_color = [0.7333 0.3333 0.4];  % Cort:LWD(0.2 0.7 0.3),NWD(0.345 0.580 0.863) 
                                        % Subc:LWD(0.5333 0.4 0.6392),NWD(0.7333 0.3333 0.4)
custom_x_label = 'Feature Weight (Score)';
custom_y_label = 'Data Index';

draw_s_curve_filled_plot(X_weights, custom_color, custom_x_label, custom_y_label);

