function [SParcelExpression_N,roi_locatedin_atls] = SUB_geneExpression(SampleGeneExpression,sample_roi_path,SUB_atls_file)
% this function is to perform correlation analysis between between 
% gene expression and group difference of gradient
% INPUT:
        % SampleGeneExpression: gene expression profile gernerated from AHBA
        % sample_roi_path: sample roi path (nifti file)
        % STR_atls_file: striatum atlas (Left or right putamen or caudate)

disp('   start performing STR_geneExpression')

sample_roi_dir = dir(fullfile(sample_roi_path,'*sphere*.nii.gz'));

% find each voxel location of each sample roi
rois_location = cell(length(sample_roi_dir),1);
rois_num = zeros(length(sample_roi_dir),1);
for i = 1:length(sample_roi_dir)
    roi_nii = load_untouch_nii(fullfile(sample_roi_path,sample_roi_dir(i).name));
    roi_img = roi_nii.img;
    roi_ind = find(roi_img);
    rois_location{i} = roi_ind;
    
    roi_name = split(sample_roi_dir(i).name,{'e','.'});
    rois_num(i) = str2double(roi_name{3});
end

% identify each sample roi corresponding to which part of the subcor_atls
SUB_atls_nii = load_untouch_nii(SUB_atls_file);
SUB_atls_img = SUB_atls_nii.img;
roi_locatedin_atls = [];
for i = 1:length(rois_location)
    subregion_num = SUB_atls_img(rois_location{i}); % each voxel of sample roi located in which part of subcor_atls
    nonzero_ind = find(subregion_num);
    if ~isempty(nonzero_ind)
        subregion_num1 = subregion_num(nonzero_ind); % excluded zero from sample roi voxel
        num_frequency = tabulate(subregion_num1);
        [maxvalue,max_frequency] = max(num_frequency(:,2));
        tmp = [rois_num(i),num_frequency(max_frequency,1)];
        roi_locatedin_atls = [roi_locatedin_atls;tmp];
    else
        continue 
    end
end

% STR gene exression
subcor_atls_num = unique(roi_locatedin_atls(:,2));
SParcelExpression = zeros(length(subcor_atls_num),size(SampleGeneExpression,2));
for i = 1:length(subcor_atls_num)
    roi_location_ind = find(roi_locatedin_atls(:,2) == subcor_atls_num(i));
    roi_location = roi_locatedin_atls(:,1);
    sample_ind  = roi_location(roi_location_ind);
    subregion_sampleGeneExpression = SampleGeneExpression(sample_ind,2:end);
    SParcelExpression(i,1) = subcor_atls_num(i);
    SParcelExpression(i,2:end) = mean(subregion_sampleGeneExpression,1);
end

subreg_idx = unique(SUB_atls_img);
subreg_idx = subreg_idx(2:end);
SParcelExpression_N = zeros(length(subreg_idx),size(SParcelExpression,2));
sge_num = SParcelExpression(:,1);
for i = 1:length(subreg_idx)
    s_ind = find(sge_num == subreg_idx(i));
    SParcelExpression_N(i,1) = subreg_idx(i);
    if ~isempty(s_ind)
        SParcelExpression_N(i,2:end) = SParcelExpression(s_ind,2:end);
    else
        SParcelExpression_N(i,2:end) = NaN;
    end
end