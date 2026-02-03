function map_to_mmp_surf(mmp_surf,Values,func_label,output_DIR)
% This function is to map the values to MMP surface
% Input:
%       mmp_surf: MMP surfer file name (.gii)
%       Values: the values that need to map on mmp surf
%       output_DIR: the file output dir path

MMP_surf_label = gifti(mmp_surf);
MMP_surf_data = MMP_surf_label.cdata;
MMP_surf_ind = unique(MMP_surf_data);
MMP_surf_ind = MMP_surf_ind(2:end);

MMP_surf_value_data = double(MMP_surf_data);

for i = 1:length(MMP_surf_ind)
    ROI_surf_ind = find(MMP_surf_data == MMP_surf_ind(i));
    if contains(mmp_surf,'lh')
        MMP_surf_value_data(ROI_surf_ind) = Values(MMP_surf_ind(i));
    else
        MMP_surf_value_data(ROI_surf_ind) = Values(MMP_surf_ind(i)+180);
    end
end
surf_name = split(mmp_surf,'.');
save(gifti(MMP_surf_value_data),fullfile(output_DIR,[surf_name{1},func_label,'.func.gii']))