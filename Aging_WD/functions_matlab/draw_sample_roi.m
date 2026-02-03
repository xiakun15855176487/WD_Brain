function draw_sample_roi(sample_id,sample_coordinates,output_dir)
% this function is to draw roi based on sample gene coordinate
% sample_id: the index of sample
% sample_coordinates: sample gene MNI coordinate
% output_dir: where the sample roi to save

disp('   start drawing sample gene roi')

% get MNI template
copyfile(fullfile('masks','MNI152_T1_2mm_brain.nii.gz'),output_dir)

% extract sample coordinate of the sample_id
voxel_cord = [];
Coord_IND = [];
for i = 1:length(sample_id)
    coord_ind  = find(sample_coordinates(:,1)==sample_id(i));
    voxel_cord = [voxel_cord;sample_coordinates(coord_ind,2:end)];
    Coord_IND = [Coord_IND;coord_ind];
end

% transform coordinate from MNI space to voxel space
voxel_cord = mni2cor(voxel_cord);

cd(output_dir)
for i = 1:size(voxel_cord,1)
    system(['fslmaths MNI152_T1_2mm_brain.nii.gz -mul 0 -add 1 -roi ',num2str(voxel_cord(i,1)),' 1 ',...
        num2str(voxel_cord(i,2)),' 1 ',num2str(voxel_cord(i,3)),' 1 0 1 roi',num2str(Coord_IND(i)) ' -odt float']);
    system(['fslmaths roi',num2str(Coord_IND(i)),' -kernel sphere 3 -fmean roi_sphere',num2str(Coord_IND(i)),' -odt float']);
end