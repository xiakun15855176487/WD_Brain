function tonii = mat_to_nii(insFile,x)
% This function dose transform matrix to nii
% INPUT:
% insFile is the mask you used to do gradient
% x is a matrix, the dimension is the same as insFile
ins_msk_nii = load_untouch_nii(insFile);
ins_msk_nii.hdr.dime.datatype=16;
ins_msk_nii.hdr.dime.bitpix=32;
ins_msk = ins_msk_nii.img;
ind_msk = find(ins_msk); % find the index in of the subcortex in the mask
[r,c,z] = size(ins_msk);
imge = zeros(r,c,z);

for i = 1:length(ind_msk)
    imge(ind_msk(i)) = x(i);
end
ins_msk_nii.img = imge;
tonii = ins_msk_nii;
end

