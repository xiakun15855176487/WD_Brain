function ind = whereind(ind1,ind2)
% locate ind1 in ind2

ind = zeros(length(ind1),1);
for i = 1:length(ind1)
    tmp_id = find(ind2 == ind1(i));
    ind(i) = tmp_id;
end
