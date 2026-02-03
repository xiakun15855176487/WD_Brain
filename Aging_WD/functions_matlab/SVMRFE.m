function features_idx = SVMRFE(X,Y)
% this function is to selecte features using SVMRFE method
% INPUTs:
%       X: feature matrix
%       Y: labels
% OUTPUTs:
%       features_idx: the selected feature's index

nFeatures = size(X,2);
remove_ratio = 0.01;
num_features_to_select = round(nFeatures * 0.2); % select top 20% features
features_idx = 1:nFeatures;

% Initialize parameters
its = 1;
fprintf('  Starting SVM-RFE ...\n');

while length(features_idx) > num_features_to_select
    % --- 1. 训练线性 SVM ---
    svmModel = fitcsvm(X(:,features_idx), Y, ...
        'KernelFunction', 'linear', 'Standardize', true, ...
        'Solver', 'ISDA', 'Verbose', 0);

    % --- 2. 计算权重平方 ---
    w = svmModel.Beta;
    w2 = w.^2;

    % --- 3. 确定要去除的特征 ---
    n_remove = max(1, round(remove_ratio * length(features_idx)));
    [~, sort_idx] = sort(w2, 'ascend');
    remove_idx = sort_idx(1:n_remove);
    features_idx(remove_idx) = [];

    fprintf('Round %d: 剩余特征数 = %d\n', its, length(features_idx));
    its = its + 1;
end
