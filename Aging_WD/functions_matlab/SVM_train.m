function [accu_mean,W_mean,cv,dist_total] = SVM_train(X,Y)
% this function is to perform SVM classification analysis
% INPUTs:
%       X: train data
%       Y: labels
% OUTPUTs:
%       mean_accu : mean accuracy across k-fold across validation
%       mean_W    : mean weight of SVM across k-fold across validation
%       cv        : cvpartition object (KFold = 5)
%       dist_total: [nSamples × 1] signed distance of each sample to SVM hyperplane

cv = cvpartition(Y, 'KFold', 5);
acc = zeros(cv.NumTestSets,1);
W_all = zeros(size(X,2), cv.NumTestSets);
dist_total = zeros(size(Y));  % Consistent with the sample order

for i = 1:cv.NumTestSets
    tr = training(cv, i);
    te = test(cv, i);

    % train SVM
    model = fitcsvm(X(tr,:), Y(tr), ...
        'KernelFunction', 'linear', ...
        'Standardize', true, ...
        'BoxConstraint', 1);

    % prediction in test data set
    ypred = predict(model, X(te,:));
    acc(i) = mean(ypred == Y(te));

    % extract linear SVM weight
    W = model.Beta;   % the weight of each feature（direction）
    W_all(:,i) = W;
    b = model.Bias;

    % compute signed distance to hyperplane
    % f(x) = W' * x + b
    dist_total(te) = (X(te,:) * W + b) ./ norm(W);
end
% mean weight
W_mean = mean(W_all, 2);
accu_mean = mean(acc);