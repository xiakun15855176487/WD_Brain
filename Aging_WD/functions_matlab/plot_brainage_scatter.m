function plot_brainage_scatter(chronological_age, predicted_age, varargin)
% plot_brainage_scatter - Draw a scatter plot for brain age prediction (including regression line, r-value, MAE).
%
% usage：
%   plot_brainage_scatter(chronological_age, predicted_age)
%   plot_brainage_scatter(chronological_age, predicted_age, 'Title', 'Hold-out set')
%
% Input：
%   chronological_age 
%   predicted_age     
%
% varargin：
%   'Title'           - tile（defalt 'Hold-out set'）
%
% Output：
%   Plot a scatter plot showing/indicating the r-value and MAE.

% ----------------------------
% set parameters
% ----------------------------
p = inputParser;
addParameter(p, 'Title', 'Hold-out set', @ischar);
parse(p, varargin{:});
titleText = p.Results.Title;

% ----------------------------
% check data
% ----------------------------
if length(chronological_age) ~= length(predicted_age)
    error('chronological_age 和 predicted_age must be same');
end

% ----------------------------
% calculate r and MAE
% ----------------------------
r = corr(chronological_age(:), predicted_age(:), 'rows', 'complete');
MAE = mean(abs(predicted_age - chronological_age), 'omitnan');

% ----------------------------
% plot
% ----------------------------
figure('Color','w');
hold on;

% scater
scatter(chronological_age, predicted_age, 30, [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.6);

% regression line
pfit = polyfit(chronological_age, predicted_age, 1);
x_fit = linspace(min(chronological_age), max(chronological_age), 100);
y_fit = polyval(pfit, x_fit);
plot(x_fit, y_fit, 'Color', [0 0.3 0.6], 'LineWidth', 2.5);

% axis
axis square;
box off;
set(gca, 'FontSize', 15, 'LineWidth', 1.2, 'FontName', 'Arial');
xlabel('Chronological age (yr)', 'FontSize', 15);
ylabel('Brain age (yr)', 'FontSize', 15);
grid on;

% title
title(titleText, 'FontSize', 15, 'FontWeight', 'normal');

% text
xText = min(chronological_age) + 0.3*(max(chronological_age)-min(chronological_age));
yText = min(predicted_age) + 0.15*(max(predicted_age)-min(predicted_age));

yRange = max(predicted_age) - min(predicted_age);
lineSpacing = 0.10 * yRange;   
text(xText, yText, sprintf('r = %.3f *', r), 'FontSize', 15, 'Color', 'k', 'FontName', 'Arial');
%text(xText, yText - lineSpacing, sprintf('MAE = %.3f', MAE), 'FontSize', 15, 'Color', 'k', 'FontName', 'Arial');


set(gca, 'XColor', 'k', 'YColor', 'k');
xlim([min(chronological_age)-2, max(chronological_age)+2]);
ylim([min(predicted_age)-2, max(predicted_age)+2]);

end
