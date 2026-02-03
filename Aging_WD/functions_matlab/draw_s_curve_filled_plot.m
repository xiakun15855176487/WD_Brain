function draw_s_curve_filled_plot(X_weights, curve_color, x_label_text, y_label_text)
% DRAW_S_CURVE_FROM_CSV
%
%   Inputs:
%       filename (string): CSV file name。
%       curve_color (vector 1x3): RGB color，such as [0.4 0.6 0.9]。
%       x_label_text (string): X label
%       y_label_text (string): Y label
%
%   Example:
%       draw_s_curve_from_csv('data.xlsx - Sheet1.csv', [0.4 0.6 0.9], 'Score Magnitude', 'Data Index');

if nargin < 4
    y_label_text = ' '; 
end
if nargin < 3
    x_label_text = 'Weight magnitude';
end
if nargin < 2
    curve_color = [0.1 0.4 0.7]; % default color
end

N_data_points = length(X_weights);
Y_index = 1:N_data_points;

% rank data
[X_weights_sorted, ~] = sort(X_weights, 'ascend'); 


%% 2.plot set
figure('Color', 'w', 'Position', [200, 200, 400, 600]); 
hold on;

x_zero_line = zeros(size(X_weights_sorted));

fill_x = [X_weights_sorted', fliplr(x_zero_line')];
fill_y = [Y_index, fliplr(Y_index)];

h_fill = fill(fill_x, fill_y, curve_color, 'EdgeColor', 'none', 'FaceAlpha', 0.8);

plot(x_zero_line, Y_index, 'k--', 'LineWidth', 0.5); 

plot(X_weights_sorted, Y_index, 'Color', curve_color * 0.7, 'LineWidth', 1.5);


%% 3. Set the axes and styles

% --- X ---
x_max_abs = max(abs([min(X_weights_sorted), max(X_weights_sorted)]));
x_lim_max = ceil(x_max_abs * 1 / 5);
xlim([-x_lim_max x_lim_max]); 
xlabel(x_label_text, 'FontSize', 14);

x_ticks = unique([-fliplr(0:5:x_lim_max), 0:5:x_lim_max]);
xticks(x_ticks);


% --- Y ---
ylim([0 N_data_points + 1]); 
ylabel(y_label_text, 'FontSize', 14);

y_step = max(1, floor(N_data_points / 20)); 
yticks(1:y_step:N_data_points); 
yticklabels({});
set(gca, 'YAxisLocation', 'left'); 


% --- style ---
ax = gca;
ax.Box = 'off';
ax.TickDir = 'out';
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
ax.FontSize = 12;
ax.LineWidth = 1;

grid off; 
hold off;

end