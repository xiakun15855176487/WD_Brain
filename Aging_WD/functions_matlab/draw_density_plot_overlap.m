function draw_density_plot_overlap(data_oa, data_hc, label_x, label_y, color_oa, color_hc)
% DRAW_DENSITY_PLOT_OVERLAP
% Generate overlapping kernel density plots to compare between the two groups
%   Input:
%       data_oa (vector): 
%       data_hc (vector): 
%       label_x (string): 
%       label_y (string): 
%       color_oa (vector 1x3)
%       color_hc (vector 1x3)
%       title_text (string)
%       effect_size (double)
%       p_value (double)
%
%   Example:
%       % 1. simulated data
%       rng(1); % Set the random seed
%       mu_hc = -1; sigma_hc = 4;
%       mu_oa = 1; sigma_oa = 6;
%       data_hc_sim = normrnd(mu_hc, sigma_hc, 500, 1);
%       data_oa_sim = normrnd(mu_oa, sigma_oa, 500, 1);
%       
%       % 2. Call function
%       draw_density_plot_overlap(data_oa_sim, data_hc_sim, ...
%           'Brain predicted age difference', ' ', ...
%           [0.1 0.4 0.7], [0.6 0.6 0.6], 'a', 0.454, 0.020);


%% 1. Data Preparation and Kernel Density Estimation

x_min = min([data_oa; data_hc]) - 5;
x_max = max([data_oa; data_hc]) + 5;
x_plot = linspace(floor(x_min/5)*5, ceil(x_max/5)*5, 500);

% calculate Kernel Density
[f_oa, xi_oa] = ksdensity(data_oa, x_plot, 'Bandwidth', 2);
[f_hc, xi_hc] = ksdensity(data_hc, x_plot, 'Bandwidth', 2);

zeros_vector = zeros(size(xi_oa));


%% 2. Plot set

figure('Color', 'w', 'Position', [200, 200, 600, 500]); % 设置图窗大小
hold on;

y_lim_max = max(max(f_oa), max(f_hc)) * 1.1; 

fill_x_hc = [xi_hc, fliplr(xi_hc)];
fill_y_hc = [f_hc, fliplr(zeros_vector)];
h_hc = fill(fill_x_hc, fill_y_hc, color_hc, 'EdgeColor', 'none', 'FaceAlpha', 0.8);

fill_x_oa = [xi_oa, fliplr(xi_oa)];
fill_y_oa = [f_oa, fliplr(zeros_vector)];
h_oa = fill(fill_x_oa, fill_y_oa, color_oa, 'EdgeColor', 'none', 'FaceAlpha', 0.6); 

plot(xi_oa, f_oa, 'Color', color_oa*0.8, 'LineWidth', 1.5);
plot(xi_hc, f_hc, 'Color', color_hc*0.7, 'LineWidth', 1.5);
%% 4. Set axis and figure legend

legend([h_oa, h_hc], {'Wilson disease', 'Healthy control'}, ...
  'Location', 'northwest', 'Box', 'off', 'FontSize', 11);

% --- axis ---
xlim([x_plot(1), x_plot(end)]);
ylim([0, y_lim_max]);
xlabel(label_x, 'FontSize', 15);
ylabel(label_y, 'FontSize', 15);

set(gca, 'YTick', []); 
grid off;
box off;
set(gca, 'FontSize', 15, 'TickDir', 'out');

ax = gca;
ax.YAxis.Axle.Visible = 'off';

hold off;

end