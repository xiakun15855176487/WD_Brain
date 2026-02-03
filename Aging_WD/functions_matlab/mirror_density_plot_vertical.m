function mirror_density_plot_vertical(male_data, female_data, colorMale, colorFemale, yLabelText)
% mirror_density_plot_count
%
% Usage：
%   mirror_density_plot_count(male_data, female_data)
%   mirror_density_plot_count(male_data, female_data, [0.8 0.2 0.2], [0.2 0.4 0.8], 'Count density')
%
% Inputs：
%   male_data    - male data（vector）
%   female_data  - female data（vector）
%   colorMale    - male color (option)
%   colorFemale  - female color (option)
%   yLabelText   - y label (defalt 'Count density')

    if nargin < 3, colorMale = [0.8 0.2 0.2]; end
    if nargin < 4, colorFemale = [0.2 0.4 0.8]; end
    if nargin < 5, yLabelText = 'Count density'; end

    %% Kernel Density estimation
    [f_m, x_m] = ksdensity(male_data, 'Bandwidth', 2);  
    [f_f, x_f] = ksdensity(female_data, 'Bandwidth', 2);
    
    % translate to count density
    f_m = f_m * numel(male_data);
    f_f = f_f * numel(female_data);

    %% Plot
    figure; hold on;
    fill(x_m,  f_m, colorMale, 'FaceAlpha', 0.35, 'EdgeColor', colorMale, 'LineWidth', 1.5);
    fill(x_f, -f_f, colorFemale, 'FaceAlpha', 0.35, 'EdgeColor', colorFemale, 'LineWidth', 1.5);
    grid on

    %% axis and label
    yline(0, 'k', 'LineWidth', 0.8);
    xlabel('Chronological age (yr)');
    ylabel(yLabelText);
    % title('Age Distribution by Sex');

    %% set axis range
    data = [male_data;female_data]; min_v = min(data); max_v = max(data);
    xlim([floor(min_v)-5 ceil(max_v)+5]);

    y_max = max([f_m,f_f]);
    ylim([-ceil(y_max) ceil(y_max)]);
    
    % y_labels = -ceil(y_max):step:ceil(y_max);
    % set(gca, 'YTick', y_labels, 'YTickLabel', abs(y_labels));
    set(gca, 'Box', 'off', 'FontSize', 15);
    set(gcf, 'Color', 'w');

    ax = gca;
    
    y_ticks = ax.YTick;
    set(gca, 'YTick', y_ticks, 'YTickLabel', abs(y_ticks));

    ax.XAxis.Axle.Visible = 'off';
    ax.YAxis.Axle.Visible = 'off';
    % ax.TickLength = [0 0];  
    % ax.Color = 'none';
    set(gca, 'XColorMode','manual', 'YColorMode','manual');

end
