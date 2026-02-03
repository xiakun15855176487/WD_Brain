function plot_density_with_p(data, threshold)
% 绘制带填充的 Nature 风格核密度图，并标注竖线与P值
% 修复了 fill 的向量长度和透明度问题
%
% 用法:
%   plot_density_with_p_fill(data, threshold, pval)

    % ---------- 样式设置 ----------
    color_line   = [0.27 0.49 0.74];  % 柔和蓝 (line & fill color)
    face_alpha   = 0.25;              % 填充透明度
    color_thresh = [0.8 0.2 0.2];     % 红线
    lw = 1.8;

    % ---------- 计算核密度 ----------
    [f, xi] = ksdensity(data);

    % 强制为行向量以避免尺寸/方向不匹配
    xi = xi(:)';   % 1 x N
    f  = f(:)';    % 1 x N

    % ---------- 绘制填充区域 ----------
    x_fill = [xi, fliplr(xi)];          % 1 x (2N)
    y_fill = [f, zeros(1, numel(f))];   % 1 x (2N)
    hfill = fill(x_fill, y_fill, color_line, 'EdgeColor', 'none');
    set(hfill, 'FaceAlpha', face_alpha);

    hold on;

    % ---------- 绘制密度曲线 ----------
    plot(xi, f, 'Color', color_line, 'LineWidth', lw);

    % ---------- 调整 ylim 后绘制竖线 ----------
    % 确保 ylim 覆盖到填充区域（fill 可能已改变 ylim）
    yl = ylim;
    % 如果上界太小（例如0），适当扩展一点
    if yl(2) <= 0
        yl(2) = max(f) * 1.1;
        ylim(yl);
    end
    plot([threshold threshold], yl, 'Color', color_thresh, 'LineWidth', 1.5);

    % ---------- 标注 P 值 ----------
    % text_x = threshold + (max(xi)-min(xi))*0.02; % 稍微偏右一点
    % text_y = yl(2) * 0.85;
    % text(text_x, text_y, sprintf('P = %.3f', pval), ...
    %     'FontAngle', 'italic', 'FontSize', 11, 'Color', 'k');

    % ---------- 坐标轴与标签 ----------
    xlabel('Accuracy', 'FontSize', 15);
    ylabel('Density', 'FontSize', 15);
    set(gca, 'Box', 'off', 'LineWidth', 1, ...
        'FontSize', 15, 'TickDir', 'out');
    xlim([min(xi) max(xi)]);

    % 去除顶部和右侧边框（更 Nature 风）
    % ax = gca;
    % ax.XRuler.Axle.Visible = 'off';
    % ax.YRuler.Axle.Visible = 'off';

    hold off;
end
