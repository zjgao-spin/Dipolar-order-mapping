% -------------------------------------------------------------------------
% Script:       RATIO_dosl_Simulation_Validation
% Authors:      Zijian Gao, Qianxue Shan, Ziqiang Yu, Weitian Chen
% Email:        zijian.gao@link.cuhk.edu.hk
% Date:         2025-12-19
% Version:      1.0
%
% Copyright (c) 2025  Zijian Gao, Qianxue Shan, Weitian Chen. All rights reserved.
%
% License:
%   Strictly for academic/research use only. Commercial use, redistribution,
%   or integration into proprietary software is prohibited without written
%   permission. Modification of this header is forbidden.
%
% Description:
%   Validates the analytical approximation of the dipoalr specific rate RATIO_dosl
%   against the exact numerical solution using White Matter 
%   tissue parameters.
%
% Outputs:
%   - Fig 1: RATIO_dosl sensitivity to Frequency Offset (dw).
%   - Fig 2: RATIO_dosl sensitivity to RF Amplitude (w1).
%   - Fig 3: Relative Error Heatmap (Exact vs. Approx).
% -------------------------------------------------------------------------

clc; clear all; close all;

addpath('..\..\Function\')

%% 1. Define Tissue Parameters (White Matter Model)
% Physiological parameters for the binary spin-bath model (Pool A: Free, Pool B: Bound).
T1a = 1840e-3;      % T1 of free pool (s)
T1b = 340e-3;       % T1 of bound pool (s)
T2a = 69e-3;        % T2 of free pool (s)
T2b = 10e-6;        % T2 of bound pool (s)

R1a = 1/T1a;
R2a = 1/T2a;
R1b = 1/T1b;

MPF = 0.15;         % Macromolecular Proton Fraction
R   = 23;           % Exchange rate constant (1/s)

%% 2. Define Simulation Space
% Range of parameters for sensitivity analysis.
T1d_plot_list = [2, 4, 6, 8] * 1e-3; % Dipolar relaxation times to compare (s)

dw_all = linspace(2000, 7000, 40);   % Frequency offset range (Hz)
w1_all = linspace(200, 800, 40);     % RF amplitude range (Hz)

n_dw = length(dw_all);
n_w1 = length(w1_all);
n_t1d_plot = length(T1d_plot_list);

% RF pulse shape parameters
B1 = 1; 
B0 = 0;

%% 3. Numerical Simulations

% --- Experiment A: Frequency Offset Dependence ---
% Fix RF amplitude (w1), vary offset (dw).
w1_fixed_val = 500; 
w1_fixed = w1_fixed_val * 2 * pi; % Convert to rad/s

RATIO_dosl_exact_dw  = zeros(n_t1d_plot, n_dw);
RATIO_dosl_approx_dw = zeros(n_t1d_plot, n_dw);

fprintf('Simulating Fig 1: Fixed w1 = %d Hz, Varying dw...\n', w1_fixed_val);
for k = 1:n_t1d_plot
    curr_T1d = T1d_plot_list(k);
    for i_dw = 1:n_dw
        dw = dw_all(i_dw) * 2 * pi;
        % Compare exact matrix solution vs. analytical approximation
        RATIO_dosl_exact_dw(k, i_dw)  = cal_RATIO_dosl_exact(R1a, R2a, R1b, MPF, R, T2b, curr_T1d, dw,  w1_fixed);
        RATIO_dosl_approx_dw(k, i_dw) = cal_RATIO_dosl_approx(R, T2b, curr_T1d, dw,  w1_fixed, B1, B0);
    end
end

% --- Experiment B: RF Amplitude Dependence ---
% Fix offset (dw), vary RF amplitude (w1).
dw_fixed_val = 5000; 
dw_fixed = dw_fixed_val * 2 * pi; % Convert to rad/s

RATIO_dosl_exact_w1  = zeros(n_t1d_plot, n_w1);
RATIO_dosl_approx_w1 = zeros(n_t1d_plot, n_w1);

fprintf('Simulating Fig 2: Fixed dw = %d Hz, Varying w1...\n', dw_fixed_val);
for k = 1:n_t1d_plot
    curr_T1d = T1d_plot_list(k);
    for i_w1 = 1:n_w1
        w1 = w1_all(i_w1) * 2 * pi;
        RATIO_dosl_exact_w1(k, i_w1)  = cal_RATIO_dosl_exact(R1a, R2a, R1b, MPF, R, T2b, curr_T1d, dw_fixed,  w1);
        RATIO_dosl_approx_w1(k, i_w1) = cal_RATIO_dosl_approx(R, T2b, curr_T1d, dw_fixed,  w1, B1, B0);
    end
end

% --- Experiment C: Error Mapping ---
% Compute relative error across the full 2D parameter space (dw, w1).
T1d_map = 6e-3; % Fixed T1d for error map
RATIO_dosl_approx_map = zeros(n_dw, n_w1);
RATIO_dosl_exact_map  = zeros(n_dw, n_w1);

fprintf('Simulating Fig 3: Computing Error Map...\n');
for i_dw = 1:n_dw
    dw = dw_all(i_dw) * 2 * pi;
    for i_w1 = 1:n_w1
        w1 = w1_all(i_w1) * 2 * pi;
        RATIO_dosl_exact_map(i_dw, i_w1)  = cal_RATIO_dosl_exact(R1a, R2a, R1b, MPF, R, T2b,T1d_map, dw,  w1);
        RATIO_dosl_approx_map(i_dw, i_w1) = cal_RATIO_dosl_approx(R, T2b,T1d_map, dw,  w1, B1, B0);
    end
end


%% 4. Visualization
fs_label = 18;  % x, y 轴标签字体大小
fs_title = 18;  % 标题字体大小
fs_tick  = 18;  % 坐标轴刻度字体大小
fs_leg   = 12;  % 图例字体大小
% General Plotting Settings
set(0, 'DefaultAxesFontSize', fs_tick);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultLineLineWidth', 1.5);
colors = lines(n_t1d_plot); % Use distinct colors

% -------------------------------------------------------------------------
% Figure 1: Varying Offset (dw)
% -------------------------------------------------------------------------
figure('Color','w', 'Position', [100, 100, 600, 500]);
hold on; box on; grid on;

legend_str = cell(1, n_t1d_plot);
p_handles = zeros(1, n_t1d_plot);

% Plot lines
for k = 1:n_t1d_plot
    % Exact Solution (Solid Line)
    p_handles(k) = plot(dw_all/1000, RATIO_dosl_exact_dw(k,:), '-', ...
        'Color', colors(k,:), 'LineWidth', 2);
    
    % Approx Solution (Markers, subsampled for clarity)
    % Sample every 3rd point to avoid overcrowding
    idx_sample = 1:3:length(dw_all); 
    plot(dw_all(idx_sample)/1000, RATIO_dosl_approx_dw(k,idx_sample), 'o', ...
        'Color', colors(k,:), 'MarkerFaceColor', 'w', 'MarkerSize', 6);
    
    legend_str{k} = ['$T_{1D}=' num2str(T1d_plot_list(k)*1000) '\,\mathrm{ms}$'];
end

% Labels and Title
% [FIXED]: Added closing $ and cleaned up syntax
xlabel('Frequency Offset $\Delta\omega^{d(1)}/2\pi$ (kHz)', 'Interpreter', 'latex', 'FontSize', fs_label);
ylabel('$RATIO_{dosl}$', 'Interpreter', 'latex', 'FontSize', fs_label);
title({['$RATIO_{dosl}$ vs. Frequency Offset'], ...
       ['($\omega_1^{d(1)}/2\pi = ' num2str(w1_fixed_val) '\,\mathrm{Hz}$)']}, ...
       'Interpreter', 'latex', 'FontSize', fs_title);

% --- Dual Legend Construction ---
h_line = plot(nan, nan, 'k-', 'LineWidth', 2);
h_circ = plot(nan, nan, 'ko', 'MarkerFaceColor', 'w');

leg1 = legend(p_handles, legend_str, 'Interpreter', 'latex', 'Location', 'northwest');
leg1_pos = leg1.Position;

ax2 = axes('Position', get(gca, 'Position'), 'Visible', 'off');
leg2 = legend(ax2, [h_line, h_circ], {'Exact', 'Approx'}, ...
              'Interpreter', 'latex');
leg2.Position = [leg1_pos(1) + leg1_pos(3) + 0.02, ... 
                 leg1_pos(2) + leg1_pos(4) - leg2.Position(4), ... 
                 leg2.Position(3), ...
                 leg2.Position(4)];

% -------------------------------------------------------------------------
% Figure 2: Varying RF Power (w1)
% -------------------------------------------------------------------------
figure('Color','w', 'Position', [750, 100, 600, 500]);
hold on; box on; grid on;

p_handles_2 = zeros(1, n_t1d_plot);

for k = 1:n_t1d_plot
    % Exact Solution
    p_handles_2(k) = plot(w1_all, RATIO_dosl_exact_w1(k,:), '-', ...
        'Color', colors(k,:), 'LineWidth', 2);
    
    % Approx Solution (Subsampled)
    idx_sample = 1:3:length(w1_all);
    plot(w1_all(idx_sample), RATIO_dosl_approx_w1(k,idx_sample), 'o', ...
        'Color', colors(k,:), 'MarkerFaceColor', 'w', 'MarkerSize', 6);
end

% Labels and Title
% [FIXED]: Removed erroneous %d and added $ signs
xlabel('RF Amplitude $\omega_1^{d(1)}/2\pi$ (Hz)', 'Interpreter', 'latex', 'FontSize', fs_label);
ylabel('$RATIO_{dosl}$', 'Interpreter', 'latex', 'FontSize', fs_label);
title({['$RATIO_{dosl}$ vs. RF Amplitude'], ...
       ['($\Delta\omega^{d(1)}/2\pi = ' num2str(dw_fixed_val/1000) '\,\mathrm{kHz}$)']}, ...
       'Interpreter', 'latex', 'FontSize', fs_title);

% --- Dual Legend Construction ---
h_line = plot(nan, nan, 'k-', 'LineWidth', 2);
h_circ = plot(nan, nan, 'ko', 'MarkerFaceColor', 'w');

leg1 = legend(p_handles_2, legend_str, 'Interpreter', 'latex', 'Location', 'northwest');
ax2 = axes('Position', get(gca, 'Position'), 'Visible', 'off');
leg2 = legend(ax2, [h_line, h_circ], {'Exact', 'Approx'}, 'Interpreter', 'latex');
leg2.Position = [leg1_pos(1) + leg1_pos(3) + 0.02, ... 
                 leg1_pos(2) + leg1_pos(4) - leg2.Position(4), ... 
                 leg2.Position(3), ...
                 leg2.Position(4)];

% -------------------------------------------------------------------------
% Figure 3: Error Map with Contours
% -------------------------------------------------------------------------
RelError = abs(RATIO_dosl_approx_map - RATIO_dosl_exact_map) ./ RATIO_dosl_exact_map * 100;
RelError(RelError > 20) = 20; 

figure('Color','w', 'Position', [1400, 100, 650, 500]);
[X, Y] = meshgrid(w1_all, dw_all/1000);

% Plot Heatmap
h = pcolor(X, Y, RelError);
set(h, 'EdgeColor', 'none'); 
shading interp;

% Colorbar settings
c = colorbar;
c.Label.String = 'Relative Error (%)';
c.Label.FontSize = 12;
caxis([0 5]); 

hold on;

% Add Contour Line (1%)
[C1, h1] = contour(X, Y, RelError, [1 1], 'w--', 'LineWidth', 2);

% Label the Contour
h_labels = clabel(C1, h1, 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', ...
    'LabelSpacing', 300);
for i = 1:length(h_labels)
    current_str = h_labels(i).String;
    % Ensure contour label has % sign
    if ~contains(current_str, '%')
        h_labels(i).String = [current_str '%'];
    end
end

% Labels and Title
% [FIXED]: Standardized double backslashes to single for consistency (optional but cleaner)
xlabel('RF Amplitude $\omega_1^{d(1)}/2\pi$ (Hz)', 'Interpreter', 'latex', 'FontSize', fs_label);
ylabel('Frequency Offset $\Delta\omega^{d(1)}/2\pi$ (kHz)', 'Interpreter', 'latex', 'FontSize', fs_label);

title(['Relative Error Map ($T_{1D} = ' num2str(T1d_map*1000) '\,\mathrm{ms}$)'], 'Interpreter', 'latex', 'FontSize', fs_title);
%subtitle('White dashed line indicates 1% error boundary','FontSize', 8);
axis tight;

%%
% -------------------------------------------------------------------------
% General Plotting & Font Settings
% -------------------------------------------------------------------------
% 定义统一的字体大小变量，方便全局调整
fs_label = 16;  % x, y 轴标签字体大小
fs_title = 16;  % 标题字体大小
fs_tick  = 14;  % 坐标轴刻度字体大小
fs_leg   = 12;  % 图例字体大小

set(0, 'DefaultAxesFontSize', fs_tick);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultTextInterpreter', 'latex'); % 全局开启 Latex
set(0, 'DefaultLegendInterpreter', 'latex');

colors = lines(n_t1d_plot); % Use distinct colors

% 创建一个宽图 (1行3列)
figure('Color','w', 'Position', [100, 100, 1800, 450]);

% -------------------------------------------------------------------------
% Subplot 1: Varying Offset (dw)
% -------------------------------------------------------------------------
subplot(1, 3, 1);
hold on; box on; grid on;

legend_str = cell(1, n_t1d_plot);
p_handles = zeros(1, n_t1d_plot);

% Plot lines
for k = 1:n_t1d_plot
    % Exact Solution (Solid Line)
    p_handles(k) = plot(dw_all/1000, RATIO_dosl_exact_dw(k,:), '-', ...
        'Color', colors(k,:), 'LineWidth', 2);
    
    % Approx Solution (Markers, subsampled)
    idx_sample = 1:5:length(dw_all); % 稍微稀疏一点，避免太拥挤
    plot(dw_all(idx_sample)/1000, RATIO_dosl_approx_dw(k,idx_sample), 'o', ...
        'Color', colors(k,:), 'MarkerFaceColor', 'w', 'MarkerSize', 6);
    
    legend_str{k} = ['$T_{1D}=' num2str(T1d_plot_list(k)*1000) '\,\mathrm{ms}$'];
end

% Labels and Title
xlabel('Frequency Offset $\Delta\omega/2\pi$ (kHz)', 'FontSize', fs_label);
ylabel('$RATIO_{dosl}$ (Hz)', 'FontSize', fs_label);
title({['$RATIO_{dosl}$ vs. Offset'], ...
       ['($\omega_1/2\pi = ' num2str(w1_fixed_val) '\,\mathrm{Hz}$)']}, ...
       'FontSize', fs_title);

% --- Combined Legend Construction ---
% 创建虚拟句柄用于显示 Exact 和 Approx 的样式
h_exact = plot(nan, nan, 'k-', 'LineWidth', 2);
h_approx = plot(nan, nan, 'ko', 'MarkerFaceColor', 'w');

% 合并图例：先显示颜色代表的 T1D，再显示线型代表的方法
legend([p_handles, h_exact, h_approx], [legend_str, {'Exact', 'Approx'}], ...
       'Location', 'northwest', 'FontSize', fs_leg);

% -------------------------------------------------------------------------
% Subplot 2: Varying RF Power (w1)
% -------------------------------------------------------------------------
subplot(1, 3, 2);
hold on; box on; grid on;

p_handles_2 = zeros(1, n_t1d_plot);

for k = 1:n_t1d_plot
    % Exact Solution
    p_handles_2(k) = plot(w1_all, RATIO_dosl_exact_w1(k,:), '-', ...
        'Color', colors(k,:), 'LineWidth', 2);
    
    % Approx Solution
    idx_sample = 1:5:length(w1_all);
    plot(w1_all(idx_sample), RATIO_dosl_approx_w1(k,idx_sample), 'o', ...
        'Color', colors(k,:), 'MarkerFaceColor', 'w', 'MarkerSize', 6);
end

% Labels and Title
xlabel('RF Amplitude $\omega_1/2\pi$ (Hz)', 'FontSize', fs_label);
ylabel('$RATIO_{dosl}$ (Hz)', 'FontSize', fs_label);
title({['$RATIO_{dosl}$ vs. RF Amplitude'], ...
       ['($\Delta\omega/2\pi = ' num2str(dw_fixed_val/1000) '\,\mathrm{kHz}$)']}, ...
       'FontSize', fs_title);

% --- Combined Legend Construction ---
% 复用上面的逻辑
h_exact = plot(nan, nan, 'k-', 'LineWidth', 2);
h_approx = plot(nan, nan, 'ko', 'MarkerFaceColor', 'w');
legend([p_handles_2, h_exact, h_approx], [legend_str, {'Exact', 'Approx'}], ...
       'Location', 'northwest', 'FontSize', fs_leg);

% -------------------------------------------------------------------------
% Subplot 3: Error Map with Contours
% -------------------------------------------------------------------------
subplot(1, 3, 3);

RelError = abs(RATIO_dosl_approx_map - RATIO_dosl_exact_map) ./ RATIO_dosl_exact_map * 100;
RelError(RelError > 20) = 20; 

[X, Y] = meshgrid(w1_all, dw_all/1000);

% Plot Heatmap
h = pcolor(X, Y, RelError);
set(h, 'EdgeColor', 'none'); 
shading interp;

% Colorbar settings
c = colorbar;
c.Label.String = 'Relative Error (%)';
c.Label.FontSize = fs_label;
c.Label.Interpreter = 'latex';
caxis([0 5]); 

hold on;

% Add Contour Line (1%)
[C1, h1] = contour(X, Y, RelError, [1 1], 'w--', 'LineWidth', 2);

% Label the Contour
h_labels = clabel(C1, h1, 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', ...
    'LabelSpacing', 300);

% Fix contour labels to include % sign
for i = 1:length(h_labels)
    current_str = h_labels(i).String;
    if ~contains(current_str, '%')
        h_labels(i).String = [current_str '\%']; % Latex 需要转义 %
    end
end

% Labels and Title
xlabel('RF Amplitude $\omega_1/2\pi$ (Hz)', 'FontSize', fs_label);
ylabel('Offset $\Delta\omega/2\pi$ (kHz)', 'FontSize', fs_label);
title(['Relative Error Map ($T_{1D} = ' num2str(T1d_map*1000) '\,\mathrm{ms}$)'], ...
      'FontSize', fs_title);

axis tight;
box on;

% 调整整体布局，防止标签被切掉
set(gcf, 'PaperPositionMode', 'auto');