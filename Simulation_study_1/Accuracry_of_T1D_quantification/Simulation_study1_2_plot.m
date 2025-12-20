%% =========================================================================
% 4. Plotting Configuration
% =========================================================================

% --- Selection Parameters (Adjust these values to slice the data) ---
% 1. Reference TSL for frequency sweeps (ms)
TSL_ref_value = 80;

% 2. w1 values to show when plotting against Delta Omega (Hz)
w1_for_dw_plot = [300, 500, 700];

% 3. Delta Omega values to show when plotting against w1 (Hz)
dw_for_w1_plot = [5000, 6000, 7000];

% 4. w1 value to show when plotting against TSL (Hz)
w1_for_tsl_plot = 500;

% 5. Delta Omega values to show when plotting against TSL (Hz)
dw_for_tsl_plot = [5000, 6000, 7000];

% --- Automatic Index Finding ---
fprintf('\n=== Plotting Configuration ===\n');

% Find TSL index
[~, tsl_ref_idx] = min(abs(TSL_vals_ms - TSL_ref_value));

% Find w1 indices (for dw plot)
w1_for_dw_idx = zeros(1, length(w1_for_dw_plot));
for i = 1:length(w1_for_dw_plot)
    [~, w1_for_dw_idx(i)] = min(abs(w1_vals_Hz - w1_for_dw_plot(i)));
end

% Find dw indices (for w1 plot)
dw_for_w1_idx = zeros(1, length(dw_for_w1_plot));
for i = 1:length(dw_for_w1_plot)
    [~, dw_for_w1_idx(i)] = min(abs(dw_vals_Hz - dw_for_w1_plot(i)));
end

% Find w1 index (for TSL plot)
[~, w1_for_tsl_idx] = min(abs(w1_vals_Hz - w1_for_tsl_plot));

% Find dw indices (for TSL plot)
dw_for_tsl_idx = zeros(1, length(dw_for_tsl_plot));
for i = 1:length(dw_for_tsl_plot)
    [~, dw_for_tsl_idx(i)] = min(abs(dw_vals_Hz - dw_for_tsl_plot(i)));
end

%% =========================================================================
% 5. Plotting Execution 
% =========================================================================

% Prepare data for plotting
T1D_est_ms  = T1D_est * 1000;
T1D_true_ms = T1D_true * 1000;

% --- 1. Global Plotting Parameters (参数化设置) ---
% Define Colors and Markers
colors = {[0, 0, 1], [1, 0, 0], [0, 0.7, 0]}; % Blue, Red, Green
markers = {'o', 's', '^'};

% Define Figure Dimensions (统一大小)
fig_width = 600;
fig_height = 500;
screen_size = get(0, 'ScreenSize');
left_pos = (screen_size(3) - fig_width) / 2;
bottom_pos = (screen_size(4) - fig_height) / 2;
fig_pos = [left_pos, bottom_pos, fig_width, fig_height];

% Define Font Sizes (字体大小参数化)
fs_title  = 26;
fs_label  = 30;
fs_legend = 14;
fs_tick   = 16;

% Define Line Properties
lw_plot = 2;      % Line width for data curves
lw_axis = 1;    % Line width for the box/axes
ms_size = 8;      % Marker size

% Set default interpreter to latex for professional look
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

% --- Figure 1: Estimated T1D vs Delta Omega (Fixed TSL, Varying w1) ---
figure(1); clf; 
set(gcf, 'Position', fig_pos, 'Color', 'w'); % Set size and white background
hold on;

for i = 1:length(w1_for_dw_idx)
    y_data = squeeze(T1D_est_ms(tsl_ref_idx, w1_for_dw_idx(i), :));
    plot(dw_vals_Hz, y_data, ['-' markers{i}], ...
        'LineWidth', lw_plot, 'Color', colors{i}, 'MarkerSize', ms_size);
end
plot(dw_vals_Hz, repmat(T1D_true_ms, size(dw_vals_Hz)), '--k', 'LineWidth', lw_plot);

xlabel('Frequency Offset $\Delta\omega/2\pi$ (Hz)', 'FontSize', fs_label);
ylabel('Estimated $T_{1D}$ (ms)', 'FontSize', fs_label);
legend_str = arrayfun(@(x) sprintf('$\\omega_1/2\\pi = %d\\,\\mathrm{Hz}$', x), w1_vals_Hz(w1_for_dw_idx), 'UniformOutput', false);
legend([legend_str, {'Ground Truth'}], 'FontSize', fs_legend, 'Location', 'best');

title(sprintf('Estimated $T_{1D}$ vs $\\Delta\\omega$ (TSL = %d ms)', round(TSL_vals_ms(tsl_ref_idx))), 'FontSize', fs_title);
xlim([min(dw_vals_Hz) max(dw_vals_Hz)]); ylim([0 10]);

% Style adjustments (Grid, Box, Ticks)
grid on; box on;
set(gca, 'FontSize', fs_tick, 'LineWidth', lw_axis); % LineWidth affects the Box thickness
saveas(gcf, 'Fig1_T1D_vs_dw.png');


% --- Figure 2: Estimated T1D vs w1 (Fixed TSL, Varying Delta Omega) ---
figure(2); clf; 
set(gcf, 'Position', fig_pos, 'Color', 'w');
hold on;

for i = 1:length(dw_for_w1_idx)
    y_data = squeeze(T1D_est_ms(tsl_ref_idx, :, dw_for_w1_idx(i)));
    plot(w1_vals_Hz, y_data, ['-' markers{i}], ...
        'LineWidth', lw_plot, 'Color', colors{i}, 'MarkerSize', ms_size);
end
plot(w1_vals_Hz, repmat(T1D_true_ms, size(w1_vals_Hz)), '--k', 'LineWidth', lw_plot);

xlabel('RF Amplitude $\omega_1/2\pi$ (Hz)', 'FontSize', fs_label);
ylabel('Estimated $T_{1D}$ (ms)', 'FontSize', fs_label);
legend_str = arrayfun(@(x) sprintf('$\\Delta\\omega/2\\pi = %d\\,\\mathrm{Hz}$', x), dw_vals_Hz(dw_for_w1_idx), 'UniformOutput', false);
legend([legend_str, {'Ground Truth'}], 'FontSize', fs_legend, 'Location', 'best');

title(sprintf('Estimated $T_{1D}$ vs $\\omega_1$ (TSL = %d ms)', round(TSL_vals_ms(tsl_ref_idx))), 'FontSize', fs_title);
xlim([min(w1_vals_Hz) max(w1_vals_Hz)]); ylim([0 10]);

grid on; box on;
set(gca, 'FontSize', fs_tick, 'LineWidth', lw_axis);
saveas(gcf, 'Fig2_T1D_vs_w1.png');


% --- Figure 3: Relative Error vs Delta Omega (Fixed TSL) ---
figure(3); clf; 
set(gcf, 'Position', fig_pos, 'Color', 'w');
hold on;

for i = 1:length(w1_for_dw_idx)
    y_data = squeeze(relative_error(tsl_ref_idx, w1_for_dw_idx(i), :));
    plot(dw_vals_Hz, y_data, ['-' markers{i}], ...
        'LineWidth', lw_plot, 'Color', colors{i}, 'MarkerSize', ms_size);
end

xlabel('Frequency Offset $\Delta\omega/2\pi$ (Hz)', 'FontSize', fs_label);
ylabel('$T_{1D}$ Relative Error (\%)', 'FontSize', fs_label);
legend_str = arrayfun(@(x) sprintf('$\\omega_1/2\\pi = %d\\,\\mathrm{Hz}$', x), w1_vals_Hz(w1_for_dw_idx), 'UniformOutput', false);
legend(legend_str, 'FontSize', fs_legend, 'Location', 'best');

title(sprintf('Relative Error vs $\\Delta\\omega$ (TSL = %d ms)', round(TSL_vals_ms(tsl_ref_idx))), 'FontSize', fs_title);
xlim([min(dw_vals_Hz) max(dw_vals_Hz)]); ylim([0 30]);

grid on; box on;
set(gca, 'FontSize', fs_tick, 'LineWidth', lw_axis);
saveas(gcf, 'Fig3_Error_vs_dw.png');


% --- Figure 4: Relative Error vs w1 (Fixed TSL) ---
figure(4); clf; 
set(gcf, 'Position', fig_pos, 'Color', 'w');
hold on;

for i = 1:length(dw_for_w1_idx)
    y_data = squeeze(relative_error(tsl_ref_idx, :, dw_for_w1_idx(i)));
    plot(w1_vals_Hz, y_data, ['-' markers{i}], ...
        'LineWidth', lw_plot, 'Color', colors{i}, 'MarkerSize', ms_size);
end

xlabel('RF Amplitude $\omega_1/2\pi$ (Hz)', 'FontSize', fs_label);
ylabel('$T_{1D}$ Relative Error (\%)', 'FontSize', fs_label);
legend_str = arrayfun(@(x) sprintf('$\\Delta\\omega/2\\pi = %d\\,\\mathrm{Hz}$', x), dw_vals_Hz(dw_for_w1_idx), 'UniformOutput', false);
legend(legend_str, 'FontSize', fs_legend, 'Location', 'best');

title(sprintf('Relative Error vs $\\omega_1$ (TSL = %d ms)', round(TSL_vals_ms(tsl_ref_idx))), 'FontSize', fs_title);
xlim([min(w1_vals_Hz) max(w1_vals_Hz)]); ylim([0 30]);

grid on; box on;
set(gca, 'FontSize', fs_tick, 'LineWidth', lw_axis);
saveas(gcf, 'Fig4_Error_vs_w1.png');


% --- Figure 5: Estimated T1D vs TSL (Fixed w1, Varying Delta Omega) ---
figure(5); clf; 
set(gcf, 'Position', fig_pos, 'Color', 'w');
hold on;

for i = 1:length(dw_for_tsl_idx)
    y_data = squeeze(T1D_est_ms(:, w1_for_tsl_idx, dw_for_tsl_idx(i)));
    plot(TSL_vals_ms, y_data, ['-' markers{i}], ...
        'LineWidth', lw_plot, 'Color', colors{i}, 'MarkerSize', ms_size);
end
plot(TSL_vals_ms, repmat(T1D_true_ms, size(TSL_vals_ms)), '--k', 'LineWidth', lw_plot);

xlabel('Time of Spin-Lock (ms)', 'FontSize', fs_label);
ylabel('Estimated $T_{1D}$ (ms)', 'FontSize', fs_label);
legend_str = arrayfun(@(x) sprintf('$\\Delta\\omega/2\\pi = %d\\,\\mathrm{Hz}$', x), dw_vals_Hz(dw_for_tsl_idx), 'UniformOutput', false);
legend([legend_str, {'Ground Truth'}], 'FontSize', fs_legend, 'Location', 'best');

title(sprintf('Estimated $T_{1D}$ vs TSL ($\\omega_1/2\\pi = %d$ Hz)', round(w1_vals_Hz(w1_for_tsl_idx))), 'FontSize', fs_title);
xlim([50 120]); ylim([0 10]);

grid on; box on;
set(gca, 'FontSize', fs_tick, 'LineWidth', lw_axis);
saveas(gcf, 'Fig5_T1D_vs_TSL.png');


% --- Figure 6: Relative Error vs TSL (Fixed w1) ---
figure(6); clf; 
set(gcf, 'Position', fig_pos, 'Color', 'w');
hold on;

for i = 1:length(dw_for_tsl_idx)
    y_data = squeeze(relative_error(:, w1_for_tsl_idx, dw_for_tsl_idx(i)));
    plot(TSL_vals_ms, y_data, ['-' markers{i}], ...
        'LineWidth', lw_plot, 'Color', colors{i}, 'MarkerSize', ms_size);
end

xlabel('Time of Spin-Lock (ms)', 'FontSize', fs_label);
ylabel('$T_{1D}$ Relative Error (\%)', 'FontSize', fs_label);
legend_str = arrayfun(@(x) sprintf('$\\Delta\\omega/2\\pi = %d\\,\\mathrm{Hz}$', x), dw_vals_Hz(dw_for_tsl_idx), 'UniformOutput', false);
legend(legend_str, 'FontSize', fs_legend, 'Location', 'best');

title(sprintf('Relative Error vs TSL ($\\omega_1/2\\pi = %d$ Hz)', round(w1_vals_Hz(w1_for_tsl_idx))), 'FontSize', fs_title);
xlim([50 120]); ylim([0 30]);

grid on; box on;
set(gca, 'FontSize', fs_tick, 'LineWidth', lw_axis);
saveas(gcf, 'Fig6_Error_vs_TSL.png');

fprintf('\nAll plots generated successfully with consistent styling.\n');