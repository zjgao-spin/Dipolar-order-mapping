% =========================================================================
% Script:       Simulation_study2_SNR.m
% Authors:      Zijian Gao,Qianxue Shan,Weitian Chen
% Email:        zijian.gao@link.cuhk.edu.hk
% Date:         2025-12-19
% Version:      1.0
%
% Copyright (c) 2025  ZZijian Gao, Qianxue Shan, Weitian Chen. All rights reserved.
%
% License:
%   Strictly for academic/research use only. Commercial use, redistribution,
%   or integration into proprietary software is prohibited without written
%   permission. Modification of this header is forbidden.
%
% Description:
%   This script performs a Monte Carlo simulation to analyze the sensitivity 
%   of T1D estimation methods (Analytical vs. Dictionary-based) under 
%   varying Signal-to-Noise Ratios (SNR).
% =========================================================================

clc; clear all; close all;

% Add path to helper functions
addpath('..\..\Function\')

% Load Dictionary Data
% Ensure 'all_jtmt' and 'var_T1d' are in this .mat file
if exist('T1d_dictionary_2025-12-13.mat', 'file')
    load('T1d_dictionary_2025-12-13.mat');
else
    warning('Dictionary file not found. Ensure the path is correct.');
end

%% 1. Define Base Parameters
% =========================================================================
% Basic Tissue Parameters (White Matter Model)
% =========================================================================
T1a_base = 1840e-3;     % T1 of free water pool (s)
T1b_base = 340e-3;      % T1 of bound pool (s)
T2a_base = 69e-3;       % T2 of free water pool (s)
T2b_base = 10e-6;       % T2 of bound pool (s)
T1d_base = 6.2e-3;      % Dipolar relaxation time (s)
MPF_base = 0.15;        % Macromolecular Proton Fraction
R_base = 23;            % Exchange rate (Hz)

% Convert Time constants to Rates (Hz)
R1a_base = 1/T1a_base;
R2a_base = 1/T2a_base;
R1b_base = 1/T1b_base;

% =========================================================================
% General Parameter Settings
% =========================================================================
% RF pulse shape parameters
B1 = 1; 
B0 = 0;
TSL = 80*1e-3;          % Spin-Lock Duration (s)

w1_base = 500 * 2 * pi; % Spin-lock frequency (rad/s)
dw_base = 5000 * 2 * pi;% Off-resonance frequency (rad/s)

% Simulation settings
SNR_all = (20:1:100);  % Range of SNR values to test
N_size = 100000;         % Number of Monte Carlo iterations
method_names = {'Analyt', 'Dic'};

%% 2. Initialize Results Structure
Results = struct();
base_values = [R1a_base, R2a_base, R1b_base, MPF_base, R_base, T2b_base, T1d_base];

%% 3. Main Calculation Loop
for i = 1:length(SNR_all)
    % Current parameter set
    current_params = base_values;
    R1a = current_params(1);
    R2a = current_params(2);
    R1b = current_params(3);
    MPF = current_params(4);
    R   = current_params(5);
    T2b = current_params(6);
    T1d = current_params(7);
    
    SNR = SNR_all(i);
    
    % Initialize temporary storage for noisy Rdosl generation
    Rdosl = zeros(1, N_size);
    
    % ====== Key Modification: Generate Rdosl only once per SNR level ======
    % This ensures both methods are tested against the exact same noise realization
    parfor i_N = 1:N_size
        Rdosl(i_N) = cal_Rdosl_acquired(R1a, R2a, R1b, MPF, R, T2b, T1d, ...
                                        dw_base, w1_base, TSL, B1, B0, SNR);
    end
    
    % Iterate through each estimation method
    for method_idx = 1:length(method_names)
        method_name = method_names{method_idx};
        
        % Initialize storage (only on the first iteration of the SNR loop)
        if i == 1
            Results.(method_name).SNR = SNR_all;
            Results.(method_name).Rdosl = zeros(length(SNR_all), N_size);
            Results.(method_name).T1d_est = zeros(length(SNR_all), N_size);
            Results.(method_name).median_est = zeros(1, length(SNR_all)); % Renamed from 'mean'
            Results.(method_name).iqr = zeros(1, length(SNR_all));         % Renamed from 'std' to match calculation
            Results.(method_name).bias_pct = zeros(1, length(SNR_all));
            Results.(method_name).cv = zeros(1, length(SNR_all));
        end
        
        % Store the generated Rdosl data for record-keeping
        Results.(method_name).Rdosl(i, :) = Rdosl;
        
        % Temporary vector for current SNR iteration results
        T1d_temp_vec = zeros(1, N_size);
        
        % Fit T1d for each realization
        parfor j = 1:N_size
            if method_idx == 1
                % Method 1: Analytical Solution (Non-Linear Least Squares usually)
                [T1d_val, ~] = cal_T1d_analytical(Rdosl(j), B1, B0);
                T1d_temp_vec(j) = T1d_val;
            else
                % Method 2: Dictionary Matching
                % Note: 'all_jtmt' and 'var_T1d' must be loaded from the .mat file
                yy_check = abs(Rdosl(j) - all_jtmt(:));
                num_T1d = find(yy_check == min(yy_check), 1); % Take only the first minimum index
                T1d_val = var_T1d(num_T1d);
                T1d_temp_vec(j) = T1d_val;
            end
        end
        
        % Store all estimates for this SNR step
        Results.(method_name).T1d_est(i, :) = T1d_temp_vec;
        
        % --- Calculate Statistics ---
        
        % Calculate Median
        mu = median(T1d_temp_vec);
        Results.(method_name).median_est(i) = mu; 
        
        % Calculate Interquartile Range (IQR)
        iqr_val = iqr(T1d_temp_vec);
        Results.(method_name).iqr(i) = iqr_val;
        
        % Calculate Relative Bias (%)
        Results.(method_name).bias_pct(i) = abs(mu - T1d_base) / T1d_base * 100;
        
        % Calculate Coefficient of Variation (CV) based on IQR and Median
        Results.(method_name).cv(i) = iqr_val / mu;
    end
    
    fprintf('Completed SNR = %d\n', SNR);
end

fprintf('Calculation complete.\n');


%% 4. Visualization
% Font size settings (Adjust globally here)
fontSize_axis = 15;      % Axis tick font size
fontSize_label = 20;     % Axis label font size
fontSize_title = 16;     % Title font size
fontSize_legend = 12;    % Legend font size
LineW = 2.5;             % Line width

fig = figure('Name', 'Sensitivity Analysis', 'Position', [100, 100, 1800, 400]);

% Define colors
color_NLS = [0, 0.4470, 0.7410];      % Blue (Analytical)
color_Dic = [0.8500, 0.3250, 0.0980]; % Red (Dictionary)

% Plot 1: Median Estimate
subplot(1,4,1); hold on; box on;
plot(SNR_all, Results.Analyt.median_est*1000, '-', 'Color', color_NLS, 'LineWidth', LineW);
plot(SNR_all, Results.Dic.median_est*1000, '-', 'Color', color_Dic, 'LineWidth', LineW);
yline(T1d_base*1000, 'k--', 'LineWidth', LineW);
ylim([5.5 8]);
xlim([min(SNR_all), max(SNR_all)]);
xlabel('SNR', 'FontSize', fontSize_label); 
ylabel('Estimated T_{1D} (ms)', 'FontSize', fontSize_label);
title('Median Estimate', 'FontSize', fontSize_title);
legend('T_{1D}-Analyt', 'T_{1D}-Dic', 'Ground Truth', 'Location', 'best', 'FontSize', fontSize_legend);
set(gca, 'FontSize', fontSize_axis);
grid on

% Plot 2: Interquartile Range 
subplot(1,4,2); hold on; box on;
plot(SNR_all, Results.Analyt.iqr*100, '-', 'Color', color_NLS, 'LineWidth', LineW);
plot(SNR_all, Results.Dic.iqr*100, '-', 'Color', color_Dic, 'LineWidth', LineW);
xlabel('SNR', 'FontSize', fontSize_label); 
ylabel('Interquartile Range (ms)', 'FontSize', fontSize_label);
xlim([min(SNR_all), max(SNR_all)])
title('Interquartile Range', 'FontSize', fontSize_title);
set(gca, 'FontSize', fontSize_axis);
grid on

% Plot 3: Relative Bias
subplot(1,4,3); hold on; box on;
plot(SNR_all, Results.Analyt.bias_pct, '-', 'Color', color_NLS, 'LineWidth', LineW);
plot(SNR_all, Results.Dic.bias_pct, '-', 'Color', color_Dic, 'LineWidth', LineW);
xlabel('SNR', 'FontSize', fontSize_label); 
ylabel('Bias (%)', 'FontSize', fontSize_label);
title('Relative Bias', 'FontSize', fontSize_title);
set(gca, 'FontSize', fontSize_axis);
grid on
xlim([min(SNR_all), max(SNR_all)])
ylim([0 30])

% Plot 4: Coefficient of Variation
subplot(1,4,4); hold on; box on;
plot(SNR_all, Results.Analyt.cv, '-', 'Color', color_NLS, 'LineWidth', LineW);
plot(SNR_all, Results.Dic.cv, '-', 'Color', color_Dic, 'LineWidth', LineW);
xlabel('SNR', 'FontSize', fontSize_label); 
ylabel('CV (IQR / Median)', 'FontSize', fontSize_label);
title('Coefficient of Variation', 'FontSize', fontSize_title);
set(gca, 'FontSize', fontSize_axis);
grid on
xlim([min(SNR_all), max(SNR_all)])

% Save figure (Optional)
% saveas(fig, 'Sensitivity_Analysis.png');
% fprintf('Figure saved as Sensitivity_Analysis.png and .fig\n');