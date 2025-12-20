% -------------------------------------------------------------------------
% Script:       Rdosl_Sensitivity_Analysis.m
% Authors:      Zijian Gao; Qianxue Shan, Weitian Chen
% Email:        zijian.gao@link.cuhk.edu.hk
% Date:         2025-12-19
% Version:      1.0
%
% Copyright (c) 2025  Zijian Gao; Qianxue Shan,Weitian Chen. All rights reserved.
%
% License:
%   Strictly for academic/research use only. Commercial use, redistribution,
%   or integration into proprietary software is prohibited without written
%   permission. Modification of this header is forbidden.
%
% Description:
%   Analyzes the sensitivity of the dipolar order saturation rate (Rdosl)
%   to various tissue parameters (T1D, R1a, R2a, R1b, MPF, R, T2b).
%   Compares the exact numerical solution against the analytical approximation
%   under different frequency offsets.
% -------------------------------------------------------------------------

clc
clear all
close all

addpath('..\..\Function\')
% =========================================================================
% Basic Tissue Parameters (White Matter Model)
% =========================================================================
T1a_base = 1840e-3;
T1b_base = 340e-3;
T2a_base = 69e-3;
T2b_base = 10e-6;
T1d_base = 6.2e-3;
MPF_base = 0.15;
R_base = 23;

% Convert to Rates (Hz)
R1a_base = 1/T1a_base;
R2a_base = 1/T2a_base;
R1b_base = 1/T1b_base;

% =========================================================================
% General Parameter Settings
% =========================================================================
% RF pulse shape parameters
B1 = 1; 
B0 = 0;

% Frequency Offsets (Hz -> rad/s)
dw_vals = [5000, 6000, 7000];
dw_all = dw_vals * 2 * pi; 

% Fixed RF Field Amplitude (Hz -> rad/s)
w1_val = 500;
w1 = w1_val * 2 * pi;


% =========================================================================
% Part 1: Rdosl vs. T1D Analysis
% =========================================================================
fprintf('========================================\n');
fprintf('Part 1: Calculating Rdosl vs T1D\n');
fprintf('========================================\n');

% T1D Range: 1ms to 10ms
T1d_all = (1:1:10) * 1e-3; 
n_t1d = length(T1d_all);
n_dw = length(dw_all);

% Pre-allocate result matrices
Rdosl_exact_t1d = zeros(n_dw, n_t1d);
Rdosl_approx_t1d = zeros(n_dw, n_t1d);

fprintf('Fixed w1=%d Hz, calculating curves for different dw...\n', w1_val);

for i = 1:n_dw
    curr_dw = dw_all(i);
    fprintf('  Processing dw = %d Hz...\n', dw_vals(i));
    
    for j = 1:n_t1d
        curr_T1d = T1d_all(j);
        
        Rdosl_exact_t1d(i, j) = cal_Rdosl_exact(R1a_base, R2a_base, R1b_base, ...
            MPF_base, R_base, T2b_base,curr_T1d,curr_dw, w1);
        Rdosl_approx_t1d(i, j) = cal_Rdosl_approx(R_base, T2b_base, curr_T1d, curr_dw,  w1, B1, B0);
    end
end

% =========================================================================
% Part 2: Sensitivity Analysis
% =========================================================================
fprintf('\n========================================\n');
fprintf('Part 2: Sensitivity Analysis\n');
fprintf('========================================\n');

steps = 25; % Number of sample points

% Define parameter variation ranges
R1a_range = linspace(0.38, 1.5, steps);
R2a_range = linspace(10, 80, steps);
R1b_range = linspace(2, 5, steps);
MPF_range = linspace(0.1, 0.2, steps);
R_range   = linspace(18, 25, steps);
T2b_range = linspace(9, 11, steps) * 1e-6;

% Define parameter structure for iteration
param_struct = struct();
param_struct(1).name = 'R1a'; param_struct(1).range = R1a_range; param_struct(1).label = '$R_{1a}$ (Hz)';
param_struct(2).name = 'R2a'; param_struct(2).range = R2a_range; param_struct(2).label = '$R_{2a}$ (Hz)';
param_struct(3).name = 'R1b'; param_struct(3).range = R1b_range; param_struct(3).label = '$R_{1b}$ (Hz)';
param_struct(4).name = 'MPF'; param_struct(4).range = MPF_range; param_struct(4).label = 'MPF';
param_struct(5).name = 'R';   param_struct(5).range = R_range;   param_struct(5).label = '$R$ (Hz)';
param_struct(6).name = 'T2b'; param_struct(6).range = T2b_range; param_struct(6).label = '$T_{2b}$';

% Pre-allocate storage for results
Results_Exact = cell(1, 6);
Results_Approx = cell(1, 6);

fprintf('Starting calculation (w1 = %d Hz)...\n', w1_val);

for p = 1:6
    p_name = param_struct(p).name;
    p_range = param_struct(p).range;
    n_vals = length(p_range);
    
    temp_exact = zeros(n_dw, n_vals);
    temp_approx = zeros(n_dw, n_vals);
    
    fprintf('  Calculating parameter: %s ...\n', p_name);
    
    for i = 1:n_vals
        % Reset all parameters to base values
        c_R1a = R1a_base;
        c_R2a = R2a_base;
        c_R1b = R1b_base;
        c_MPF = MPF_base;
        c_R   = R_base;
        c_T2b = T2b_base;
        c_T1d = T1d_base;
        
        % Update the current variable parameter
        val = p_range(i);
        switch p_name
            case 'R1a', c_R1a = val;
            case 'R2a', c_R2a = val;
            case 'R1b', c_R1b = val;
            case 'MPF', c_MPF = val;
            case 'R',   c_R = val;
            case 'T2b', c_T2b = val;
        end
        
        % Calculate for different dw
        for d = 1:n_dw
            dw = dw_all(d);
            temp_exact(d, i) = cal_Rdosl_exact(c_R1a, c_R2a, c_R1b, c_MPF, ...
                c_R, c_T2b,c_T1d, dw, w1);
            temp_approx(d, i) = cal_Rdosl_approx(c_R, c_T2b, c_T1d, dw, w1, B1,B0);
        end
    end
    Results_Exact{p} = temp_exact;
    Results_Approx{p} = temp_approx;
end

% =========================================================================
% Visualization
% =========================================================================

% General Plot Settings
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultLineLineWidth', 1.5);
colors = lines(n_dw);

% -------------------------------------------------------------------------
% Figure 1: Rdosl vs T1D
% -------------------------------------------------------------------------
figure('Color','w', 'Position', [100, 100, 500, 500]);
hold on; box on; grid on;

legend_str = cell(1, n_dw);
p_handles = zeros(1, n_dw);

for i = 1:n_dw
    % Exact Solution (Solid Line)
    p_handles(i) = plot(T1d_all*1000, Rdosl_exact_t1d(i,:), '-', ...
        'Color', colors(i,:), 'LineWidth', 2);
    
    % Approx Solution (Circle Markers)
    idx_sample = 1:1:length(T1d_all);
    plot(T1d_all(idx_sample)*1000, Rdosl_approx_t1d(i,idx_sample), 'o', ...
        'Color', colors(i,:), 'MarkerFaceColor', 'w', 'MarkerSize', 6);
    
    legend_str{i} = ['$\Delta\omega=' num2str(dw_vals(i)/1000) '\,\mathrm{kHz}$'];
end

xlim([1 10]);
ylim([0 2]);

xlabel('$T_{1D}$ (ms)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$R_{dosl}$ (Hz)', 'Interpreter', 'latex', 'FontSize', 18);
title({['$R_{dosl}$ vs. $T_{1D}$'], ...
       ['($\omega_1/2\pi = ' num2str(w1_val) '\,\mathrm{Hz}$)']}, ...
       'Interpreter', 'latex', 'FontSize', 17);

% Dual Legend Setup
h_line = plot(nan, nan, 'k-', 'LineWidth', 2);
h_circ = plot(nan, nan, 'ko', 'MarkerFaceColor', 'w');

leg1 = legend(p_handles, legend_str, 'Interpreter', 'latex', 'Location', 'northeast');
leg1_pos = leg1.Position;

ax2 = axes('Position', get(gca, 'Position'), 'Visible', 'off');
leg2 = legend(ax2, [h_line, h_circ], {'Exact', 'Approx'}, ...
              'Interpreter', 'latex', 'Location', 'northeast');

leg2.Position = [leg1_pos(1), ...
                 leg1_pos(2) - leg2.Position(4) - 0.01, ... 
                 leg2.Position(3), ...
                 leg2.Position(4)];

% -------------------------------------------------------------------------
% Figure 2: Sensitivity Analysis (6 Subplots)
% -------------------------------------------------------------------------
figure('Color','w', 'Position', [150, 50, 1200, 700]);
sgtitle(['Sensitivity Analysis of $R_{dosl}$ ($\omega_1/2\pi = ' num2str(w1_val) '\,\mathrm{Hz}$)'], ...
    'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

lw = 2;

for p = 1:6
    subplot(2, 3, p);
    hold on; box on; grid on;
    
    x_vals = param_struct(p).range;
    
    % Handle unit display scaling
    x_scale = 1;
    if strcmp(param_struct(p).name, 'T2b')
        x_scale = 1e6;
        xlabel_str = '$T_{2b}$ ($\mu$s)';
    elseif strcmp(param_struct(p).name, 'MPF') % Fixed string comparison
        x_scale = 1e2;
        xlabel_str = '$MPF$ (\%)';
    else
        xlabel_str = param_struct(p).label;
    end
    
    % Plot curves
    for d = 1:n_dw
        % Exact (Solid Line)
        plot(x_vals * x_scale, Results_Exact{p}(d, :), '-', ...
            'Color', colors(d,:), 'LineWidth', lw);
        
        % Approx (Circle Markers, Downsampled)
        idx_sample = 1:3:length(x_vals);
        plot(x_vals(idx_sample) * x_scale, Results_Approx{p}(d, idx_sample), 'o', ...
            'Color', colors(d,:), 'MarkerFaceColor', 'w', 'MarkerSize', 6);
    end
    
    xlim([min(x_vals) * x_scale, max(x_vals) * x_scale]);
    ylim([0 1.3]);
    
    set(gca, 'FontSize', 14, 'FontWeight','normal');
    
    xlabel(xlabel_str, 'Interpreter', 'latex', 'FontSize', 22, 'FontWeight','bold');
    ylabel('$R_{dosl}$ (Hz)', 'Interpreter', 'latex', 'FontSize', 22, 'FontWeight','bold');
end

fprintf('\n========================================\n');
fprintf('Calculation and plotting completed!\n');
fprintf('========================================\n');