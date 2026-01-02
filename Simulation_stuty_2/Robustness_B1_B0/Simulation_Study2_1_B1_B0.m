% =========================================================================
% Script:       Simulation_study2_1_B1_B0.m
% Authors:      Zijian Gao,Qianxue Shan,Weitian Chen
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
%   This script performs a simulation to analyze the sensitivity 
%   of RATIO_dosl with B1 and B0 inhomogeneity.
% =========================================================================

clc; clear all; close all;

% Add path to helper functions
addpath('..\..\Function\')


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
B1_base = 1; 
B0_base = 0;

B1_all = linspace(0.7, 1.3, 100);
n_b1 = length(B1_all);
B0_all = linspace(-100, 100, 400); %Hz
n_b0 = length(B0_all);

TSL = 80*1e-3;          % Spin-Lock Duration (s)

w1_base = 500 * 2 * pi; % Spin-lock frequency (rad/s)
dw_base = 5000 * 2 * pi;% Off-resonance frequency (rad/s)

%% 2.Simulation Loop
RATIO_dosl_base =  cal_RATIO_dosl_acquired(R1a_base, R2a_base, R1b_base, MPF_base, R_base, T2b_base, T1d_base, ...
                                        dw_base, w1_base, TSL, B1_base, B0_base);
for i_b1 = 1:n_b1
    for i_b0 = 1:n_b0
        B1=B1_all(i_b1);
        B0=B0_all(i_b0)*2*pi;
        RATIO_dosl_all(i_b1,i_b0)=cal_RATIO_dosl_acquired(R1a_base, R2a_base, R1b_base, MPF_base, R_base, T2b_base, T1d_base, ...
            dw_base, w1_base, TSL, B1, B0);
    end
end


R_dosl_BMP=RATIO_dosl_all;

%% 4. Visualization
fontSize_axis=18;
b0_Hz = B0_all;

%b1_nominal_idx = find(abs(b1_all - 1.0) < 0.01, 1);
b1_nominal_idx = n_b1/2;
b0_nominal_idx = n_b0/2; % 
R_dosl_reference = RATIO_dosl_base;

Relative_Error = abs((R_dosl_BMP - R_dosl_reference) / R_dosl_reference) * 100;

figure('Position', [100 100 1400 500]);

% --- Subplot 1: 2D Plot vs B0 ---
subplot(1,2,1);
plot(b0_Hz, R_dosl_BMP(b1_nominal_idx,:), 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410]);
grid on;
xlabel('$B_0$ Offset (Hz)', 'FontSize', 18);
ylabel('$RATIO_{dosl}$', 'FontSize', 18);
title('$RATIO_{dosl}$ vs $B_0$ ($B_1$ = 1.0)', 'FontSize', 18, 'FontWeight', 'bold');
% Note: You might want to remove the hardcoded ylim if data changes, 
% but I kept it as per your original code.
ylim_vals = [0.2 0.6]; 
ylim(ylim_vals); 
xlim([-100 100])
set(gca, 'FontSize', fontSize_axis);
hold on;
plot([0 0], ylim_vals, 'r--', 'LineWidth', 1.5);
hold off;

% --- Subplot 2: 2D Plot vs B1 ---
subplot(1,2,2);
plot(B1_all, R_dosl_BMP(:,b0_nominal_idx), 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980]);
grid on;
xlabel('$B_1$ Scale (n.u.)', 'FontSize', 18);
ylabel('$RATIO_{dosl}$', 'FontSize', 18);
title('$RATIO_{dosl}$ vs $B_1$ ($B_0$ = 0 Hz)', 'FontSize', 18, 'FontWeight', 'bold');

ylim_vals = [0.2 0.6]; 
ylim(ylim_vals); 

xlim([0.7 1.3])
set(gca, 'FontSize', fontSize_axis);
hold on;
plot([1.0 1.0], ylim_vals, 'r--', 'LineWidth', 1.5);
hold off;

