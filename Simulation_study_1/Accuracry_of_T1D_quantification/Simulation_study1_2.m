% -------------------------------------------------------------------------
% Script:       T1D_Estimation.m
% Authors:      Zijian Gao,Weitian Chen
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
%   This script performs a numerical simulation to evaluate the accuracy of 
%   T1D (Dipolar Relaxation Time) estimation.
% 
% Dependencies:
%   - cal_Rdosl_acquired.m
%   - cal_T1d_analytical.m
%   - Simulation_study1_2_plot.m
% -------------------------------------------------------------------------

clc; clear all; close all;

% Add function path (Adjust as necessary to point to your library)
addpath('..\..\Function\')


%% 1. Parameter Initialization


% --- Basic Tissue Parameters (White Matter Model) ---
% These parameters define the biological properties of the simulated tissue.
T1a_base = 1840e-3;   % T1 of the free water pool (s)
T1b_base = 340e-3;    % T1 of the bound (macromolecular) pool (s)
T2a_base = 69e-3;     % T2 of the free water pool (s)
T2b_base = 10e-6;     % T2 of the bound pool (s) - very short (semi-solid)
T1d_base = 6.2e-3;    % Dipolar relaxation time (s) - The target parameter
MPF_base = 0.15;      % Macromolecular Proton Fraction (pool size ratio)
R_base   = 23;        % Magnetization exchange rate between pools (Hz)

% Define the Ground Truth for T1D
T1D_true = T1d_base;

% --- Convert Time Constants to Rates (Hz) ---
R1a_base = 1/T1a_base;
R2a_base = 1/T2a_base;
R1b_base = 1/T1b_base;

% --- Simulation Grid Settings ---
% Define the range of experimental conditions to simulate.

% RF pulse shape parameters
B1 = 1; % B1 scaling factor (1 = ideal B1 field)
B0 = 0; % B0 offset (Hz) (0 = on resonance)

% Frequency Offsets (Delta Omega)
% Range: 2000 to 7000 Hz off-resonance
dw_vals_Hz = 2000:1000:7000;
dw_vals_rad = dw_vals_Hz * 2 * pi; % Convert to radians/s

% RF Field Amplitude (Omega 1)
% Range: 200 to 800 Hz power
w1_vals_Hz = 200:100:800;
w1_vals_rad = w1_vals_Hz * 2 * pi; % Convert to radians/s

% Spin Lock Time (TSL)
% Range: 50 to 150 ms duration
TSL_vals_ms = 50:10:150;
TSL_vals_s  = TSL_vals_ms * 1e-3; % Convert to seconds

% Get dimensions for the loops
n_dw  = length(dw_vals_Hz);
n_w1  = length(w1_vals_Hz);
n_tsl = length(TSL_vals_ms);

% Pre-allocate arrays for speed
% Rdosl_acquired: Stores the synthetic observed relaxation rate
% T1D_est: Stores the recovered T1D parameter
Rdosl_acquired = zeros(n_tsl, n_w1, n_dw);
T1D_est        = zeros(n_tsl, n_w1, n_dw); 

%% 2. Simulation and Estimation Loop

fprintf('Starting simulation and estimation...\n');

% Iterate over all Spin-Lock Times
for i_tsl = 1:n_tsl
    curr_TSL = TSL_vals_s(i_tsl);
    
    % Iterate over all RF Amplitudes
    for i_w1 = 1:n_w1
        curr_w1 = w1_vals_rad(i_w1);
        
        % Iterate over all Frequency Offsets
        for i_dw = 1:n_dw
            curr_dw = dw_vals_rad(i_dw);
            
            % ---------------------------------------------------------
            % Step 1: Generate Synthetic Data (Forward Model)
            % ---------------------------------------------------------
            % Calculate the theoretical R1rho dispersion (Rdosl) based on 
            % the defined tissue parameters and current experimental setup.
            Rdosl = cal_Rdosl_acquired(R1a_base, R2a_base, R1b_base, MPF_base, ...
                                       R_base, T2b_base, T1d_base, ...
                                       curr_dw, curr_w1, curr_TSL, B1, B0);
                                   
            % Store the "acquired" data point
            Rdosl_acquired(i_tsl, i_w1, i_dw) = Rdosl;
            
            % ---------------------------------------------------------
            % Step 2: Estimate Parameter (Inverse Model)
            % ---------------------------------------------------------
            % Attempt to estim T1D from the synthetic Rdosl value using 
            % the analytical formula. This tests the validity of the 
            % analytical approximation.
            [T1d_val, ~] = cal_T1d_analytical(Rdosl, B1, B0, 'w1', curr_w1, 'dw', curr_dw);
            
            % Store the estimated parameter
            T1D_est(i_tsl, i_w1, i_dw) = T1d_val;
        end
    end
end


%% 3. Error Analysis

% Calculate Relative Error (%): |Est - True| / True * 100
relative_error = abs(T1D_est - T1D_true) / T1D_true * 100;

% Calculate Absolute Error (ms): |Est - True|
absolute_error = abs(T1D_est - T1D_true) * 1000;

% Display Statistics to the Command Window
% 'omitnan' is used to safely calculate means even if some estimations failed (NaN)
fprintf('\n=== T1D Estimation Error Statistics ===\n');
fprintf('Ground Truth T1D: %.2f ms\n', T1D_true*1000);
fprintf('Mean Est. T1D:    %.2f ± %.2f ms\n', mean(T1D_est(:),'omitnan')*1000, std(T1D_est(:),'omitnan')*1000);
fprintf('Mean Rel. Error:  %.2f ± %.2f %%\n', mean(relative_error(:),'omitnan'), std(relative_error(:),'omitnan'));
fprintf('Mean Abs. Error:  %.2f ± %.2f ms\n', mean(absolute_error(:),'omitnan'), std(absolute_error(:),'omitnan'));


%% 4. Visualization

% Call external script to plot the results

Simulation_study1_2_plot
