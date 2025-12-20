% -------------------------------------------------------------------------
% Script:       T1D_Estimation.m
% Authors:      Zijian Gao,Weitian Chen
% Email:        zijian.gao@link.cuhk.edu.hk
% Date:         2025-12-20
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
%   This script performs a Monte Carlo simulation to evaluate the robustness 
%   of T1d (Dipolar Relaxation Time) estimation methods against tissue
%   parameters from MT pool
% 
% Dependencies:
%   - T1d_dictionary.mat (Pre-calculated dictionary)
%   - cal_Rdosl_acquired.m      (Forward model function)
%   - cal_T1d_analytical.m      (Inverse analytical function)
%   - cal_T1d_dictionary.m      (Inverse dictionary function)
% -------------------------------------------------------------------------

clc; clear; close all;

%% 0. Initialization
% Load the pre-calculated dictionary for the inverse model
% Ensure this file exists in your path
try
    dictionary = load('T1d_dictionary_B1_2025-11-27.mat');
    fprintf('Dictionary loaded successfully.\n');
catch
    error('Dictionary file not found. Please check the filename and path.');
end

%% 1. Define Simulation Parameters
% =========================================================================
% 1.1 Basic Tissue Parameters (White Matter Model - Ground Truth)
% =========================================================================
T1a_base = 1840e-3;     % T1 of free water pool (s)
T1b_base = 340e-3;      % T1 of bound pool (s)
T2a_base = 69e-3;       % T2 of free water pool (s)
T2b_base = 10e-6;       % T2 of bound pool (s)
T1d_base = 6.2e-3;      % Dipolar relaxation time (s) [Target Ground Truth]
MPF_base = 0.15;        % Macromolecular Proton Fraction (Nominal)
R_base   = 23;          % Exchange rate (Hz) (Nominal)

% Convert Time constants to Rates (Hz) for calculation
R1a_base = 1/T1a_base;
R2a_base = 1/T2a_base;
R1b_base = 1/T1b_base;

% =========================================================================
% 1.2 Scan Protocol Parameters
% =========================================================================
B1_base = 1;            % Nominal B1 (normalized, 1 = 100%)
B0_base = 0;            % Nominal B0 (Hz, 0 = On-resonance)

% Spin-Lock Parameters
TSL     = 80 * 1e-3;    % Spin-Lock Duration (s)
w1_base = 500 * 2 * pi; % Spin-lock frequency (rad/s)
dw_base = 5000 * 2* pi; % Off-resonance frequency (rad/s)

% =========================================================================
% 1.3 Monte Carlo Simulation Ranges
% =========================================================================
% The simulation will vary these parameters uniformly within ranges
range_R1b = [2, 5];     % R1 of bound pool (Hz)
range_MPF = [0.1, 0.2]; % Macromolecular Proton Fraction
range_R   = [15, 30];   % Exchange Rate (Hz)

num_sims  = 10000;        % Number of Monte Carlo iterations

%% 2. Simulation Setup
% =========================================================================
% Generate Random Distributions
% =========================================================================
fprintf('Generating %d random parameter sets...\n', num_sims);

% Generate uniform random distributions: r = a + (b-a).*rand(N,1)
R1b_rand = range_R1b(1) + (range_R1b(2)-range_R1b(1)) * rand(num_sims, 1);
MPF_rand = range_MPF(1) + (range_MPF(2)-range_MPF(1)) * rand(num_sims, 1);
R_rand   = range_R(1)   + (range_R(2)-range_R(1))   * rand(num_sims, 1);

% Pre-allocate memory for results (Optimization)
Results_MC = struct();
Results_MC.Rdosl      = zeros(1, num_sims); % Observed R1rho
Results_MC.T1d_Analyt = zeros(1, num_sims); % Estimated T1d (Analytical)
Results_MC.T1d_Dic    = zeros(1, num_sims); % Estimated T1d (Dictionary)

%% 3. Monte Carlo Simulation Loop
% =========================================================================
% Forward Model (Generate Signal) -> Inverse Model (Fit T1d)
% =========================================================================
fprintf('Starting simulation loop...\n');
tic; % Start timer

for i = 1:num_sims
    
    % 3.1 Get current random parameters
    p_R1b = R1b_rand(i);
    p_MPF = MPF_rand(i);
    p_R   = R_rand(i);
    
    % ---------------------------------------------------------------------
    % A. Forward Model: Generate "Observed" R1rho (Rdosl)
    % ---------------------------------------------------------------------
    % Note: T1d is fixed at T1d_base (Ground Truth) for this study
    Rdosl = cal_Rdosl_acquired(R1a_base, R2a_base, p_R1b, p_MPF, p_R, ...
                               T2b_base, T1d_base, dw_base, w1_base, ...
                               TSL, B1_base, B0_base);    

    % Store observed Rdosl
    Results_MC.Rdosl(i) = Rdosl;

    % ---------------------------------------------------------------------
    % B. Inverse Model: Estimate T1d
    % ---------------------------------------------------------------------
    
    % Method 1: Analytical Solution
   % [T1d_val, ~] = cal_T1d_analytical_direct(Rdosl, B1_base, B0_base);
    [T1d_val, ~] = cal_T1d_analytical(Rdosl, B1_base, B0_base);
    Results_MC.T1d_Analyt(i) = T1d_val;   

    % Method 2: Dictionary Matching
    T1D_dic_val = cal_T1d_dictionary(Rdosl, B1_base, dictionary);
    Results_MC.T1d_Dic(i) = T1D_dic_val; 
    
end

sim_time = toc;
fprintf('Simulation finished in %.2f seconds.\n', sim_time);

%% 4. Statistical Analysis & Save Results
% =========================================================================
% Calculate Error Metrics and Statistics, then Save to TXT
% =========================================================================

% --- 4.1 Observed Signal Statistics ---
Vals_Rdosl = Results_MC.Rdosl;
Stats.Rdosl.mean = mean(Vals_Rdosl);
Stats.Rdosl.std  = std(Vals_Rdosl);

% --- 4.2 T1d Estimation Performance ---
% Ground Truth
GT = T1d_base; 

% Extract Results
Est_Analyt = Results_MC.T1d_Analyt;
Est_Dic    = Results_MC.T1d_Dic;

% 1. Calculate Absolute Errors (Estimate - GT)
Err_Analyt = Est_Analyt - GT;
Err_Dic    = Est_Dic - GT;

% 2. Calculate Relative Errors (%)
RelErr_Analyt = (Err_Analyt ./ GT) * 100;
RelErr_Dic    = (Err_Dic    ./ GT) * 100;

% 3. Calculate Statistical Metrics
% Helper function: [Mean, Std, RMSE, MAPE]
calc_metrics = @(est, err) [mean(est), std(est), sqrt(mean(err.^2)), mean(abs(err)./GT)*100];

% Calculate metrics for both methods
M_Analyt = calc_metrics(Est_Analyt, Err_Analyt);
M_Dic    = calc_metrics(Est_Dic, Err_Dic);

% --- 4.3 Display and Save to File ---
output_filename = 'T1d_Simulation_Results.txt';
fileID = fopen(output_filename, 'w'); % Open file for writing

% Define a helper function to print to both screen and file
print_both = @(fid, fmt, varargin) fprintf(fid, fmt, varargin{:});

% Loop twice: once for screen (1), once for file (fileID)
targets = [1, fileID]; 

for fid = targets
    print_both(fid, '\n--- Observed Signal Statistics ---\n');
    print_both(fid, 'Rdosl (R1rho): Mean = %.4f Hz, Std = %.4f Hz\n', Stats.Rdosl.mean, Stats.Rdosl.std);
    
    print_both(fid, '\n--- T1d Estimation Performance Analysis ---\n');
    print_both(fid, '\n================================================================\n');
    print_both(fid, '  T1d Estimation Summary (N = %d Simulations)\n', num_sims);
    print_both(fid, '================================================================\n');
    print_both(fid, 'Ground Truth T1d: %.2f ms\n', GT * 1000);
    print_both(fid, '----------------------------------------------------------------\n');
    print_both(fid, '%-15s | %-10s | %-10s | %-10s | %-10s\n', 'Method', 'Mean(ms)', 'Std(ms)', 'RMSE(ms)', 'MAPE(%%)');
    print_both(fid, '----------------------------------------------------------------\n');
    print_both(fid, '%-15s | %-10.2f | %-10.2f | %-10.2f | %-10.2f\n', ...
        'Analytical', M_Analyt(1)*1000, M_Analyt(2)*1000, M_Analyt(3)*1000, M_Analyt(4));
    print_both(fid, '%-15s | %-10.2f | %-10.2f | %-10.2f | %-10.2f\n', ...
        'Dictionary', M_Dic(1)*1000, M_Dic(2)*1000, M_Dic(3)*1000, M_Dic(4));
    print_both(fid, '================================================================\n');
end

fclose(fileID); % Close the file
fprintf('\nResults saved to: %s\n', output_filename);
