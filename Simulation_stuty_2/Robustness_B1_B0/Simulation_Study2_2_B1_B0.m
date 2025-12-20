% -------------------------------------------------------------------------
% Script:       Simulation_study2_2_B1_B0.m
% Authors:      Zijian Gao,Qianxue Shan,Weitian Chen
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
%   of T1d (Dipolar Relaxation Time) estimation methods against field 
%   inhomogeneities
% 
% Dependencies:
%   - T1d_dictionary_B1_2025-11-27.mat (Pre-calculated dictionary)
%   - cal_Rdosl_acquired.m      (Forward model function)
%   - cal_T1d_analytical.m      (Inverse analytical function)
%   - cal_T1d_dictionary.m      (Inverse dictionary function)
% -------------------------------------------------------------------------
clc
clear all
close all

% Load the pre-calculated dictionary
dictionary = load('T1d_dictionary_B1_2025-11-27.mat');

%% 1. Define Base Parameters
% =========================================================================
% 1.1 Basic Tissue Parameters (White Matter Model)
% =========================================================================
T1a_base = 1840e-3;     % T1 of free water pool (s)
T1b_base = 340e-3;      % T1 of bound pool (s)
T2a_base = 69e-3;       % T2 of free water pool (s)
T2b_base = 10e-6;       % T2 of bound pool (s)
T1d_base = 6.2e-3;      % Dipolar relaxation time (s) [Ground Truth]
MPF_base = 0.15;        % Macromolecular Proton Fraction
R_base   = 23;          % Exchange rate (Hz)

% Convert Time constants to Rates (Hz)
R1a_base = 1/T1a_base;
R2a_base = 1/T2a_base;
R1b_base = 1/T1b_base;

% =========================================================================
% 1.2 General Scan Parameters
% =========================================================================
% RF pulse shape parameters
B1_base = 1;            % Nominal B1 (normalized)
B0_base = 0;            % Nominal B0 (Hz)

% Simulation Ranges
range_B1 = [0.8, 1.3];  % B1 inhomogeneity range (normalized units)
range_B0 = [-100, 100]; % B0 off-resonance range (Hz)

% Spin-Lock Parameters
TSL     = 80 * 1e-3;    % Spin-Lock Duration (s)
w1_base = 500 * 2 * pi; % Spin-lock frequency (rad/s)
dw_base = 5000 * 2* pi; % Off-resonance frequency (rad/s)

%% 2. Simulation Setup
% =========================================================================
% Generate Random Distributions for Monte Carlo
% =========================================================================
num_sims = 10000; 

% Generate uniform random distributions: r = a + (b-a).*rand(N,1)
B1_rand = range_B1(1) + (range_B1(2)-range_B1(1)) * rand(num_sims, 1);
B0_rand = range_B0(1) + (range_B0(2)-range_B0(1)) * rand(num_sims, 1);

% Pre-allocate memory for results to improve speed
Results_MC = struct();

% Intermediate R1rho (Rdosl) values
Results_MC.Rdosl.B1ihn = zeros(1, num_sims);
Results_MC.Rdosl.B0ihn = zeros(1, num_sims);

% T1d Estimation Results
Results_MC.B1ih.T1d_Analyt.B1cor   = zeros(1, num_sims);
Results_MC.B1ih.T1d_Analyt.B1uncor = zeros(1, num_sims);
Results_MC.B1ih.T1d_Dic.B1cor      = zeros(1, num_sims);
Results_MC.B1ih.T1d_Dic.B1uncor    = zeros(1, num_sims);

Results_MC.B0ih.T1d_Analyt         = zeros(1, num_sims);
Results_MC.B0ih.T1d_Dic            = zeros(1, num_sims);

%% 3. Monte Carlo Simulation Loop
% =========================================================================
% Forward Model (Generate Signal) -> Inverse Model (Fit T1d)
% =========================================================================
tic; % Start timer

for i = 1:num_sims
    
    % Get current random parameters
    p_B1 = B1_rand(i);
    p_B0 = B0_rand(i);

    % ---------------------------------------------------------------------
    % A. Forward Model: Generate "Observed" R1rho (Rdosl)
    % ---------------------------------------------------------------------
    % Calculate Rdosl with B1 inhomogeneity (B0 is ideal)
    Rdosl_B1ihn = cal_Rdosl_acquired(R1a_base, R2a_base, R1b_base, MPF_base, R_base, ...
                                     T2b_base, T1d_base, dw_base, w1_base, TSL, p_B1, B0_base);
    
    % Calculate Rdosl with B0 inhomogeneity (B1 is ideal)
    Rdosl_B0ihn = cal_Rdosl_acquired(R1a_base, R2a_base, R1b_base, MPF_base, R_base, ...
                                     T2b_base, T1d_base, dw_base, w1_base, TSL, B1_base, p_B0);
    
    % Store observed rates
    Results_MC.Rdosl.B1ihn(i) = Rdosl_B1ihn;
    Results_MC.Rdosl.B0ihn(i) = Rdosl_B0ihn;

    % ---------------------------------------------------------------------
    % B. Inverse Model: Estimate T1d under B1 Inhomogeneity
    % ---------------------------------------------------------------------
    
    % Scenario 1: B1 Corrected (We know the true B1)
    B1_known = p_B1;
    [T1d_val, ~] = cal_T1d_analytical(Rdosl_B1ihn, B1_known, B0_base);
    Results_MC.B1ih.T1d_Analyt.B1cor(i) = T1d_val;   

    T1D_dic_val = cal_T1d_dictionary(Rdosl_B1ihn, B1_known, dictionary);
    Results_MC.B1ih.T1d_Dic.B1cor(i) = T1D_dic_val; 

    % Scenario 2: B1 Uncorrected (We assume B1 is nominal 1.0)
    B1_assumed = B1_base; 
    [T1d_val, ~] = cal_T1d_analytical(Rdosl_B1ihn, B1_assumed, B0_base);
    Results_MC.B1ih.T1d_Analyt.B1uncor(i) = T1d_val;   

    T1D_dic_val = cal_T1d_dictionary(Rdosl_B1ihn, B1_assumed, dictionary);
    Results_MC.B1ih.T1d_Dic.B1uncor(i) = T1D_dic_val; 

    % ---------------------------------------------------------------------
    % C. Inverse Model: Estimate T1d under B0 Inhomogeneity
    % ---------------------------------------------------------------------
    % We assume B0 is 0 (Ignored), but signal has off-resonance effects
    B0_assumed = B0_base; 
    B1_nominal = B1_base;
    
    [T1d_val, ~] = cal_T1d_analytical(Rdosl_B0ihn, B1_nominal, B0_assumed);
    Results_MC.B0ih.T1d_Analyt(i) = T1d_val;   

    T1D_dic_val = cal_T1d_dictionary(Rdosl_B0ihn, B1_nominal, dictionary);
    Results_MC.B0ih.T1d_Dic(i) = T1D_dic_val; 
    
end

toc; % End timer

%% 4. Statistical Analysis & Error Calculation
% =========================================================================
% Calculate Mean, Std, and Relative Error
% =========================================================================

% --- 4.1 Basic Statistics for Observed Rates ---
Vals_Rdosl_B1ihn = Results_MC.Rdosl.B1ihn;
Stats.Rdosl_B1ihn.mean = mean(Vals_Rdosl_B1ihn);
Stats.Rdosl_B1ihn.std  = std(Vals_Rdosl_B1ihn);

Vals_Rdosl_B0ihn = Results_MC.Rdosl.B0ihn;
Stats.Rdosl_B0ihn.mean = mean(Vals_Rdosl_B0ihn);
Stats.Rdosl_B0ihn.std  = std(Vals_Rdosl_B0ihn);

fprintf('\n--- Observed R1rho Statistics ---\n');
fprintf('Rdosl (B1 Inhomogeneity):\n  Mean: %.4f Hz, Std: %.4f Hz\n', Stats.Rdosl_B1ihn.mean, Stats.Rdosl_B1ihn.std);
fprintf('Rdosl (B0 Inhomogeneity):\n  Mean: %.4f Hz, Std: %.4f Hz\n', Stats.Rdosl_B0ihn.mean, Stats.Rdosl_B0ihn.std);

% --- 4.2 T1d Error Analysis ---
fprintf('\n================================================\n');
fprintf('Monte Carlo Error Analysis Results (N=%d)\n', num_sims);
fprintf('Ground Truth T1d: %.4f ms\n', T1d_base * 1e3);
fprintf('================================================\n');

% Helper function for stats
calc_stats = @(x, true_val) struct('mean', mean(x), ...
                                   'std', std(x), ...
                                   'error_pct', (mean(x) - true_val)/true_val * 100);

% --- Analyze B1 Effects ---
fprintf('\n[4.1] B1 Inhomogeneity Effects (Range: %.1f - %.1f)\n', range_B1(1), range_B1(2));

stats_B1cor_ana   = calc_stats(Results_MC.B1ih.T1d_Analyt.B1cor, T1d_base);
stats_B1uncor_ana = calc_stats(Results_MC.B1ih.T1d_Analyt.B1uncor, T1d_base);
stats_B1cor_dic   = calc_stats(Results_MC.B1ih.T1d_Dic.B1cor, T1d_base);
stats_B1uncor_dic = calc_stats(Results_MC.B1ih.T1d_Dic.B1uncor, T1d_base);

fprintf('%-20s | %-12s | %-12s | %-12s\n', 'Scenario', 'Mean (ms)', 'Std (ms)', 'Error (%)');
fprintf('----------------------------------------------------------------\n');
fprintf('%-20s | %-12.4f | %-12.4f | %-12.2f%%\n', 'Analyt (B1 Cor)', stats_B1cor_ana.mean*1e3, stats_B1cor_ana.std*1e3, stats_B1cor_ana.error_pct);
fprintf('%-20s | %-12.4f | %-12.4f | %-12.2f%%\n', 'Analyt (B1 Uncor)', stats_B1uncor_ana.mean*1e3, stats_B1uncor_ana.std*1e3, stats_B1uncor_ana.error_pct);
fprintf('%-20s | %-12.4f | %-12.4f | %-12.2f%%\n', 'Dict (B1 Cor)', stats_B1cor_dic.mean*1e3, stats_B1cor_dic.std*1e3, stats_B1cor_dic.error_pct);
fprintf('%-20s | %-12.4f | %-12.4f | %-12.2f%%\n', 'Dict (B1 Uncor)', stats_B1uncor_dic.mean*1e3, stats_B1uncor_dic.std*1e3, stats_B1uncor_dic.error_pct);

% --- Analyze B0 Effects ---
fprintf('\n[4.2] B0 Inhomogeneity Effects (Range: %.1f - %.1f Hz)\n', range_B0(1), range_B0(2));

stats_B0_ana = calc_stats(Results_MC.B0ih.T1d_Analyt, T1d_base);
stats_B0_dic = calc_stats(Results_MC.B0ih.T1d_Dic, T1d_base);

fprintf('%-20s | %-12s | %-12s | %-12s\n', 'Scenario', 'Mean (ms)', 'Std (ms)', 'Error (%)');
fprintf('----------------------------------------------------------------\n');
fprintf('%-20s | %-12.4f | %-12.4f | %-12.2f%%\n', 'Analyt (B0 Ignored)', stats_B0_ana.mean*1e3, stats_B0_ana.std*1e3, stats_B0_ana.error_pct);

fprintf('%-20s | %-12.4f | %-12.4f | %-12.2f%%\n', 'Dict (B0 Ignored)', stats_B0_dic.mean*1e3, stats_B0_dic.std*1e3, stats_B0_dic.error_pct);
