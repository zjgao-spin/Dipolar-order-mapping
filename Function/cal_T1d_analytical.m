function [T1d_val, fit_info] = cal_T1d_analytical(R_dosl, B1, B0, varargin)
% FIT_T1D_SINGLE_POINT - Single-point T1d parameter fitting
%
% Input Parameters (Scalar):
%   R_dosl   - Observed R_dosl value (numeric)
%   B1       - B1 correction factor (numeric, 1.0 = nominal)
%   B0       - B0 offset (numeric, Hz, usually converted to rad/s inside)
%   varargin - Optional parameters (Name-Value pairs)
%
% Output Parameters:
%   T1d_val  - Fitted T1d value (Unit: s)
%   fit_info - Structure containing residuals, exit flags, etc.

    % 1. Parse input parameters
    p = inputParser;
    addRequired(p, 'R_dosl', @isnumeric);
    addRequired(p, 'B1', @isnumeric);
    addRequired(p, 'B0', @isnumeric);
    
    % Default physical parameters (Consistent with original code)
    addParameter(p, 'dw', 5000*2*pi, @isnumeric);
    addParameter(p, 'w1', 500*2*pi, @isnumeric);
    addParameter(p, 'T2b', 10e-6, @isnumeric);
    addParameter(p, 'R', 23, @isnumeric);
    addParameter(p, 'T1d_init', 3.5e-3, @isnumeric);
    addParameter(p, 'T1d_bounds', [0, inf], @isnumeric);
    
    parse(p, R_dosl, B1, B0, varargin{:});
    
    % Extract parameters
    dw = p.Results.dw;
    w1 = p.Results.w1;
    T2b = p.Results.T2b;
    R = p.Results.R;
    T1d_init = p.Results.T1d_init;
    T1d_bounds = p.Results.T1d_bounds;
    
    % R1b = 1/(340*1e-3);
    
    % 2. Pre-calculate physical constants
    % Calculate actual RF frequency (rad/s)
    w1 = B1 * w1;
    dw = dw + B0*2*pi;
    
    % Calculate RF saturation rate (Note: Requires RF_MT function in your path)
    % If RF_MT is missing, you need to provide the code or replace this logic.
    try
        % This call is just to check if RF_MT works or to pre-calc if needed, 
        % though it's recalculated inside the approx function usually.
        % Here it seems unused in the main scope but good for validation.
        RF_val = RF_MT(T2b, w1, dw, 'SuperLorentzian');
    catch ME
        error('Cannot call RF_MT function. Please ensure it is in the path. Error: %s', ME.message);
    end

    % 3. Define optimization problem
    % Objective function: Difference between theoretical R and observed R_dosl
    objective_func = @(T1d) cal_Rdosl_approx(R, T2b, dw, T1d, w1) - R_dosl;
    
    % Optimization options
    % Note: 'levenberg-marquardt' does not handle bounds. 
    % If bounds are provided, lsqnonlin automatically switches to 'trust-region-reflective'.
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                          'Algorithm', 'levenberg-marquardt', ...
                          'FunctionTolerance', 1e-3, ...
                          'StepTolerance', 1e-2, ...
                          'MaxIterations', 10000);

    % options = optimoptions('lsqnonlin', 'Display', 'off', ...
    %     'FunctionTolerance', 1e-3, ...
    %     'StepTolerance', 1e-2, ...
    %     'MaxIterations', 500);

    % 4. Execute fitting
    try
        [T1d_fit, resnorm, residual, exitflag, output] = lsqnonlin(objective_func, ...
            T1d_init, T1d_bounds(1), T1d_bounds(2), options);
        
        % Process results
        T1d_val = T1d_fit; 
        
        % Assemble fitting information
        fit_info.residual = residual;
        fit_info.resnorm = resnorm;
        fit_info.exitflag = exitflag;
        fit_info.iterations = output.iterations;
        fit_info.R_theory = cal_Rdosl_approx(R, T2b, dw, T1d_fit, w1); % Recalculate with fitted T1d
        
    catch ME
        warning('Fitting failed: %s', ME.message);
        T1d_val = NaN;
        fit_info = struct();
    end
end

function Rdosl_approx = cal_Rdosl_approx(R, T2b, dw, T1d, w1)
    % Dipolar order term
    D = 1/sqrt(15)/T2b;
    N_scale = 5; % Scaling factor for the low-power condition
    
    % Note: Assuming RF_MT function exists in path
    % Rrfb1: High power saturation rate
    Rrfb1 = RF_MT(T2b, w1, dw, 'SuperLorentzian');
    % Rrfb2: Low power saturation rate (Scaled down)
    Rrfb2 = RF_MT(T2b, w1/N_scale, dw/N_scale, 'SuperLorentzian');
    
    % Analytical approximation formula for R_DOSL
    numerator = (R + Rrfb2) * Rrfb1^2 * T1d * (dw/D)^2;
    denominator = (R + Rrfb1) * (Rrfb2 * (Rrfb1 * T1d * (dw/D)^2 + 1) - Rrfb1);

    Rdosl_approx = abs(numerator / denominator);
end