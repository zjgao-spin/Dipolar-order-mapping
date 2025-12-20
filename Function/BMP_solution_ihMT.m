function M_t = BMP_solution_ihMT(t, M, R1a, R2a, dw, w1, w2, kab, kba, M0a, M0b, T2b, R1b, T1d, rfmt)
% BMP_SOLUTION_IHMT Solves the Bloch-McConnell-Provotorov equations for ihMT.
%
%   dM/dt = A * M + C
%
% Inputs:
%   t    - Evolution time duration (s)
%   M    - Initial magnetization state vector [Mxa; Mya; Mza; Mzb; Md]
%   R1a  - Longitudinal relaxation rate of free pool (1/s)
%   R2a  - Transverse relaxation rate of free pool (1/s)
%   dw   - Frequency offset (rad/s)
%   w1   - RF amplitude x-component (rad/s)
%   w2   - RF amplitude y-component (rad/s)
%   kab  - Exchange rate from pool A to B (1/s)
%   kba  - Exchange rate from pool B to A (1/s)
%   M0a  - Equilibrium magnetization of free pool
%   M0b  - Equilibrium magnetization of bound pool
%   T2b  - T2 of bound pool, the lineshape prameter
%   R1b  - Longitudinal relaxation rate of bound pool (1/s)
%   T1d  - Dipolar relaxation time (s)
%   rfmt - RF saturation rate for the bound pool (1/s)
%
% Output:
%   M_t  - Magnetization vector at time t

    %% 1. Define Physical Constants
    % Local dipolar field strength (D)
    D = 1/sqrt(15)/T2b; 
    
    % Initial state
    M_init = M; 

    %% 2. System Matrix A
    % Represents relaxation, precession, exchange, and RF saturation interactions.
    % State vector: [Mxa, Mya, Mza, Mzb, beta]'
    A = [-R2a,   dw,          -w2,        0,               0;                % Mxa
        -dw,   -R2a,          w1,        0,               0;                % Mya
         w2,   -w1, -(R1a + kab),      kba,               0;                % Mza
          0,     0,          kab, -(R1b + kba + rfmt),  rfmt * dw;          % Mzb
          0,     0,            0,   rfmt * dw / (D^2), -(1/T1d + rfmt * (dw / D)^2)]; % beta

    %% 3. Recovery Vector C
    % Drives the system towards thermal equilibrium.
    C = [0; 0; R1a * M0a; R1b * M0b; 0];

    %% 4. Analytical Solution
    % Solves the linear differential equation dM/dt = A*M + C
    % Solution: M(t) = exp(A*t) * (M_init - M_ss) + M_ss, where M_ss = -inv(A)*C
    M_t = expm(A * t) * (M_init + A \ C) - A \ C;

end