function RATIO = cal_RATIO_dosl_acquired(R1a, R2a, R1b, MPF, R, T2b, T1d,dw, w1, TSL, B1, B0, SNR)
% COMPUTE_ACQUIRED_RDOSL Calculates the Rdosl.
% Simulates magnetization evolution under ihMT pulse sequences.

%% 1. Model Initialization (Two-pool model)
M0b = MPF;          % Macromolecular pool fraction
M0a = 1-M0b;        % Free water pool fraction
fb = M0b/M0a;       % Pool ratio

% Exchange rates (Detailed balance: kba*M0a = kab*M0b)
kba = R*M0a;        % A -> B rate
kab = kba*fb;       % B -> A rate

base_freq_offset = dw;
base_SL_val = w1;   % Base spin-lock amplitude

b1_scale = B1;      % B1 inhomogeneity scaling
b0 = B0;             % Main field offset

%% 3. Spin Lock Protocol Settings
Ts_long = TSL/2;
Ts_short = 0.5*1e-3;
eTSL_dual = Ts_short.*2;           % Duration of one SL segment = 2* Ts
N_eTSL_dual = round(TSL / eTSL_dual);  % Calculate number of segments based on TSL
eTSL_single = Ts_long.*2;           % Duration of one SL segment = 2* Ts
N_eTSL_single = round(TSL / eTSL_single);  % Calculate number of segments based on TSL


N_scale = 5;           % Scaling factor for the second condition


%% 4. RF Hardware Parameters
gyr = 42.58e6;         % Gyromagnetic ratio (Hz/T)
Bexp = 13.5e-6;        % Expected B1 (T)
flip_angl_180 = pi;    % 180 degree flip angle
a_rf_exp = gyr * Bexp * 2 * pi;
pw_rf_180 = abs(flip_angl_180 / a_rf_exp); % Pulse width for 180
wait_time = 2e-3;      % Crusher/Wait time

%% 5. Simulation Loop
% temp_MzA stores final Z-magnetization of Pool A
% Dimensions: [Scale Condition, Dual/Single Mode, Toggle State]
temp_MzA = zeros(2,2,2);

for i_n = 1:2
    % Define RF parameters based on scaling condition (High vs Low power/offset)
    if i_n == 2
        curr_freq_offset = base_freq_offset / N_scale;
        curr_SL_val = base_SL_val / N_scale;
    else
        curr_freq_offset = base_freq_offset;
        curr_SL_val = base_SL_val;
    end

    % --- Pre-calculate RF parameters with B1 scaling ---
    b1_fx = 400; % Fixed RF amplitude basis

    a_rf_inh = gyr * (Bexp * b1_scale) * 2 * pi;  % Actual RF strength
    a_rf_fx_inh = (b1_fx * b1_scale) * 2 * pi;    % Scaled fixed RF
    a_rf_fx_exp = b1_fx * 2 * pi;                 % Nominal fixed RF

    % Calculate Tip-down/up pulse parameters
    flip_ang1 = atan(curr_SL_val / curr_freq_offset);
    pw_rf = abs(flip_ang1 / a_rf_fx_exp);

    % --- Iterate over experiment conditions ---
    % i_dual: 1=Alternating (ihMT, preserves dipolar), 2=Continuous (MT, saturates dipolar)
    % i_tog:  1=Toggle ON (Inversion prep), 2=Toggle OFF

    for i_dual = 1:2
        isdual = (i_dual == 1);

        for i_tog = 1:2
            isToggle = (i_tog == 1);

            % Initialize Magnetization [Mx; My; Mza; Mzb; Mzd]
            Mall = [0; 0; M0a; M0b; 0];

            % --- Step 1: Toggle (180 Inversion Pulse) ---
            if isToggle && pw_rf_180 > 0
                t = pw_rf_180;
                w1_pulse = a_rf_inh;
                rfmt = RF_MT(T2b, w1_pulse, 0, 'SuperLorentzian');
                % Solve Bloch-McConnell-Provotorov equations
                Mall = BMP_solution_ihMT(t, Mall, R1a, R2a, b0, w1_pulse, 0, kab, kba, M0a, M0b, T2b, R1b, T1d, rfmt);
            end

            % --- Step 2: Crusher / Wait Time ---
            if wait_time > 0
                rfmt_wait = RF_MT(T2b, 0, 0, 'SuperLorentzian');
                Mall(1) = 0; % Hard crusher (destroy transverse mag)
                Mall(2) = 0;
                Mall = BMP_solution_ihMT(wait_time, Mall, R1a, R2a, b0, 0, 0, kab, kba, M0a, M0b, T2b, R1b, T1d, rfmt_wait);
            end

            % --- Step 3: Tip Down (Align to effective field) ---
            if pw_rf > 0
                t = pw_rf;
                w1_pulse = a_rf_fx_inh;
                rfmt = RF_MT(T2b, w1_pulse, 0, 'SuperLorentzian');
                Mall = BMP_solution_ihMT(t, Mall, R1a, R2a, b0, w1_pulse, 0, kab, kba, M0a, M0b, T2b, R1b, T1d, rfmt);
            end

            % --- Step 4: Spin Lock (Core Sequence) ---
            if eTSL_dual ~= 0
                w2_pos = curr_SL_val * b1_scale;
                w2_neg = -curr_SL_val * b1_scale;
                dw_pos = curr_freq_offset + b0;
                dw_neg = -curr_freq_offset + b0;

                % Pre-calculate absorption lineshapes
                rfmt_sl_pos = RF_MT(T2b, abs(w2_pos), dw_pos, 'SuperLorentzian');
                rfmt_sl_neg = RF_MT(T2b, abs(w2_neg), dw_neg, 'SuperLorentzian');

                if isdual
                    % Dual frequency: Alternating pos/neg offset preserves Dipolar Order
                    t_half = eTSL_dual./2;
                    for k = 1:N_eTSL_dual
                        Mall = BMP_solution_ihMT(t_half, Mall, R1a, R2a, dw_pos, 0, w2_pos, kab, kba, M0a, M0b, T2b, R1b, T1d, rfmt_sl_pos);
                        Mall = BMP_solution_ihMT(t_half, Mall, R1a, R2a, dw_neg, 0, w2_neg, kab, kba, M0a, M0b, T2b, R1b, T1d, rfmt_sl_neg);
                    end
                else
                    % Continuous/Single: Long duration saturates Dipolar Order
                    t_half = eTSL_single./2;
                    for k = 1:N_eTSL_single
                        Mall = BMP_solution_ihMT(t_half, Mall, R1a, R2a, dw_pos, 0, w2_pos, kab, kba, M0a, M0b, T2b, R1b, T1d, rfmt_sl_pos);
                        Mall = BMP_solution_ihMT(t_half, Mall, R1a, R2a, dw_neg, 0, w2_neg, kab, kba, M0a, M0b, T2b, R1b, T1d, rfmt_sl_neg);
                    end
                end
            end

            % --- Step 5: Tip Up (Return to Z-axis) ---
            if pw_rf > 0
                t = pw_rf;
                w1_pulse = -a_rf_fx_inh; % Negative amplitude to tip back
                rfmt = RF_MT(T2b, abs(w1_pulse), 0, 'SuperLorentzian');
                Mall = BMP_solution_ihMT(t, Mall, R1a, R2a, b0, w1_pulse, 0, kab, kba, M0a, M0b, T2b, R1b, T1d, rfmt);
            end

            % Store Result (Z-magnetization of pool A)
            temp_MzA(i_n, i_dual, i_tog) = Mall(3);

        end % End Toggle Loop
    end % End Dual Loop
end % End Scale Loop

MzA = temp_MzA;

%% 6. Calculate R_DOSL

% Extract magnetizations for calculation
Md_1 = MzA(1, 1, 2); % Dual 
Ms_1 = MzA(1, 2, 2); % Single
Md_2 = MzA(2, 1, 2); % Dual

if exist('SNR', 'var') && ~isempty(SNR) && SNR > 0 && ~isinf(SNR)
sigma = abs(Md_1 ./ SNR);
% Generate independent Gaussian noise for each measurement
noise1 = sigma * randn(1,1);
noise2 = sigma * randn(1,1);
noise3 = sigma * randn(1,1);

Md_1 = MzA(1, 1, 2) + noise1 ; % Dual
Ms_1 = MzA(1, 2, 2) +noise2; % Single
Md_2 = MzA(2, 1, 2) +noise3; % Dual

end

% Calculate R1rho difference (isolating dipolar contribution)
Rdosl1 = abs(-log(Md_1 ./ Ms_1) ./ TSL);
Rdosl2 = abs(-log(Md_2 ./ Ms_1) ./ TSL);

% Final Ratio
RATIO = squeeze(Rdosl1 ./ Rdosl2);

end