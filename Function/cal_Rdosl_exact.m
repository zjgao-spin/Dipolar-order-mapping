function Rdosl_exact = cal_Rdosl_exact(R1a, R2a, R1b, MPF, R, T2b, T1d, dw, w1)
    M0b = MPF;
    M0a = 1-M0b;
    fb = M0b/M0a;
    kba = R*M0a;
    kab = kba*fb;
    D = 1/sqrt(15)/T2b;
    N_scale = 5;
   
    R1rho = zeros(2,2);
    
    for i = 1:2
        if i == 1
            curr_dw = dw;
            curr_w1 = w1;
        else
            curr_dw = dw/N_scale;
            curr_w1 = w1/N_scale;
        end
        
        % Note: Assuming RF_MT function exists in path
        Rrfb = RF_MT(T2b, curr_w1, curr_dw, 'SuperLorentzian');
        
        % 4x4 Matrix (No Dipolar Order)
        A = [-R2a, curr_dw, 0, 0;
             -curr_dw, -R2a, curr_w1, 0;
             0, -curr_w1, -R1a-kab, kba;
             0, 0, kab, -R1b-Rrfb-kba];
        
        % 5x5 Matrix (With Dipolar Order)
        AD = [-R2a, curr_dw, 0, 0, 0;
              -curr_dw, -R2a, curr_w1, 0, 0;
              0, -curr_w1, -R1a-kab, kba, 0;
              0, 0, kab, -R1b-Rrfb-kba, Rrfb*curr_dw;
              0, 0, 0, Rrfb*(curr_dw/D/D), -(1/T1d+Rrfb*(curr_dw/D)^2)];
          
        R1rho(i,1) = min(abs(real(eig(A))));
        R1rho(i,2) = min(abs(real(eig(AD))));
    end
    
    DeltaR1rho1 = R1rho(1,1) - R1rho(1,2);
    DeltaR1rho2 = R1rho(2,1) - R1rho(1,2);
    
    if abs(DeltaR1rho2) < 1e-10
        Rdosl_exact = NaN;
    else
        Rdosl_exact = abs(DeltaR1rho1/DeltaR1rho2);
    end
end