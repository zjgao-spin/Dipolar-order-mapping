function RATIO_dosl_approx = cal_RATIO_dosl_approx(R, T2b, T1d, dw, w1, B1, B0)

    w1 = w1*B1;
    dw = dw+B0.*2*pi;
    D = 1/sqrt(15)/T2b;
    N_scale = 5;

    Rrfb1 = RF_MT(T2b, w1, dw, 'SuperLorentzian');
    Rrfb2 = RF_MT(T2b, w1/N_scale, dw/N_scale, 'SuperLorentzian');
    
    numerator = (R + Rrfb2) * Rrfb1^2 * T1d * (dw/D)^2;
    denominator = (R + Rrfb1) * (Rrfb2 * (Rrfb1 * T1d * (dw/D)^2 + 1) - Rrfb1);
    
    RATIO_dosl_approx = abs(numerator / denominator);

end