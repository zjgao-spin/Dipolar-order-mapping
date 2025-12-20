function rfmt=RF_MT(T2c,w1,dw,lineshape)

if strcmp(lineshape,'SuperLorentzian') %SuperLorentzian
%     step=1/10000;
%     sup_theta = 0:step:pi/2;
    
    cutoff=10;
    G=zeros(1,length(dw));
    for j=1:length(dw)
        if abs(dw(j)) >= cutoff  %% the superlorentz has a pole: this avoids infinities
            % see Morrsion and Henkelman 1995.  units = s.  Seems weird to
            % me.  Need to multiply by w1^2 * pi to get saturation rate.
            du=.0001;
            u=0:du:1;
            integrand2=sqrt(2/pi)*T2c./abs(3*u.^2-1) .* exp(-2*(dw(j)*T2c./abs(3*u.^2-1)).^2);
            G(j)=w1.^2.*pi.*sum(integrand2)*du;
        else
            X = [-1.1*cutoff -cutoff cutoff 1.1*cutoff];
            Y = RF_MT(T2c,w1,X,lineshape);
            G(j)=interp1(X,Y,dw(j),'spline');
        end
        rfmt=G;
        
    end

elseif strcmp(lineshape,'Gaussian')
    rfmt=w1.^2*T2c.*sqrt(pi/2).*exp(-(dw.*T2c).^2./2); %Gaussian
else
        rfmt=w1.^2*T2c./(1+(dw.*T2c).^2); %Lorentzian
end

end
