%GO Field Calculation
function [Vgoth0, Vgoph0, Egoth0, Egoph0] = GOField(E0, th, ph, fL, k0)
    %Calculating the GOField
    constant_term1 = -2./(1+cos(th)).*E0;
    % Y-pol
    %     Egoth0 = constant_term1.*sin(ph);
    %     Egoph0 = constant_term1.*cos(ph);
    % X pol
    Egoth0 = constant_term1.*cos(ph);
    Egoph0 = - constant_term1.*sin(ph);

    %Calculating GO Voltage
    constant_term2 = exp(1j*k0*fL)./fL; %R = fL
    Vgoth0 = Egoth0./constant_term2;
    Vgoph0 = Egoph0./constant_term2;
end