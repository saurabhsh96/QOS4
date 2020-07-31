function [Prx, etaA] = PowerRec(Vath0, Vaph0, Vgoth0, Vgoph0, Prad2, zeta0 ...
    ,th, ph, Rr, E0)
    %Elementary 
    dth = th(1,2) - th(1,1);
    dph = ph(2,1) - ph(1,1);

    %Vatx0Mag.*Vgo0Mag
    Vgx = (2/zeta0).*sum((Vath0.*Vgoth0 + Vaph0.*Vgoph0).*sin(th), 'all').*dth.*dph;
    %Prad1 = (1/(2*zeta0))*sum(((abs(Vatx0Mag).^2).*sin(th)), 'all').*dth.*dph;

    %Calculating Recieved power
    Prx = ((abs(Vgx)).^2)./(16*Prad2);

    %Calculating aperture efficiency
    etaA = Prx./(((abs(E0)).^2)./(2*zeta0).*pi.*Rr.^2);
end 