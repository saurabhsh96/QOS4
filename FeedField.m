%Function to calculate the feedfield
function [Vath, Vaph, Efth, Efph, Prad] = FeedField(th, ph, fL, k0, Jf, Rf)
    %Get the airy pattern
    %[th, ph] = meshgrid(eps:pi/180:pi/2, eps:pi/180:2*pi);
    J = JFT(k0, Rf, th, Jf);
    
    %Constant term, phase and amp
    constant_term = exp(-1j*k0*fL)/fL;
    
    %Voltage calculations
    Vath = squeeze(J(2,:,:)).*cos(th).*sin(ph);
    Vaph = squeeze(J(2,:,:)).*cos(ph);
    
    %Field calculations
    Efth = Vath.*constant_term;
    Efph = Vaph.*constant_term;
    EfMag = sqrt((abs(Efth)).^2 + (abs(Efph)).^2);
    
    %Prad calculations
    Intensity = fL^2.*abs(EfMag).^2./(2*377); %Assuming Z0 = 377

    %Calculate prad
    Prad = (sum(Intensity.*sin(th), 'all'))*(th(1,2)-th(1,1))*(ph(2,1)-ph(1,1));
end