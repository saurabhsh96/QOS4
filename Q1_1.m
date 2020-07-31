%Fourier optics Q1_1
clear;
%% Defining inputs

%Freq of operation
freq = 100e9;

%Speed of light in vaccum
c = 3e8;

%Wavelength
lam = c/freq;

%Wavenumber
k0 = 2*pi/lam;

%Feed diam
Df = 6*lam;
%Rad of feed
Rf = Df / 2;

%Antenna Diam
Dr = 50*lam;
%Antenna Rad
Rr = Dr/2;

%Antenna F number
fNum = 3;
%Focal length
fL = fNum*Dr;

%Amplitude of incident light
E0 = 1;

%Incident angle
%Theta0 = 0;

%Defining current specifications
Jf = [0; 1; 0]; %Orientation along Y

%Defining free space impedance
zeta0 = 377;

%% Field calculations Q1.1 and Q1.2
%Defining meshgrid
drad = pi/180;
% dph = drad/2;
% dth = drad/2;

%Defining rho
% dRho = 0.00001;
% rho_arr1 = eps:dRho:Rr;
theta0 = 2*atan(Dr/(4*fL));
% [rho, ph] = meshgrid(rho_arr1, eps:dph:2*pi);
% th = 2*atan(rho./(2*fL));
%r_dash = fL.*(1+(tan(th/2)).^2);
theta = linspace(eps, theta0, 500);
phi_req = linspace(eps, 2*pi, 500);
[th, ph] = meshgrid(theta, phi_req);
dth = th(1,2) - th(1,1);
dph = ph(2,1) - ph(1,1);

%Calculating Vgo
[Vgoth0, Vgoph0, Egoth0, Egoph0] = GOField(E0, th, ph, fL, k0);
EgoMag = sqrt(abs(Egoth0).^2 + abs(Egoph0).^2);
EgoMax = max(max(EgoMag));

%Calculating Vtx
[Vath0, Vaph0, Efth0, Efph0, Prad] = FeedField(th, ph, fL, k0, Jf, Rf);
EfoMag = sqrt(abs(Efth0).^2 + abs(Efph0).^2);
EfoMax = max(max(EfoMag));

%Just to get the proper Prad calculations
thetaF = linspace(eps, pi/2-drad, 360);
phi_reqF = linspace(eps, 2*pi, 360);
[thF, phF] = meshgrid(thetaF, phi_reqF);
[Vath0_, Vaph0_, Efth0_, Efph0_, Prad2] = FeedField(thF, phF, fL, k0, Jf, Rf);


% Plotting GO Fields Q1_1 and Q1_2
%Plotting Emag and Theta
figure(1);
plot([-th(1,size(th, 2):-1:1)*180/pi, th(1,:)*180/pi], [mag2db(EgoMag(1,size(th, 2):-1:1)/EgoMax), mag2db(EgoMag(1,:)/EgoMax)], 'LineWidth', 2); hold on;
%plot(th(1,:)*180/pi, mag2db(EfoMag(1,:)/EfoMax), 'LineWidth', 2);
plot([-th(1,size(th, 2):-1:1)*180/pi, th(1,:)*180/pi], [mag2db(EfoMag(1,size(th, 2):-1:1)/EfoMax), mag2db(EfoMag(1,:)/EfoMax)], 'LineWidth', 2);
title('Magnitude of E-field (\phi = 0 [deg])');
xlabel('\theta(in deg)');
ylabel('Normalized E-field [dBV]');
legend('Normalized GO Field', 'Normalized Feed Far Field');
ylim([-6, 0]);
grid on;

%Phi = 90
figure(2);
plot([-th(1,size(th, 2):-1:1)*180/pi, th(1,:)*180/pi], [mag2db(EgoMag(91,size(th, 2):-1:1)/EgoMax), mag2db(EgoMag(91,:)/EgoMax)], 'LineWidth', 1.5); hold on;
%plot(th(1,:)*180/pi, mag2db(EfoMag(1,:)/EfoMax), 'LineWidth', 2);
plot([-th(1,size(th, 2):-1:1)*180/pi, th(1,:)*180/pi], [mag2db(EfoMag(91,size(th, 2):-1:1)/EfoMax), mag2db(EfoMag(91,:)/EfoMax)], 'LineWidth', 1.5);
title('Magnitude of E-field \phi = 90 [deg]');
xlabel('\theta (in deg)');
ylabel('Normalized E-field [dBV]');
ylim([-6, 0]);
legend('Normalized GO Field', 'Normalized Feed Far Field');
grid on;

%% Q1.3 
%Calculating Prx and Vgx
[Prx, etaA] = PowerRec(Vath0, Vaph0, Vgoth0, Vgoph0, Prad2, zeta0, th, ph, Rr, E0);

%Voltages
%Vatx0Mag = sqrt((abs(Vath0)).^2 + (abs(Vaph0)).^2);
%Vgo0Mag = sqrt((abs(Vgoth0)).^2 + (abs(Vgoph0)).^2);
%% Q1.4 Changing the Theta(PW) range 0:6 deg

%Defining meshgrid for incident angles
thi_arr = linspace(-6*drad, 6*drad, 160);
phi_arr = linspace(eps, 2*pi, 4);
[thi, phi] = meshgrid(thi_arr, phi_arr);
dthi = thi(1,2) - thi(1,1);
dphi = phi(2,1) - phi(1,1);
% dthi = drad/4;
% dphi = drad/4;
% [thi, phi] = meshgrid(eps:dthi:6*drad, eps:dphi:2*pi);

%Loop over fields
delTh = (1 - cos(th))./(1+cos(th));
%Defining electric fields for different incident angles (also defining
%voltages)

EgothDel = zeros([size(thi), size(Egoth0)]);
EgophDel = zeros([size(thi), size(Egoph0)]);
VgothDel = zeros([size(thi), size(Vgoth0)]);
VgophDel = zeros([size(thi), size(Vgoph0)]);


%Recieved power calculations
PrxR = zeros(size(thi));
etaAR = zeros(size(thi));

for indPh = 1:size(phi, 1)
    for ind = 1:size(thi, 2)
        %Defining Krho
        kRhox = k0.*sin(th).*cos(ph);
        kRhoy = k0.*sin(th).*sin(ph);
        
        %Defining DelK
        delKx = k0.*sin(thi(indPh, ind)).*cos(phi(indPh, ind));
        delKy = k0.*sin(thi(indPh, ind)).*sin(phi(indPh, ind));
        
        %Define DelRho
        delRhox = (fL/k0).*delKx;
        delRhoy = (fL/k0).*delKy;
        
        %Dot products
        kDel = delRhox.*kRhox + delRhoy.*kRhoy;
        
        %Exponent Term
        expTerm1 = exp(-1j.*kDel);
        expTerm2 = exp(-1j.*kDel.*delTh);
        
        %Calculation of the fields
        EgothDel(indPh, ind, :, :) = Egoth0.*expTerm1.*expTerm2;
        EgophDel(indPh, ind, :, :) = Egoph0.*expTerm1.*expTerm2;
        VgothDel(indPh, ind, :, :) = Vgoth0.*expTerm1.*expTerm2;
        VgophDel(indPh, ind, :, :) = Vgoph0.*expTerm1.*expTerm2;
        
        [PrxR(indPh, ind), etaAR(indPh, ind)] = PowerRec(Vath0, Vaph0, ...
            squeeze(VgothDel(indPh,ind,:,:)), ...
            squeeze(VgophDel(indPh,ind,:,:)), Prad2, zeta0, th, ph, Rr, E0);
    end
end

%% Q1.5; Power patterns
% using the same Prad as above

%Tx calculations
EFRMag = Q2(zeta0, freq, Df, Jf, fL, Dr, thi_arr, thi, phi);
Ptx = (abs(EFRMag).^2)./(2*zeta0);
PtxMax = max(max(Ptx));

% for indPh = 1:size(thi, 1)
%     for ind = 1:size(thi, 2)
%         [PrxR(indPh, ind), etaAR(indPh, ind)] = PowerRec(Vath0, Vaph0, ...
%             squeeze(VgothDel(indPh,ind,:,:)), ...
%             squeeze(VgophDel(indPh,ind,:,:)), Prad2, zeta0, th, ph, Rr, E0);
%     end
% end

%Plotting 
%phi = 0
PrxRMax0 = (max(PrxR(1,:)));
figure();
plot([-thi(1,size(thi, 2):-1:1)/drad, thi(1,:)/drad], [pow2db(PrxR(1, size(thi, 2):-1:1)./PrxRMax0), pow2db(PrxR(1,:)./PrxRMax0)], 'LineWidth', 2.0);
hold on;
plot([-thi(1,size(thi, 2):-1:1)/drad, thi(1,:)/drad], [pow2db(Ptx(1, size(thi, 2):-1:1)./PtxMax), pow2db(Ptx(1,:)./PtxMax)], 'LineWidth', 1.5);
title('Recieved Power at feed and Radiation Pattern with Tx approach(\phi_i = 0)');
xlabel('Incident Angle \theta_i [deg]');
ylabel('Normalized Power [dB]');
legend('FO-Rx', 'Tx');
ylim([-42, 0]);
grid on;

%phi = 90
PrxRMax90 = (max(PrxR(2,:)));
figure();
plot([-thi(1,size(thi, 2):-1:1)/drad, thi(1,:)/drad], [pow2db(PrxR(2, size(thi, 2):-1:1)./PrxRMax90), pow2db(PrxR(2,:)./PrxRMax90)], 'LineWidth', 1.5);
hold on;
plot([-thi(1,size(thi, 2):-1:1)/drad, thi(1,:)/drad], [pow2db(Ptx(2, size(thi, 2):-1:1)./PtxMax), pow2db(Ptx(2,:)./PtxMax)], 'LineWidth', 1.5);
title('Recieved Power at feed and Radiation Pattern with Tx approach(\phi_i = 90)');
xlabel('Incident Angle \theta_i [deg]');
ylabel('Normalized Power [dB]');
legend('FO-Rx', 'Tx');
ylim([-42, 0])
grid on;

%% Q1.6 Directivity comparisons

%Rx Approach
%Directivity
denom = sum(PrxR.*sin(thi), 'all').*dthi.*dphi;
Dir = 4.*pi.*PrxR./(denom);

Dmax = max(max(Dir));

Gain = etaAR.*Dir;
Gmax = max(max(Gain));

dRho = 0.0001;
rho_dash1R = eps:dRho:Rr;
theta0R = 2*atan(Rr*2/(4*fL));
[rho_dashR, phi_dashR] = meshgrid(rho_dash1R, eps:drad:2*pi);
th_dashR = 2*atan(rho_dashR./(2*fL));
r_dashR = fL.*(1+(tan(th_dashR/2)).^2);

[taperS, etaS, Area] = Spillover(freq, 1, Jf, Rf, 1000*lam, ...
            th_dashR, phi_dashR, thF, phF, rho_dashR, fL, Rr);
etaATx = taperS.*etaS;

%Maximum directiviy
DirTxM = (4*pi/lam^2).*Area;

%Directivity
DirTx = DirTxM.*taperS;

%Gain
GTx = DirTxM.*etaATx;

%% Part 2: Q2.1
dx = Df;

kx = k0.*sin(th).*cos(ph);

VthDx =  Vath0.*exp(1j*kx*dx);
VphDx =  Vaph0.*exp(1j*kx*dx);

%Recieved power calculations
PrxRDx = zeros(size(thi));
etaARDx = zeros(size(thi));

for indPh = 1:size(thi, 1)
    for ind = 1:size(thi, 2)
        [PrxRDx(indPh, ind), etaARDx(indPh, ind)] = PowerRec(VthDx, VphDx, squeeze(VgothDel(indPh,ind,:,:)), squeeze(VgophDel(indPh,ind,:,:)), Prad2, zeta0, th, ph, Rr, E0);
    end
end

PrxRDMax0 = (max(PrxRDx(1,:)));

figure();
%Uses theta from -6 to 6
plot(thi(1,:)/drad, pow2db(PrxR(1,:)./PrxRMax0), 'LineWidth', 1.5); hold on;
plot(thi(1,:)/drad, pow2db(PrxRDx(1,:)./PrxRDMax0), 'LineWidth', 1.5);

title('Recieved Power (\phi_i = 0 [deg] plane)');
xlabel('Incident Angle \theta_i [deg]');
ylabel('Normalized Power [dB]');
legend('No Displacement', '6 \lambda Displacement');
ylim([-45, 0])
grid on;

%% Q2.2 Effect on efficiency
Df_dash = [2*lam*fNum, lam*fNum, 0.5*lam*fNum];
dx_dash = Df_dash;

Prx_dash = zeros([size(Df_dash, 1), size(thi)]);
etaA_dash = zeros([size(Df_dash, 1), size(thi)]);
Prad_dash = zeros(size(Df_dash, 1));

names = ["Df=dx=2\lambda_0f#", "Df=dx=\lambda_0f#", ...
    "Df=dx=0.5\lambda_0f#"];
figure();
for indU = 1:size(dx_dash, 2)
    [VathTemp, VaphTemp, EfthTemp, EfphTemp, PradTemp] = FeedField(th, ph, ...
        fL, k0, Jf, Df_dash(indU)/2);
    [Vath0_, Vaph0_, Efth0_, Efph0_, Prad_dash(indU)] = FeedField(thF, phF, fL, k0, Jf, Df_dash(indU)/2);
    VthDxT =  VathTemp.*exp(1j*kx*dx_dash(indU));
    VphDxT =  VaphTemp.*exp(1j*kx*dx_dash(indU));
    for indPh = 1:size(thi, 1)
        for ind = 1:size(thi, 2)
            [Prx_dash(indU, indPh, ind), etaA_dash(indU, indPh, ind)] = ...
                PowerRec(VthDxT, VphDxT, squeeze(VgothDel(indPh,ind,:,:)), ...
                squeeze(VgophDel(indPh,ind,:,:)), Prad_dash(indU), zeta0, th, ph, Rr, E0);
        end
    end
    plot(thi(1,:)/drad, pow2db(squeeze(Prx_dash(indU,1,:)) ...
        ./max(max(Prx_dash(indU,:,:)))), 'DisplayName', names(indU), ...
        'LineWidth', 1.5);
    hold on;
    
    title('Recieved Power (\phi_i = 0 [deg] plane)');
    xlabel('Incident Angle \theta_i [deg]');
    ylabel('Normalized Power [dB]');
    ylim([-45, 0])
    grid on;
end
hold off;
legend show;