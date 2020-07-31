clear;

%% Defining Input parameters
%Zeta
zeta = 377;

%Freq of operation
freq = 100e9;

%Speed of light
c = 3e8;

%Wavelength
lam = c/freq;

%Wavenumber
k0 = 2*pi/lam;

%Diameter of the feed
Df = 6*lam;

%Radius of feed
af = Df/2;

%Defining the current distribution of feed
jf = [0 1 0]'; %Only across X why? Because X polarized 

%Jf i.e. Magnitude of current
Jf = 1;

%Diameter of Parabola
Dr = 50*lam;

%F number
fhash = 3;

%Focal length
fref = fhash*Dr;

%Plane wave amplitude
Epw = 1; %V/m

%Plane wave incident angle
theta_pw = 0;

%% 1.1 and 1.2 Routine to evaluate the GO voltage and far field from feed to circular sphere
%Defining the meshgrid
drad = pi/360;
dth = 0.0001;

R = fref; %Taking the radius of equivalent sphere to be equal to focal length
theta0 = 2*atan(Dr/(4*R)); %Max angle theta0

%rho_var = eps:drho:Dr/2;
phi_var = eps:drad:2*pi;
th_var = eps:dth:theta0;
[th,ph] = meshgrid(th_var,phi_var);
%th = 2*atan(rho./(2*fref));

%Calculating GO field at incident angle of theta = 0
Spar = -2./(1+cos(th));

erGO_th_0 = Spar.*Epw.*sin(ph); 
erGO_ph_0 = Spar.*Epw.*cos(ph);
erGO_mag_0 = sqrt(erGO_th_0.^2 + erGO_ph_0.^2);
erGO_mag_max = max(max(erGO_mag_0));

%Calculating Eff of the circular feed over same FO Sphere:
Af = pi*(af)^2; %Area of the feed 
x = k0*sin(th)*af;
J1 = besselj(1,x); %Bessel function
Airy = (Af*J1)./x;

Vatx_ph = Airy.*cos(ph); %Tx Volatge
Vatx_th = Airy.*cos(th).*sin(ph);
Vatx_mag = sqrt(Vatx_ph.^2 + Vatx_th.^2);

Efar_ph = Vatx_th.*exp(-1j*k0*R)./R; %Far-field computation
Efar_th = Vatx_ph.*exp(-1j*k0*R)./R;
Efar_mag = sqrt(abs(Efar_ph).^2 + abs(Efar_th).^2); %Magnitude
Efar_mag_max = max(max(Efar_mag)); % Max value

%Plotting GO field and Eff wrt Theta at Phi = 0
figure(1);
plot([-th(1,1663:-1:1)*180/pi,th(1,1:1:1663)*180/pi], [mag2db(erGO_mag_0(1,1663:-1:1)/erGO_mag_max),mag2db(erGO_mag_0(1,1:1:1663)/erGO_mag_max)],'LineWidth',2); 
hold on;
plot([-th(1,1663:-1:1)*180/pi,th(1,1:1:1663)*180/pi], [mag2db(Efar_mag(1,1663:-1:1)/Efar_mag_max),mag2db(Efar_mag(1,1:1:1663)/Efar_mag_max)],'LineWidth',2); 
title('GO Electric Field and EFar field at feed at \phi = 0');
xlabel('\theta in [degrees]');
ylabel('normalized |E_F_a_r|, normalized |Er_G_O| in [dB]');
legend('Er_G_O field','E_F_a_r Field');
xlim([-9,9]);
ylim([-5,0]);
set(gca,'ytick',[-5:0.5:0]);

%Plotting GO field and Eff wrt Theta at Phi = 90
figure(2);
plot([-th(1,1663:-1:1)*180/pi,th(1,1:1:1663)*180/pi], [mag2db(erGO_mag_0(91,1663:-1:1)/erGO_mag_max),mag2db(erGO_mag_0(91,1:1:1663)/erGO_mag_max)],'LineWidth',2); 
hold on;
plot([-th(1,1663:-1:1)*180/pi,th(1,1:1:1663)*180/pi], [mag2db(Efar_mag(91,1663:-1:1)/Efar_mag_max),mag2db(Efar_mag(91,1:1:1663)/Efar_mag_max)],'LineWidth',2); 
title('GO Electric Field and EFar field at feed at \phi = 90');
xlabel('\theta in [degrees]');
ylabel('normalized |E_F_a_r|, normalized |Er_G_O| in [dB]');
legend('Er_G_O field','E_F_a_r Field');
xlim([-9,9]);
ylim([-5,0]);
set(gca,'ytick',[-5:0.5:0]);

%% 1.3 power at input terminal of feed - impedance match. Get aperture efficiency.

VrGO_mag = erGO_mag_0*R./exp(1j*k0*R);

%Computing Vg
var1 = (sum(Vatx_mag.*VrGO_mag.*sin(th), 'all'))*drad*(th(1,2)-th(1,1));
Vg = 2*var1/zeta;

%Computing P radiated
drad_feed = pi/360;
theta_feed = eps:drad_feed:pi/2;
phi_feed = eps:drad_feed:2*pi;
[th_feed,ph_feed] = meshgrid(theta_feed,phi_feed);

%Vatx for Prad computation
x1 = k0*sin(th_feed)*af;
J1_feed = besselj(1,x1); %Bessel function
Airy_feed = (Af*J1_feed)./x1;
Vatx_ph_feed = Airy_feed.*cos(ph_feed); %Tx Volatge
Vatx_th_feed = Airy_feed.*cos(th_feed).*sin(ph_feed);
Vatx_mag_feed = sqrt(Vatx_ph_feed.^2 + Vatx_th_feed.^2);

var2 = (sum((abs(Vatx_mag_feed).^2).*sin(th_feed), 'all'))*drad_feed*(th_feed(1,2)-th_feed(1,1));
Prad = var2/(2*zeta);

%Computing Prx
Prx = abs(Vg)^2/(16*Prad);

%Computing aperture efficiency
Ar = pi*(Dr/2)^2;
Eta_app = (Prx*2*zeta)/((abs(Epw))^2*Ar);

%% Incident plane wave angle theta = [0:6] degrees
dth_pw = 0.001;
drad_2 = pi/4;
theta_pw1 = eps:dth_pw:6*pi/180; %Incidentplane wave in theta
phi_pw = eps:drad_2:2*pi; %Incident plane wave in phi

[th_pw,ph_pw] = meshgrid(theta_pw1,phi_pw); %th_pw goes in col ph_pw goes in rows

deln_th = (1-cos(th))./(1+cos(th));

var1_loop = zeros(size(th_pw));
Vg_loop = zeros(size(th_pw));
Prx_loop = zeros(size(th_pw));


for m = 1: size(ph_pw,1) % looping over rows
    for k = 1:size(th_pw,2) % looping over cols
        krho_x = k0*sin(th).*cos(ph); % krho vector in x
        krho_y = k0*sin(th).*sin(ph); % krho vector in y

        del_krhoi_x = k0*sin(th_pw(m,k)).*cos(ph_pw(m,k));
        del_krhoi_y = k0*sin(th_pw(m,k)).*sin(ph_pw(m,k));
        %del_krhoi = sqrt(del_krhoi_x.^2 + del_krhoi_y.^2);

        del_rhofp_x = (R.*del_krhoi_x)./k0;
        del_rhofp_y = (R.*del_krhoi_y)./k0;

        krho_delrho_x = krho_x.*del_rhofp_x;
        krho_delrho_y = krho_y.*del_rhofp_y;

        e1 = exp(-1j*(krho_delrho_x + krho_delrho_y));
        e2 = exp(-1j*(krho_delrho_x + krho_delrho_y).*deln_th);
        erGO_th = erGO_th_0.*e1.*e2;
        erGO_ph = erGO_ph_0.*e1.*e2;

        vrGo_th = erGO_th.*R./exp(1j*k0*R);
        vrGO_ph = erGO_ph.*R./exp(1j*k0*R);

        erGO_mag_loop = sqrt(abs(erGO_th).^2 + abs(erGO_ph).^2);

        VrGO_mag_loop = erGO_mag_loop*R./exp(1j*k0*R);

        varX =  (vrGo_th.*Vatx_th + vrGO_ph.*Vatx_ph);
        %Computing Vg
        var1_loop(m,k) = (sum(varX.*sin(th), 'all')).*drad.*(th(1,2)-th(1,1));
        Vg_loop(m,k)= 2*var1_loop(m,k)/zeta;

        %Computing Prx
        Prx_loop(m,k) = (abs(Vg_loop(m,k)))^2/(16*Prad);
    end
end
erGO_mag_loop_max = max(max(erGO_mag_loop));
Prx_loop_max = max(Prx_loop(1,:));

% figure(3)
% plot([-th(1,1663:-1:1)*180/pi,th(1,1:1:1663)*180/pi], [mag2db(erGO_mag_loop(1,1663:-1:1)/erGO_mag_loop_max),mag2db(erGO_mag_loop(1,1:1:1663)/erGO_mag_loop_max)],'LineWidth',2); 
% title('GO Electric Field at \phi = 0');
% xlabel('\theta in [degrees]');
% ylabel('normalized |Er_G_O| in [dB]');
% legend('Er_G_O field');
% set(gca,'ytick',[-5:0.5:0]);

%% Question 1.5 : power received in dB
% Computing for Antenna in Tx:
%[th_obs,EFRMag,EFRMagMax] = antenna_tx(zeta,freq,Df,fref,Dr);

figure(4);
plot([-th_pw(1,size(th_pw, 2):-1:1)*180/pi,th_pw(1,:)*180/pi], [pow2db(Prx_loop(1,size(th_pw, 2):-1:1)/Prx_loop_max),pow2db(Prx_loop(1,:)/Prx_loop_max)],'LineWidth',2); 
title('Received Power in dB by reflector feed as a function of incident angle');
xlabel('\theta_i in [degrees]');
ylabel('Received Power P_r_x in [dB]');

%figure(5);
%plot(th_obs(1,:)*180/pi, mag2db(EFRMag(1,:)/EFRMagMax)); 
% title('Normalized E-field phi = 0');
% xlabel('Normalized E-field [dBV]');
% ylabel('Observation angle, Theta (in deg)');
% ylim([-35, 0]);
% legend('Far-field of antenna in Tx');
% 
% figure();
% plot(th_obs(1,:)*180/pi, mag2db(EFRMag(3,:)/EFRMagMax)); 
% title('Normalized E-field phi = 90');
% xlabel('Normalized E-field [dBV]');
% ylabel('Observation angle, Theta (in deg)');
% ylim([-35, 0]);
% legend('Far-field of antenna in Tx');