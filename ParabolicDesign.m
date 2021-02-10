close all
clear
clc

%% Variables

%Frequency
ftx = [1020 1040]; %MHz
frx = [1080 1100]; %MHz

%BandWidth -3dB in Degrees 
BWh = 7; 
BWv = [-5 15]; 

%% Preliminariy calculations

c = 3*(10^8); %m/s^2    Speed of light in vacuum 
wl = c/(fmin*10^6); %m  Wavelength                
k = ((2*pi)/wl); %rad/m Wavenumber

% Waveguide 0.97 to 1.45 GHz WR770

b = 0.9779; % m     Waveguide height
a = 0.19558;% m     Waveguide width

%% Horn dimesnions

A = a; % m      Horn dimension in the x axis
B = 0.9;% m     Horn dimension in the y axis
dx = A;% m      Distance between the two horns

Beta = (2*pi/wl)*sqrt(1-(wl/(2*A))^2);
s = 1/4;                                
R2 = (B^2)/(8*wl*s);

%% Parabolic reflector dimensions

% Widht of the reflector
widht_2W = 1.2*wl*180/(BWh*pi); %m 

% Height of the reflector 
height_2H = 1.17*wl*180/((BWv(2)-BWv(1)).*pi); %m 

% Focal distance
% It is stablish a relation F/D = 0.4
% We take D as the major dimension, in this case, Widht

F_D = 0.4  ; 
F = F_D*widht_2W; %m    Focal distance  

% Offset h
% We estimate twice the wavelength

h = 2*wl; %m

% View angles from the focal point
% Angles in degrees
detheta0 = (180/pi)*atan(4*F*((height_2H/2)+h)/(4*F^2-((height_2H/2)+h)^2)); %Offset angle 
dethetamin = (180/pi)*atan((4*F*h)/(4*(F^2) - (h^2))); %Minimum angle
dethetamax = (180/pi)*atan((4*F*(height_2H+h))/((4*(F^2))-(height_2H+h)^2)); %Maximum angle
dealpha = (180/pi)*atan((widht_2W/2*4*F)/(((4*F*(h+height_2H/2))^(2)+(4*F^2-(height_2H/2+h)^(2)-(widht_2W/2)^2)^(2))^(1/2))); %Widht Angle 

% Angles in radian
theta0 = atan(4*F*((height_2H/2)+h)/(4*F^2-((height_2H/2)+h)^2)); 
thetamin = atan((4*F*h)/(4*(F^2) - (h^2))); 
thetamax = atan((4*F*(height_2H+h))/((4*(F^2))-(height_2H+h)^2)); 
alpha = atan((widht_2W/2*4*F)/(((4*F*(h+height_2H/2))^(2)+(4*F^2-(height_2H/2+h)^(2)-(widht_2W/2)^2)^(2))^(1/2))); 

%% Calculation of the radiated field from the Horns

%Vectors 

[THETA,PHI] = meshgrid(linspace(-pi/2,pi/2,49),linspace(0,2*pi,49));
u=sin(THETA).*cos(PHI);
v=sin(THETA).*sin(PHI);

theta = linspace(-pi/2,pi/2,1000);
phi = linspace(0,2*pi,1000);

% Array factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FAsuma = exp(-1i.*((2*pi)/wl).*(dx/2).*sin(THETA).*cos(PHI))+exp(1i.*((2*pi)/wl).*(dx/2).*sin(THETA).*cos(PHI));
FAresta = exp(-1i.*((2*pi)/wl).*(dx/2).*sin(THETA).*cos(PHI))-exp(1i.*((2*pi)/wl).*(dx/2).*sin(THETA).*cos(PHI));

HEFAsuma = exp(-1i.*((2*pi)/wl).*(dx/2).*sin(theta))+exp(1i.*((2*pi)/wl).*(dx/2).*sin(theta));
HEFAresta = exp(-1i.*((2*pi)/wl).*(dx/2).*sin(theta))-exp(1i.*((2*pi)/wl).*(dx/2).*sin(theta));


figure('Name','Array Factor','NumberTitle','off')
set(gcf, 'WindowState', 'maximized');

% ----- Sum --------------------------------------------------------------

subplot(1,2,1)
plot(linspace(-90,90,1000), 20.*log10(abs(HEFAsuma)) - max(20.*log10(abs(HEFAsuma))), '-r');
grid on
set(gcf, 'PaperType', 'A4'); 
title('Sum diagram Array Factor')
xlabel('\theta [º]');
ylabel('Normalize Power [dB]');
axis([-90 90 -50 10]);

% ----- Difference --------------------------------------------------------

subplot(1,2,2)
plot(linspace(-90,90,1000), 20.*log10(abs(HEFAresta)) - max(20.*log10(abs(HEFAresta))), '-r');
grid on
set(gcf, 'PaperType', 'A4');
title('Difference diagram array factor')
xlabel('\theta [º]');
ylabel('Normalize Power [dB]');
axis([-90 90 -50 10]);

% Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ex = @(x,y) 0;
Ey = @(x,y) E0.*cos((pi.*x)/A).*exp(-1i*(Beta/2).*((y.^2)/R2));

[Px, Py, HEtheta, HEphi, EEtheta, EEphi] = funcion_DIMmejorada(A, B, Ex, Ey, k);  

Etheta = (1i.*k)/(2.*pi).*(Px.*cos(PHI)+Py.*sin(PHI)); 
Ephi = ((-1i.*k)/(2.*pi)).*(cos(THETA).*(Px.*sin(PHI)-Py.*cos(PHI))); 
Erad = sqrt(Etheta.^2+Ephi.^2); 
Esuma = Erad.*abs(FAsuma);
Eresta = Erad.*abs(FAresta);

%% Representations
% Sum diagram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Diagrama Suma','NumberTitle','off')
set(gcf, 'WindowState', 'maximized');

% ----- 3D ----------------------------------------------------------------

x1= u;
y1= v;
z1= 20*log10(abs(Esuma));
z1 = z1 - max(max(z1));
z1(z1 <= -30) = -30;
subplot(1,2,1)
surf(x1,y1,z1);
title('Sum diagram)');
xlabel('u [rad]');
ylabel('v [rad]');
zlabel('Total field');
axis([-1 1 -1 1 -30 0]);

c = colorbar;
c.Location = 'southoutside';
colormap parula
caxis([-30 0])


% ----- H plane -----------------------------------------------------------
% Phi = 0º

subplot(2,2,2)
plot(linspace(-pi/2,pi/2,1000)*(180/pi), 20.*log10(abs(HEFAsuma).*abs(sqrt(HEtheta.^2+HEphi.^2)))-max(20.*log10(abs(HEFAsuma).*abs(sqrt(HEtheta.^2+HEphi.^2)))), '-r');
set(gcf, 'PaperType', 'A4'); 
xlabel('[º]');
ylabel('Erad Sum [dB]');
axis([-90 90 -60 10]);
xline(-35.23, 'b');
xline(35.23, 'b');
hline=refline(0, -12); hline.Color='r'; hline.LineStyle='--'; hline.DisplayName;
legend('H plane', 'BW_{-12dB}');

% ----- E plane -----------------------------------------------------------
% Phi = 90º

subplot(2,2,4)
plot(linspace(-pi/2,pi/2,1000)*(180/pi), 20.*log10(abs(HEFAsuma).*abs(sqrt(EEtheta.^2+EEphi.^2)))-max(20.*log10(abs(HEFAsuma).*abs(sqrt(EEtheta.^2+EEphi.^2)))), '-r');
set(gcf, 'PaperType', 'A4'); 
xlabel('[º]');
ylabel('Erad Suma [dB]');
axis([-90 90 -60 10]);
xline(-13.96, 'b');
xline(13.96, 'b');
hline=refline(0, -10); hline.Color='r'; hline.LineStyle='--'; hline.DisplayName;
legend('E plane', 'BW_{-10dB}');


% Difference diagram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Difference diagram','NumberTitle','off')
set(gcf, 'WindowState', 'maximized');

% ----- 3D ----------------------------------------------------------------

subplot(1,2,1)
x1= sin(THETA).*cos(PHI);
y1= sin(THETA).*sin(PHI);
z1= 20*log10(abs(Eresta));
surf(x1,y1,z1 - max(max(z1)));
title('Difference diagram');
xlabel('u');
ylabel('v');
zlabel('Total field');
axis([-1 1 -1 1 -40 0]);

c = colorbar;
c.Location = 'southoutside';
colormap parula
caxis([-60 10])

% ----- Plano H -----------------------------------------------------------
% Phi = 0º

subplot(2,2,2)
plot(linspace(-pi/2,pi/2,1000)*(180/pi), 20.*log10(abs(HEFAresta).*abs(sqrt(HEtheta.^2+HEphi.^2)))-max( 20.*log10(abs(HEFAresta).*abs(sqrt(HEtheta.^2+HEphi.^2)))), '-r');
set(gcf, 'PaperType', 'A4'); 
xlabel('[º]');
ylabel('Dield [dB]');
legend('H plane');
axis([-90 90 -60 0]);

% ----- Plano E -----------------------------------------------------------
% Phi = 90º

subplot(2,2,4)
plot(linspace(-pi/2,pi/2,1000)*(180/pi), 20.*log10(abs(HEFAresta).*abs(sqrt(EEtheta.^2+EEphi.^2)))-max( 20.*log10(abs(HEFAresta).*abs(sqrt(EEtheta.^2+EEphi.^2)))), '-r');
set(gcf, 'PaperType', 'A4'); 
xlabel('[º]');
ylabel('Field [dB]');
legend('E plane');
axis([-90 90 -60 0]);