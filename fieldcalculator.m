function [Px, Py, HEtheta, HEphi, EEtheta, EEphi] = fieldcalculator(Lx,Ly, Ex, Ey, k)
% obtener el campo radiado de una apertura dado el campo de apertura por
% las componentes Ex y Ey.
% Debemos aplicar el 2 principio de equivalencia teniendo en cuenta 
% que conocemos los campos electricos segun theta y phi.

% Precisión de los cálculos
res = 49;
Hres = 1000;
Eres = 1000;

%% Declaraciones de vectores Theta y Phi

% Plano H

Htheta = linspace(-pi/2,pi/2,Hres);
Hphi = 0;
Hu = sin(Htheta)*cos(Hphi);
Hv = sin(Htheta)*sin(Hphi);

HPx = zeros(1,Hres); 
HPy = zeros(1,Hres); 

% Plano E

Etheta = linspace(-pi/2,pi/2,Eres);
Ephi = pi/2;
Eu = sin(Etheta)*cos(Ephi);
Ev = sin(Etheta)*sin(Ephi);

EPx = zeros(1,Eres); 
EPy = zeros(1,Eres); 

% 3D

[THETA,PHI] = meshgrid(linspace(-pi/2,pi/2,res), linspace(0,2*pi,res));
u = sin(THETA).*cos(PHI);
v = sin(THETA).*sin(PHI);

Px = zeros(res,res);
Py = zeros(res,res);

%% Operaciones

% Plano H
% Realizamos la operación para Px
for n=1:Hres
    %Primero defino campo
    P = @(x,y) Ex(x,y).*exp(1i*(k*Hu(n).*x)).*exp(1i*(k*Hv(n).*y));
    %Luego integro campo
    HPx(n) = integral2(P,-Lx/2,Lx/2,-Ly/2,Ly/2);
end
% Realizamos la integral para Py
for n=1:Hres
    %Primero defino campo
    P = @(x,y) Ey(x,y).*exp(1i*(k*Hu(n).*x)).*exp(1i*(k*Hv(n).*y));
    %Luego integro campo
    HPy(n) = integral2(P,-Lx/2,Lx/2,-Ly/2,Ly/2); 
end
HEtheta = (1i.*k)/(2.*pi).*(HPx.*cos(Hphi)+HPy.*sin(Hphi));
HEphi = ((-1i.*k)/(2.*pi)).*(cos(Htheta).*(HPx.*sin(Hphi)-HPy.*cos(Hphi)));


% Plano E
% Realizamos la operación para Px
for n=1:Eres
    %Primero defino campo
    P = @(x,y) Ex(x,y).*exp(1i*(k*Eu(n).*x)).*exp(1i*(k*Ev(n).*y));
    %Luego integro campo
    EPx(n) = integral2(P,-Lx/2,Lx/2,-Ly/2,Ly/2);
end
% Realizamos la integral para Py
for n=1:Eres
    %Primero defino campo
    P = @(x,y) Ey(x,y).*exp(1i*(k*Eu(n).*x)).*exp(1i*(k*Ev(n).*y));
    %Luego integro campo
    EPy(n) = integral2(P,-Lx/2,Lx/2,-Ly/2,Ly/2); 
end
EEtheta = (1i.*k)/(2.*pi).*(EPx.*cos(Ephi)+EPy.*sin(Ephi));
EEphi = ((-1i.*k)/(2.*pi)).*(cos(Etheta).*(EPx.*sin(Ephi)-EPy.*cos(Ephi)));


% 3D
% Realizamos la operación para Px
for n1=1:res
    for n2=1:res
        %Primero defino campo
        P = @(x,y) Ex(x,y).*exp(1i*(k*u(n1,n2).*x)).*exp(1i*(k*v(n1,n2).*y));
        %Luego integro campo
        Px(n1,n2) = integral2(P,-Lx/2,Lx/2,-Ly/2,Ly/2); 
    end
end
% Realizamos la integral para Py
for n1=1:res
    for n2=1:res
        %Primero defino campo
        P = @(x,y) Ey(x,y).*exp(1i*(k*u(n1,n2).*x)).*exp(1i*(k*v(n1,n2).*y));
        %Luego integro campo
        Py(n1,n2) = integral2(P,-Lx/2,Lx/2,-Ly/2,Ly/2); 
    end
end

end
