clear
clc


M=361; 
theta=linspace(0,2*pi,M);
phi=linspace(0,2*pi,M);


Nx=16;
Mz=12;
lamda=1;
k=2*pi/lamda; 

dx=input('Distance between dipoles (multiples of lamda): '); 
dz=dx; 

phi_input=input('Degrees of maximum radiation in radians: ');

deltax = -k*dx*cos(phi_input)*sin(pi/2);  
deltaz = -k*dz*cos(pi/2);   


psix= k*dx*cos(phi)*sin(pi/2) + deltax;
psiz= k*dz*cos(pi/2) + deltaz;
AFx=0;
AFz=0;


for n=1:Nx
    AFx = AFx + exp(1j*(n-1)*psix);
end

for m=1:Mz
    AFz = AFz + exp(1j*(m-1)*psiz);
end


AF=AFx.*AFz;
AFmag=abs(AF);

polarplot(phi, AFmag);
title('Radiation pattern of E-plane');
