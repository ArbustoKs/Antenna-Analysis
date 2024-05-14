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



AF=zeros(length(theta), length(phi));

for i=1:length(theta)
    for j=1:length(phi)
        AF(i,j) = sum(exp(1j*(0:Nx-1)*k*dx*(sin(theta(i))*cos(phi(j)) - cos(phi_input)))) * sum(exp(1j*(0:Mz-1)*k*dz*cos(theta(i))));

    end
end




E=abs(AF)./abs(max(max(AF)));



dtheta=2*pi/100;
dphi=pi/100;

[PHI,THETA]=meshgrid(phi, theta);

Int=sum(sum((E.^2).*sin(THETA)*dtheta*dphi));
D=4*pi/Int;
DdB=10.*log10(D)
















