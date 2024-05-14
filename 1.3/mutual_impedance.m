clear
clc


ETA = 120*pi;
d=0:0.01:3;
lamda=1;
k=2*pi/lamda; 


u0=k*d;
u1=k*(sqrt(d.^2+(1/2)^2)+1/2);
u2=k*(sqrt(d.^2+(1/2)^2)-1/2);
R21m=(ETA/(4*pi))*(2*cosint(u0)-cosint(u1)-cosint(u2));
X21m=-(ETA/(4*pi))*(2*sinint(u0)-sinint(u1)-sinint(u2));


R=plot(d/lamda,R21m);

hold on 

X=plot(d/lamda,X21m);

hold off


xlabel('d/λ')
ylabel('Mutual Impedance Z21m')
legend([R,X], 'R21m', 'X21m')