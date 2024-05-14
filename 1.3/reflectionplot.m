clear
clc

ETA = 120*pi;
d=0:0.01:1;
lamda=1;
k=2*pi/lamda; 


u0=k*d;
u1=k*(sqrt(d.^2+(1/2)^2)+1/2);
u2=k*(sqrt(d.^2+(1/2)^2)-1/2);
R1m=(ETA/(4*pi))*(2*cosint(u0)-cosint(u1)-cosint(u2));
X1m=-(ETA/(4*pi))*(2*sinint(u0)-sinint(u1)-sinint(u2));

Z1m=R1m+1j*X1m;



u3=k*2*d;
u4=k*(sqrt((2*d).^2+(1/2)^2)+1/2);
u5=k*(sqrt((2*d).^2+(1/2)^2)-1/2);
R2m=(ETA/(4*pi))*(2*cosint(u3)-cosint(u4)-cosint(u5));
X2m=-(ETA/(4*pi))*(2*sinint(u3)-sinint(u4)-sinint(u5));

Z2m=R2m+1j*X2m;


Zs=73.1+1j*42.5;


Zin = Zs - 2*((Z1m.^2)./(Zs+Z2m));


S11=(Zin-50)./(Zin+50);
disp(S11);


plot(d,abs(S11));
xlabel('d')
title('Reflection coefficient')

