clear
clc

ETA = 120*pi;
lamda=1;
M=100;  
d=linspace(0, 1, M);
h=linspace(0, 1, M); 


k=2*pi/lamda; 


%b (Mutual impedance for distance=d)

u0=zeros(length(d), length(h));
u1=zeros(length(d), length(h));
u2=zeros(length(d), length(h));

R1m=zeros(length(d), length(h));
X1m=zeros(length(d), length(h));
b=zeros(length(d), length(h));

for i=1:length(d)
    u0(i)=k*d(i);
    u1(i)=k*(sqrt(d(i).^2+(1/2)^2)+1/2);
    u2(i)=k*(sqrt(d(i).^2+(1/2)^2)-1/2);
    R1m(i)=(ETA/(4*pi))*(2*cosint(u0(i))-cosint(u1(i))-cosint(u2(i)));
    X1m(i)=-(ETA/(4*pi))*(2*sinint(u0(i))-sinint(u1(i))-sinint(u2(i)));

    b(i,1)=R1m(i)+1j*X1m(i);


end



%c (Mutual impedance for distance=2d)

u3=zeros(length(d), length(h));
u4=zeros(length(d), length(h));
u5=zeros(length(d), length(h));

R2m=zeros(length(d), length(h));
X2m=zeros(length(d), length(h));
c=zeros(length(d), length(h));

for i=1:length(d)
    u3(i)=k*2*d(i);
    u4(i)=k*(sqrt((2*d(i)).^2+(1/2)^2)+1/2);
    u5(i)=k*(sqrt((2*d(i)).^2+(1/2)^2)-1/2);
    R2m(i)=(ETA/(4*pi))*(2*cosint(u3(i))-cosint(u4(i))-cosint(u5(i)));
    X2m(i)=-(ETA/(4*pi))*(2*sinint(u3(i))-sinint(u4(i))-sinint(u5(i)));

    c(i,1)=R2m(i)+1j*X2m(i);
    


end



%e (Mutual impedance for distance=2h)

u9=zeros(length(d), length(h));
u10=zeros(length(d), length(h));
u11=zeros(length(d), length(h));

R4m=zeros(length(d), length(h));
X4m=zeros(length(d), length(h));
e=zeros(length(d), length(h));


for i=1:length(h)
    u9(i)=k*2*h(i);
    u10(i)=k*(sqrt((2*h(i)).^2+(1/2)^2)+1/2);
    u11(i)=k*(sqrt((2*h(i)).^2+(1/2)^2)-1/2);
    R4m(i)=(ETA/(4*pi))*(2*cosint(u9(i))-cosint(u10(i))-cosint(u11(i)));
    X4m(i)=-(ETA/(4*pi))*(2*sinint(u9(i))-sinint(u10(i))-sinint(u11(i)));

    e(1,i)=R4m(i)+1j*X4m(i);
    


end




%t  (Mutual impedance for the smaller diagonal distance)

m=zeros(length(d), length(h));

u6=zeros(length(d), length(h));
u7=zeros(length(d), length(h));
u8=zeros(length(d), length(h));

R3m=zeros(length(d), length(h));
X3m=zeros(length(d), length(h));
t=zeros(length(d), length(h));


for i=1:length(d)
    for j=1:length(h)
        
        m(i,j)=sqrt((2*h(j)).^2+(d(i).^2) );
        
        u6(i,j)=k*m(i,j);
        u7(i,j)=k*(sqrt(m(i,j).^2+(1/2)^2)+1/2);
        u8(i,j)=k*(sqrt(m(i,j).^2+(1/2)^2)-1/2);
        
        R3m(i,j)=(ETA/(4*pi))*(2*cosint(u6(i,j))-cosint(u7(i,j))-cosint(u8(i,j)));
        X3m(i,j)=-(ETA/(4*pi))*(2*sinint(u6(i,j))-sinint(u7(i,j))-sinint(u8(i,j)));
        
        t(i,j) = R3m(i,j)+1j*X3m(i,j);
        
    end
    
end




%f (Mutual impedance for the bigger diagonal distance)

n=zeros(length(d), length(h));

u12=zeros(length(d), length(h));
u13=zeros(length(d), length(h));
u14=zeros(length(d), length(h));

R5m=zeros(length(d), length(h));
X5m=zeros(length(d), length(h));
f=zeros(length(d), length(h));


for i=1:length(d)
    for j=1:length(h)
        
        n(i,j)=sqrt((2*h(j)).^2+(2*d(i)).^2);
        
        u12(i,j)=k*n(i,j);
        u13(i,j)=k*(sqrt(n(i,j).^2+(1/2)^2)+1/2);
        u14(i,j)=k*(sqrt(n(i,j).^2+(1/2)^2)-1/2);
        
        R5m(i,j)=(ETA/(4*pi))*(2*cosint(u12(i,j))-cosint(u13(i,j))-cosint(u14(i,j)));
        X5m(i,j)=-(ETA/(4*pi))*(2*sinint(u12(i,j))-sinint(u13(i,j))-sinint(u14(i,j)));
        
        f(i,j)=R5m(i,j)+1j*X5m(i,j);

        
        
    end
    
end




Zs=zeros(length(d), length(h));

for i=1:length(d)
    for j=1:length(h)
       Zs(i,j)=73.1+1j*42.1;
    end
    
end


Zin=zeros(length(d), length(h));
S11=zeros(length(d), length(h));

for i=1:length(d)
    for j=1:length(h)
        Zin(i,j) = (Zs(i,j)-e(i,j)) + (2*b(i,j)-2*t(i,j)).*((t(i,j)-b(i,j))./(Zs(i,j)+c(i,j)-f(i,j)-e(i,j)));
        S11(i,j) =(Zin(i,j)-50)./(Zin(i,j)+50);
    end
    
end

[d,h]=meshgrid(d,h);


contour(d,h,abs(S11), 0:0.01:0.3)
xlabel('d');
ylabel('h');
zlabel('S11');
title('Reflection coefficient less than 0.3');
colorbar


