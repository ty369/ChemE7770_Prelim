function f=ACDC2(t,Y,S)
% global S
% Y(1)=X, Y(2)=Y, Y(3)=Z f(1)=dx/dt, f(2)=dy/dt, f(3)=dz/dt
% Parameters from Suppl. Material S-1 Table
alphax=3.9*10e-2; alphay=4.3*10e-3;
betax=6.1; betay=5.7;
zx=1.3*10e-5; yz=11*10e-3; xz=12*10e-2; xy=7.9*10e-4;
nzx=2.32; nxz=2; nxy=2; nyz=2;
deltay=1.05; deltaz=1.04;

f(1,1)=(alphax+betax.*S)./(1+S+((Y(3)./zx).^nzx))-Y(1);
f(2,1)=(alphay+betay.*S)./(1+S+(Y(1)./xy).^nxy)-deltay.*Y(2);
f(3,1)=1./(1+((Y(1)./xz).^nxz)+(Y(2)./yz).^nyz)-deltaz.*Y(3);
