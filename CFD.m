% EGRM321-002
% Fall 2014
% Numerical Methods
% Project 2
% Peter Donaldson

clc
clear 
close all

%The pleat dimesions
l=0.0275;
h=0.0011;

%initial y positions
xp=0;
yp=[0.05 0.25 0.5 0.75 0.95];
%yp(1)=0.05h
%yp(2)=0.25h
%yp(3)=0.5h
%yp(4)=0.75h
%yp(5)=0.95h

%Number of trajectories
[lol,Num_Y]=size(yp);

%particle diameter (m)
%d_p
part_dia=10.^-6;
%Particle diameter (kg/m^3)
%P_p
part_dens=1000;
%Air Viscosity
%mu
air_visc=1.85*10.^-5;
%Average Air Velocity At Enterence
Ui_avg=1;



for iter1=1:Num_Y
  %Stokes number
Stokes=(part_dens*part_dia.^2*Ui_avg)/(18*air_visc*2*yp(iter1));

%Initial Angle
Ang_i=75*(0.78*(yp(iter1)/h).^2+0.16*(yp(iter1)/h))*exp(-1.61*Stokes);


end
%set up initial heights
y_p(iter1)=yp(iter1)*h;


%Iteration Size
H=0.01;

x_p=0;
 Xsize=l/H;
for y=1:5

i3=1;

%for x=0:H:l

x=0;
u_avg=Ui_avg*(1-(x_p/(l+h)))
v_avg=Ui_avg*(y_p(y)/(l+y_p(y)))
u1=(3/2)*u_avg*(1-(y_p(y)/h).^2) * H
v1=(v_avg*sin((pi/2)*(y_p(y)/h))) * H
u_avg=Ui_avg*(1-((x_p+H/2)/(l+h)));
u2=((3/2)*u_avg*(1-((y_p(y)+u1/2)/h)).^2)*H
v2=(v_avg*sin((pi/2)*((y_p(y)+v1/2)/h))) * H
u_avg=Ui_avg*(1-((x_p+H/2)/(l+h)));
u3=((3/2)*u_avg*((1-(y_p(y)+u2/2))/h).^2)*H
v3=(v_avg*sin((pi/2)*((y_p(y)+v2/2)/h))) * H
u_avg=Ui_avg*(1-((x_p+H)/(l+h)));
u4=((3/2)*u_avg*((1-(y_p(y)+u3/2))/h).^2)*H
v4=(v_avg*sin((pi/2)*((y_p(y)+u3/2)/h))) * H
y_p(y)=y_p(y)+((v1+2*v2+v3*2+v4)/6);
x_p=x_p+((u1+2*u2+2*u3+u4)/6)
i3=i3+1;
Y_p(i3)=y_p(y);
%Y_p(i3,y)=y_p(y);
X_p(i3)=x_p;
%end
i3;
figure

plot(X_p,Y_p,'+')
grid on;


end



part_relax=(part_dia.^2*part_dens*air_visc.^-1)/18;

%Average air velocity in pleats in x direction;


%Average air velocity in pleats in y direction;






% Velocity grid size
% Xsize=l/H;
% Ysize=Num_Y;
% 
% u_matr=zeros(Xsize,Ysize);
% v_matr=zeros(Xsize,Ysize);




% for n=1
% 
% k1=delta*dT;
% k2=delta*(Ti)/(rho*Vsys1*Cp));
% k2=delta*dT(Tc+.5*k1);
% k3=delta*dT(Tc+.5*k2);
% k3=delta*((-h*A*((Tc+delta*k2*.5)-Ti))/(rho*Vsys1*Cp));
% k4=delta*dT(Tc+k3);
% k4=delta*((-h*A*((Tc+delta*k3)-Ti))/(rho*Vsys1*Cp));
% 
% 
% 
% 
% 
% for n=1:F
% Ra=(g*dB*(L.^3)*(Tc-Ti))/(v*alp);
% Nu=0.54*Ra.^(1/4);
% h=(Nu*kair)/L;
% 
% dT=(-h*A*(Tc-Ti))/(rho*Vsys1*Cp);
% k1=delta*dT;
% k2=delta*((-h*A*((Tc+delta*k1*.5)-Ti))/(rho*Vsys1*Cp));
% k2=delta*dT(Tc+.5*k1);
% k3=delta*dT(Tc+.5*k2);
% k3=delta*((-h*A*((Tc+delta*k2*.5)-Ti))/(rho*Vsys1*Cp));
% k4=delta*dT(Tc+k3);
% k4=delta*((-h*A*((Tc+delta*k3)-Ti))/(rho*Vsys1*Cp));
% 
% 
% 
% Tc=Tc+k1/6+k2/3+k3/3+k4/6;
% Tct2(n)=Tc;
% end