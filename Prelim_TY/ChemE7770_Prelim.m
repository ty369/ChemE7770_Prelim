%% ChemE 7770 Prelim
clear all; close all; clc;

%% 1a).
n=[19 21 41 67 86 93 93];
m=3*1e-13;%gDW/cell
Nc=1e8;%cells/ml
V=1; % assume 1ml
B=m.*Nc.*V;
nmol=n.*(6.02.*1e23).^(-1).*1e9;
convert=nmol./B;
display(convert)

%% 2a 2b See hand-written part
%% 2c
clear all; 
tspan=linspace(0,30,50);
Y0=[10e-4 10e-4];
Sarray=linspace(10e-2,3.5);
% Increase S to a large number but still have the same trend, so I pick 3.5
% to enable zoom in at the region.
for i=1:length(Sarray)
    S=Sarray(i);
[t,Y]=ode45(@(t,Y)ACDC(t,Y,S),tspan,Y0);
A(:,i)=Y(:,1);
B(:,i)=Y(:,2);
end

% figure (1)
% subplot(2,1,1)
% for j=1:size(A,1)
% plot(t,A(:,j));
% hold on
% title('2c. X vs Time');
% xlabel('Time, t');
% ylabel('X');
% end
% subplot(2,1,2)
% for j=1:size(A,1)
% plot(Sarray,A(j,:),'.');
% hold on
% title('2c. X vs S');
% xlabel('S');
% ylabel('X');
% end

% figure(2) 
% plot(Sarray, A(size(A,1),:),'.');
% title('2c. X vs Time at S=3.5');
% xlabel('Time, t');
% ylabel('X');

fprintf('Matlab seems cannot show the swist part in between, I then use Excel to plot the middle parts.');

%% 2d
clear all; 
tspan=linspace(0,50);
Y0=[0 0 0];
Sarray=[0.02 10 10e5];
for i=1:length(Sarray)
    S=Sarray(i);
[t,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan,Y0);
A(:,i)=Y(:,1);
B(:,i)=Y(:,2);
C(:,i)=Y(:,3);
end
% For X Y Z at different S, we can plot:
% plot(t,A(:,1),t,B(:,1),t,C(:,1))--> S=0.02;
% plot(t,A(:,2),t,B(:,2),t,C(:,2))--> S=10;
% plot(t,A(:,3),t,B(:,3),t,C(:,3))--> S=10e5;

figure (3)
subplot(2,2,1)
plot(t,A(:,1))
title('2d. X vs Time, with S=0.02');
xlabel('Time, t');
ylabel('Gene Expression');
subplot(2,2,2)
plot(t,A(:,2));
title('2d. X vs Time, with S=10');
xlabel('Time, t');
ylabel('Gene Expression');
subplot(2,2,3)
plot(t,A(:,3));
title('2d. X vs Time, with S=10e5');
xlabel('Time, t');
ylabel('Gene Expression');
disp('From the plots, it seems like the values approx. matches the value shown in Fig.2');

%% 2e) I. S: 0.7-->100
clear all; 
tspan=linspace(0,100);
tspan2=linspace(100,150);
Y0=[0 0 0];
S=0.7;
[t,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan,Y0);
A=Y;
S=100;
Y01=[A(end,1) A(end,2) A(end,3)];
[t2,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan2,Y01);
B=Y;
Y02=1.25.*[A(end,1) A(end,2) A(end,3)];
[t2,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan2,Y02);
C=Y;
Y03=(1-0.25).*[A(end,1) A(end,2) A(end,3)];
[t2,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan2,Y03);
D=Y;

figure (4) % Plot of Z vs time with S change to 100
subplot(2,1,1)
plot(t,A(:,3),t2,B(:,3));
title('2e I. Z vs Time, Cell 1 with S changes from 0.7 to 100');
xlabel('Time, t');
ylabel('Gene Expression');
legend('S=0.7','S=100');

subplot(2,1,2)
plot(t2,B(:,3),t2,C(:,3),t2,D(:,3));
title('2e I. Z vs Time, with S changes from 0.7 to 100');
xlabel('Time, t');
ylabel('Gene Expression');
legend('Cell 1','Cell 2','Cell 3');


%% 2e II. S: 1010-->100
clear all; 
tspan=linspace(0,100);
tspan2=linspace(100,150);
Y0=[0 0 0];
S=1010;
[t,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan,Y0);
A=Y;
S=100;
Y01=[A(end,1) A(end,2) A(end,3)];
[t2,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan2,Y01);
B=Y;
Y02=1.25.*[A(end,1) A(end,2) A(end,3)];
[t2,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan2,Y02);
C=Y;
Y03=(1-0.25).*[A(end,1) A(end,2) A(end,3)];
[t2,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan2,Y03);
D=Y;

figure (5)
subplot(2,1,1)
plot(t,A(:,3),t2,B(:,3));
title('2e II. Z vs Time, Cell 1 with S changes from 1010 to 100');
xlabel('Time, t');
ylabel('Gene Expression');
legend('S=1010','S=100');

subplot(2,1,2)
plot(t,B(:,3),t,C(:,3),t,D(:,3));
title('2e II. Z vs Time, with S changes from 1010 to 100');
xlabel('Time, t');
ylabel('Gene Expression');
legend('Cell 1','Cell 2','Cell 3');

% Explanation from the paper: "This difference in behavior is a consequence of the different
% initial gene expression states in relation to oscillatory spiral center. In a Hopf bifurcation,
% oscillations arise through an attracting spiral losing its stability and becoming a repulsing spiral.
% Hence, oscillations originating from a Hopf bifurcation start their transient close to the unstable 
% spiral center, and a small variation in the initial condition can lead to a substantial difference in the
% final oscillation phase. Small initial differences are amplified, resulting in lack of coherence of 
% oscillations for a population of cells undergoing the bifurcation. By contrast, cells passing
% through the saddle-node bifurcation toward the limit cycle do so at expression levels that are far from 
% those associated with the attracting oscillatory regime. Consequently, they have the
% same initial phase, and stochastic trajectories are canalized together toward the oscillatory state."

%% 2f
clear all; 
tspan=linspace(0,100);
Y0=[0 0 0];
%Sarray=[0.7 100 1010];
tspan2=linspace(100,150);
S=105;
[t,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan,Y0);
A=Y;
S=100;
Y01=[A(end,1) A(end,2) A(end,3)];
[t2,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan2,Y01);
B=Y;
Y02=1.25.*[A(end,1) A(end,2) A(end,3)];
[t2,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan2,Y02);
C=Y;
Y03=(1-0.25).*[A(end,1) A(end,2) A(end,3)];
[t2,Y]=ode45(@(t,Y)ACDC2(t,Y,S),tspan2,Y03);
D=Y;

figure (6)
subplot(2,1,1)
plot(t,A(:,3),t2,B(:,3));
title('Z vs Time, Cell 1 with S changes from 105 to 100');
xlabel('Time, t');
ylabel('Gene Expression');
legend('S=105','S=100');

subplot(2,1,2)
plot(t2,B(:,3),t2,C(:,3),t2,D(:,3));
title('2f. Z vs Time, with S changes from 105 to 100');
xlabel('Time, t');
ylabel('Gene Expression');
legend('Cell 1','Cell 2','Cell 3');
