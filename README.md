# thesis
clc;         % Clear Command window 
clearvars;   % removes all variables from the currently active workspace
% The Main program
global rou EI l mu
l = 50;        % Bridge Length (m)
Nele = 500;    % Number of element
mode = Nele-1; % Input('how many mode shape')
Nom = 2; %No. of modes used
N = Nom+4;    % 4 dof from vehicle
% Input of mass matrix parameters (Half Car Model) 
ms = 17000;    % Half car body mass (kg)
J = 9e4;       % Rotary Inertia (kg.m^2)
mt1 = 600;     % Front axle mass (kg) 
mt2 = 1000;    % Rear axle mass (kg)
% Input of stiffness matrix parameters (Half Car Model)
ks1 = 232000;  % Front suspension stifness (N/m)
ks2 = 746000;  % Rear suspension stifness (N/m)
kt1 = 1570000; % Front wheel stifness (N/m)
kt2 = 3140000; % Rear Wheel stiffness (N/m)
% Input of Damping matrix parameters (Half Car Model)
cs1 = 50000;   % Front suspension Damping (N.s/m)
cs2 = 70000;   % Rear suspension Damping (N.s/m)
ct1 = 0;       % Front wheel Damping (N.s/m)
ct2 = 0;       % Rear wheel Damping (N.s/m)
% Input of the Distance From CG
a1 = 3.8;      % Distance of CG from front wheel (m)
a2 = 1.5;      % Distance of CG from rear wheel (m)
a = a1+a2;     % Distance between front & rear wheel (m)
% Input of the Velocity of vehicle
vel = 40;      % the speed of the car (km/hr)
velocity = vel*1000/3600; % the speed of the car (m/s)
global rou EI l mu
rou = 6818.5;       % Bridge mass per unit Length (kg/m) 
EI = 6.961766e10;   % Flexural rigidity of Bridge ()
mu = 0.0;           % Damping coefficient of Bridge per unit Length ()
T = (l+a)/velocity; % Total time rquired to pass the bridge by Vehicle (sec)
dL = l/Nele;        % Length of each element (m)
Nada = round(a/dL);
dt = T/(Nele+Nada);  %0.000001*T;
t = 0:dt:T;
t = t';
x = t*velocity;
L = l+a;


%% Roughness Calculation
k = 3;           % Values For ISO Road A-B Roughness Classification, from 3 to 9
LT  = length(t);  % Number of data points
e  = velocity*dt;     % Sampling Interval (m)
dn = 1/L;        % Frequency Band 
n0 = 0.1;        % Spatial Frequency (cycles/m)
n  = dn:dn:LT*dn; % Spatial Frequency Band
phi =  2*pi*rand(size(n)); % Random Phase Angle
Amp1 = sqrt(dn)*(2^k)*(1e-3)*(n0./n); % Amplitude for Road  Class A-B
x = 0:e:L; % Abscissa Variable from 0 to L
r = zeros(size(x));
for i=1:length(x)
    r(i) = sum(Amp1.*cos(2*pi*n*x(i)+ phi));
end
r=r';
% fc = 0.0328;%explain%
% % fc = 0.1;
% fn = 1/e/2;
% order = 6;
% [b14,a14] = butter(order,(fc/fn), 'high');
% r = filtfilt (b14, a14, r1);
% save('r40','r');
%% Modal Parameters undamaged
C = EI/dL^3;
Nd = (Nele+1);
DOF = 2*Nd;%each node two DOF, y, rotation

%Element Stiffness Matrix
ke = C*[12 6*dL -12 6*dL;
        6*dL 4*dL^2 -6*dL 2*dL^2;
        -12 -6*dL 12 -6*dL;
        6*dL 2*dL^2 -6*dL 4*dL^2];
    
%Element Mass Matrix
mass = rou*dL; %Element mass
me = diag([mass/2 0 mass/2 0]); %????????????????

%Global Stiffness
K = zeros(DOF);
ML = zeros(DOF);
for j = 1:Nele
    i = 2*j-1;
    K(i+(0:3),i+(0:3))=K(i+(0:3),i+(0:3))+ke;
    ML(i+(0:3),i+(0:3))=ML(i+(0:3),i+(0:3))+me;
end

%Static Condensation
Ktt = K(1:2:end,1:2:end);
Kto = K(1:2:end,2:2:end);
Kot = K(2:2:end,1:2:end);
Koo = K(2:2:end,2:2:end);
Kt = Ktt - Kot'*inv(Koo)*Kot;
M = ML(1:2:end,1:2:end);

%Delete support DOF
rDOF = [1 Nd];
Kt(rDOF,:)=[];
Kt(:,rDOF)=[];
M(rDOF,:)=[];
M(:,rDOF)=[];

%Eigen Solution
[ph,Lambda] = eig(Kt,M);
[Lambda, I] = sort(diag(Lambda));
f = sqrt(Lambda)/(2*pi);
ph = ph(:,I);

Mn = diag(ph'*M*ph);
x1=zeros(1,length(ph));
phi=[x1;ph;x1];
phi1=zeros(Nele+Nada+1,Nele-1);
phi1(1:Nele+1,1:Nele-1)=phi;
phi1_dot=zeros(Nele+Nada+1,Nele-1);
phi2=zeros((Nele+Nada+1),Nele-1);
phi2(Nada+1:Nele+Nada+1,1:Nele-1)=phi;
phi2_dot=zeros((Nele+Nada+1), Nele-1);
save('Mode_un.mat','phi1','phi2','phi1_dot','phi2_dot','Lambda');
%% MCK and Force (undamaged)
load('r40.mat');
r_dot=(r(2:length(r))-r(1:length(r)-1))/dt;
r_dot= [r_dot; 0];
A=zeros(N,length(t));
V=zeros(N,length(t));
X=zeros(N,length(t));
Q=zeros(N,length(t));
A0=zeros(N,1);
V0=zeros(N,1);
X0=zeros(N,1);

for jjj=1:length(t)
 if t(jjj)*velocity<l%length l
 delta1=1;

 irr1=r(jjj);
 irr_dot1=r_dot(jjj);
 else delta1=0;
 irr1=0;
 irr_dot1=0;
 end
 
 if t(jjj)*velocity<=(l+a) && t(jjj)*velocity>=a
 delta2=1;
 irr2=r(jjj-2);
 irr_dot2=r_dot(jjj-2);
 else delta2=0;
 irr2=0;
 irr_dot2=0;
 end

[M,C,K] = MCK_un(ms,mt1,mt2,J,ks1,ks2,cs1,cs2,kt1,kt2,ct1,ct2,a1,a2,N,delta1,delta2,t(jjj),jjj,velocity);
Q = force_un(ct1,ct2,kt1,kt2,ms,a1,a2,mt1,mt2,irr1,irr2,irr_dot1,irr_dot2,velocity,t(jjj),jjj,delta1,delta2,N);
 
 A0=M\(Q-C*V0-K*X0);
 Kc = K+(2*C)/dt+(4*M)/dt^2;
 Qc = Q+C*((2*X0)/dt+V0)+M*((4*X0)/dt^2+(4*V0)/dt+A0);
 X1 = Kc\Qc;
 V1 = 2*(X1-X0)/dt-V0;
 A1 = M\(Q-C*V1-K*X1);

 X0=X1;
 V0=V1;
 A0=A1;
 X(:,jjj)=X0;
 V(:,jjj)=V0;
 A(:,jjj)=A0;
 
end
%%
Mv=M(1:4,1:4);
Kv=K(1:4,1:4);
[ph,Lambda] = eig(Kv,Mv);
[Lambda, I] = sort(diag(Lambda));
f = sqrt(Lambda)/(2*pi);
%% Response (undamaged)
yb1 = sqrt(2/rou/L)*phi1(:,1:Nom)*X(5:4+Nom,:); % yb1 = sqrt(2/rou/L)*phi1(:,1:mode)*X(5:end,:);
Vb1 = sqrt(2/rou/L)*phi1(:,1:Nom)*V(5:4+Nom,:); % Vb1 = sqrt(2/rou/L)*phi1(:,1:mode)*V(5:end,:);
Ab1 = sqrt(2/rou/L)*phi1(:,1:Nom)*A(5:4+Nom,:); % Ab1 = sqrt(2/rou/L)*phi1(:,1:mode)*A(5:end,:);

yb1_contact = diag(yb1)';
yb2=zeros(Nele+Nada+1,Nele+Nada+1);
yb2(Nada+1:Nele+Nada+1,Nada+1:Nele+Nada+1)=yb1(1:Nele+1,Nada+1:end);
yb2_contact = diag(yb2)';
%% save file (undamaged)
yb_un = yb1(2:Nele,:)'; % yb = yb1(2:Nele,:);
Vb_un = Vb1(2:Nele,:)'; % Vb = Vb1(2:Nele,:);
Ab_un = Ab1(2:Nele,:)'; % Ab = Ab1(2:Nele,:);
save('yb_un','yb_un');
save('Vb_un','Vb_un');
save('Ab_un','Ab_un');
t1=t';
plot(t1,(yb_un(:,250)),'-')
xlabel('Time (sec)');
ylabel('Displacement (m)');
figure;
plot(t1,(Vb_un(:,250)),'-')
xlabel('Time (sec)');
ylabel('Velocity (m/s)');
figure;
plot(t1,(Ab_un(:,250)),'-')
xlabel('Time (sec)');
ylabel('Accleration (m/s^2)');
%% Modal Parameters Damaged (D13)
C = EI/dL^3;
Nd = (Nele+1);
DOF = 2*Nd;

% Element Stiffness Matrix
ke = C*[12 6*dL -12 6*dL;
        6*dL 4*dL^2 -6*dL 2*dL^2;
        -12 -6*dL 12 -6*dL;
        6*dL 2*dL^2 -6*dL 4*dL^2];
    
% Element Mass Matrix
mass = rou*dL; %Element mass
me = diag([mass/2 0 mass/2 0]);

% Global Stiffness
K = zeros(DOF);
ML = zeros(DOF);

% Damage introduction
beta = 0.4;          % Damage factor
ked = (1-beta)*ke;   % Damaged element stiffness
for j = 1:240
    i = 2*j-1;
    K(i+(0:3),i+(0:3))=K(i+(0:3),i+(0:3))+ke;
end

for j = 241:260 % Element 13 to 14 are damaged
    i = 2*j-1;
    K(i+(0:3),i+(0:3))=K(i+(0:3),i+(0:3))+ked;
end

for j = 261:Nele
    i = 2*j-1;
    K(i+(0:3),i+(0:3))=K(i+(0:3),i+(0:3))+ke;
end

for j = 1:Nele
    i = 2*j-1;
    ML(i+(0:3),i+(0:3))=ML(i+(0:3),i+(0:3))+me;
end


%Static Condensation
Ktt = K(1:2:end,1:2:end);
Kto = K(1:2:end,2:2:end);
Kot = K(2:2:end,1:2:end);
Koo = K(2:2:end,2:2:end);
Kt = Ktt - Kot'*inv(Koo)*Kot;
M = ML(1:2:end,1:2:end);

%Delete support DOF
rDOF = [1 Nd];
Kt(rDOF,:)=[];
Kt(:,rDOF)=[];
M(rDOF,:)=[];
M(:,rDOF)=[];
%Eigen Solution
[ph,Lambda] = eig(Kt,M);
[Lambda, I] = sort(diag(Lambda));
f = sqrt(Lambda)/(2*pi);
ph = ph(:,I);

Mn = diag(ph'*M*ph);
x1=zeros(1,length(ph));
phi=[x1;ph;x1];
phi1=zeros(Nele+Nada+1,Nele-1);
phi1(1:Nele+1,1:Nele-1)=phi;
phi1_dot=zeros(Nele+Nada+1,Nele-1);
phi2=zeros((Nele+Nada+1),Nele-1);
phi2(Nada+1:Nele+Nada+1,1:Nele-1)=phi;
phi2_dot=zeros((Nele+Nada+1), Nele-1);
save('Mode_D13.mat','phi1','phi2','phi1_dot','phi2_dot','Lambda');
%% MCK and Force (D13)
% r_dot=(r(2:length(r))-r(1:length(r)-1))/dt;
% r_dot= [r_dot; 0];
A=zeros(N,length(t));
V=zeros(N,length(t));
X=zeros(N,length(t));
Q=zeros(N,length(t));
A0=zeros(N,1);
V0=zeros(N,1);
X0=zeros(N,1);

for jjj=1:length(t)
 if t(jjj)*velocity<l%length l
 delta1=1;

 irr1=r(jjj);
 irr_dot1=r_dot(jjj);
 else delta1=0;
 irr1=0;
 irr_dot1=0;
 end
 
 if t(jjj)*velocity<=(l+a) && t(jjj)*velocity>=a
 delta2=1;
 irr2=r(jjj-2);
 irr_dot2=r_dot(jjj-2);
 else delta2=0;
 irr2=0;
 irr_dot2=0;
 end

[M,C,K] = MCK_D13(ms,mt1,mt2,J,ks1,ks2,cs1,cs2,kt1,kt2,ct1,ct2,a1,a2,N,delta1,delta2,t(jjj),jjj,velocity);
Q = force_D13(ct1,ct2,kt1,kt2,ms,a1,a2,mt1,mt2,irr1,irr2,irr_dot1,irr_dot2,velocity,t(jjj),jjj,delta1,delta2,N);
 
 A0=M\(Q-C*V0-K*X0);
 Kc = K+(2*C)/dt+(4*M)/dt^2;
 Qc = Q+C*((2*X0)/dt+V0)+M*((4*X0)/dt^2+(4*V0)/dt+A0);
 X1 = Kc\Qc;
 V1 = 2*(X1-X0)/dt-V0;
 A1 = M\(Q-C*V1-K*X1);

 X0=X1;
 V0=V1;
 A0=A1;
 X(:,jjj)=X0;
 V(:,jjj)=V0;
 A(:,jjj)=A0;
 
end

%% Response (D13)
yb1 = sqrt(2/rou/L)*phi1(:,1:Nom)*X(5:4+Nom,:); % yb1 = sqrt(2/rou/L)*phi1(:,1:mode)*X(5:end,:);
Vb1 = sqrt(2/rou/L)*phi1(:,1:Nom)*V(5:4+Nom,:); % Vb1 = sqrt(2/rou/L)*phi1(:,1:mode)*V(5:end,:);
Ab1 = sqrt(2/rou/L)*phi1(:,1:Nom)*A(5:4+Nom,:); % Ab1 = sqrt(2/rou/L)*phi1(:,1:mode)*A(5:end,:);

yb1_contact = diag(yb1)';
yb2=zeros(Nele+Nada+1,Nele+Nada+1);
yb2(Nada+1:Nele+Nada+1,Nada+1:Nele+Nada+1)=yb1(1:Nele+1,Nada+1:end);
yb2_contact = diag(yb2)';
%% Save file (D13)
yb_D13 = yb1(2:Nele,:)'; % yb = yb1(2:Nele,:);
Vb_D13 = Vb1(2:Nele,:)'; % Vb = Vb1(2:Nele,:);
Ab_D13 = Ab1(2:Nele,:)'; % Ab = Ab1(2:Nele,:);
save('yb_D13','yb_D13');
save('Vb_D13','Vb_D13');
save('Ab_D13','Ab_D13');
plot(t1,(yb_D13(:,250)),'-')
xlabel('Time (sec)');
ylabel('Displacement (m)');
figure;
plot(t1,(Vb_D13(:,250)),'-')
xlabel('Time (sec)');
ylabel('Velocity (m/s)');
figure;
plot(t1,(Ab_D13(:,250)),'-')
xlabel('Time (sec)');
ylabel('Accleration (m/s^2)');
%% Modal Parameter Damaged (D7)
C = EI/dL^3;
Nd = (Nele+1);
DOF = 2*Nd;

% Element Stiffness Matrix
ke = C*[12 6*dL -12 6*dL;
        6*dL 4*dL^2 -6*dL 2*dL^2;
        -12 -6*dL 12 -6*dL;
        6*dL 2*dL^2 -6*dL 4*dL^2];
    
% Element Mass Matrix
mass = rou*dL; %Element mass
me = diag([mass/2 0 mass/2 0]);

% Global Stiffness
K = zeros(DOF);
ML = zeros(DOF);

% Damage introduction
beta = 0.4;          % Damage factor
ked = (1-beta)*ke;   % Damaged element stiffness
for j = 1:120
    i = 2*j-1;
    K(i+(0:3),i+(0:3))=K(i+(0:3),i+(0:3))+ke;
end

for j = 121:140 % Element 7 
    i = 2*j-1;
    K(i+(0:3),i+(0:3))=K(i+(0:3),i+(0:3))+ked;
end

for j = 141:Nele
    i = 2*j-1;
    K(i+(0:3),i+(0:3))=K(i+(0:3),i+(0:3))+ke;
end

for j = 1:Nele
    i = 2*j-1;
    ML(i+(0:3),i+(0:3))=ML(i+(0:3),i+(0:3))+me;
end


%Static Condensation
Ktt = K(1:2:end,1:2:end);
Kto = K(1:2:end,2:2:end);
Kot = K(2:2:end,1:2:end);
Koo = K(2:2:end,2:2:end);
Kt = Ktt - Kot'*inv(Koo)*Kot;
M = ML(1:2:end,1:2:end);

%Delete support DOF
rDOF = [1 Nd];
Kt(rDOF,:)=[];
Kt(:,rDOF)=[];
M(rDOF,:)=[];
M(:,rDOF)=[];
%Eigen Solution
[ph,Lambda] = eig(Kt,M);
[Lambda, I] = sort(diag(Lambda));
f = sqrt(Lambda)/(2*pi);
ph = ph(:,I);

Mn = diag(ph'*M*ph);
x1=zeros(1,length(ph));
phi=[x1;ph;x1];
phi1=zeros(Nele+Nada+1,Nele-1);
phi1(1:Nele+1,1:Nele-1)=phi;
phi1_dot=zeros(Nele+Nada+1,Nele-1);
phi2=zeros((Nele+Nada+1),Nele-1);
phi2(Nada+1:Nele+Nada+1,1:Nele-1)=phi;
phi2_dot=zeros((Nele+Nada+1), Nele-1);
save('Mode_D7.mat','phi1','phi2','phi1_dot','phi2_dot','Lambda');
%% MCK and Force (D7)
% r_dot=(r(2:length(r))-r(1:length(r)-1))/dt;
% r_dot= [r_dot; 0];
A=zeros(N,length(t));
V=zeros(N,length(t));
X=zeros(N,length(t));
Q=zeros(N,length(t));
A0=zeros(N,1);
V0=zeros(N,1);
X0=zeros(N,1);

for jjj=1:length(t)
 if t(jjj)*velocity<l%length l
 delta1=1;

 irr1=r(jjj);
 irr_dot1=r_dot(jjj);
 else delta1=0;
 irr1=0;
 irr_dot1=0;
 end
 
 if t(jjj)*velocity<=(l+a) && t(jjj)*velocity>=a
 delta2=1;
 irr2=r(jjj-2);
 irr_dot2=r_dot(jjj-2);
 else delta2=0;
 irr2=0;
 irr_dot2=0;
 end

[M,C,K] = MCK_D7(ms,mt1,mt2,J,ks1,ks2,cs1,cs2,kt1,kt2,ct1,ct2,a1,a2,N,delta1,delta2,t(jjj),jjj,velocity);
Q = force_D7(ct1,ct2,kt1,kt2,ms,a1,a2,mt1,mt2,irr1,irr2,irr_dot1,irr_dot2,velocity,t(jjj),jjj,delta1,delta2,N);
 
 A0=M\(Q-C*V0-K*X0);
 Kc = K+(2*C)/dt+(4*M)/dt^2;
 Qc = Q+C*((2*X0)/dt+V0)+M*((4*X0)/dt^2+(4*V0)/dt+A0);
 X1 = Kc\Qc;
 V1 = 2*(X1-X0)/dt-V0;
 A1 = M\(Q-C*V1-K*X1);

 X0=X1;
 V0=V1;
 A0=A1;
 X(:,jjj)=X0;
 V(:,jjj)=V0;
 A(:,jjj)=A0;
 
end

%% Response (D7)
yb1 = sqrt(2/rou/L)*phi1(:,1:Nom)*X(5:4+Nom,:); % yb1 = sqrt(2/rou/L)*phi1(:,1:mode)*X(5:end,:);
Vb1 = sqrt(2/rou/L)*phi1(:,1:Nom)*V(5:4+Nom,:); % Vb1 = sqrt(2/rou/L)*phi1(:,1:mode)*V(5:end,:);
Ab1 = sqrt(2/rou/L)*phi1(:,1:Nom)*A(5:4+Nom,:); % Ab1 = sqrt(2/rou/L)*phi1(:,1:mode)*A(5:end,:);

yb1_contact = diag(yb1)';
yb2=zeros(Nele+Nada+1,Nele+Nada+1);
yb2(Nada+1:Nele+Nada+1,Nada+1:Nele+Nada+1)=yb1(1:Nele+1,Nada+1:end);
yb2_contact = diag(yb2)';
%% Save file (D7)
yb_D7 = yb1(2:Nele,:)'; % yb = yb1(2:Nele,:);
Vb_D7 = Vb1(2:Nele,:)'; % Vb = Vb1(2:Nele,:);
Ab_D7 = Ab1(2:Nele,:)'; % Ab = Ab1(2:Nele,:);
save('yb_D7','yb_D7');
save('Vb_D7','Vb_D7');
save('Ab_D7','Ab_D7');
plot(t1,(yb_D7(:,250)),'-')
xlabel('Time (sec)');
ylabel('Displacement (m)');
figure;
plot(t1,(Vb_D7(:,250)),'-')
xlabel('Time (sec)');
ylabel('Velocity (m/s)');
%figure;
plot(t1,(Ab_D7(:,250)),'-')
xlabel('Time (sec)');
ylabel('Accleration (m/s^2)');
