% ViscoelasticCycles_dense.m

% ViscoelasticCycles physical problem but with dense matrices
% Physical problem:
% - fault goes 35km deep
% - transition depth is 40 km
% - entire problem is 200km x 200km
%
% x1 (x) is in/out of page
% x2 (y) is left/right of page
% x3 (z) is up/down

% val's x2c is like my shearY_chat (50x1)
% whereas xx2c is grid of centers (50x50)

% ---     Kernels & functions      ---

G = 30e3;

% Boxcar function
boxc=@(x) (x+0.5>=0)-(x-0.5>=0);

% Shear Stress kernels for distributed deformation
s1312=@(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    +log((x2-L/2).^2+(x3+D+W).^2) - log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    -log((x2-L/2).^2+(x3+D).^2) + log((x2+L/2).^2+(x3+D).^2));

s1212=@(D,L,W,x2,x3) G/pi*( ...
    atan((x3-D)./(x2+L/2))-atan((x3-D)./(x2-L/2)) ...
    +atan((x3-D-W)./(x2-L/2))-atan((x3-D-W)./(x2+L/2)) ...
    -atan((x3+D+W)./(x2-L/2))-atan((x3+D)./(x2+L/2)) ...
    +atan((x3+D)./(x2-L/2))+atan((x3+D+W)./(x2+L/2)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

s1213=@(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    -log((x2-L/2).^2+(x3+D+W).^2) + log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    +log((x2-L/2).^2+(x3+D).^2) - log((x2+L/2).^2+(x3+D).^2));

s1313=@(D,L,W,x2,x3) G/pi*( ...
    atan((x2+L/2)./(x3-D))  -atan((x2-L/2)./(x3-D)) ...
    -atan((x2+L/2)./(x3-D-W))+atan((x2-L/2)./(x3-D-W)) ...
    +atan((x2+L/2)./(x3+D))  -atan((x2-L/2)./(x3+D)) ...
    -atan((x2+L/2)./(x3+D+W))+atan((x2-L/2)./(x3+D+W)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

% Displacement Kernels for distributed deformation
uk12=@(D,L,W,x2,x3) 1/(2*pi)*( ...
    (x3-D-W).*log((x2-L/2).^2+(x3-D-W).^2) ...
    -(x3-D-W).*log((x2+L/2).^2+(x3-D-W).^2) ...
    -(x3-D).*log((x2-L/2).^2+(x3-D).^2) ...
    +(x3-D).*log((x2+L/2).^2+(x3-D).^2) ...
    +2*(x2-L/2).*(atan((x3-D-W)./(x2-L/2))-atan((x3-D)./(x2-L/2))) ...
    +2*(x2+L/2).*(atan((x3-D)./(x2+L/2))-atan((x3-D-W)./(x2+L/2))) ...
    +(x3+D+W).*log((x2+L/2).^2+(x3+D+W).^2) ...
    -(x3+D+W).*log((x2-L/2).^2+(x3+D+W).^2) ...
    -(x3+D).*log((x2+L/2).^2+(x3+D).^2) ...
    +(x3+D).*log((x2-L/2).^2+(x3+D).^2) ...
    +2*(x2+L/2).*(atan((x3+D+W)./(x2+L/2))-atan((x3+D)./(x2+L/2))) ...
    +2*(x2-L/2).*(atan((x3+D)./(x2-L/2))-atan((x3+D+W)./(x2-L/2))) ...
    );

uk13=@(D,L,W,x2,x3) 1/(2*pi)*( ...
    (x2-L/2).*log((x2-L/2).^2+(x3-D-W).^2) ...
    -(x2+L/2).*log((x2+L/2).^2+(x3-D-W).^2) ...
    -(x2-L/2).*log((x2-L/2).^2+(x3-D).^2) ...
    +(x2+L/2).*log((x2+L/2).^2+(x3-D).^2) ...
    +2*(x3-W-D).*(atan((x2-L/2)./(x3-D-W))-atan((x2+L/2)./(x3-D-W))) ...
    +2*(x3-D).*(atan((x2+L/2)./(x3-D))-atan((x2-L/2)./(x3-D))) ...
    +(x2-L/2).*log((x2-L/2).^2+(x3+D+W).^2) ...
    -(x2+L/2).*log((x2+L/2).^2+(x3+D+W).^2) ...
    -(x2-L/2).*log((x2-L/2).^2+(x3+D).^2) ...
    +(x2+L/2).*log((x2+L/2).^2+(x3+D).^2) ...
    +2*(x3+W+D).*(atan((x2-L/2)./(x3+D+W))-atan((x2+L/2)./(x3+D+W))) ...
    +2*(x3+D).*(atan((x2+L/2)./(x3+D))-atan((x2-L/2)./(x3+D))) ...
    );

% Displacement kernels for fault slip
u1h=@(x2,x3,y2,y3,W) ...
    (+atan((x3-y3)./(x2-y2))-atan((x3+y3)./(x2-y2)) ...
     -atan((x3-y3-W)./(x2-y2))+atan((x3+y3+W)./(x2-y2)) ...
    )/2/pi;

% --      General Params      --
probL = 200e3;
probW = 200e3;
lambdaZ = 35e3; %fault depth
transition = 40e3; %where shear zone starts
Transition = transition;

ss.Ny = 51; % num elems for the shear mesh
ss.Nz = 51;

% shear mesh
eps = 1e-12;

% shear patch edges
nc = (-ss.Nz/2:ss.Nz/2);
shearZhat = transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*transition;
shearYhat = tan(nc*pi/(2.5*max(nc)))*32e3;
shearX = zeros(1,ss.Ny*ss.Nz);

% shear patch centers
shearX_c = shearX;
shearY_chat = zeros(1,ss.Ny);
shearZ_chat = zeros(1,ss.Nz);
for idx=(1:length(shearZhat)-1)
  shearZ_chat(idx) = shearZhat(idx) + abs(shearZhat(idx+1) - shearZhat(idx))/2;
  shearY_chat(idx) = shearYhat(idx) + abs(shearYhat(idx+1) - shearYhat(idx))/2;
end

% grid and flatten
shearZhat(end)=[]; shearYhat(end)=[];
[shearZ shearY] = ndgrid(shearZhat, shearYhat);
shearY = shearY(:)';
shearZ = shearZ(:)';

[shearZ_c shearY_c] = ndgrid(shearZ_chat, shearY_chat);
shearZ_c = shearZ_c(:)';
shearY_c = shearY_c(:)';

c.X = [shearX_c; shearY_c; shearZ_c];
c.Y = [shearX; shearY; shearZ];

L = shearYhat(2:end) - shearYhat(1:end-1);
W = shearZhat(2:end) - shearZhat(1:end-1);
L(length(L)+1) = L(1);
W(length(W)+1) = W(1);
c.eta=3;

% ---       Stress Kernels from Shear Zones       ---

ss.k1212 = zeros(length(shearY_chat)*length(shearZ_chat));
ss.k1213 = zeros(length(shearY_chat)*length(shearZ_chat));
ss.k1312 = zeros(length(shearY_chat)*length(shearZ_chat));
ss.k1313 = zeros(length(shearY_chat)*length(shearZ_chat));

% fields from shear zones
for ky=1:length(shearY_chat)
  for kz=1:length(shearZ_chat)
    ss.k1212(:,(kz-1)*ss.Ny+ky) = s1212(shearZ_c(kz)+transition, L(ky), W(kz), shearY_c'-shearYhat(ky)', shearZ_c');
    ss.k1213(:,(kz-1)*ss.Ny+ky) = s1213(shearZ_c(kz)+transition, L(ky), W(kz), shearY_c'-shearYhat(ky)', shearZ_c');
    ss.k1312(:,(kz-1)*ss.Ny+ky) = s1312(shearZ_c(kz)+transition, L(ky), W(kz), shearY_c'-shearYhat(ky)', shearZ_c');
    ss.k1313(:,(kz-1)*ss.Ny+ky) = s1313(shearZ_c(kz)+transition, L(ky), W(kz), shearY_c'-shearYhat(ky)', shearZ_c');
  end
end

ss.x3c = shearZ_chat;
ss.x2c = shearY_chat;

% ---     Rheology      ---
% Confining pressure (MPa) and Temperature (K)
k  = 3.138; % thermal conductivity (W/m/K)
Cp = 1171 ; % specific heat (J/kg/K)
Rm = 3330 ; % mantle density (kg/m^3)

Pconf       = Rm*9.8*ss.x3c/1e6;  % Shear zones
Pconf_fault = Rm*9.8; % Faults

Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
Age_plate = 2e15; % seconds
ss.Tprof  = 300+1380*erf(ss.x3c/(sqrt(4* Kappa * Age_plate)));  % Kelvin

% default friction properties (velocity-weakening friction)
% effective confining pressure (MPa)
ss.sigmab = 1000;

% frictional parameters
ss.aW = 1e-3;
ss.bW = ss.aW+2.1e-4;

ss.aE = 1e-3;
ss.bE = ss.aE;

% static friction coefficient
ss.mu0 = 0.2;

% characteristic weakening distance (m)
ss.L = 0.012;

% plate velocity (m/s)
ss.V_plate = 1e-9;

% reference slip rate (m/s)
ss.Vo = 1e-6;

% shear wave speed (m/s)
ss.Vs = 3e3;

% Velocity-strengthening at edges
% West fault: 5-15km velocity-weakening
%topW    = floor(5e3/(Transition/ss.M));
%bottomW = ceil(15e3/(Transition/ss.M));
%ss.bW(1:topW)      = ss.aW(1:topW)-2.1e-4*ones(topW,1);
%ss.bW(bottomW:end) = ss.aW(bottomW:end)-2.1e-4*ones(length(ss.aW(bottomW:end)),1);

% East fault: 3-20km velocity-weakening
%topE    = floor(3e3/(Transition/ss.M));
%bottomE = ceil(20e3/(Transition/ss.M));
%ss.bE(1:topE)      = ss.aE(1:topE)-2.1e-4*ones(topE,1);
%ss.bE(bottomE:end) = ss.aE(bottomE:end)-2.1e-4*ones(length(ss.aE(bottomE:end)),1);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                   R H E O L O G Y                    %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% Values taken for wet olivine - Hirth, G. and D. Kohlstedt (2003)

% Driving strain rate (1/s)
ss.e12p_plate = 1e-14*ones(length(ss.x2c)*length(ss.x3c),1);
ss.e13p_plate =      zeros(length(ss.x2c)*length(ss.x3c),1);

% Rheological Parameters
% Reference Strain Rate (for stress in MPa)
ss.Adif = 1e6*ones(length(ss.x3c)*length(ss.x2c),1);
ss.Adis = 90 *ones(length(ss.x3c)*length(ss.x2c),1);

% Power-Law Exponent
ss.n = 3.5*ones(length(ss.x3c)*length(ss.x2c),1);

% Activation Energy Wet Oliving (J/mol)
ss.Qdif = 335e3*ones(length(ss.x3c)*length(ss.x2c),1);
ss.Qdis = 480e3*ones(length(ss.x3c)*length(ss.x2c),1);

% Activation Volume (m^3/mol)
ss.Voldif = 4e-6*ones(length(ss.x3c)*length(ss.x2c),1);
ss.Voldis = 11e-6*ones(length(ss.x3c)*length(ss.x2c),1);

% Grain size (m)
ss.d    = 1e-2*ones(length(ss.x3c)*length(ss.x2c),1);
ss.pexp = 3*ones(length(ss.x3c)*length(ss.x2c),1);

% Water fugacity (H/10^6 Si)
ss.COH = 1000*ones(length(ss.x3c)*length(ss.x2c),1);
ss.r   = 1.2*ones(length(ss.x3c)*length(ss.x2c),1);

% Pressure (Pa)
ss.P = repmat(1e6*Pconf',length(ss.x2c),1);
ss.P = reshape(ss.P,[length(ss.x2c)*length(ss.x3c),1]);

% Temperature (K)
Te0 = repmat(ss.Tprof',length(ss.x2c),1);
Te0 = reshape(Te0,[length(ss.x2c)*length(ss.x3c),1]);

% Coefficients for dislocation and diffusion creep
ss.Const_dis = ss.Adis.*exp(-(ss.Qdis+ss.P.*ss.Voldis)./(8.314.*Te0)).*ss.COH.^(ss.r);
ss.Const_diff = ss.Adif.*exp(-(ss.Qdif+ss.P.*ss.Voldif)./(8.314.*Te0)).*ss.COH.^(ss.r).*ss.d.^(-ss.pexp);

% Strengh profile
s120 = (ss.e12p_plate./ss.Const_dis).^(1./ss.n);
s130 = zeros(size(s120));
e120 = zeros(size(s120));
e130 = zeros(size(s120));

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%         N U M E R I C A L   S O L U T I O N           %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% state parameters
ss.dgfF=3;
ss.dgfS=4;
ss.M=0;
%% Initialize State Vector
Y0=zeros(ss.M*ss.dgfF+ss.M*ss.dgfF+length(ss.x2c)*length(ss.x3c)*ss.dgfS,1);

% Fault patches
Y0(1:ss.dgfF:ss.M*ss.dgfF)=0;
Y0(2:ss.dgfF:ss.M*ss.dgfF)=ss.strength_W;
Y0(3:ss.dgfF:ss.M*ss.dgfF)=log(ss.Vo./ss.V_plate);

Y0(ss.M*ss.dgfF+1:ss.dgfF:2*ss.M*ss.dgfF)=0;
Y0(ss.M*ss.dgfF+2:ss.dgfF:2*ss.M*ss.dgfF)=ss.strength_E;
Y0(ss.M*ss.dgfF+3:ss.dgfF:2*ss.M*ss.dgfF)=log(ss.Vo./ss.V_plate);

% Shear zones
Y0(2*ss.M*ss.dgfF+1:ss.dgfS:end)=s120;
Y0(2*ss.M*ss.dgfF+2:ss.dgfS:end)=s130;
Y0(2*ss.M*ss.dgfF+3:ss.dgfS:end)=e120;
Y0(2*ss.M*ss.dgfF+4:ss.dgfS:end)=e130;

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) odeViscoelastic_so_d(t,y,ss);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6);
[t,Y]=ode45(yp,[0 1e4],Y0,options);
toc
%%
% Compute the instantaneous derivative
Yp=zeros(length(t)-1,size(Y,2));
for k=1:length(t)-1
    Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
end
