% BP1_visco_m.m
% dense matrix version of BP1_visco
% ---     Kernels & functions      ---

G = 30e3;

% Boxcar function
boxc=@(x) (x+0.5>=0)-(x-0.5>=0);

% Heaviside function
HS=@(x) 0+x>=0;

% ramp function
Ramp=@(x) x.*BC(x-1/2)+HS(x-1);

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

    % Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
  -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
  +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
  )/2/pi;

s13h=@(x2,x3,y2,y3,Wf) G*( ...
  (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
  -(x2-y2)./((x2-y2).^2+(x3-y3-Wf).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+Wf).^2) ...
  )/2/pi;

% Displacement kernels for fault slip
u1h=@(x2,x3,y2,y3,W) ...
    (+atan((x3-y3)./(x2-y2))-atan((x3+y3)./(x2-y2)) ...
     -atan((x3-y3-W)./(x2-y2))+atan((x3+y3+W)./(x2-y2)) ...
    )/2/pi;

% ---       General params      ---
probL = 200e3;
probW = 200e3;

ss.lambdaZ = 40e3; % fault depth extent
ss.M = 400; %number of fault cells
ss.dz = ss.lambdaZ/ss.M;

ss.transition = 40e3;
ss.Ny = 51;
ss.Nz = 51;
ss.Nx = ss.Nz;

% FAULT
% width of fault patches
Wf = ones(ss.M,1)*ss.dz;
% fault patch edges (top left)
faultX = zeros(1,ss.M);
faultY = zeros(1,ss.M);
faultZ = linspace(0, ss.lambdaZ-ss.dz, ss.M);
% tops of fault patches
ss.fpTops = faultZ';

% fault patch centers
faultX_c = faultX;
faultY_c = faultY;
faultZ_c = faultZ+(ss.dz/2);

% SHEAR
eps = 1e-12;
nc = (-ss.Nz/2:ss.Nz/2);
shearZhat = ss.transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*ss.transition;
shearYhat = tan(nc*pi/(2.5*max(nc)))*32e3;
shearX = zeros(1,ss.Ny*ss.Nz);

% shear patch centers
shearX_c = shearX;
ss.shearY_chat = zeros(1,ss.Ny);
ss.shearZ_chat = zeros(1,ss.Nz);
for idx=(1:length(shearZhat)-1)
  ss.shearZ_chat(idx) = shearZhat(idx) + abs(shearZhat(idx+1) - shearZhat(idx))/2;
  ss.shearY_chat(idx) = shearYhat(idx) + abs(shearYhat(idx+1) - shearYhat(idx))/2;
end

% grid and flatten
shearZhat(end)=[]; shearYhat(end)=[];
[shearZ shearY] = ndgrid(shearZhat, shearYhat);
shearY = shearY(:)';
shearZ = shearZ(:)';

[shearZ_c shearY_c] = ndgrid(ss.shearZ_chat, ss.shearY_chat);
shearZ_c = shearZ_c(:)';
shearY_c = shearY_c(:)';

% combo mesh
comboX = [faultX shearX];
comboY = [faultY shearY];
comboZ = [faultZ shearZ];

comboX_c = [faultX_c shearX_c];
comboY_c = [faultY_c shearY_c];
comboZ_c = [faultZ_c shearZ_c];

disp('mesh created')

% ---         Stress Kernels        ---

ss.k12 = zeros(length(shearY_c),ss.M);
ss.k13 = zeros(length(shearY_c),ss.M);

ss.k12f= zeros(ss.M,ss.M);

ss.k1212 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));
ss.k1213 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));
ss.k1312 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));
ss.k1313 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));

ss.k1212f= zeros(length(ss.fpTops), length(shearY_c)*length(faultZ_c));
ss.k1312f= zeros(length(ss.fpTops), length(shearY_c)*length(faultZ_c));
ss.k1213f= zeros(length(ss.fpTops), length(shearY_c)*length(faultZ_c));
ss.k1313f= zeros(length(ss.fpTops), length(shearY_c)*length(faultZ_c));

disp('beginning kernels')

% fields from faults
for k=1:ss.M
  ss.k12(1:ss.M,k)=s12h(faultY_c(:), faultZ_c(:), 0, ss.fpTops(k), Wf(k));
  ss.k13(1:ss.M,k)=s13h(faultY_c(:), faultZ_c(:), 0, ss.fpTops(k), Wf(k));

  ss.k12f(:,k)=s12h(faultY,ss.fpTops+ss.dz/2,faultY,ss.fpTops(k),Wf(k));
end

% fields from shear zones
for ky=1:length(shearY_chat)
  for kz=1:length(shearZ_chat)
    ss.k1212(:,(kz-1)*ss.Ny+ky) = s1212(shearZ_c(kz)+transition, L(ky), W(kz), ...
      shearY_c'-shearYhat(ky)', shearZ_c');
    ss.k1213(:,(kz-1)*ss.Ny+ky) = s1213(shearZ_c(kz)+transition, L(ky), W(kz), ...
      shearY_c'-shearYhat(ky)', shearZ_c');
    ss.k1312(:,(kz-1)*ss.Ny+ky) = s1312(shearZ_c(kz)+transition, L(ky), W(kz), ...
      shearY_c'-shearYhat(ky)', shearZ_c');
    ss.k1313(:,(kz-1)*ss.Ny+ky) = s1313(shearZ_c(kz)+transition, L(ky), W(kz), ...
      shearY_c'-shearYhat(ky)', shearZ_c');

    ss.k1212f(:,(kz-1)*ss.Nx+kx)=s1212(shearZ_c(kz))

  end
end

disp('kernels created')

% ---     Rheology      ---
% Confining pressure (MPa) and Temperature (K)
k  = 3.138; % thermal conductivity (W/m/K)
Cp = 1171 ; % specific heat (J/kg/K)
Rm = 3330 ; % mantle density (kg/m^3)

Pconf       = Rm*9.8*shearZ_chat/1e6;  % Shear zones
Pconf_fault = Rm*9.8; % Faults

Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
Age_plate = 2e15; % seconds
ss.Tprof  = 300+1380*erf(shearZ_chat/(sqrt(4* Kappa * Age_plate)));  % Kelvin

% ---       fault properties      ---
%  -      FRICTIONAL PARAMETERS    -

% density (kg/m^3)
rho = 2670;
% shear wave speed (m/s)
Vs = 3464;
% shear modulus (MPa)
G = rho*Vs^2/1e6;

% reference friction coefficient
ss.fo=0.6*ones(size(ss.fpTops));

% Dieterich-Ruina R+S frictional parameters (velocity-weakening friction)
ss.a=1e-2+Ramp((ss.fpTops-15e3)/3e3)*(0.025-0.01);
ss.b=0.015*ones(size(ss.fpTops));

% TODO: is the sigma the same as sigmab?
% effective normal stress (MPa)
ss.sigma=50.0*ones(size(ss.fpTops));

% characteristic weakening distance (m)
ss.Drs=8e-3*ones(size(ss.fpTops));

% plate rate (m/s)
ss.Vpl=1e-9*ones(size(ss.fpTops));

% reference slip rate (m/s)
ss.Vo=1e-6*ones(size(ss.fpTops));

% Radiation damping coefficient
ss.eta = G./(2*Vs);

% static friction coefficient
ss.mu0 = 0.2*ones(size(ss.fpTops));

% fault strength
ss.strength = ss.sigma.*(ss.mu0+(ss.a-ss.b).*log(ss.Vpl./ss.Vo))+...
  G*ss.Vpl./(2*Vs);

% Estimates of some key parameters
VWp = find(ss.a < ss.b); % VW region
% Critical nucleation size ( h* = pi/2 G b D_rs / (b-a)^2 / sigma )
hstar=min(pi/2*G*ss.Drs(VWp).*ss.b(VWp)./(ss.b(VWp)-ss.a(VWp)).^2./ss.sigma(VWp));

% Quasi-static cohesive zone ( coh0 = 9/32 G D_rs /(b*sigma) )
% Note that for this QD simulation the cohesive zone will not change,
% which would not be the case for a fully dynamic simulation
coh = min(9/32*pi*G*ss.Drs(VWp)./ss.b(VWp)./ss.sigma(VWp));

% Estimate of recurrence time ( T ~ 5(b-a)*sigma / G * R/Vpl )
Ti = 5*mean((ss.b(VWp)-ss.a(VWp)).*ss.sigma(VWp)).*0.5.*(ss.fpTops(VWp(end))...
  -ss.fpTops(VWp(1)))./(G*mean(ss.Vpl(VWp)));

% Print information about discretization
fprintf('\nGrid size = %.2f (m)\n', ss.dz);
fprintf('VW zone = %.2f (km)\n', (ss.fpTops(VWp(end))-ss.fpTops(VWp(1)))/1e3);
fprintf('Critical nucleation size = %.2f (m)\n',hstar);
fprintf('QS Cohesive zone = %.2f (m)\n',coh);
fprintf('Est. Recurrence time = %.2f (yr)\n', Ti/3.15e7);

% ---         Visco properties      ---
%  -              RHEOLOGY           -

% Driving strain rate (1/s)
ss.e12p_plate = 1e-14*ones(length(ss.shearY_chat)*length(ss.shearZ_chat),1);
ss.e13p_plate =      zeros(length(ss.shearY_chat)*length(ss.shearZ_chat),1);

% Confining pressure (MPa) and Temperature (K)
k  = 3.138; % thermal conductivity (W/m/K)
Cp = 1171 ; % specific heat (J/kg/K)
Rm = 3330 ; % mantle density (kg/m^3)

Pconf       = Rm*9.8*ss.shearZ_chat/1e6;  % Shear zones
Pconf_fault = Rm*9.8*(faultZ'+ss.dz); % Faults

Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
Age_plate = 2e15; % seconds
ss.Tprof  = 300+1380*erf(ss.shearZ_chat/(sqrt(4* Kappa * Age_plate)));  % Kelvin

% Rheological Parameters
% Reference Strain Rate (for stress in MPa)
ss.Adif = 1e6*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);
ss.Adis = 90 *ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);

% Power-Law Exponent
ss.n = 3.5*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);

% Activation Energy Wet Oliving (J/mol)
ss.Qdif = 335e3*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);
ss.Qdis = 480e3*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);

% Activation Volume (m^3/mol)
ss.Voldif = 4e-6*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);
ss.Voldis = 11e-6*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);

% Grain size (m)
ss.d    = 1e-2*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);
ss.pexp = 3*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);

% Water fugacity (H/10^6 Si)
ss.COH = 1000*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);
ss.r   = 1.2*ones(length(ss.shearZ_chat)*length(ss.shearY_chat),1);

% Pressure (Pa)
ss.P = repmat(1e6*Pconf',length(ss.shearY_chat),1);
ss.P = reshape(ss.P,[length(ss.shearY_chat)*length(ss.shearZ_chat),1]);

% Temperature (K)
Te0 = repmat(ss.Tprof',length(ss.shearY_chat),1);
Te0 = reshape(Te0,[length(ss.shearY_chat)*length(ss.shearZ_chat),1]);

% Coefficients for dislocation and diffusion creep
ss.Const_dis = ss.Adis.*exp(-(ss.Qdis+ss.P.*ss.Voldis)./(8.314.*Te0)).*ss.COH.^(ss.r);
ss.Const_diff = ss.Adif.*exp(-(ss.Qdif+ss.P.*ss.Voldif)./(8.314.*Te0)).* ...
  ss.COH.^(ss.r).*ss.d.^(-ss.pexp);

% Strengh profile
ss.s120 = (ss.e12p_plate./ss.Const_dis).^(1./ss.n);
ss.s130 = zeros(size(ss.s120));
ss.e120 = zeros(size(ss.s120));
ss.e130 = zeros(size(ss.s120));

% characteristic weakening distance (m)
ss.L = 0.012;

ss.dgfF = 4;
ss.dgfS = 4;

% state vector init
Y0=zeros(r.ss.M*r.ss.dgfF+length(r.ss.shearY_chat)*length(r.ss.shearZ_chat)*r.ss.dgfS,1);

% Fault patches
Y0(1:r.ss.dgfF:r.ss.M*r.ss.dgfF)=zeros(size(r.ss.fpTops));
Y0(2:r.ss.dgfF:r.ss.M*r.ss.dgfF)=r.ss.strength;
Y0(3:r.ss.dgfF:r.ss.M*r.ss.dgfF)=r.ss.a./r.ss.b.*log(2*r.ss.Vo./r.ss.Vpl.*sinh((Y0(2:r.ss.dgfF:r.ss.M*r.ss.dgfF)- ...
  r.ss.eta.*r.ss.Vpl)./r.ss.a./r.ss.sigma))-r.ss.fo./r.ss.b;
Y0(4:r.ss.dgfF:r.ss.M*r.ss.dgfF)=log(r.ss.Vo./r.ss.Vpl);

% Shear zones
Y0(r.ss.M*r.ss.dgfF+1:r.ss.dgfS:end)=r.ss.s120;
Y0(r.ss.M*r.ss.dgfF+2:r.ss.dgfS:end)=r.ss.s130;
Y0(r.ss.M*r.ss.dgfF+3:r.ss.dgfS:end)=r.ss.e120;
Y0(r.ss.M*r.ss.dgfF+4:r.ss.dgfS:end)=r.ss.e130;

% initialize the function handle with set constitutive parameters
yp=@(t,y) odeBP1v(t,y,r.ss, hm);

% ODE45 Settings
% Initial step of 1e-5 seconds
% Relative tolerance of 3e-8
% [0 3e10] = simulate 3e10 seconds, 3.15e7 seconds / year
disp('begin solving')
tic
options=odeset('Refine',1,'RelTol',1e-8,'InitialStep',1e-5);
[t,Y]=ode45_2(yp,[0 1*3.15e7],Y0,options);
disp('Done solving');
toc

Y=Y';

% ---       Figures        ---
disp(size(r.ss.Vo))
disp(size(Y(:,3:r.ss.dgfF:r.ss.M*r.ss.dgfF)'))
y.V = r.ss.Vo.*exp(-Y(:,3:r.ss.dgfF:r.ss.M*r.ss.dgfF)'); % Slip rate (m/s)
y.tau = Y(:,2:r.ss.dgfF:r.ss.M*r.ss.dgfF);            % Shear stress (MPa)
y.Vmax = zeros(length(t),1);          % Maximum slip rate (m/s)
y.Vcenter = y.V(floor(r.ss.M/2),:);          % Slip rate at center of VW region
for ti = 1:length(t)
  y.Vmax(ti) = max(y.V(:,ti));
end

clf;
imagesc(log10(y.V)); colorbar;
title('Slip Rate')
xlabel('time steps')
ylabel('fault mesh block')
saveas(gcf, 'figures/BP1v_slip.png')
