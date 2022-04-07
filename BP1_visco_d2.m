clear all;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%            Unified earthquake cycles of               %
%          fault slip and viscoelastic strain           %
%                                                       %
%                  Valere Lambert, 2016                 %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Simulates unified earthquake cycles including slip on
% two parallel strike-slip faults and viscoelastic strain
% beneath the brittle-ductile transition in 2D antiplane.

% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%            P H Y S I C A L   M O D E L                %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% rigidity (MPa)
G=30e3;

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

% Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

s13h=@(x2,x3,y2,y3,Wf) G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
    -(x2-y2)./((x2-y2).^2+(x3-y3-Wf).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        M E S H                       %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %

%            -25 km       0       25 km
% _______________________________________________ 0 km    --------> x2
%               |                   |                    |
%               |                   |                    |
%               |      Brittle      |                    |
%               |                   |                    V
%               |                   |                    x3
%          West | fault        East | fault
% ----------------------------------------------- 35 km
%
%              Viscoelastic Shear Zones
%

disp('begin mesh...')

% Fault Meshes
probL = 200e3;
probW = 200e3;

ss.lambdaZ = 40e3; % fault depth extent
ss.M = 400; %number of fault cells
ss.dz = ss.lambdaZ/ss.M; dz = ss.dz;

ss.transition = 35e3; Transition = ss.transition;
ss.Ny = 50;
ss.Nz = 60;
ss.Nx = ss.Nz;

% FAULT
% fault patch edges (top left)
yf = 0;
faultX = zeros(1,ss.M);
faultY = zeros(1,ss.M);
faultZ = linspace(0, ss.lambdaZ-ss.dz, ss.M);
% tops of fault patches
ss.fpTops = faultZ';

% fault patch centers
faultX_c = faultX;
faultY_c = faultY;
faultZ_c = faultZ+(ss.dz/2);

% width of fault patches
Lf = zeros(ss.M,1);
Wf = ones(ss.M,1)*ss.dz;

% SHEAR
eps = 1e-12;
nc = (-ss.Ny/2:ss.Ny/2);
shearZhat = ss.transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*ss.transition;
shearYhat = tan(nc*pi/(2.5*max(nc)))*32e3;
shearX = zeros(1,ss.Ny*ss.Nz);

% shear patch centers
shearX_c = shearX;
ss.shearY_chat = zeros(ss.Ny,1);
ss.shearZ_chat = zeros(ss.Nz,1);

% sizing
L = zeros(1,ss.Ny);
W = zeros(1,ss.Nz);

for idx=(1:length(shearZhat)-1)
  ss.shearZ_chat(idx) = shearZhat(idx) + abs(shearZhat(idx+1) - shearZhat(idx))/2;
  W(idx) = abs(shearZhat(idx) - shearZhat(idx+1));
end
for idx=(1:length(shearYhat)-1)
  ss.shearY_chat(idx) = shearYhat(idx) + abs(shearYhat(idx+1) - shearYhat(idx))/2;
  L(idx) = abs(shearYhat(idx) - shearYhat(idx+1));
end

L(end) = L(1);
W(end) = abs(shearZhat(end-1) - shearZhat(end));

% grid and flatten
shearZhat(end)=[]; shearYhat(end)=[];
[shearY shearZ] = ndgrid(shearYhat, shearZhat);

[shearY_c shearZ_c] = ndgrid(ss.shearY_chat, ss.shearZ_chat);

% convert between naming conventions
% fault:
ss.y3f = ss.fpTops;
xx2c = shearY_c;
xx3c = shearZ_c;
ss.x2c = ss.shearY_chat;
ss.x3c = ss.shearZ_chat;
ss.polesz = shearZhat;

disp('mesh done.')
disp('beginning kernels...')

%%
% Stress kernels from fault
ss.k12W=zeros(length(xx2c(:)),ss.M);
ss.k13W=zeros(length(xx2c(:)),ss.M);

ss.KWW=zeros(ss.M,ss.M);

% Fields from Faults
for k=1:ss.M
    % we evaluate the stress at the center of the shear zones
    ss.k12W(:,k)=s12h(xx2c(:), xx3c(:), yf, ss.y3f(k), Wf(k));
    ss.k13W(:,k)=s13h(xx2c(:), xx3c(:), yf, ss.y3f(k), Wf(k));

    % we evaluate the stress at the center of the fault patches
    ss.KWW(:,k)=s12h(yf, ss.y3f+dz/2, yf, ss.y3f(k), Wf(k));
end

% Stress kernels from shear zones
ss.k1312=zeros(length(ss.x2c)*length(ss.x3c));
ss.k1313=zeros(length(ss.x2c)*length(ss.x3c));
ss.k1212=zeros(length(ss.x2c)*length(ss.x3c));
ss.k1213=zeros(length(ss.x2c)*length(ss.x3c));

ss.k1212fW=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));
ss.k1312fW=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));

% Fields from Shear zones
for kx=1:length(ss.x2c)
    for kz=1:length(ss.x3c)
        % we evaluate the stress at the center of the shear zones
        ss.k1212(:,(kz-1)*ss.Nx+kx)=s1212(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));
        ss.k1312(:,(kz-1)*ss.Nx+kx)=s1312(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));
        ss.k1213(:,(kz-1)*ss.Nx+kx)=s1213(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));
        ss.k1313(:,(kz-1)*ss.Nx+kx)=s1313(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));

        % we evaluate stress at the center of the fault patches
        ss.k1212fW(:,(kz-1)*ss.Nx+kx)=s1212(ss.polesz(kz),L(kx),W(kz),yf-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1312fW(:,(kz-1)*ss.Nx+kx)=s1312(ss.polesz(kz),L(kx),W(kz),yf-ss.x2c(kx),ss.y3f(:)+dz/2);

    end
end

disp('kernels done.')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Confining pressure (MPa) and Temperature (K)
k  = 3.138; % thermal conductivity (W/m/K)
Cp = 1171 ; % specific heat (J/kg/K)
Rm = 3330 ; % mantle density (kg/m^3)

Pconf       = Rm*9.8*ss.x3c/1e6;  % Shear zones
Pconf_fault = Rm*9.8*(ss.y3f+dz); % Faults

Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
Age_plate = 2e15; % seconds
ss.Tprof  = 300+1380*erf(ss.x3c/(sqrt(4* Kappa * Age_plate)));  % Kelvin

% default friction properties (velocity-weakening friction)
% effective confining pressure (MPa)
ss.sigmab = 1000;

% frictional parameters
ss.aW = 1e-3*ones(size(ss.y3f));
ss.bW = ss.aW+2.1e-4*ones(size(ss.y3f));

ss.aE = 1e-3*ones(size(ss.y3f));
ss.bE = ss.aE+2.1e-4*ones(size(ss.y3f));

% static friction coefficient
ss.mu0 = 0.2*ones(size(ss.y3f));

% characteristic weakening distance (m)
ss.L = 0.012*ones(size(ss.y3f));

% plate velocity (m/s)
ss.V_plate = 1e-9*ones(size(ss.y3f));

% reference slip rate (m/s)
ss.Vo = 1e-6*ones(size(ss.y3f));

% shear wave speed (m/s)
ss.Vs = 3e3*ones(size(ss.y3f));

% Velocity-strengthening at edges
% fault: 5-15km velocity-weakening
topW    = floor(5e3/(Transition/ss.M));
bottomW = ceil(15e3/(Transition/ss.M));
ss.bW(1:topW)      = ss.aW(1:topW)-2.1e-4*ones(topW,1);
ss.bW(bottomW:end) = ss.aW(bottomW:end)-2.1e-4*ones(length(ss.aW(bottomW:end)),1);

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

% Fault Strength
ss.strength_E = ss.sigmab.*(ss.mu0+(ss.aE-ss.bE).*log(ss.V_plate./ss.Vo))+G*ss.V_plate./(2*ss.Vs);
ss.strength_W = ss.sigmab.*(ss.mu0+(ss.aW-ss.bW).*log(ss.V_plate./ss.Vo))+G*ss.V_plate./(2*ss.Vs);

% Plot strength profiles
figure(1);clf;
subplot(2,1,1);
plot(ss.y3f/1e3,ss.strength_E,ss.y3f/1e3,ss.strength_W)
xlabel('Depth (km)')
ylabel('Strength (MPa)');
subplot(2,1,2);
plot(ss.x3c/1e3,log10(s120(1:length(ss.x2c):end)))
xlim([ss.x3c(1)/1e3,ss.x3c(end)/1e3]);
xlabel('Depth (km)')
ylabel('Strength (MPa) log10');


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%         N U M E R I C A L   S O L U T I O N           %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% state parameters
ss.dgfF=3;
ss.dgfS=4;
%% Initialize State Vector
Y0=zeros(ss.M*ss.dgfF+length(ss.x2c)*length(ss.x3c)*ss.dgfS,1);

% Fault patches
Y0(1:ss.dgfF:ss.M*ss.dgfF)=zeros(size(ss.y3f));
Y0(2:ss.dgfF:ss.M*ss.dgfF)=ss.strength_W;
Y0(3:ss.dgfF:ss.M*ss.dgfF)=log(ss.Vo./ss.V_plate);

% Shear zones
Y0(ss.M*ss.dgfF+1:ss.dgfS:end)=s120;
Y0(ss.M*ss.dgfF+2:ss.dgfS:end)=s130;
Y0(ss.M*ss.dgfF+3:ss.dgfS:end)=e120;
Y0(ss.M*ss.dgfF+4:ss.dgfS:end)=e130;

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) odeBP1v_d(t,y,ss);
disp('begin solving...')
tic
% Solve the system
options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6);
[t,Y]=ode45_2(yp,[0 1*3.15e7],Y0,options);
disp('done solving.')
toc
%%
% Compute the instantaneous derivative
Yp=zeros(length(t)-1,size(Y,2));
for k=1:length(t)-1
    Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
end

Y = Y';

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Strain rate at center
Ep=sqrt(Yp(:,2*ss.M*ss.dgfF+floor(length(ss.x2c)/2)*ss.dgfS+3:ss.dgfS*length(ss.x2c):end)'.^2 +...
        Yp(:,2*ss.M*ss.dgfF+floor(length(ss.x2c)/2)*ss.dgfS+4:ss.dgfS*length(ss.x2c):end)'.^2);

% strain rate over whole ductile area
Epall = sqrt( Yp(:,ss.M*ss.dgfF+3:ss.dgfS:end)'.^2 +...
             Yp(:,ss.M*ss.dgfF+4:ss.dgfS:end)'.^2);

% Velocity
y.V = ss.Vo.*exp(Y(:,4:ss.dgfF:ss.M*ss.dgfF)'); % Slip rate (m/s)
y.tau = Y(:,2:ss.dgfF:ss.M*ss.dgfF);            % Shear stress (MPa)
y.Vmax = zeros(length(t),1);          % Maximum slip rate (m/s)
y.Vcenter = y.V(floor(ss.M/2),:);          % Slip rate at center of VW region
for ti = 1:length(t)
  y.Vmax(ti) = max(y.V(:,ti));
end

% fault slip figure
clf;
imagesc(log10(y.V)); colorbar;
title('Slip Rate')
xlabel('time steps')
ylabel('fault mesh block')
saveas(gcf, 'figures/BP1vD_slip.png')

% ---         Movies        ---
Smovie=true;
if Smovie
  disp('begin shear movie')
  clf;
  fig = figure;
  fname = 'figures/BP1vD_strain.gif';
  for idx = 1:size(Epall, 2)
    oneE = Epall(:,idx);
    oneEsq = reshape(oneE, [ss.Ny, ss.Nz]);
    imagesc(oneEsq'); colorbar;
    title(idx)
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);

    [A,map] = rgb2ind(im{idx},256);
    if idx==1
      imwrite(A,map,fname,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
      imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',0.1);
    end

  end
disp('shear movie done')
end
