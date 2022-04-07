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
dz = ss.lambdaZ/ss.M;

Tranisition = 35e3;
ss.Ny = 50;
ss.Nz = 60;

% FAULT
% fault patch edges (top left)
yf = 0;
faultX = zeros(1,ss.M);
faultY = zeros(1,ss.M);
faultZ = linspace(0, ss.lambdaZ-dz, ss.M);
% tops of fault patches
ss.fpTops = faultZ';

% fault patch centers
faultX_c = faultX;
faultY_c = faultY;
faultZ_c = faultZ+(dz/2);

% width of fault patches
Lf = zeros(ss.M,1);
Wf = ones(ss.M,1)*dz;

% SHEAR
eps = 1e-12;
nc = (-ss.Ny/2:ss.Ny/2);
shearZhat = Transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*Transition;
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
%ss.fpTops = ss.fpTops;
%xx2c = shearY_c;
%xx3c = shearZ_c;
%ss.x2c = ss.shearY_chat;
%ss.x3c = ss.shearZ_chat;
%ss.polesz = shearZhat;

disp('mesh done.')
disp('beginning kernels...')

%%

% stress on fault by fault
ss.ff12=zeros(ss.M,ss.M);

% Stress kernels on shear zones from fault
ss.sf12 = zeros(length(shearY_c(:)),ss.M);
ss.sf12 = zeros(length(shearY_c(:)),ss.M);

% Fields from Faults
for k=1:ss.M
    % we evaluate the stress at the center of the fault patches
    ss.ff12(:,k)=s12h(yf, ss.fpTops+dz/2, yf, ss.fpTops(k), Wf(k));

    % we evaluate the stress at the center of the shear zones
    ss.sf12(:,k) = s12h(shearY_c(:), shearZ_c(:), yf, ss.fpTops(k), Wf(k));
    ss.sf13(:,k) = s13h(shearY_c(:), shearZ_c(:), yf, ss.fpTops(k), Wf(k));
end

% Stress kernels on shear zones from shear zones
ss.ss1312 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));
ss.ss1313 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));
ss.ss1212 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));
ss.ss1213 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));

% stress kernels on fault from shear zones
ss.fs1212 = zeros(length(ss.fpTops), length(ss.shearY_chat)*length(ss.shearZ_chat));
ss.fs1312 = zeros(length(ss.fpTops), length(ss.shearY_chat)*length(ss.shearZ_chat));

% Fields from Shear zones
for ky=1:length(ss.shearY_chat)
    for kz=1:length(ss.shearZ_chat)
        % we evaluate the stress at the center of the shear zones
        ss.k1212(:,(kz-1)*ss.Ny+ky) = s1212(shearZhat(kz), L(ky), W(kz), ...
          shearY_c(:)-ss.shearY_chat(ky), shearZ_c(:));
        ss.k1312(:,(kz-1)*ss.Ny+ky) = s1312(shearZhat(kz), L(ky), W(kz), ...
          shearY_c(:)-ss.shearY_chat(ky), shearZ_c(:));
        ss.k1213(:,(kz-1)*ss.Ny+ky) = s1213(shearZhat(kz), L(ky), W(kz), ...
          shearY_c(:)-ss.shearY_chat(ky),shearZ_c(:));
        ss.k1313(:,(kz-1)*ss.Ny+ky) = s1313(shearZhat(kz), L(ky), W(kz), ...
          shearY_c(:)-ss.shearY_chat(ky),shearZ_c(:));

        % we evaluate stress at the center of the fault patches
        ss.k1212fW(:,(kz-1)*ss.Ny+ky) = s1212(shearZhat(kz), L(ky), W(kz), ...
          yf-ss.shearY_chat(ky), ss.fpTops(:)+dz/2);
        ss.k1312fW(:,(kz-1)*ss.Ny+ky) = s1312(shearZhat(kz), L(ky), W(kz), ...
          yf-ss.shearY_chat(ky), ss.fpTops(:)+dz/2);

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

Pconf       = Rm*9.8*ss.shearZ_chat/1e6;  % Shear zones
Pconf_fault = Rm*9.8*(ss.fpTops+dz); % Faults

Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
Age_plate = 2e15; % seconds
ss.Tprof  = 300+1380*erf(ss.shearZ_chat/(sqrt(4* Kappa * Age_plate)));  % Kelvin

% default friction properties (velocity-weakening friction)
% effective confining pressure (MPa)
ss.sigmab = 1000;

% frictional parameters
ss.a = 1e-3*ones(size(ss.fpTops));
ss.b = ss.a+2.1e-4*ones(size(ss.fpTops));

% static friction coefficient
ss.mu0 = 0.2*ones(size(ss.fpTops));

% characteristic weakening distance (m)
ss.L = 0.012*ones(size(ss.fpTops));

% plate velocity (m/s)
ss.V_plate = 1e-9*ones(size(ss.fpTops));

% reference slip rate (m/s)
ss.Vo = 1e-6*ones(size(ss.fpTops));

% shear wave speed (m/s)
ss.Vs = 3e3*ones(size(ss.fpTops));

% Velocity-strengthening at edges
% fault: 5-15km velocity-weakening
top    = floor(5e3/(Transition/ss.M));
bottom = ceil(15e3/(Transition/ss.M));
ss.b(1:top)      = ss.a(1:top)-2.1e-4*ones(top,1);
ss.b(bottom:end) = ss.a(bottom:end)-2.1e-4*ones(length(ss.a(bottom:end)),1);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                   R H E O L O G Y                    %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% Values taken for wet olivine - Hirth, G. and D. Kohlstedt (2003)

% Driving strain rate (1/s)
ss.e12p_plate = 1e-14*ones(length(ss.shearY_chat)*length(ss.shearZ_chat),1);
ss.e13p_plate =      zeros(length(ss.shearY_chat)*length(ss.shearZ_chat),1);

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
ss.Const_dis = ss.Adis.*exp(-(ss.Qdis+ss.P.*ss.Voldis)./(8.314.*Te0)) ...
  .*ss.COH.^(ss.r);
ss.Const_diff = ss.Adif.*exp(-(ss.Qdif+ss.P.*ss.Voldif)./(8.314.*Te0)) ...
  .*ss.COH.^(ss.r).*ss.d.^(-ss.pexp);

% Strengh profile
s120 = (ss.e12p_plate./ss.Const_dis).^(1./ss.n);
s130 = zeros(size(s120));
e120 = zeros(size(s120));
e130 = zeros(size(s120));

% Fault Strength
ss.strength_W = ss.sigmab.*(ss.mu0+(ss.a-ss.b).*log(ss.V_plate./ss.Vo))+ ...
  G*ss.V_plate./(2*ss.Vs);

% Plot strength profiles
figure(1);clf;
subplot(2,1,1);
plot(ss.fpTops/1e3,ss.strength_E,ss.fpTops/1e3,ss.strength_W)
xlabel('Depth (km)')
ylabel('Strength (MPa)');
subplot(2,1,2);
plot(ss.shearZ_chat/1e3,log10(s120(1:length(ss.shearY_chat):end)))
xlim([ss.shearZ_chat(1)/1e3,ss.shearZ_chat(end)/1e3]);
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
Y0=zeros(ss.M*ss.dgfF+length(ss.shearY_chat)*length(ss.shearZ_chat)*ss.dgfS,1);

% Fault patches
Y0(1:ss.dgfF:ss.M*ss.dgfF)=zeros(size(ss.fpTops));
Y0(2:ss.dgfF:ss.M*ss.dgfF)=ss.strength_W;
Y0(3:ss.dgfF:ss.M*ss.dgfF)=log(ss.Vo./ss.V_plate);

% Shear zones
Y0(ss.M*ss.dgfF+1:ss.dgfS:end)=s120;
Y0(ss.M*ss.dgfF+2:ss.dgfS:end)=s130;
Y0(ss.M*ss.dgfF+3:ss.dgfS:end)=e120;
Y0(ss.M*ss.dgfF+4:ss.dgfS:end)=e130;

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) ode2Faults2(t,y,ss);
disp('begin solving...')
tic
% Solve the system
options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6);
[t,Y]=ode45_2(yp,[0 1*3.15e7],Y0,options);
disp('done solving.')
toc
%%

Y = Y';

% Compute the instantaneous derivative
Yp=zeros(length(t)-1,size(Y,2));
for k=1:length(t)-1
    Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Strain rate at center
Ep=sqrt(Yp(:,ss.M*ss.dgfF+floor(length(ss.shearY_chat)/2)* ...
          ss.dgfS+3:ss.dgfS*length(ss.shearY_chat):end)'.^2 +...
        Yp(:,ss.M*ss.dgfF+floor(length(ss.shearY_chat)/2)* ...
          ss.dgfS+4:ss.dgfS*length(ss.shearY_chat):end)'.^2);

% strain rate over whole ductile area
Epall = sqrt( Yp(:,ss.M*ss.dgfF+3:ss.dgfS:end)'.^2 +...
             Yp(:,ss.M*ss.dgfF+4:ss.dgfS:end)'.^2);

% Velocity
y.V = Yp(:,1:ss.dgfF:ss.M*ss.dgfF); % Slip rate (m/s)
y.V = y.V';
y.tau = Y(:,2:ss.dgfF:ss.M*ss.dgfF);            % Shear stress (MPa)

% fault slip figure
clf;
imagesc(log10(y.V)); colorbar;
title('Slip Rate')
xlabel('time steps')
ylabel('fault mesh block')
saveas(gcf, 'figures/BP1vD2_slip.png')

% strain rate at center of ductile region
clf;
imagesc(log10(Ep)); colorbar;
title('Strain Rate of Center of Ductile Region')
xlabel('Time Steps')
ylabel('Block')
saveas(gcf, 'figures/BP1vD2_strainCenter.png')

% ---         Movies        ---
Smovie=false;
if Smovie
  disp('begin shear movie')
  clf;
  fig = figure;
  fname = 'figures/BP1vD2_strain.gif';
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
