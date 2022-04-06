% BP1_visco_m.m
% dense matrix version of BP1_visco

function varargout = bp1vD (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ----------------------------------- Public ---------------------------------
function out = run()

  % ---     Kernels & functions      ---

  addpath('../hmmvp-okada/matlab')
  addpath('ODEsolving')

  G = 30e3;

  % Boxcar function
  boxc=@(x) (x+0.5>=0)-(x-0.5>=0);

  % Heaviside function
  HS=@(x) 0+x>=0;

  % ramp function
  Ramp=@(x) x.*boxc(x-1/2)+HS(x-1);

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

  ss.transition = 35e3;
  ss.Ny = 50;
  ss.Nz = 60;
  ss.Nx = ss.Nz;

  % FAULT
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
  ss.shearY_chat = zeros(1,ss.Ny);
  ss.shearZ_chat = zeros(1,ss.Nz);

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
  shearY = shearY(:)';
  shearZ = shearZ(:)';

  [shearY_c shearZ_c] = ndgrid(ss.shearY_chat, ss.shearZ_chat);
  shearZ_c = shearZ_c(:)';
  shearY_c = shearY_c(:)';

  % TWOFAULTS MESH BEGIN

  y2_W=-25e3;
  y3=0e3;

  % East Fault (m)
  y2_E=25e3;

  % Brittle-Ductile Tranisition Depth (m)
  Transition=35e3;

  % Number of patches
  ss.M=400;

  dz=Transition/ss.M;
  fpoles=y3+(0:ss.M)'*dz;
  % Tops of fault patches
  ss.y3f=fpoles(1:end-1);
  ss.fpTops = ss.y3f; % adapt for my naming scheme
  % Width of fault patches
  Wf=ones(ss.M,1)*dz;
  %%
  % Shear Zone Mesh
  ss.Nx=50;
  ss.Nz=50;
  eps=1e-12;

  % Shear zone Grid
  % edges along x3
  ss.polesz=Transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*Transition;

  % center of shear zone (x3)
  ss.x3c=(ss.polesz(2:end)+ss.polesz(1:end-1))/2;

  % shear zone width
  W=ss.polesz(2:end)-ss.polesz(1:end-1);

  xx3t=repmat(ss.polesz(1:end-1)',ss.Nx,1);  % tops
  xx3b=repmat(ss.polesz(2:end)',  ss.Nx,1);  % bottoms
  xx3c=repmat(ss.x3c',ss.Nx,1);              % centers

  % edges along x2
  ss.polesxc=(2*y2_W/1e3:(2*y2_E-2*y2_W)/(26e3):2*y2_E/1e3)'*1e3;
  edges= floor(ss.Nx-length(ss.polesxc)+1)/2;
  ss.polesxl=flipud(min(ss.polesxc)-tan((0:edges)'*pi/(2.2*(edges)+eps))*Transition);
  ss.polesxr=max(ss.polesxc)+tan((0:edges)'*pi/(2.2*(edges)+eps))*Transition;
  ss.polesx=[ss.polesxl(1:end-1);ss.polesxc;ss.polesxr(2:end)];

  % center of shear zone (x2)
  ss.x2c=(ss.polesx(2:end)+ss.polesx(1:end-1))/2;

  % shear zone length
  L=ss.polesx(2:end)-ss.polesx(1:end-1);
  xx2l=repmat(ss.polesx(1:end-1),1,length(ss.x3c));  % left edge
  xx2r=repmat(ss.polesx(2:end),  1,length(ss.x3c));  % right edge
  xx2c=repmat(ss.x2c,1,length(ss.x3c));              % center

  % adapt naming convention
  shearZhat = ss.polesz;
  shearY_c = xx2c(:);
  shearZ_c = xx3c(:)';
  ss.shearY_chat = ss.x2c';
  ss.shearZ_chat = ss.x3c';


  % TWOFAULTS MESH END

  % plot mesh
  clf;
  hold on;
  scatter(shearY, -1*shearZ, 0.25, 'red');
  scatter(faultY, -1*faultZ, 0.25, 'red');
  scatter(shearY_c, -1*shearZ_c, 0.25, 'blue');
  scatter(faultY_c, -1*faultZ_c, 0.25, 'blue');
  hold off;
  saveas(gcf, 'figures/BP1vD_mesh.png');

  disp('mesh created')

  % ---         Stress Kernels        ---

  % fault - fault
  ss.ff12= zeros(ss.M,ss.M);

  % shear - shear
  ss.ss1212 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));
  ss.ss1213 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));
  ss.ss1312 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));
  ss.ss1313 = zeros(length(ss.shearY_chat)*length(ss.shearZ_chat));

  % on shear by fault
  ss.sf12 = zeros(length(shearY_c),ss.M);
  ss.sf13 = zeros(length(shearY_c),ss.M);

  % on fault by shear
  ss.fs1212 = zeros(length(ss.fpTops), ss.Ny*ss.Nz);
  ss.fs1312 = zeros(length(ss.fpTops), ss.Ny*ss.Nz);

  disp('beginning kernels')
  % fields from faults
  for k=1:ss.M
    % stress at center of fault patches
    ss.ff12(:,k)=s12h(0, ss.fpTops+ss.dz/2, 0, ss.fpTops(k), Wf(k));

    % stress at center of shear zones
    ss.sf12(:,k)=s12h(shearY_c(:), shearZ_c', 0, ss.fpTops(k), Wf(k));
    ss.sf13(:,k)=s13h(shearY_c(:), shearZ_c', 0, ss.fpTops(k), Wf(k));

  end

  %disp(size(ss.fs1212))

  % fields from shear zones
  % fills entire columns at a time. ie all strain from a given source
  for ky=1:length(ss.shearY_chat)
    for kz=1:length(ss.shearZ_chat)

      % stress at center of shear zones
      ss.ss1212(:,(kz-1)*ss.Ny+ky) = s1212(shearZhat(kz), L(ky), W(kz), ...
        shearY_c(:)-ss.shearY_chat(ky)', shearZ_c');
      ss.ss1213(:,(kz-1)*ss.Ny+ky) = s1213(shearZhat(kz), L(ky), W(kz), ...
        shearY_c(:)-ss.shearY_chat(ky)', shearZ_c');
      ss.ss1312(:,(kz-1)*ss.Ny+ky) = s1312(shearZhat(kz), L(ky), W(kz), ...
        shearY_c(:)-ss.shearY_chat(ky)', shearZ_c');
      ss.ss1313(:,(kz-1)*ss.Ny+ky) = s1313(shearZhat(kz), L(ky), W(kz), ...
        shearY_c(:)-ss.shearY_chat(ky)', shearZ_c');

      % stress at center of fault patches
      ss.fs1212(:,(kz-1)*ss.Ny+ky)=s1212(shearZhat(kz), L(ky), W(kz), ...
        y2_W-ss.shearY_chat(ky)', ss.fpTops+ss.dz/2);
      ss.fs1312(:,(kz-1)*ss.Ny+ky)=s1312(shearZhat(kz), L(ky), W(kz), ...
        y2_W-ss.shearY_chat(ky)', ss.fpTops+ss.dz/2);

        %if (ky == 51) & (kz == 51)
        %  disp( (kz-1)*ss.Ny+ky )
        %  disp(L(ky))
        %  disp(W(kz))
        %  disp(shearZhat(kz))
        %  tmp = shearY_c'-shearYhat(ky)';
        %  disp( tmp(42) )
        %  disp( shearZ_c(42))
        %  disp(ss.ss1212(42,(kz-1)*ss.Ny+ky))
        %end

    end
  end

  %disp(size(ss.fs1212))

  out.ff12 = ss.ff12;
  out.ss1212 = ss.ss1212;
  out.ss1213 = ss.ss1213;
  out.ss1312 = ss.ss1312;
  out.ss1313 = ss.ss1313;
  out.sf12 = ss.sf12;
  out.sf13 = ss.sf13;
  out.fs1212 = ss.fs1212;
  out.fs1312 = ss.fs1312;

  out.p.shearYc = shearY_c';
  out.p.shearZc = shearZ_c';
  out.p.fpTops = ss.fpTops;
  out.p.w = Wf;

  % figure for ff12 kernel
  %clf;
  %imagesc(log10(ss.ff12)); colorbar;
  %saveas(gcf, 'figures/BP1vD_ff12.png')

  disp('kernels created')


  % ---     Rheology      ---
  % Confining pressure (MPa) and Temperature (K)
  k  = 3.138; % thermal conductivity (W/m/K)
  Cp = 1171 ; % specific heat (J/kg/K)
  Rm = 3330 ; % mantle density (kg/m^3)

  Pconf       = Rm*9.8*ss.shearZ_chat/1e6;  % Shear zones
  Pconf_fault = Rm*9.8; % Faults

  Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
  Age_plate = 2e15; % seconds
  ss.Tprof  = 300+1380*erf(ss.shearZ_chat/(sqrt(4* Kappa * Age_plate)));  % Kelvin

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
  fprintf('Est. Recurrence time = %.2f (yr)\n\n', Ti/3.15e7);

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

  % artificial strain test
  % size(ss.e120) -> (Nx x Ny) x 1
  ast=false;
  if ast
    idx = 1000;
    ss.e120(idx) = ss.e120(idx)+1e-13;
  end

  ss.dgfF = 4;
  ss.dgfS = 4;

  % state vector init
  Y0=zeros(ss.M*ss.dgfF + ss.Ny*ss.Nz*ss.dgfS,1);

  % Fault patches
  Y0(1:ss.dgfF:ss.M*ss.dgfF)=zeros(ss.M,1);
  Y0(2:ss.dgfF:ss.M*ss.dgfF)=max(ss.a).*ss.sigma.*asinh(ss.Vpl./ss.Vo/2.* ...
    exp((ss.fo+ss.b.*log(ss.Vo./ss.Vpl))./max(ss.a))) + ss.eta.*ss.Vpl;
  Y0(3:ss.dgfF:ss.M*ss.dgfF)=ss.a./ss.b.*log(2*ss.Vo./ss.Vpl.* ...
    sinh((Y0(2:ss.dgfF:ss.M*ss.dgfF)-ss.eta.*ss.Vpl)./ss.a./ss.sigma))-ss.fo./ss.b;
  Y0(4:ss.dgfF:ss.M*ss.dgfF)=log(ss.Vpl./ss.Vo);

  % Shear zones
  Y0(ss.M*ss.dgfF+1:ss.dgfS:end)=ss.s120;
  Y0(ss.M*ss.dgfF+2:ss.dgfS:end)=ss.s130;
  Y0(ss.M*ss.dgfF+3:ss.dgfS:end)=ss.e120;
  Y0(ss.M*ss.dgfF+4:ss.dgfS:end)=ss.e130;

  % initialize the function handle with set constitutive parameters
  yp=@(t,y) odeBP1v_d(t,y,ss);

  % ODE45 Settings
  % Initial step of 1e-5 seconds
  % Relative tolerance of 3e-8
  % [0 3e10] = simulate 3e10 seconds, 3.15e7 seconds / year
  disp('begin solving')
  tic
  options=odeset('Refine',1,'RelTol',1e-8,'InitialStep',1e-5,'MaxStep',3e6);
  [t,Y]=ode45_2(yp,[0 1*3.15e7],Y0,options);
  disp('Done solving');
  toc

  Y=Y'; %necessary for ode45_2

  % ---       Figures        ---
  % instantaneous derivative
  Yp=zeros(length(t)-1,size(Y,2));
  for k=1:length(t)-1
    Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
  end

  % strain rate at center
  Ep=sqrt(Yp(:,ss.M*ss.dgfF+floor(length(ss.shearY_chat)/2)*ss.dgfS+3:ss.dgfS* ...
     length(ss.shearY_chat):end)'.^2 + Yp(:,ss.M*ss.dgfF+floor(length(ss.shearY_chat)/2)*...
     ss.dgfS+4:ss.dgfS*length(ss.shearY_chat):end)'.^2);

  % strain rate over whole ductile area
  Epall = sqrt( Yp(:,ss.M*ss.dgfF+3:ss.dgfS:end)'.^2 +...
               Yp(:,ss.M*ss.dgfF+4:ss.dgfS:end)'.^2);

  % stress
  %Epall = Y(:,ss.M*ss.dgfF+1:ss.dgfS:end)' - ss.s120;

  % ---       Figures        ---
  %disp(size(r.ss.Vo
  disp(size(Y))
  %disp(size(Y(:,4:r.ss.dgfF:r.ss.M*r.ss.dgfF)'))
  y.V = ss.Vo.*exp(Y(:,4:ss.dgfF:ss.M*ss.dgfF)'); % Slip rate (m/s)
  y.tau = Y(:,2:ss.dgfF:ss.M*ss.dgfF);            % Shear stress (MPa)
  y.Vmax = zeros(length(t),1);          % Maximum slip rate (m/s)
  y.Vcenter = y.V(floor(ss.M/2),:);          % Slip rate at center of VW region
  for ti = 1:length(t)
    y.Vmax(ti) = max(y.V(:,ti));
  end

  clf;
  plot(ss.a - ss.b);
  title('A minus B')
  saveas(gcf, 'figures/BP1vD_aminusb.png')

  clf;
  imagesc(log10(y.V)); colorbar;
  title('Slip Rate')
  xlabel('time steps')
  ylabel('fault mesh block')
  saveas(gcf, 'figures/BP1vD_slip.png')

  clf;
  imagesc(log10(Ep)); colorbar;
  title('Strain Rate of Center of Ductile Region')
  xlabel('Time Steps')
  ylabel('Block')
  saveas(gcf, 'figures/BP1vD_strainCenter.png')

  % ---         Movies        ---
  Smovie=false;
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

  Fmovie=false;
  if Fmovie
    % fault movie
    disp('begin fault movie')
    clf;
    fig = figure;
    fname='figures/BP1vD_slip.gif';
    for idx = 1:size(Y,1)
      oneV = y.V(:,idx);
      caxis([-15, -13]);
      imagesc(oneV); colorbar;
      title(idx)
      drawnow
      frame = getframe(fig);
      im{idx} = frame2im(frame);

      [A, map] = rgb2ind(im{idx}, 256);
      if idx==1
        imwrite(A,map,fname,'gif','LoopCount',Inf,'DelayTime',0.1);
      else
        imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',0.1);
      end
    end
  disp('fault movie done')
  end


end
