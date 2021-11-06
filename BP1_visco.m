% BP1_visco.m

function varargout = bp1v (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ---------------------------- Public ----------------------------

% set up h-matrices/kernels
% there will be 9 Kernels

function r = build()
  addpaths();

  % 9 kernels:
  % - 1 for fault-fault interaction (S12)
  % - 4 for shear-shear interaction (1212, 1213, 1312, 1313)
  % - 4 for fault-shear interaction (1212, 1213, 1312, 1313))
  % - 2 for shear-fault interaction (s12, s13)

  % problem specifications
  % - fault goes 40 km deep
  %  - shear TBD

  % ---       meshes      ---

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

  % combo mesh
  comboX = [faultX shearX];
  comboY = [faultY shearY];
  comboZ = [faultZ shearZ];

  comboX_c = [faultX_c shearX_c];
  comboY_c = [faultY_c shearY_c];
  comboZ_c = [faultZ_c shearZ_c];

  % KERNELS
  % ---       kvf params        ---
  c.command = 'compress';
  c.lambdaZ = ss.lambdaZ;
  c.tol = 1e-8;
  c.G = 30e3;
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  % ---       s12 kernel for fault-fault interaction        ---
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/BP1v_ff-s12';
  c.write_hd_filename = './tmp/BP1v_ff-s12-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.Y = [faultX; faultY; faultZ];
  c.X = [faultX_c; faultY_c; faultZ_c];

  kvf('Write', c.kvf, c, 4);
  disp('run these in a shell:')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.s12 = c.write_hmat_filename;

  % ---       shear1212 kernel for shear-shear interaction      ---
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/BP1v_ss-shear1212';
  c.write_hd_filename = './tmp/BP1v_ss-shear1212-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.transition = ss.transition;
  c.Y = [shearX; shearY; shearZ];
  c.X = [shearX_c; shearY_c; shearZ_c];

  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1212 = c.write_hmat_filename;

  % ---       shear 1213 kernel for shear-shear interaction     ---
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/BP1v_ss-shear1213';
  c.write_hd_filename = './tmp/BP1v_ss-shear1213-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is the same for all shear-shear kernels
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1213 = c.write_hmat_filename;

  % ---       shear 1312 kernel for shear-shear interaction       ---
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/BP1v_ss-shear1312';
  c.write_hd_filename = '/tmp/BP1v_ss-shear1312-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp/okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1312 = c.write_hmat_filename;

  % ---       shear 1313 kernel for shear-shear interaction       ---
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/BP1v_ss-shear1313';
  c.write_hd_filename = './tmp/BP1v_ss-shear1313-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1313 = c.write_hmat_filename;

  % ---       shear 1212 kernels for fault-shear interaction      ---
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/BP1v_fs-shear1212';
  c.write_hd_filename = './tmp/BP1v_fs-shear1212-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.X = [comboX_c, comboY_c, comboZ_c];
  c.Y = [comboX, comboY, comboZ];

  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1212 = c.write_hmat_filename;

  % ---       shear1312 kernels for fault-shear interaction     ---
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/BP1v_fs-shear1312';
  c.write_hd_filename = './tmp/BP1v_fs-shear1312-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1312 = c.write_hmat_filename;

  % ---         s12 kernel for shear-fault interaction      ---
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/BP1v_sf-s12';
  c.write_hd_filename = './tmp/BP1v_sf-s12-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  cmd = ['    ../hmmvo-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf12 = c.write_hmat_filename;

  % ---       s13 kernel for shear-fault interaction      ---
  c.greens_fn = 'okadaS13';
  c.write_hmat_filename = './tmp/BP1v_sf-s13';
  c.write_hd_filename = './tmp/BP1v_sf-s13-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf13 = c.write_hmat_filename;

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
  ss.e12p_plate = 1e-14*ones(length(shearY_chat)*length(ss.x3c),1);
  ss.e13p_plate =      zeros(length(shearY_chat)*length(ss.x3c),1);

  % Rheological Parameters
  % Reference Strain Rate (for stress in MPa)
  ss.Adif = 1e6*ones(length(ss.x3c)*length(shearY_chat),1);
  ss.Adis = 90 *ones(length(ss.x3c)*length(shearY_chat),1);

  % Power-Law Exponent
  ss.n = 3.5*ones(length(ss.x3c)*length(shearY_chat),1);

  % Activation Energy Wet Oliving (J/mol)
  ss.Qdif = 335e3*ones(length(ss.x3c)*length(shearY_chat),1);
  ss.Qdis = 480e3*ones(length(ss.x3c)*length(shearY_chat),1);

  % Activation Volume (m^3/mol)
  ss.Voldif = 4e-6*ones(length(ss.x3c)*length(shearY_chat),1);
  ss.Voldis = 11e-6*ones(length(ss.x3c)*length(shearY_chat),1);

  % Grain size (m)
  ss.d    = 1e-2*ones(length(ss.x3c)*length(shearY_chat),1);
  ss.pexp = 3*ones(length(ss.x3c)*length(shearY_chat),1);

  % Water fugacity (H/10^6 Si)
  ss.COH = 1000*ones(length(ss.x3c)*length(shearY_chat),1);
  ss.r   = 1.2*ones(length(ss.x3c)*length(shearY_chat),1);

  % Pressure (Pa)
  ss.P = repmat(1e6*Pconf',length(shearY_chat),1);
  ss.P = reshape(ss.P,[length(shearY_chat)*length(ss.x3c),1]);

  % Temperature (K)
  Te0 = repmat(ss.Tprof',length(shearY_chat),1);
  Te0 = reshape(Te0,[length(shearY_chat)*length(ss.x3c),1]);

  % Coefficients for dislocation and diffusion creep
  ss.Const_dis = ss.Adis.*exp(-(ss.Qdis+ss.P.*ss.Voldis)./(8.314.*Te0)).*ss.COH.^(ss.r);
  ss.Const_diff = ss.Adif.*exp(-(ss.Qdif+ss.P.*ss.Voldif)./(8.314.*Te0)).* ...
    ss.COH.^(ss.r).*ss.d.^(-ss.pexp);

  % Strengh profile
  ss.s120 = (ss.e12p_plate./ss.Const_dis).^(1./ss.n);
  ss.s130 = zeros(size(s120));
  ss.e120 = zeros(size(s120));
  ss.e130 = zeros(size(s120));

  r.ss = ss;

end

function y = solve(r)
  addpaths();

  % Use ode45 (Runge-Kutta 4th / 5th order accurate integration) to
  % solve the ODE time integration problem with adaptive time-steps
  % yp = f(t,y)
  % Y = [slip; stress; state variable; log10(slip rate / ref slip rate)]
  % Degrees of Freedom
  r.ss.dgfF=3;
  r.ss.dgfS=4;

  % state vector init
  Y0=zeros(r.ss.M*ss.dgfF+length(r.shearY_chat)*length(r.ss.shearZ_chat)*ss.dgfS,1);

  % Fault patches
  Y0(ss.M*ss.dgfF+1:ss.dgfF:2*ss.M*ss.dgfF)=zeros(size(ss.faultZ));
  Y0(ss.M*ss.dgfF+2:ss.dgfF:2*ss.M*ss.dgfF)=ss.strength;
  Y0(ss.M*ss.dgfF+3:ss.dgfF:2*ss.M*ss.dgfF)=log(ss.Vo./ss.b.Vpl);

  % Shear zones
  Y0(ss.M*ss.dgfF+1:ss.dgfS:end)=s120;
  Y0(ss.M*ss.dgfF+2:ss.dgfS:end)=s130;
  Y0(ss.M*ss.dgfF+3:ss.dgfS:end)=e120;
  Y0(ss.M*ss.dgfF+4:ss.dgfS:end)=e130;

  % initial conditions (starts at steady state w zero slip)
  Y0=zeros(r.ss.M*r.ss.dgf,1);
  Y0(1:r.ss.dgf:end)=zeros(r.ss.M,1);
  Y0(2:r.ss.dgf:end)=max(r.ss.a).*r.ss.sigma.*asinh(r.ss.Vpl./r.ss.Vo/2.* ...
    exp((r.ss.fo+r.ss.b.*log(r.ss.Vo./r.ss.Vpl))./max(r.ss.a))) + r.ss.eta.*r.ss.Vpl;
  Y0(3:r.ss.dgf:end)=r.ss.a./r.ss.b.*log(2*r.ss.Vo./r.ss.Vpl.*sinh((Y0(2:r.ss.dgf:end)-...
    r.ss.eta.*r.ss.Vpl)./r.ss.a./r.ss.sigma))-r.ss.fo./r.ss.b;
  Y0(4:r.ss.dgf:end)=log(r.ss.Vpl./r.ss.Vo);

  % load HM kernels
  hm.s12    = hmmvp('init', r.s12,     4);
  hm.ss1212 = hmmvp('init', r.ss1212, 32);
  hm.ss1213 = hmmvp('init', r.ss1213, 32);
  hm.ss1312 = hmmvp('init', r.ss1312, 32);
  hm.ss1313 = hmmvp('init', r.ss1313, 32);
  hm.sf12   = hmmvp('init', r.sf12,   32);
  hm.sf13   = hmmvp('init', r.sf13,   32);

  % initialize the function handle with set constitutive parameters
  yp=@(t,y) DieterichRuinaRegAging_BP1v(t,y,r.ss, hm);

  % ODE45 Settings
  % Initial step of 1e-5 seconds
  % Relative tolerance of 3e-8
  % [0 3e10] = simulate 3e10 seconds, 3.15e7 seconds / year
  disp('begin solving')
  tic
  options=odeset('Refine',1,'RelTol',1e-8,'InitialStep',1e-5);
  [t,Y]=ode45(yp,[0 500*3.15e7],Y0,options);
  disp('Done solving');
  toc

  % ---       Figures        ---
  y.V = r.ss.Vo.*exp(Y(:,4:r.ss.dgf:end)'); % Slip rate (m/s)
  y.tau = Y(:,2:r.ss.dgf:end);            % Shear stress (MPa)
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

end

% ---------------------------- Private ---------------------------

% add paths to hmmvp
function addpaths()
  addpath('../hmmvp-okada/matlab')
end

% boxcar function
function bc = BC(x)
  boxc=@(x) (x+0.5>=0)-(x-0.5>=0);
  bc = boxc(x);
end

% ramp function
function rmp = Ramp(x)
  rmp = x.*BC(x-1/2)+HS(x-1);
end

function h = HS(x)
  h = 0+x>=0;
end
