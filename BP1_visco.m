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

  ss.lambdaZ = 40e3; % fault depth extent
  ss.M = 400; %number of fault cells
  ss.dz = ss.lambdaZ/ss.M;

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
  % TODO

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

  % TODO: other kernels

  % ---       fault properties      ---

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
  fprintf('Est. Recurrence time = %.2f (yr)\n', Ti/3.15e7);

  r.ss = ss;

end

function y = solve(r)
  addpaths();

  % Use ode45 (Runge-Kutta 4th / 5th order accurate integration) to
  % solve the ODE time integration problem with adaptive time-steps
  % yp = f(t,y)
  % Y = [slip; stress; state variable; log10(slip rate / ref slip rate)]
  % Degrees of Freedom
  r.ss.dgf=4;

  % initial conditions (starts at steady state w zero slip)
  Y0=zeros(r.ss.M*r.ss.dgf,1);
  Y0(1:r.ss.dgf:end)=zeros(r.ss.M,1);
  Y0(2:r.ss.dgf:end)=max(r.ss.a).*r.ss.sigma.*asinh(r.ss.Vpl./r.ss.Vo/2.* ...
    exp((r.ss.fo+r.ss.b.*log(r.ss.Vo./r.ss.Vpl))./max(r.ss.a))) + r.ss.eta.*r.ss.Vpl;
  Y0(3:r.ss.dgf:end)=r.ss.a./r.ss.b.*log(2*r.ss.Vo./r.ss.Vpl.*sinh((Y0(2:r.ss.dgf:end)-...
    r.ss.eta.*r.ss.Vpl)./r.ss.a./r.ss.sigma))-r.ss.fo./r.ss.b;
  Y0(4:r.ss.dgf:end)=log(r.ss.Vpl./r.ss.Vo);

  % load HM kernels
  hm.s12 = hmmvp('init', r.s12, 4);

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
  imagesc(y.V); colorbar;
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