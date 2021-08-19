% ViscoelasticCycles.m

function varargout = vc (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% --------------------- Public ---------------------------

% set up h-matrices
% there will be 11 kernels
function r = build()
  addpaths();

  % we'll have 9 kernels:
  % - 1 for fault-fault interaction (s12)
  % - 4 for shear-shear interaction (1212, 1213, 1312, 1313)
  % - 2 for fault-shear interaction (1212, 1213)
  % - 2 for shear-fault interaction (1212, 1213)

  % fault mesh
  lambdaZ = 40e3; %(m)

  % visco mesh
  transition = 35e3;
  vL = 25e3;
  vW = 15e3;
  ss.transition = transition;

  dz = 100;
  ss.dz = 100;
  % fault has 400 elems
  % visco is 250x150

  % set up kvf general propts
  c.lambdaZ = lambdaZ; ss.lambdaZ = lambdaZ
  c.dz = dz;
  c.tol = 1e-3;
  c.G = 30e3;
  c.command = 'compress';
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  % - s12 kernel for fault-fault interaction; 1D -
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/VC_ff-s12';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.W = lambdaZ;
  c.L = 0;

  n = lambdaZ/dz;
  x = zeros(1,n);
  y = zeros(1,n) + 15000;
  z = linspace(c.dz,lambdaZ,n);
  c.X = [x;y;z];

  y2 = 15000;
  ss.y2 = y2;
  % tops of fault patches
  ss.y3f = z - dz;
  % width of fault patches
  ss.Wf = ones(n,1)*dz;
  ss.M = n;

  kvf('Write', c.kvf, c, 4);
  disp('run these in a shell:')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.s12 = c.write_hmat_filename;

  % - shear 1212 kernel for shear-shear interaction -
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_ss-shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.transition=transition;
  c.L = vL;
  c.W = vW;
  n = vL/dz;
  m = vW/dz;
  x = zeros(1,n*m);
  y = linspace(c.dz,vL,n);
  z = linspace(c.dz,vW,m);
  [Y Z] = ndgrid(y,z);
  c.X = [x;Y(:)';Z(:)'];

  % shear mesh
  ss.Nx = n;
  ss.Nz = m;
  ss.polesz = transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*transition;
  % center of shear zone (x3)
  ss.x3c=(ss.polesz(2:end)+ss.polesz(1:end-1))/2;
  W=ss.polesz(2:end)-ss.polesz(1:end-1);

  ss.polesxc=(2*y2/1e3:(2*y2)/(26e3):2*y2/1e3)'*1e3;
  edges= floor(ss.Nx-length(ss.polesxc)+1)/2;
  ss.polesxl=flipud(min(ss.polesxc)-tan((0:edges)'*pi/(2.2*(edges)+eps))*transition);
  ss.polesxr=max(ss.polesxc)+tan((0:edges)'*pi/(2.2*(edges)+eps))*transition;
  ss.polesx=[ss.polesxl(1:end-1);ss.polesxc;ss.polesxr(2:end)];

  % center of shear zone (x2)
  ss.x2c=(ss.polesx(2:end)+ss.polesx(1:end-1))/2;

  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1212 = c.write_hmat_filename;

  % - shear 1213 kernel for shear-shear interaction -
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/VC_ss-shear1213';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is the same for all shear-shear interactions
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1213 = c.write_hmat_filename;

  % - shear 1312 kernel for shear-shear interaction -
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/VC_ss-shear1312';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1312 = c.write_hmat_filename;

  % - shear 1313 kernel for shear-shear interaction -
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/VC_ss-shear1313';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1313 = c.write_hmat_filename;

  % - shear 1212 kernel for fault-shear interaction -
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_fs-shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % update geometry
  fn = lambdaZ/dz;                %fault
  fx = zeros(1,fn);
  fy = zeros(1,fn) + 15000;
  fz = linspace(c.dz,lambdaZ,fn);
  sn = vL/dz;                     %shear
  sm = vW/dz;
  sx = zeros(1,sn*sm);
  sy = linspace(c.dz,vL,sn);
  sz = linspace(c.dz,vW,sm);
  [Ys Zs] = ndgrid(sy,sz);
  bx = [fx sx];                   %together
  by = [fy Ys(:)'];
  bz = [fz Zs(:)'];
  c.X = [bx;by;bz];
  % write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1212 = c.write_hmat_filename;

  % - shear 1213 kernel for fault-shear interaction -
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/VC_fs-shear1213';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is same as fs1212
  % %write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1213 = c.write_hmat_filename;

  % - shear 1212 kernel for shear-fault interaction -
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_sf-shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % update geometry - need to adjust fault locs to be receiver
  fn = lambdaZ/dz;                %fault
  fx = zeros(1,fn);
  fy = zeros(1,fn) + 15000 + 0.5*dz;
  fz = linspace(c.dz,lambdaZ,fn);
  bx = [fx sx];                   %together
  by = [fy Ys(:)'];
  bz = [fz Zs(:)'];
  c.X = [bx;by;bz];
  % write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf1212 = c.write_hmat_filename;

  % - shear 1313 kernel for shear-fault interaction -
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/Vc_sf-shear1313';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is same as sf1212
  % write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf1313 = c.write_hmat_filename;

  r.ss = ss;

end

% runs numerical solution
function solve(r)
  addpaths();

  % - frictional parameters -

  % Confining pressure (MPa) and Temperature (K)
  k = 3.138; % thermal conductivity (W/m/K)
  Cp = 1171 ; % specific heat (J/kg/K)
  Rm = 3330 ; % mantle density (kg/m^3)

  Pconf       = Rm*9.8*r.ss.x3c/1e6;  % Shear zones
  Pconf_fault = Rm*9.8*(r.ss.y3f+r.ss.dz); % Faults

  Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
  Age_plate = 2e15; % seconds
  ss.Tprof  = 300+1380*erf(r.ss.x3c/(sqrt(4* Kappa * Age_plate)));  % Kelvin

  % default friction properties (velocity-weakening friction)
  % effective confining pressure (MPa)
  ss.sigmab = 1000;

  % frictional parameters
  r.ss.a = 1e-3*ones(size(r.ss.y3f));
  r.ss.b = r.ss.a+2.1e-4*ones(size(r.ss.y3f));

  r.ss.mu0 = 0.2*ones(size(r.ss.y3f));
  % characteristic weakening distance (m)
  r.ss.L = 0.012*ones(size(r.ss.y3f));
  % plate velocity (m/s)
  r.ss.V_plate = 1e-9*ones(size(r.ss.y3f));
  % reference slip rate (m/s)
  r.ss.Vo = 1e-6*ones(size(r.ss.y3f));
  % shear wave speed (m/s)
  r.ss.Vs = 3e3*ones(size(r.ss.y3f));

  % Velocity-strengthening at edges
  top    = floor(5e3/(r.ss.lambdaZ/r.ss.M));
  bottom = ceil(30e3/(r.ss.lambdaZ/r.ss.M));
  r.ss.b(1:top)      = r.ss.a(1:top)-2.1e-4*ones(1,top);
  r.ss.b(bottom:end) = r.ss.a(bottom:end)-2.1e-4*ones(length(r.ss.a(bottom:end)),1);

  % - Rheology -
  r.ss.e12p_plate = 0.0 %1e-14*ones(length(ss.x2c)*length(ss.x3c),1);
  r.ss.e13p_plate = 0.0     %zeros(length(ss.x2c)*length(ss.x3c),1);
  % Rheological Parameters
  % Reference Strain Rate (for stress in MPa)
  r.ss.Adif = 1e6*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.Adis = 90 *ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Power-Law Exponent
  r.ss.n = 3.5*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Activation Energy Wet Oliving (J/mol)
  r.ss.Qdif = 335e3*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.Qdis = 480e3*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Activation Volume (m^3/mol)
  r.ss.Voldif = 4e-6*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.Voldis = 11e-6*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Grain size (m)
  r.ss.d    = 1e-2*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.pexp = 3*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Water fugacity (H/10^6 Si)
  r.ss.COH = 1000*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.R   = 1.2*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Pressure (Pa)
  r.ss.P = repmat(1e6*Pconf',length(r.ss.x2c),1);
  r.ss.P = reshape(r.ss.P,[length(r.ss.x2c)*length(r.ss.x3c),1]);
  % Temperature (K)
  Te0 = repmat(r.ss.Tprof',length(r.ss.x2c),1);
  Te0 = reshape(Te0,[length(r.ss.x2c)*length(r.ss.x3c),1]);
  % Coefficients for dislocation and diffusion creep
  r.ss.Const_dis = r.ss.Adis.*exp(-(r.ss.Qdis+r.ss.P.*r.ss.Voldis)./(8.314.*Te0)).*r.ss.COH.^(r.ss.R);
  r.ss.Const_diff =r.ss.Adif.*exp(-(r.ss.Qdif+r.ss.P.*r.ss.Voldif)./(8.314.*Te0)).*r.ss.COH.^(r.ss.R).*r.ss.d.^(-r.ss.pexp);
  % Strengh profile
  s120 = (r.ss.e12p_plate./r.ss.Const_dis).^(1./r.ss.n);
  s130 = zeros(size(s120));
  e120 = zeros(size(s120));
  e130 = zeros(size(s120));
  % Fault Strength
  r.ss.strength = r.ss.sigmab.*(r.ss.mu0+(r.ss.a-r.ss.b).*log(r.ss.V_plate./r.ss.Vo))+G*r.ss.V_plate./(2*r.ss.Vs);

  % - Numerical Solution -

  % state parameters
  r.ss.dgfF=3;
  r.ss.dgfS=4;
  %% Initialize State Vector
  Y0=zeros(r.ss.M*r.ss.dgfF+length(r.ss.x2c)*length(r.ss.x3c)*r.ss.dgfS,1);

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
  yp=@(t,y) odeViscoelastic(t,y,ss);
  tic
  % Solve the system
  options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6);
  [t,Y]=ode45(yp,[0 1e10],Y0,options);
  toc

end

% -------------------- Private ---------------------------

function addpaths()
  addpath('../hmmvp-okada/matlab')
end

function f = getFname(p)
  f = sprintf('./tmp/VC_%s_tol%f_lz%d_n%d', p.greens_fn, p.tol, p.lambdaZ, p.n);
end

% boxcar function
function bc = BC(x)
  boxc=@(x) (x+0.5>=0)-(x-0.5>=0);
  bc = boxc(x);
end
