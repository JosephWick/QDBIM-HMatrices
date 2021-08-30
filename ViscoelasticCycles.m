% ViscoelasticCycles.m

function varargout = vc (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% --------------------- Public ---------------------------

% set up h-matrices
% there will be 9 kernels
function r = build()
  addpaths();

  % we'll have 9 kernels:
  % - 1 for fault-fault interaction (s12)
  % - 4 for shear-shear interaction (1212, 1213, 1312, 1313)
  % - 4 for fault-shear interaction (1212, 1213, 1312, 1313))
  % - 2 for shear-fault interaction (s12, s13)

  % problem specifications
  % 200km long (L), 200km deep (W)
  % transition depth is 40 km
  % fault goes 40km deep

  % x1 (x) is in/out of page
  % x2 (y) is left/right of page
  % x3 (z) is up/down

  % ---     general params      ---
  probL = 200e3;
  probW = 200e3;
  lambdaZ = 40e3; % fault depth
  transition = 40e3; %where shear zone starts

  % fault mesh
  ss.M = 400;
  ss.dz = lambdaZ/ss.M;
  % fault patch edges (top left)
  faultX = zeros(1,ss.M);
  faultY = ones(1,ss.M)*probL/2;
  faultZ = linspace(0,lambdaZ-ss.dz, ss.M);
  %tops of fault patches
  ss.y3f = faultZ;
  % fault patch centers
  faultX_c = faultX;
  faultY_c = faultY;
  faultZ_c = faultZ+ss.dz/2;

  % shear mesh
  eps = 1e-12;
  shearYsize = 100;
  ss.Ny = probL/shearYsize;
  ss.Nz = ss.Ny;
  shearZhat = transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*transition;
  % shear patch edges
  shearX = zeros(1,ss.Nz*ss.Ny);
  shearYhat = linspace(0,probL-shearYsize, ss.Ny);
  [shearZ shearY] = ndgrid(shearYhat, shearZhat);
  shearY = shearY(:)';
  shearZ = shearZ(:)';
  % shear patch centers
  shearX_c = shearX;
  shearY_chat = shearYhat+shearYsize/2;
  shearZ_chat = zeros(1,length(shearZhat));
  for idx=(1:length(shearZhat)-1)
    shearZ_c(idx) = (shearZhat(idx+1) - shearZhat(idx))/2;
  end
  shearZ_c(length(shearZ_c)) = shearZ_c(length(shearZ_c)-1);
  [shearZ_c shearY_c] = ndgrid(shearY_chat, shearZ_chat);

  y2 = probL/2;
  ss.Nx = length(shearX);
  ss.Nz = length(shearZ);
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

  % combo mesh
  disp(size(faultY_c))
  disp(size(shearY_c))
  comboX = [faultX shearX];
  comboY = [faultY shearY];
  comboZ = [faultZ shearZ];

  comboX_c = [faultX_c shearX_c];
  comboY_c = [faultY_c shearY_c];
  comboZ_c = [faultZ_c shearZ_c];

  % ---       kvf params      ---
  c.command = 'compress';
  c.lambdaZ = lambdaZ; ss.lambdaZ = lambdaZ;
  c.dz = dz;
  c.tol = 1e-3;
  c.G = 30e3;
  c.command = 'compress';
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  % ---       s12 kernel for fault-fault interaction      ---
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/VC_ff-s12';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.W = lambdaZ;
  c.L = 0;
  ss.Wf = 0;
  c.Y = [faultX; faultY; faultZ];
  c.X = [faultX_c; faultY_c; faultZ_c];

  kvf('Write', c.kvf, c, 4);
  disp('run these in a shell:')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.s12 = c.write_hmat_filename;

  % ---       shear1212 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_ss-shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.transition = transition;
  c.L = probL;
  c.W = probW;
  c.Y = [shearX; shearY; shearZ];
  c.X = [shearX_c; shearY_c; shearZ_c];

  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1212 = c.write_hmat_filename;

  % ---       Shear 1213 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/VC_ss-shear1213';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is the same for all shear-shear interaction
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1213 = c.write_hmat_filename;

  % ---     Shear 1312 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/VC_ss-shear1312';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf]
  disp(cmd)
  r.ss1312 = c.write_hmat_filename;

  % ---     shear 1313 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/VC_ss-shear1313';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1313 = c.write_hmat_filename;

  % ---     shear 1212 kernels for fault-shear interaction ---
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_fs-shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.X = [comboX_c; comboY_c; comboZ_c];
  c.Y = [comboX; comboY; comboZ];

  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1212 = c.write_hmat_filename;

  % ---       shear 1312 kernel for fault-shear interaction ---
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/VC_fs-shear1312';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1312 = c.write_hmat_filename;

  % ---       s12 kernel for shear-fault interaction ---
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/VC_sf-12';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf12 = c.write_hmat_filename;

  % ---     s13 kernel for shear-fault interaction ---
  c.greens_fn = 'okadaS13';
  c.write_hmat_filename = './tmp/VC_sf-s13';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd);
  r.sf13 = c.write_hmat_filename;

  r.ss = ss;

end

% set up h-matrices
% there will be 11 kernels
function r = buildOld()
  addpaths();

  % we'll have 9 kernels:
  % - 1 for fault-fault interaction (s12)
  % - 4 for shear-shear interaction (1212, 1213, 1312, 1313)
  % - 4 for fault-shear interaction (1212, 1213, 1312, 1313))
  % - 2 for shear-fault interaction (s12, s13)

  % fault mesh
  lambdaZ = 40e3; %(m)

  % visco mesh
  transition = 40e3;
  vL = 200e3;
  vW = 200e3;
  ss.transition = transition;

  dz = 100;
  ss.dz = 100;
  % fault has 400 elems
  % visco is 250x150

  % set up kvf general propts
  c.lambdaZ = lambdaZ; ss.lambdaZ = lambdaZ;
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
  y = zeros(1,n) + 100e3;
  z = linspace(c.dz,lambdaZ,n);
  c.X = [x;y;z];

  y2 = 100e3;
  ss.y2 = y2;
  % tops of fault patches
  ss.y3f = z - dz;
  ss.y3f = ss.y3f';
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
  n = vL/dz/2;
  m = vW/dz/2;
  z = transition+tan((0:m)'*pi/(2.2*(m+eps)))*transition;
  x = zeros(1,n*m);
  y = linspace(c.dz,vL,n);
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

  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1212 = c.write_hmat_filename;

  % - shear 1213 kernel for shear-shear interaction -
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/VC_ss-shear1213';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is the same for all shear-shear interactions
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1213 = c.write_hmat_filename;

  % - shear 1312 kernel for shear-shear interaction -
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/VC_ss-shear1312';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1312 = c.write_hmat_filename;

  % - shear 1313 kernel for shear-shear interaction -
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/VC_ss-shear1313';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
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
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1212 = c.write_hmat_filename;

  % - shear 1312 kernel for fault-shear interaction -
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/VC_fs1312';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geoemtry is same as fs1212
  % write
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1312 = c.write_hmat_filename;

  % - s12 kernel for shear-fault interaction -
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/VC_sf-s12';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % update geometry - need to adjust fault locs to be receiver
  fy = zeros(1,fn) + 15000 + 0.5*dz;
  by = [fy Ys(:)'];
  c.X = [bx;by;bz];
  % write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf12 = c.write_hmat_filename;

  % - s13 kernel for shear-fault interaction -
  c.greens_fn = 'okadaS13';
  c.write_hmat_filename = './tmp/VC_sf-s13';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is same as sf12
  % write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf13 = c.write_hmat_filename;

  r.ss = ss;

end

% runs numerical solution
function out = solve(r)
  addpaths();

  % rigidity (MPa)
  G=30e3;

  % - frictional parameters -

  % Confining pressure (MPa) and Temperature (K)
  k = 3.138; % thermal conductivity (W/m/K)
  Cp = 1171 ; % specific heat (J/kg/K)
  Rm = 3330 ; % mantle density (kg/m^3)

  Pconf       = Rm*9.8*r.ss.x3c/1e6;  % Shear zones
  Pconf_fault = Rm*9.8*(r.ss.y3f+r.ss.dz); % Faults

  Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
  Age_plate = 2e15; % seconds
  r.ss.Tprof  = 300+1380*erf(r.ss.x3c/(sqrt(4* Kappa * Age_plate)));  % Kelvin

  % default friction properties (velocity-weakening friction)
  % effective confining pressure (MPa)
  r.ss.sigmab = 1000;

  % frictional parameters
  r.ss.a = 1e-3*ones(size(r.ss.y3f));
  r.ss.b = 1e-4*ones(size(r.ss.y3f)); %r.ss.a+2.1e-4*ones(size(r.ss.y3f));

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
  bottom = ceil(15e3/(r.ss.lambdaZ/r.ss.M));
  r.ss.b(1:top)      = r.ss.a(1:top)-2.1e-4*ones(top,1);
  r.ss.b(bottom:end) = r.ss.a(bottom:end)-2.1e-4*ones(length(r.ss.a(bottom:end)),1);

  % - Rheology -
  r.ss.e12p_plate = 1e-14; %1e-14*ones(length(ss.x2c)*length(ss.x3c),1);
  r.ss.e13p_plate = 0.0; %zeros(length(ss.x2c)*length(ss.x3c),1);
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
  r.ss.Const_diff=r.ss.Adif.*exp(-(r.ss.Qdif+r.ss.P.*r.ss.Voldif)./(8.314.*Te0)).*r.ss.COH.^(r.ss.R).*r.ss.d.^(-r.ss.pexp);
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
  Y0(1:r.ss.dgfF:r.ss.M*r.ss.dgfF)=zeros(size(r.ss.y3f));
  Y0(2:r.ss.dgfF:r.ss.M*r.ss.dgfF)=r.ss.strength;
  Y0(3:r.ss.dgfF:r.ss.M*r.ss.dgfF)=log(r.ss.Vo./r.ss.V_plate);

  % Shear zones
  Y0(r.ss.M*r.ss.dgfF+1:r.ss.dgfS:end)=s120;
  Y0(r.ss.M*r.ss.dgfF+2:r.ss.dgfS:end)=s130;
  Y0(r.ss.M*r.ss.dgfF+3:r.ss.dgfS:end)=e120;
  Y0(r.ss.M*r.ss.dgfF+4:r.ss.dgfS:end)=e130;

  % get the h-matrices
  hm.s12    = hmmvp('init', r.s12,     4);
  hm.ss1212 = hmmvp('init', r.ss1212, 32);
  hm.ss1213 = hmmvp('init', r.ss1213, 32);
  hm.ss1312 = hmmvp('init', r.ss1312, 32);
  hm.ss1313 = hmmvp('init', r.ss1313, 32);
  hm.fs1212 = hmmvp('init', r.fs1212, 32);
  hm.fs1312 = hmmvp('init', r.fs1312, 32);
  hm.sf12   = hmmvp('init', r.sf12,   32);
  hm.sf13   = hmmvp('init', r.sf13,   32);

  % initialize the function handle with
  % set constitutive parameters
  yp=@(t,y) odeViscoelastic(t,y,r.ss, hm);
  tic
  % Solve the system
  options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6);
  [t,Y]=ode45(yp,[0 1e10],Y0,options); %1e10
  toc
  % Compute the instantaneous derivative
  Yp=zeros(length(t)-1,size(Y,2));
  for k=1:length(t)-1
    Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
  end

  % - figures -
  % Strain rate at center
  Ep=sqrt(Yp(:,2*r.ss.M*r.ss.dgfF+floor(length(r.ss.x2c)/2)*r.ss.dgfS+3:r.ss.dgfS*length(r.ss.x2c):end)'.^2 +...
        Yp(:,2*r.ss.M*r.ss.dgfF+floor(length(r.ss.x2c)/2)*r.ss.dgfS+4:r.ss.dgfS*length(r.ss.x2c):end)'.^2);

  % Velocity
  V=Yp(:,1:r.ss.dgfF:r.ss.M*r.ss.dgfF);

  % Maximum Velocity
  Vmax=zeros(length(t)-1,1);
  for ts=1:length(t)-1
      Vmax(ts)=max(V(ts,:));
  end

  %%
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %                    Function of Time                   %
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  figure(1);clf;set(gcf,'name','Time Evolution')
  f1=subplot(5,1,1);cla;
  pcolor(t(1:end-1)/3.15e7,r.ss.y3f/1e3,log10(V')), shading flat
  set(gca,'YDir','reverse');

  h=colorbar('Location','NorthOutside');
  caxis([min(min(log10(V))) max(max(log10(V)))]);
  colormap(f1,parula);
  title(h,'Slip Rate West (m/s)')
  xlabel('Time (yr)')
  ylabel('Depth (km)');

  f2=subplot(5,1,2);cla;
  pcolor(t(1:end-1)/3.15e7,r.ss.y3f/1e3,log10(V')), shading flat
  set(gca,'YDir','reverse');

  h=colorbar('Location','NorthOutside');
  caxis([min(min(log10(V))) max(max(log10(V)))]);
  colormap(f2,parula);
  title(h,'Slip Rate East (m/s)')
  xlabel('Time (yr)')
  ylabel('Depth (km)');

  f3=subplot(5,1,3);cla;
  %pcolor(t(1:end-1)/3.15e7,r.ss.x3c/1e3,log10(Ep)), shading flat
  set(gca,'YDir','reverse');

  %caxis([log10(min(min(Ep))) log10(max(max(Ep)))]);
  h1=colorbar('Location','NorthOutside');
  colormap(f3,hot);
  title(h1,'Strain Rate (1/s)')
  xlabel('Time (Yr)')
  ylabel('Depth (km)');

  subplot(5,1,4);cla;
  plot(t(1:end-1)/3.15e7,log10(Vmax))
  xlabel('Time (Yr)')
  ylabel('Velocity (m/s) log10')
  title('Maximum slip rates on faults')

  subplot(5,1,5);cla;
  plot(t(1:end-1)/3.15e7,log10(V(:,floor((top+bottom)/2))))
  xlabel('Time (Yr)')
  ylabel('Velocity (m/s) log10')
  title('Time series at center of seismogenic zones')

  % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %                Function of Time Steps                 %
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  f1=figure(2);clf;set(gcf,'name','Time Step Evolution')

  subplot(5,1,1);cla;
  pcolor(1:length(t)-1,r.ss.y3f/1e3,log10(V')), shading flat
  set(gca,'YDir','reverse');

  h=colorbar('Location','NorthOutside');
  caxis([min(min(log10(V))) max(max(log10(V)))]);
  colormap(f1,parula);
  title(h,'Slip Rate West (m/s)')
  xlabel('Time Steps')
  ylabel('Depth (km)');

  out.t = t;
  out.V = V;
  out.E = Ep;

  f3=subplot(5,1,3);cla;
  pcolor(1:length(t)-1, r.ss.x3c(1:end-1)/1e3, log10(Ep)), shading flat
  set(gca,'YDir','reverse');

  %caxis([log10(min(min(Ep))) log10(max(max(Ep)))]);
  h1=colorbar('Location','NorthOutside');
  colormap(f3,hot);
  title(h1,'Strain Rate (1/s)')
  xlabel('Time Steps')
  ylabel('Depth (km)');

  subplot(5,1,4);cla;
  plot(1:length(t)-1,log10(Vmax))
  xlabel('Time Steps')
  ylabel('Velocity (m/s) log10')
  title('Maximum slip rates on faults')

  subplot(5,1,5);cla;
  plot(1:length(t)-1,log10(V(:,floor((top+bottom)/2))))
  xlabel('Time Steps')
  ylabel('Velocity (m/s) log10')
  title('Time series at center of seismogenic zones')

  saveas(f1, 'figures/VC_f1.png')
  saveas(f2, 'figures/VC_f2.png')
  saveas(f3, 'figures/VC_f3.png')

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
