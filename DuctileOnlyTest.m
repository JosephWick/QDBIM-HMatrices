% DuctileOnlyTest.m
% files prefixed with 'DO_'

function varargout =  dot (varargin)

  [varargout{1:nargout}] = feval(vargin{:});

end

% ------------------ Public -------------------------

% set up h-matrices
% there will be 2 kernels, ss1212, ss1213, ss1312, ss1313
function r = build()
  addpaths();

  % problem specifications
  % ductile zone starts 40 km deep

  % x1 (x) is in/out of page
  % x2 (y) is left/rifht of page
  % x3 (z) is up/down

  % ---     General params    ---
  transition = 40e3;
  probL = 200e3;
  probW = 200e3;

  ss.Ny = 51;
  ss.Nz = 51;

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
    shearZ_chat(idx) = shearZhat(idx) + (shearZhat(idx+1) - shearZhat(idx))/2;
    shearY_chat(idx) = shearYhat(idx) + (shearYhat(idx+1) - shearYhat(idx))/2;
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

  y2 = probL/2;
  ss.Nx = ss.Nz;

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

  % ---     kvf params      ---
  c.command = 'compress';
  c.lambdaZ = lamdbaZ; ss.lamdbaZ = lamdbaZ;
  c.tol = 1e-8
  c.G = 30e3;
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  % ---       shear1212 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_ss-shear1212_200';
  c.write_hd_filename = './tmp/VC_ss-shear1212-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.transition = transition;
  c.L = probL;
  c.W = probW;
  c.Y = [shearX; shearY; shearZ];
  c.X = [shearX_c; shearY_c; shearZ_c];

  clf;
  imagesc(c.X); colorbar;
  saveas(gcf, 'figures/VCX.png')

  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp('Run these in a shell in this directory')
  disp(cmd)
  r.ss1212 = c.write_hmat_filename;
  r.c=c;
  % ---       Shear 1213 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/VC_ss-shear1213_200';
  c.write_hd_filename = './tmp/VC_ss-shear1213-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is the same for all shear-shear interaction
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1213 = c.write_hmat_filename;

  % ---     Shear 1312 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/VC_ss-shear1312_200';
  c.write_hd_filename = './tmp/VC_ss-shear1312-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1312 = c.write_hmat_filename;

  % ---     shear 1313 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/VC_ss-shear1313_200';
  c.write_hd_filename = './tmp/VC_ss-shear1313-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1313 = c.write_hmat_filename;

end

% run numerical sim
function out = solve(r)
  addpaths();

  G = 30e3;

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
  % a is 'direct effect' - strength parameter
  % b is 'evolution effect' - weakening parameter
  r.ss.a = 1e-3*ones(size(r.ss.y3f));
  r.ss.b = r.ss.a+2.1e-4*ones(size(r.ss.y3f)); % make + for with vw region

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
  %r.ss.b(1:top)      = r.ss.a(1:top)-2.1e-4;
  %r.ss.b(bottom:end) = r.ss.a(bottom:end)-2.1e-4;

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


end

  % ------------------ Private -------------------------
  function addpaths()
    addpath('../hmmvp-okada/matlab')
  end

  % boxcar function
  function bc = BC(x)
    boxc=@(x) (x+0.5>=0)-(x-0.5>=0);
    bc = boxc(x);
  end
