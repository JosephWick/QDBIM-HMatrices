% ViscoelasticCycles.m

function varargout = vc (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% --------------------- Public ---------------------------

% set up h-matrices
% there will be 11 kernels
function r = build()
  addpaths();

  % fault mesh
  lambdaZ = 40e3; %(m)

  % visco mesh
  transition = 35e3;
  vL = 25e3;
  vW = 15e3;

  dz = 100;
  % fault has 400 elems
  % visco is 250x150

  % set up kvf general propts
  c.lambdaZ = lambdaZ;
  c.dz = dz;
  c.tol = 1e-3;
  c.G = 30e3;
  c.command = 'compress';
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  % s12 kernel for fault-fault interaction; 1D
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/VC_ff-s12';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.W = lambdaZ;
  c.L = 0;

  n = lambdaZ/dz;
  x = zeros(1,n);
  y = zeros(1,n);
  z = linspace(c.dz,lambdaZ,n);
  c.X = [x;y;z];

  kvf('Write', c.kvf, c, 4);
  disp('run this in a shell:')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)

  % shear 1212 kernel for shear-shear interaction
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_ff-shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];

  n = vL/dz;
  m = vW/dz;
  x = zeros(1,n);
  y = linspace(c.dz,vL,n);
  z = linspace(c.dz,vW,m);
  [Y Z] = ndgrid(y,z);
  c.X = [x;Y(:)';Z(:)'];

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
