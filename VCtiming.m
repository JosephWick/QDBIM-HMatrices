function varargout = vct (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ------------------------- Public ----------------------------------

% time single MVP action using hmmvp and with dense matrices
function b = build()
  addpaths();

  probL = 200e3;
  probW = 200e3;

  N = 100;
  size = probL/N;

  transition = 50e3;

  % grid edges
  shearX = zeros(1,N*N);
  shearYhat = linspace(0, probL-size, N);
  shearZhat = linspace(0, probW-size, N)+transition;

  [shearZ shearY] = ndgrid(shearYhat, shearZhat);
  shearY = shearY(:)';
  shearZ = shearZ(:)';

  % grid centers
  shearX_c = shearX;
  shearY_c = shearY + size/2;
  shearZ_c = shearZ + size/2;

  c.X = [shearX_c; shearY_c; shearZ_c];
  c.Y = [shearX; shearY; shearZ];

  c.transition = transition;

  c.dz = size;
  c.tol = 1e-5;
  c.G = 30e3;
  c.command = 'compress';
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  c.L = probL;
  c.W = probW;

  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VCT';
  c.write_hd_filename = './tmp/VCT-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);

  disp('run this in a shell: ')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)

  b.N = N;
  b.hm = c.write_hmat_filename;

end

function m = time(b)
  addpaths();

  %load hm
  fs1212 = hmmvp('init', b.hm);

  % make vector (boxcar slip)
  n = hmmvp('getn', fs1212);
  X = zeros(n,1);
  X(0.1*n:0.9*n) = 1.0;

  % hmmvp mvp
  disp('hmmvp:')
  tic
  p = hmmvp('mvp', fs1212, X);
  toc

  % dense mvp
  disp('dense:')
  fs1212_d = hmmvp('extract', fs1212, (1:1:n), (1:1:n));
  tic
  pd = fs1212_d*X;
  toc

  m = fs1212_d;

end

% ----------------------- Private ----------------------------------
function addpaths()
  addpath('../hmmvp-okada/matlab')
end
