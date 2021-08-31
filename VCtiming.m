function varargout = vct (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ------------------------- Public ----------------------------------

% time single MVP action using hmmvp and with dense matrices
function build()
  addpaths();

  probL = 200e3;
  probW = 200e3;

  N = 10;

  size = probL/N;

  % grid edges
  shearX = zeros(N*N);
  shearYhat = linspace(0, probL-size, N);
  shearZhat = linspace(0, probW-size, N);

  [shearZ shearY] = ndgrid(shearYhat, shearZhat);
  shearY = shearY';
  shearY = shearY(:)';
  shearZ = shearZ(:)';

  % grid centers
  shearX_c = shearX + size/2;
  shearY_c = shearY + size/2;
  shearZ_c = shearZ + size/2;

  c.X = [shearX_c; shearY_c; shearZ_c];
  c.Y = [shearX; shearY; shearZ];

  c.command = 'compress';
  c.dz = size;
  c.tol = 1e-5;
  c.G = 30e3;
  c.command = 'compress';
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VCT';
  c.write_hd_filename = './tmp/VCT-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);

  disp('run this in a shell: ')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)

end

function time(b)
  addpaths();

  % Slip velocities
  V = (2.*b.ss.Vs.*b.ss.a.*b.ss.sigmab./G).*...
       Lambert_W(G*b.ss.Vo./(2*b.ss.Vs.*b.ss.a.*b.ss.sigmab).* ...
       exp((tauF-b.ss.mu0.*b.ss.sigmab-b.ss.sigmab.*b.ss.b.*th)./ ...
       (b.ss.sigmab.*b.ss.a)));

  % slices for kernels
  % greater than ss.M
  gM = (b.ss.M+1:1: (b.ss.Nx*ss.Nz)+b.ss.M);
  % less than or equal to ss.M
  lM = (1:1:b.ss.M);

  fs1212 = hmmvp('init', b.fs1212, 32);
  m = hmmvp('getm', fs1212);)
  X = rand(m,1);

  % hmmvp mvp
  disp('hmmvp:')
  tic
  p = hmmvp('mvp', fs1212, X, lM, gM);
  toc

  % dense mvp
  fs1212_d = hmmvp('extract', b.fs1212, (1:1:m), (1:1:m));
  tic
  pd = fs1212_d(lm,gm)*X(lm);


end

% ----------------------- Private ----------------------------------
function addpaths()
  addpath('../hmmvp-okada/matlab')
end
