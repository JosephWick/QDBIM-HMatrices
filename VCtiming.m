function varargout = vct (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ------------------------- Public ----------------------------------

% time single MVP action using hmmvp and with dense matrices
function b = build()
  addpaths();

  probL = 200e3;
  probW = 200e3;

  % 20, 40, 80, 160, 320

  transition = 40e3; %where shear zone starts
  shearYsize = 625;

  % grid edges
  ss.Ny = probL/shearYsize;
  ss.Nz = ss.Ny;
  nc = (-ss.Nz/2:ss.Nz/2);
  shearZhat = transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*transition;
  shearYhat = tan(nc*pi/(2.5*max(nc)))*32e3 ;
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

  c.transition = transition;

  c.dz = shearYsize;
  c.tol = 1e-6;
  c.G = 30e3;
  c.command = 'compress';
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  c.L = probL;
  c.W = probW;

  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VCT_8';
  c.write_hd_filename = './tmp/VCT-hd_8';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 4);

  disp('run this in a shell: ')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)

  b.N = ss.Ny;
  b.hm = c.write_hmat_filename;

end

function m = time(b)
  addpaths();

  %load hm
  fs1212 = hmmvp('init', b.hm, 1);

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
  %disp('dense:')
  %fs1212_d = hmmvp('extract', fs1212, (1:1:n), (1:1:n));
  %tic
  %pd = fs1212_d*X;
  %toc

  m = 1;%fs1212_d;

end

function greensFn(b)
  addpaths();

  % load hm
  hm = hmmvp('init', b, 16);

  % extract a row
  n = hmmvp('getn', hm);)
  row = hmmvp('extract', 200, 1, 1:n);

  plot(row);
  saveas(gcf, 'figures/greensFnRow.png')

end

% ----------------------- Private ----------------------------------
function addpaths()
  addpath('../hmmvp-okada/matlab')
end
