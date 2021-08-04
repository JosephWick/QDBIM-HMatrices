function varargout = testhmmvp (varargin)

  % trying to do the benchmark circular displacement test with hmmvp

  [varargout{1:nargout}] = feval(varargin{:});

end

% ------------------------ Public ---------------------------

% create h-matrix
function r = build()
  addpaths();

  c.command = 'compress';
  c.write_hmat_filename = './tmp/hmmvpTest';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'okada';
  %c.greens_fn = 'inverse-r';
  %c.greens_fn = 'test';
  c.err_method = 'mrem-fro';
  c.tol = 1.0e-3;
  c.order = 3;
  c.delta = 1.0e-4;
  c.allow_overwrite = 1;

  c.halfspace = 1;

  c.mu = 30e9;
  c.nu = 0.25;
  c.dz = 1.0;
  c.W = 100.0;
  c.L = 100.0;
  c.dip = 90;

  c.d1 = 1;
  c.d2 = 0;
  c.d3 = 0;

  c.eta = 3;

  n = 10;
  x = linspace(0,100,n);
  y = linspace(0,100,n);
  c.X = [x; y; zeros(1,n)];

  kvf('Write', c.kvf, c, 1);

  disp('run this in a shell:')
  fprintf('    ../hmmvp-okada/bin/hmmvpbuild_omp ./tmp/hmmvpTest.kvf \n')

  r.c = c;

end

function plot_hm(r)
  addpaths();

  hm_fname = r.c.write_hmat_filename;
  hm = hmmvp('init', hm_fname, 4);
  m = hmmvp('getm', hm);
  n = hmmvp('getn', hm);
  rs = (1:1:m); cs = (1:1:n);
  fullM = hmmvp('extract', hm, rs, cs);

  clf;
  imagesc(fullM); title('h-matrix'); colorbar;
  saveas(gcf, 'figures/test_hmmvp_fullHM.png')

end

function y = test_mvp(r)
  addpaths();

  hm_fname = r.c.write_hmat_filename;
  hm = hmmvp('init', hm_fname, 4);
  m = hmmvp('getm', hm);
  n = hmmvp('getn', hm);

  % let's get the whole hm to plot it
  rs = (1:1:m); cs = (1:1:n);
  fullM = hmmvp('extract', hm, rs, cs);

  % let's do the vector multiplication
  % make a vector to multiply by
  x = linspace(-50,50,10);
  y = linspace(-50,50,10);
  [X Y] = meshgrid(x,y);
  a = 10;
  R = (X.^2 + Y.^2)/a;

  d = real(sqrt(a - R));

  x = reshape(d, [100, 1]);

  y = hmmvp('mvp', hm, x);

  clf;
  subplot(221); imagesc(fullM); title('h-matrix'); colorbar;
  subplot(222); imagesc(d); title('slip'); colorbar;
  subplot(223); imagesc(reshape(y, [10,10])); title('output'); colorbar;
  subplot(224); plot(y(50:60)); title('row');
  saveas(gcf, 'figures/test_hmmvp.png')

end

% ------------------------- Private ----------------------------

function addpaths()
  addpath('../hmmvp-okada/matlab');
end
