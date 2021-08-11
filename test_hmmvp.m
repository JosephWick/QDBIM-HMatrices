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

  len = 40000;
  n = 400;

  c.mu = 30e9;
  c.nu = 0.25;
  c.W = len;
  c.L = 1.0;
  c.dz = len/n;
  c.dip = 90;

  c.d1 = 1;
  c.d2 = 0;
  c.d3 = 0;

  c.eta = 3;

  x = zeros(1,n);
  y = linspace(0,len,n);
  z = zeros(1,n);
  c.X = [x; y; z];

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

function y = benchmark(r)
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
  d = zeros(m,n);
  s = int8((m/2)-(m/4));
  f = int8((m/2)+(m/4));
  x(s:f, s:f) = 1;
  x = reshape(d, [1,m*n]);

  y = hmmvp('mvp', hm, x);
  y = reshape(y, [m,n]);
  row = y(int8(m/2));

  clf;
  subplot(221); imagesc(fullM); title('h-matrix'); colorbar;
  subplot(222); imagesc(d); title('slip'); colorbar;
  subplot(223); imagesc(y); title('output'); colorbar;
  subplot(224); plot(row); title('row');
  saveas(gcf, 'figures/test_hmmvp.png')

end

% ------------------------- Private ----------------------------

function addpaths()
  addpath('../hmmvp-okada/matlab');
end
