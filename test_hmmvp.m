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
  c.greens_fn = 'inverse-r';
  c.err_method = 'mrem-fro';
  c.tol = 1.0e-3;
  c.order = 3;
  c.delta = 1.0e-4;
  c.allow_overwrite = 1;

  n = 10;
  x = linspace(-1,1,n);
  [X Y] = ndgrid(x,x);
  c.X = [X(:)'; Y(:)'; zeros(1,n^2)];

  kvf('Write', c.kvf, c, 1);

  disp('run this in a shell:')
  fprintf('    ../hmmvp-main/bin/hmmvpbuild_omp ./tmp/hmmvpTest.kvf \n')

  r.c = c;

end

function y = test_mvp(r)
  addpaths();

  hm_fname = r.c.write_hmat_filename;
  hm = hmmvp('init', hm_fname, 4);
  m = hmmvp('getm', hm);
  n = hmmvp('getn', hm);

  % let's get the whole hm to plot it
  rs = (1:1:m); cs = (1:1:n);
  m = hmmvp('extract', hm, rs, cs);

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
  subplot(221); imagesc(m); title('h-matrix'); colorbar;
  subplot(222); imagesc(d); title('slip'); colorbar;
  subplot(223); imagesc(reshape(y, [10,10])); title('output'); colorbar;
  subplot(224); plot(y(500:510)); title('row');
  saveas(gcf, 'figures/test_hmmvp.png')

end

% ------------------------- Private ----------------------------

function addpaths()
  addpath('../hmmvp-main/matlab');
end
