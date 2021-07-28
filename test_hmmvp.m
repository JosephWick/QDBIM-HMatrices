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

  n = 100;
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
  m = hmmvp('getm', hm)
  n = hmmvp('getn', hm)

  y=0;

end

% ------------------------- Private ----------------------------

function addpaths()
  addpath('../hmmvp-main/matlab');
end
