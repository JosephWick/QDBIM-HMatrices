% QDBIM2D_hm
% basically QDBIM2D but with hmmvp stress kernel

function varargout = exmb (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% --------------------- Public -----------------------------

function r = build()
  addpaths();

  % properties
  c.lambdaZ=40e3;
  c.n = 400;
  c.dz = c.lambdaZ/c.n;
  c.tol = 1.0e-3;

  % housekeeping
  c.command = 'compress';
  c.write_hmat_filename = getFname(c);
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.err_method  = 'mrem_fro';
  c.allow_overwrite = 1;

  c.greens_fn = 'okadaS12';

  x = zeros(1,n);
  y = linspace(0, -c.lambdaZ, c.n);
  z = zeros(1,n);

  disp('run this in a shell:')
  cmd = ['     ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf '\n'];
  disp(cmd)

  r.c = c;

end


% -------------------- Private -----------------------------

function addpaths()
  addpath('../hmmvp-okada/matlab');
end

function f=  getFname(p)
  f = sprintf('./tmp/QDBIM_tol%f_lz%d_n%d', p.tol, p.lambdaZ, p.n);
end
