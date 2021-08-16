% ViscoelasticCycles.m

function varargout = vc (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% --------------------- Public ---------------------------

% set up h-matrices
% there will be 11 kernels
function r = build()
  addpaths();

  % properties

  % fault mesh
  y3=0e3;
  y2 = -25e3;

  % Brittle-Ductile transition depth
  transDepth = 35e3;
  % number of elements
  ss.M = 120;
  dz = transDepth/r.c.N;

  fpoles = y3+(0:ss.M)'*dz;
  % tops of fault patches
  ss.y3f = fpoles(1:end-1);
  % width of fault patches
  Wf = ones(ss.M,1)*dz;

  % shear zone mesh
  ss.Nx = 50;
  ss.Nz = 50;
  eps = 1e-12;

  % hmmvp args
  c.command= 'compress';
  c.err_method = 'mrem-fro';
  c.allow_overwrite = 1;


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
