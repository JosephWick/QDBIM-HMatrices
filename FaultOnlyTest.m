% FaultOnlyTest.m
% files prefixed with 'FO_'

function varargout =  fot (varargin)

  [varargout{1:nargout}] = feval(vargin{:});

end

% ------------------ Public -------------------------

% set up h-matrices
% there will be one kernel (fault-fault)
function r = build()
  addpaths();

  % problem specifications
  % fault goes 40km deep

  % x1 (x) is in/out of page
  % x2 (y) is left/right of page
  % x3 (z) is up/down

  % --      General Params    --
  % depth  extend of fault
  lambdaZ = 40e3;

  % number of elements in mesh
  ss.M = 400;

  ss.dz = lamdbaZ/ss.M;
  % fault patch edges(top left)
  faultX = zeros(1,ss.M);
  faultY = zeros(1,ss.M);
  faultZ = linspace(0, lambdaZ-ss.dz, ss.M);
  % tops of fault patches
  ss.y3f = faultZ;
  % fault patch centers
  faultX_c = faultX;
  faultY_c = faultY;
  faultZ_c = faultZ + ss.dz/2;

  % --     kvf params     --
  c.command = 'compress';
  c.lamdbaZ = lamdbaZ; ss.lamdbaZ = lamdbaZ;
  c.dz = s.dz;
  c.tol = 1e-8;
  c.G = 30e3;
  c.allow_overwrite = 1;
  c.err_method = 'mrem-frp';

  % --      s12 kernel for fault-fault interaction      --
  c.greens_fn  = 'okadaS12';
  c.write_hmat_filename = './tmp/FO_ff-s12';
  c.write_hd_filename = './tmp/FO_ff-s12-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.Y = [faultX; faultY; faultZ];
  c.X = [faultX_c; faultX_c; faultZ_c];

  c.kvf('Write', c.kvf, c, 4);
  disp('run this in a shell:')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.s12 = c.write_hmat_filename;

  r.ss = ss;

end

% numerical solution
function out = solve(r)
  addpaths();

  % - frictional parameters -
  % reference friction coefficient
  ss.fo = 0.6*ones(size(y3));

  % Dieterich-Ruina RS frictional params (vw friction)
  ss.a = 1e-3*ones(size(r.ss.y3f));
  ss.b = r.ss.a - 2.1e-4*ones(size(r.ss.y3f)); % make + for with vw region

  % effective normal stress (MPa)
  ss.sigma=50.0*ones(size(y3));

  % characteristic weakening distance (m)
  ss.L=8e-3*ones(size(y3));

  % plate rate (m/s)
  ss.Vpl=1e-9*ones(size(y3));

  % reference slip rate (m/s)
  ss.Vo=1e-6*ones(size(y3));

  % rigidity (MPa)
  G = 30e3;

  % Radiation damping coefficient
  ss.eta = G./(2*Vs);

  r.ss.dgf = 4;

  Y0 = zeros(r.c.n*r.ss.dgf,1);
  Y0(1:r.ss.dgf:end) = zeros(r.c.n,1);
  Y0(2:r.ss.dgf:end) = max(ss.a).*ss.sigma.*asinh(ss.Vpl./ss.Vo/2.*exp((ss.fo+ss.b.*log(ss.Vo./ss.Vpl))./max(ss.a))) + ss.eta.*ss.Vpl;
  Y0(3:r.ss.dgf:end)=ss.a./ss.b.*log(2*ss.Vo./ss.Vpl.*sinh((Y0(2:ss.dgf:end)-ss.eta.*ss.Vpl)./ss.a./ss.sigma))-ss.fo./ss.b;
  Y0(4:r.ss.dgf:end)=log(ss.Vpl./ss.Vo);

  % load kernel as hm
  hm = hmmvp('init', r.c.write_hmat_filename);

  yp=@(t,y) DieterichRuinaRegAging(t,y,r,hm);

  tic % val had err at 1e-8
  options=odeset('Refine',1,'RelTol',1e-6,'InitialStep',1e-5);
  [t,Y]=ode45(yp,[0 1e10],Y0,options); %
  disp('Done solving');
  toc

end

% ------------------ Private -------------------------
function addpaths()
  addpath('../hmmvp-okada/matlab')
end

% boxcar function
function bc = BC(x)
  boxc=@(x) (x+0.5>=0)-(x-0.5>=0);
  bc = boxc(x);
end
