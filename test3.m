% test3.m
% can I make a simpler version of the ex.m problem

function varargout = t3 (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% create_kvfs
% create 3 key-value files:
%  - r.cm: mesh file
%  - r.cb: fault build
%  - r.cc: compressed elasticity matrix
% cut down version of that in ex.m
function r = create_kvfs ()
  addpaths();
  path_to_dc3dm = '../dc3dm-main/bin/dc3dm';

  o = setopts();
  p = make_props(o);

  %r.cm = write_mesh_kvf(o, p);
  %r.cb = write_build_kvf(o);
  %r.cc = write_compress_kvf(o, r.cb);

  disp('Run these in a shell in  this directory')
  cmds = 'mbc';
  for (i = 1:numel(cmds))
    fprintf('../dc3dm-main/bin/dc3dm %s.kvf\n', r.(['c' cmds(i)]).kvf);
  end

  clf;
  imagesc(p.x, p.y, p.a - p.b); title('a - b'); colorbar;
  saveas(gcf, 'figures/test3_fig1.png')

end

% setopts()
% set up some options
% pretty much ripped from ex.m
% all lengths are in [m]
function o = setopts()
  o.rfac = 2; % not sure what this does
  o.len_fac = 1; % seems to change the size of the vw region
  o.vary_fac = 2; %impacts scale of sigma
  o.want_free_surface = 1;
  o.tol = 1e-5;
  o.problem = 1;
  o.do_uniform = 0;
  o.neighborhood = 8;
  o.max_len = inf;
  o.nthreads = 4;
  o.dir = './tmp/';
end

% make_props()
% set up frictional and other properties for the fault
% mostly ripped from ex.m with some changes:
%  -
function p = make_props(o)

  dip_len = o.len_fac*1e3;
  strike_len = o.len_fac*1e3;
  n = 1001;

  p.x = linspace(-0.5*strike_len, 0.5*strike_len, n);
  p.y = linspace(-0.5*dip_len, 0.5*dip_len, n);
  [X Y] = meshgrid(p.x, p.y);

  one = ones(size(X));
  bl = 0.01; bs = 0.005;
  amb_vw = -0.005; amb_vs = 0.005;
  d_c = 1e-4;
  sigma_s = 1e6; sigma_l = o.vary_fac*sigma_s;

  w_amb = 1; % left over from sigmoid transition stuff  

  p.mu = 3e10;
  p.nu = 0.25;
  p.b = bl*(1 - w_amb) + bs*w_amb;
  p.a = p.b + amb_vw*(1 - w_amb) + amb_vs*w_amb;
  p.d_c = d_c*one;
  p.sigma = 50*one;
  p.h_star = 1.377*p.mu/(1 - p.nu)*p.d_c./(p.sigma.*p.b);

end



% ----------------- Helpers ------------------

% addpaths()
% adds path to the dc3dm/matlab directory
function addpaths()
  addpath /auto/home/jmwick/EarthquakeSimulation/dc3dm-main/matlab;
end
