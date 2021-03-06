% test3.m
% can I make a simpler version of the ex.m problem

function varargout = t3 (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ---------------- Public -------------

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

  r.cm = write_mesh_kvf(o, p);
  r.cb = write_build_kvf(o);
  r.cc = write_compress_kvf(o, r.cb);

  disp('Run these in a shell in  this directory')
  cmds = 'mbc';
  for (i = 1:numel(cmds))
    fprintf('../dc3dm-main/bin/dc3dm %s.kvf\n', r.(['c' cmds(i)]).kvf);
  end

  clf;
  imagesc(p.x, p.y, p.a - p.b); title('a - b'); colorbar;
  saveas(gcf, 'figures/test3_fig1.png')

  clf;
  imagesc(p.x, p.y, p.h_star); title('h*'); colorbar;
  saveas(gcf, 'figures/test3_hstar.png')

end

% timeMVP()
% perform a test MVP operation
function time_MVP (r)
  addpaths();

  % load H-matrix
  hm_file = r.cc.hm_write_filename;
  hm = hmmvp('init', hm_file);

  % make a test vector and multiply
  n = hmmvp('getn', hm);
  x = randn(n, 1);
  % time hmmvp multiplication
  disp('HMMVP Multiplication')
  tic
  y = hmmvp('mvp', hm, x);
  toc

  % compare against non-hmmvp multiplication
  m = randn(n,n);
  x = randn(n,1);
  disp('Non-HMMVP Multiplication')
  tic
  y2 = m*x;
  toc

  % clean memory
  hmmvp('cleanup', hm);

end

function r = test_nm()
  addpaths();

  c.X = [linspace(1,10,10); linspace(1,10,10); linspace(1,10,10)];
  c.Y = [linspace(1,5,5);   linspace(1,5,5);   linspace(1,5,5)];

  c.greens_fn = 'test';
  c.write_hmat_filename = './tmp/test_nm';
  c.write_hd_filename = './tmp/test_nm-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.allow_overwrite = 1;
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)

  r = c.write_hmat_filename;

end

function load_nm(r)
  addpaths();

  hm = hmmvp('init', r, 4);

  disp(hmmvp('getn', hm))
  disp(hmmvp('getm', hm))

end

% ---------------------- Private but still relevant ---------------

% setopts()
% set up some options
% pretty much ripped from ex.m
% all lengths are in [m]
function o = setopts()
  o.rfac = 2; % not sure what this does
  o.len_fac = 1; % seems to change the size of the vw region
  o.vary_fac = 2; % impacts scale of sigma
  o.want_free_surface = 1; % boolean for if one of the surface is free
  o.tol = 1e-5;
  o.problem = 1; % used in the write_kvf methods to define some properties
  o.do_uniform = 0; % boolean for uniform mesh
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

  radius = 0.2*o.len_fac*dip_len;

  one = ones(size(X));
  bl = 0.01; bs = 0.005;
  amb_vw = -0.005; amb_vs = 0.005;
  d_c = 1e-4;
  sigma_s = 1e6; sigma_l = o.vary_fac*sigma_s;

  w_amb = 1; % left over from sigmoid transition stuff

  p.mu = 3e10;
  p.nu = 0.25;

  % b will be constant, a  will be constant except for the circle
  p.b = 15*one;
  p.a = 2*one;
  % making circle
  circle = X.^2 + Y.^2 <= radius.^2;
  p.a = p.a+(circle*10);

  p.d_c = d_c*one;
  p.sigma = 50*one;
  p.h_star = -100*((p.a - p.b)+1);
  % here, h_star is const (bc sigma const) and as a result the resolutoin
  % func will also be constant

end

% write_mesh_kvf
function c = write_mesh_kvf(o, p)
  % get a filename
  c.mesh_write_filename = make_base_fn(o);

  c.do_uniform = o.do_uniform;

  c.min_len = 0;
  if (isinf(o.max_len))
    c.max_len = min(diff(p.x([1 end])), diff(p.y([1 end])))/8;
  else
    c.max_len = o.max_len;
  end

  % make a tensor mesh, on which the resolution function f is set
  % f has same units as x and y
  c.x = p.x;
  c.y = p.y;
  c.f = p.h_star/(o.rfac*5);

  clf;
  imagesc(p.x, p.y, c.f); title('resolution function f'); colorbar;
  saveas(gcf, 'figures/test3_resFunc.png')

  c.command='mesh';
  c.kvf = [make_base_fn(o) '_m'];

  dc3dm.WriteKvf(c.kvf, c, true);

end

% write_build_kvf
function c = write_build_kvf(o)

  bfn = make_base_fn(o);
  c.mesh_read_filename = bfn;
  c.build_write_filename = sprintf('%s_p%d', bfn, o.problem);

  % depth_min == 0 means there is a free surface at the zero-depth boundary
  c.depth_min=100;
  if (o.want_free_surface) c.depth_min = 0; end
  switch (o.problem)
    case 1  % subduction, also the only option rn
      % dipdeg>0 means that north boundary is at surface and south boundary has
      % the velocity coundition below
      c.dipdeg = 90;
      % fault is periodic along-strike
      c.ewpbc = 0;
      % this is velocity boundary condition at depth
      c.svbc = 0;
    otherwise
      error(sprintf('%d is not a valid problem number.', o.problem));
  end

  c.neightborhood = o.neighborhood;
  c.bc_periodic_nlayers = 3;

  c.command = 'build';
  c.kvf = [bfn '_b'];

  dc3dm.WriteKvf(c.kvf, c, true);

end

% write_compress_kvf
function c = write_compress_kvf(o, cb)

  switch (o.problem)
    case 1
      v = [1 2 0];
      v = v/norm(v);
      c.src_disl = v;
      c.rcv_traction = c.src_disl;
      c.component = 1;
    otherwise
      error(sprintf('%d is not a valid problem number.', o.problem));
    end

    bfn = make_base_fn(o);
    c.mesh_read_filename = bfn;
    c.build_read_filename = cb.build_write_filename;
    c.tol = o.tol;

    c.hm_write_filename = sprintf('%s_to1$1.1f, cb.build_read_filename', ...
     -log10(c.tol));

    c.allow_overwrite = 1;
    c.mu = 3e10;
    c.nu = 0.25;
    c.nthreads = o.nthreads;

    c.command = 'compress';
    c.kvf = [bfn '_c'];
    dc3dm.WriteKvf(c.kvf, c, true);

end

% ----------------- Helpers ------------------

% addpaths()
% adds path to the dc3dm/matlab directory
function addpaths()
  addpath /auto/home/jmwick/EarthquakeSimulation/dc3dm-main/matlab;
end

% make_base_fn()
% generates a file name based on the inputted options
function bfn = make_base_fn (o)
  bfn = sprintf('%sdc3t_rf%1.2flf%1.2fvf%1.2fnbr%du%d', o.dir, o.rfac, ...
                o.len_fac, o.vary_fac, o.neighborhood, o.do_uniform);
end
