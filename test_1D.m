% test_1d.m
%
% Joseph Wick
% 7/22/2021
%
% testing if dc3dm can be used on a 2D problem with a 1D fault
% Basic operations are as follows:
%  - set up h-matrix for a 1D fault
%  - multiply this fault against a simple slip vector as a sanity check

% argout
function varargout = t1d (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ----------------------- Public -------------------------------

% write_kvfs ()
% creates kvf's for hmatrix and prints shell commands to build hmatrix
function r = create_kvfs ()
  addpaths();
  path_to_dc3dm = '../dc3dm-main/bin/dc3dm';

  p = setup_problem();

  r.cm = write_mesh_kvf(p);
  r.cb = write_build_kvf(p);
  r.cc = write_compress_kvf(p, r.cb);

  clf;
  imagesc(p.x, p.y, p.f); title('res function'); colorbar;
  saveas(gcf, 'figures/test_1D_res.png')

  disp('run these in a shell in this directory:')
  cmds = 'mbc';
  for (i = 1:numel(cmds))
    fprintf('     ../dc3dm-main/bin/dc3dm %s.kvf\n', r.(['c' cmds(i)]).kvf);
  end
end

% mvp_test ()
function y = mvp_test (r)
  addpaths();

  % first load in the h-matrix
  hm_fname = r.cc.hm_write_filename;
  hm = hmmvp('init', hm_fname, 4);
  m = hmmvp('getm', hm);
  n = hmmvp('getn', hm);

  % make a vector to multiply by
  a = zeros(100,100);
  a(45:54,45:54) = 1;
  x = reshape(a, [10000, 1]);

  %x = [0,0,0,0, 0,1,1,0, 0,1,1,0, 0,0,0,0];
  %x = x';

  y = hmmvp('mvp', hm, x);

  clf;
  plot(y(5000:5100));
  saveas(gcf, 'figures/test_1D_mvp-5k51.png')

  clf;
  imagesc(a); colorbar;
  saveas(gcf, 'figures/test_1D_a.png')

  hmmvp('cleanup', hm);

end

% get_elem_sizes ()
% plot the sizes of elements in the hmatrix (should be uniform)
function get_elem_sizes (r)

  clf;
  dc3dm.mViewBuild(r.cb); axis xy;
  saveas(gcf, 'figures/test_1D_sizes.png')

end

% get_mn ()
% returns m, n of the hmatrix
function get_mn (r)
  addpaths();

  hm_fname = r.cc.hm_write_filename;
  hm = hmmvp('init', hm_fname, 4);
  m = hmmvp('getm', hm)
  n = hmmvp('getn', hm)

end

% get_whole_hm ()
% returns the entire hm
% prolly dont do this unless its pretty small
function m = get_whole_hm (r)
  addpaths();

  hm_fname = r.cc.hm_write_filename;
  hm = hmmvp('init', hm_fname, 4);
  m = hmmvp('getm', hm);
  n = hmmvp('getn', hm);

  rs = (1:1:m); cs = (1:1:n);
  m = hmmvp('extract', hm, rs, cs);

  clf;
  imagesc(m); title('hm'); colorbar;
  saveas(gcf, 'figures/test_1D_hm.png')

end

% ----------------------- Private -------------------------------

% addpaths ()
% add path to dc3dm matlab files
function addpaths()
  addpath /auto/home/jmwick/EarthquakeSimulation/dc3dm-main/matlab;
end

% make_base_fn ()
% creates name(s) for hm related files
function bfn = make_base_fn (p)
  bfn = sprintf('%stest_1d', p.dir);
end

% setup_problem ()
% sets parameters for the fault
function p = setup_problem()

  p.rfac = 2;
  p.want_free_surface = 0;
  p.do_uniform = 1;
  p.neighborhood = 8;
  p.min_len = 0;
  p.max_len = 10;
  p.nthreads = 4;
  p.dir = './tmp/';
  p.dip_len = 1000;
  p.strike_len = 1000;

  p.tol = 1e-5;

  p.problem = 0;

  n = 1001;
  p.x = linspace(-0.5*p.strike_len, 0.5*p.strike_len, n);
  p.y = linspace(-0.5*p.dip_len, 0.5*p.dip_len, n);

  p.mu = 3e10;
  p.nu = 0.25;
  %'finding' nucleation size
  h_star = 2;

  p.f = (h_star*ones(n,n))/p.rfac;

end

% write_mesh_kvf ()
function c = write_mesh_kvf(p)
  c.mesh_write_filename = make_base_fn(p);
  c.do_uniform = p.do_uniform;

  c.min_len = p.min_len;
  c.max_len = p.max_len;

  c.x = p.x;
  c.y = p.y;

  c.f = p.f;

  c.command = 'mesh';
  c.kvf = [make_base_fn(p) '_m'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

% write_build_kvf ()
function c = write_build_kvf(p)
  bfn = make_base_fn(p);
  c.mesh_read_filename = bfn;
  c.build_write_filename = sprintf('%s_p%d', bfn, p.problem);

  c.depth_min = 10;
  c.dipdeg = 90;
  c.svbc = 2;

  c.neighborhood = p.neighborhood;
  c.bc_periodic_nlayers = 3;

  c.command = 'build';
  c.kvf = [bfn '_b'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

% write_compress_kvf ()
function c = write_compress_kvf (p, cb)
  v = [1 2 0];
  v = v/norm(v);
  c.src_disl = v;
  c.rcv_traction = c.src_disl;
  c.component = 1;

  bfn = make_base_fn(p);
  c.mesh_read_filename = bfn;
  c.build_read_filename = cb.build_write_filename;
  c.tol = p.tol;
  c.hm_write_filename = sprintf( ...
    '%s_tol%1.1f', cb.build_write_filename, -log10(c.tol));

  c.allow_overwrite = 1;
  c.mu = p.mu;
  c.nu = p.nu;
  c.nthreads = p.nthreads;

  c.command = 'compress';
  c.kvf = [bfn '_c'];
  dc3dm.WriteKvf(c.kvf, c, true);
end
