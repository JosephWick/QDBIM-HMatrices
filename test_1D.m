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

  disp('run these in a shell in this directory')
  cmds = 'mbc'
  for (i = 1:numel(cmds))
    fprintf('../dc3dm-main/bin/dc3dm %s.kvf\n', r.(['c' cmds(i)]).kvf);
  end
end

% ----------------------- Private -------------------------------
function addpaths()
  addpath /auto/home/jmwick/EarthquakeSimulation/dc3dm-main/matlab;
end

function bfn = make_base_fn ()
  bfn = 'test_1d';
end

function p = setup_problem()

  p.rfac = 2;
  p.want_free_surface = 1;
  p.do_uniform = 0;
  p.neighborhood = 8;
  p.max_len = inf;
  p.nthreads = 4;
  p.dir = './tmp';
  p.dip_len = 100;

  p.tol = 1e-5;

  p.problem = 0;

  n = 101
  p.x = [1.0]
  p.y = linspace(-0.5*p.dip_len, 0.5*p.dip_len, n);

  p.mu = 3e10;
  p.nu = 0.25;
  %vw region
  p.f = rectangularPulse(0, 30, p.y);

end

% write_mesh_kvf ()
function c = write_mesh_kvf(p)
  c.mesh_write_filename = make_base_fn;
  c.do_uniform = p.do_uniform;

  c.min_len = 0;
  c.max_len = p.max_len;

  c.x = p.x;
  c.y = p.y;

  c.f = p.f;

  c.command = 'mesh';
  c.kvf = [make_base_fn() '_m'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

% write_build_kvf ()
function c = write_build_kvf(p)
  bfn = make_base_fn();
  c.mesh_read_filename = bfn;
  c.build_write_filename = sprintf('%s_p%d', bfn, p.problem);

  c.depth_min = 0;
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

  bfn = make_base_fn();
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
