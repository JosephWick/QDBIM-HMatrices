function varargout = exmb (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% addpaths()
% add path to dc3dm matlab
function addpaths()
  addpath /auto/home/jmwick/EarthquakeSimulation/dc3dm-main/matlab;
end

function r = create_kvfs ()
  addpaths();
  path_to_dc3dm = '../dc3dm-main/bin/dc3dm';

  o = setopts();
  p = make_props(o);

  clf;
  subplot(221); imagesc(p.x, p.y, p.a - p.b); title('a - b'); colorbar;
  subplot(222); imagesc(p.x, p.y, p.a./p.b); title('a/b'); colorbar;
  subplot(223); imagesc(p.x, p.y, p.sigma); title('\sigma'); colorbar;
  subplot(224); imagesc(p.x, p.y, p.h_star); title('h^*_b'); colorbar;
  saveas(gcf, 'figures/test2_fig1.png')

  r.cm = write_mesh_kvf(o, p);
  r.cb = write_build_kvf(o);
  r.cc = write_compress_kvf(o, r.cb);

  disp('Run these in a shell in this directory')
  cmds = 'mbc';
  for (i = 1:numel(cmds))
    fprintf('../dc3dm-main/bin/dc3dm %s.kvf\n', r.(['c' cmds(i)]).kvf);
  end

end

function get_hm(r)
  addpaths();

  hm_filename = r.cc.hm_write_filename;
  id = hmmvp('init', hm_filename, 4);
  m = hmmvp('getm', id)
end

% setopts()
% sets up basic values for simulation
% not entiirely sure what a lot of these do
function o = setopts ()
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
% sets the frictional and other properties of the domain
% also not entirely sure what they all do
function p = make_props (o)

  dip_len = o.len_fac*1e3;
  strike_len = o.len_fac*1e3;
  n = 10;                           %num cells

  p.x = linspace(-0.5*strike_len, 0.5*strike_len, n);
  p.y = linspace(-0.5*dip_len, 0.5*dip_len, n);
  [X Y] = meshgrid(p.x, p.y)

  %radius = 0.1*o.len_fac*dip_len; % radius of disk
  %t_width = 0.5*radius; % transition width

  sidelen = 0.1*o.len_fac*dip_len*2;
  t_width = 0.25*sidelen;

  r = sqrt(X.^2 + Y.^2) / sidelen
  alpha = calc_transition_width(t_width/sidelen, 0.99);
  w_sigma = calc_sigmoid(r, 1, alpha, 0, max(r(:)), 0, 1)
  w_amb = calc_sigmoid(r, 2, alpha, 0, max(r(:)), 0, 1)

  one = ones(size(X));
  bl = 0.01; bs = 0.005;
  amb_vw = -0.005; amb_vs = 0.005;
  d_c = 1e-4;
  sigma_s = 1e6; sigma_l = o.vary_fac*sigma_s;
  p.mu = 3e10;
  p.nu = 0.25;
  p.b = bl*(1 - w_amb) + bs*w_amb;
  %p.b=0.015*ones(size(y3));
  p.a = p.b + amb_vw*(1 - w_amb) + amb_vs*w_amb;
  %p.a = 1e-2+Ramp((y3-15e3)/3e3)*(0.025-0.01);

  p.a

  p.d_c = d_c*one;
  p.sigma = sigma_l*(1 - w_sigma) + sigma_s*w_sigma;
  p.h_star = 1.377*p.mu/(1 - p.nu)*p.d_c./(p.sigma.*p.b);

  %aLin = linspace(100, 0, n);
  %bLin = linspace(0, 100, n);
  %[z p.a] = meshgrid(aLin)
  %[z p.b] = meshgrid(bLin)

end

function p = my_props(o)
  Lambdaz=40e3;         % Depth extent of the frictional domain
  M=10;                % Number of cells (400)
  dz=Lambdaz/M;         % Cell size
  y3=(0:M-1)'*dz;       % Top of slip patch
  W=ones(M,1)*dz;       % Down-dip width of slip patch

  len = o.len_fac*1e3;

  p.x = linspace(-0.5*strike_len, 0.5*strike_len, n);
  p.y = linspace(-0.5*dip_len, 0.5*dip_len, n);
  [X Y] = meshgrid(p.x, p.y)

  one = ones(size(X));
  bl = 0.01; bs = 0.005;
  amb_vw = -0.005; amb_vs = 0.005;
  d_c = 1e-4;
  sigma_s = 1e6; sigma_l = o.vary_fac*sigma_s;
  p.mu = 3e10;
  p.nu = 0.25;
  p.b=0.015*one(size(y3)); %from val's QDBIM2D
  p.a = meshgrid(1e-2+Ramp((y3-15e3)/3e3)*(0.025-0.01)); % from val's QDBIM2D

  p.d_c = d_c*one; %from ex.m

  p.sigma = sigma_l*(1 - w_sigma) + sigma_s*w_sigma; %from ex.m
  p.h_star = 1.377*p.mu/(1 - p.nu)*p.d_c./(p.sigma.*p.b); % from ex.m

end

% calc_transition_width()
function alpha = calc_transition_width (width, at_p)
% alpha = calc_transition_width(width, at_p)
%   The transition has width 'width' in the sense that
%       abs(diff(y(x +/- width/2))) = at_p * abs(y1 - ye).
% So at_p should be something like 0.9.
%   We do this calculation for the exponent p = 1 only.
  assert(at_p > 0 && at_p < 1);
  assert(width > 0);
  alpha = -2/width*log(2/(1 + at_p) - 1);
  assert(alpha > 0);
end

% calc_sigmoid()
function y = calc_sigmoid (x, xt, a, xs, xe, ys, ye)
% y = calc_sigmoid(x, xt, alpha, p, xs, xe, ys, ye)
%   xt is the transition point.
%   ys(xs), ye(xe) are the values at reference points.
%   alpha is the constant specifying the transition width. See
% calc_transition_width for more.
  assert(xe > xs);
  fn = @(x) 1./(1 + exp(-a.*x));
  y0s = fn(xs - xt);
  y0e = fn(xe - xt);
  y = ys + (ye - ys).*(fn(x - xt) - y0s)./(y0e - y0s);
end

% ramp function
function r = Ramp(x)
  r = x.*BC(x-1/2)+HS(x-1);
end

function bc = BC(x)
  bc = (x+0.5>=0)-(x-0.5>=0);
end

% write_mesh_kvf
% creates the key-value file for dc3dm mesh
function c = write_mesh_kvf (o, p)
  c.mesh_write_filename = make_base_fn(o);
  c.do_uniform = o.do_uniform;
  % Min and max element lengths don't really matter unless they are used to
  % bound the resolution function f. Here f is well behaved so I set the min
  % to 0 and make sure there are at least 8 elements in each direction of the
  % domain. (max_len is all that matters if the do_uniform option is true.)
  c.min_len = 0;
  if (isinf(o.max_len))
    c.max_len = min(diff(p.x([1 end])), diff(p.y([1 end])))/8;
  else
    c.max_len = o.max_len;
  end
  % Create a tensor mesh ...
  c.x = p.x;
  c.y = p.y;
  % ... on which f, the resolution function, is set. f has the same units as
  % x and y. Make sure there are 5 o.rfac elements in each direction per h*.
  c.f = p.h_star/(o.rfac*5);
  c.command = 'mesh';
  c.kvf = [make_base_fn(o) '_m'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

% write_build_kvf
% creatse the key-value file for dc3dm build
function c = write_build_kvf (o)
  bfn = make_base_fn(o);
  c.mesh_read_filename = bfn;
  c.build_write_filename = sprintf('%s_p%d', bfn, o.problem);

  % Setting depth_min to 0 makes the 0-depth boundary (N or S depending on
  % the sign of dipdeg) have a free surface boundary condition.
  c.depth_min = 100;
  if (o.want_free_surface) c.depth_min = 0; end
  switch (o.problem)
    case 1 % subduction
      % Positive dip makes the north boundary be at the surface and the south
      % boundary to have the velocity BC.
      c.dipdeg = 12;
      % The fault is periodic along-strike.
      c.ewpbc = 0;
      % This is the velocity boundary condition at depth.
      c.svbc = 2;
    otherwise
      error(sprintf('%d is not a valid problem number.', o.problem));
  end

  c.neighborhood = o.neighborhood;
  c.bc_periodic_nlayers = 3;

  c.command = 'build';
  c.kvf = [bfn '_b'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

% write_compress_kvf
% creates the key-value file for dc3dm compress
function c = write_compress_kvf (o, cb)
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
  c.hm_write_filename = sprintf( ...
    '%s_tol%1.1f', cb.build_write_filename, -log10(c.tol));

  c.allow_overwrite = 1;
  c.mu = 3e10;
  c.nu = 0.25;
  c.nthreads = o.nthreads;

  c.command = 'compress';
  c.kvf = [bfn '_c'];

  dc3dm.WriteKvf(c.kvf, c, true);
end

% make_base_fn()
% generates a file name based on the inputted options
function bfn = make_base_fn (o)
  bfn = sprintf('%sdc3t_rf%1.2flf%1.2fvf%1.2fnbr%du%d', o.dir, o.rfac, ...
                o.len_fac, o.vary_fac, o.neighborhood, o.do_uniform);
end
