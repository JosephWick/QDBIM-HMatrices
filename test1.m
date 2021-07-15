function varargout = exmb (varargin)

  addpath /auto/home/jmwick/EarthquakeSimulation/dc3dm-main/matlab;

end


% create_kvfs()
% DIR: >> r = test1('create_kvfs');
% creates 3 key-value files:
%  - r.cm: mesh file
%  - r.cb: fault build
%  - r.cc: compressed elasticity matrix
function r = create_kvfs ()
  o = setopts();
  p = make_props(o);

  clf;
  subplot(221); imagesc(p.x, p.y, p.a - p.b); title('a - b'); colorbar;
  subplot(222); imagesc(p.x, p.y, p.a./p.b); title('a/b'); colorbar;
  subplot(223); imagesc(p.x, p.y, p.sigma); title('\sigma'); colorbar;
  subplot(224); imagesc(p.x, p.y, p.h_star); title('h^*_b'); colorbar;

  r.cm = write_mesh_kvf(o, p);
  r.cb = write_build_kvf(o);
  r.cc = write_compress_kvf(o, r.cb);

end

% -------------- PRIVATE ---------------------

% setopts()
% sets up basic values for simulation
% not entiirely sure what a lot of these do
function o = setopts ()
  o.rfac = 1;
  o.len_fac = 1;
  o.vary_fac = 2;
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
  n = 1001;

  p.x = linspace(-0.5*strike_len, 0.5*strike_len, n);
  p.y = linspace(-0.5*dip_len, 0.5*dip_len, n);
  [X Y] = meshgrid(p.x, p.y);

  radius = 0.1*o.len_fac*dip_len; % radius of disk
  t_width = 0.5*radius; % transition width

  r = sqrt(X.^2 + Y.^2) / radius;
  alpha = calc_transition_width(t_width/radius, 0.99);
  w_sigma = calc_sigmoid(r, 1, alpha, 0, max(r(:)), 0, 1);
  w_amb = calc_sigmoid(r, 2, alpha, 0, max(r(:)), 0, 1);

  one = ones(size(X));
  bl = 0.01; bs = 0.005;
  amb_vw = -0.005; amb_vs = 0.005;
  d_c = 1e-4;
  sigma_s = 1e6; sigma_l = o.vary_fac*sigma_s;

  p.mu = 3e10;
  p.nu = 0.25;
  p.b = bl*(1 - w_amb) + bs*w_amb;
  p.a = p.b + amb_vw*(1 - w_amb) + amb_vs*w_amb;
  p.d_c = d_c*one;
  p.sigma = sigma_l*(1 - w_sigma) + sigma_s*w_sigma;
  p.h_star = 1.377*p.mu/(1 - p.nu)*p.d_c./(p.sigma.*p.b);
end
