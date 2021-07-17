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
  saveas(gcf, 'figures/test1_fig1.png')

  r.cm = write_mesh_kvf(o, p);
  r.cb = write_build_kvf(o);
  r.cc = write_compress_kvf(o, r.cb);

  disp('Run these in a shell in this directory')
  cmds = 'mbc';
  for (i = 1:numel(cmds))
    fprintf('../dc3dm-main/bin/dc3dm %s.kvf\n', r.(['c' cmds(i)]).kvf);
  end

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
  n = 10;

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
  p.a = p.b + amb_vw*(1 - w_amb) + amb_vs*w_amb;
  p.d_c = d_c*one;
  p.sigma = sigma_l*(1 - w_sigma) + sigma_s*w_sigma;
  p.h_star = 1.377*p.mu/(1 - p.nu)*p.d_c./(p.sigma.*p.b);

  aLin = linspace(100, 0, n);
  bLin = linspace(0, 100, n);
  p.a = meshgrid(aLin)
  p.b = meshgrid(bLin)

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
