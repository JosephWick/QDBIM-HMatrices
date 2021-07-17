function varargout = exmb (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end


% create_kvfs()
% DIR: >> r = test1('create_kvfs');
% creates 3 key-value files:
%  - r.cm: mesh file
%  - r.cb: fault build
%  - r.cc: compressed elasticity matrix
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

% analyze()
function analyze (r)
  addpaths();
  % see element sizes
  clf;
  dc3dm.mViewBuild(r.cb); axis xy;
  saveas(gcf, 'figures/test1_fig2.png')

  r.cb.build_write_filename

  % read mesh produced by dc3dm build
  rid = dc3dm.mRead(r.cb.build_write_filename);
  % get the elements
  rs = dc3dm.mRects(rid);
  % element centers:
  [cx cy] = dc3dm.mCC(rs);
  dc3dm.mClear(rid);

  % Compare these data with those in the .elem file, which is a plain text
  % file formatted as comma separated values. This format is easy to parse,
  % though Matlab happens to have a built-in reader for it, csvread.
  %%ers = parse_elem_file(r.cb.build_write_filename);
  % Confirm ers contains equivalent information to rs.
  %%fprintf('Should be <= ~%1.1e: %1.1e %1.1e %1.1e %1.1e\n', ...
  %%        eps(10), relerr(cx(:), ers(1,:)'), relerr(cy(:), ers(2,:)'), ...
  %%        relerr(rs(3,:), ers(3,:)), relerr(rs(4,:), ers(4,:)));

end

% get_hMatrix
% load hMatrix into memory from r.cc elasticity .kvf
function get_hMatrix (r)
  addpaths();
  hm_filename = r.cc.hm_write_filename;

  % Load the H-matrix A into memory. id is a pointer to this H-matrix. Use 4
  % threads if they are available.
  id = hmmvp('init', hm_filename, 4);

end

function r = run_test ()
  addpaths();

  o = setopts();
  o.want_free_surface = 0;
  r.p = make_props(o);

  % AIGA-8: approximate IGA with neightborhood = 8
  o.neighborhood =8;
  r.n8.o = 0;
  r.n8.cm = write_mesh_kvf(o, r.p);
  r.n8.cb = write_build_kvf(o);
  r.n8.cc = write_compress_kvf(o, r.n8.cb);

  clf;
  subplot(221); imagesc(r.p.x, r.p.y, r.p.a - r.p.b); title('a - b'); colorbar;
  subplot(222); imagesc(r.p.x, r.p.y, r.p.a./r.p.b); title('a/b'); colorbar;
  subplot(223); imagesc(r.p.x, r.p.y, r.p.sigma); title('\sigma'); colorbar;
  subplot(224); imagesc(r.p.x, r.p.y, r.p.h_star); title('h^*_b'); colorbar;
  saveas(gcf, 'figures/test1_fig3.png')

  % AIGA - 0: approximate IG with neighborhood = 0
  o.neightborhood = 0
  r.n0.o = o
  r.n0.cm = write_mesh_kvf(o, r.p);
  r.n0.cb = write_build_kvf(o);
  r.n0.cc = write_compress_kvf(o, r.n0.cb);



end

function demo_mvp_slip (r)
%DIR Carry out typical operations for a real problem. Run
% >> test1('demo_mvp_slip', r);

  % Get element centers. (Could use the .elem file.)
  rid = dc3dm.mRead(r.cb.build_write_filename);
  rs = dc3dm.mRects(rid);
  [cx cy] = dc3dm.mCC(rs);
  md = dc3dm.mData(rs);

  % Make a slip distribution. Respect the BCs.
  slip_fn = @(x, y) cos(2*pi*(x + 0.39*diff(md.xlim))/diff(md.xlim)).* ...
            sin(2*pi*y/diff(md.ylim)) + ...
            0.*y/diff(md.ylim);
  slip = slip_fn(cx, cy);

  % Get boundary condition data.
  bc = dc3dm.ReadBoundaryConditions(r.cc.hm_write_filename);
  % Boundary values are ordered (east, north, west, south); only those for
  % non-0 velocity-BC matter. Here I fill in all values even though only the
  % S one is necessary.
  bdy_vals = [slip_fn(md.xlim(2), 0), slip_fn(0, md.ylim(2)), ...
              slip_fn(md.xlim(1), 0), slip_fn(0, md.ylim(1))];

  % Compute traction.
  id = hmmvp('init', r.cc.hm_write_filename, 4, 1);
  traction = hmmvp('mvp', id, slip) + bc*bdy_vals(:);
  hmmvp('cleanup', id);

  % Plot. Show the raw (const interp) mesh values and two different
  % interpolations.
  x = CC(linspace(md.xlim(1), md.xlim(2), round(diff(md.xlim)/md.dx) + 1));
  y = CC(linspace(md.ylim(1), md.ylim(2), round(diff(md.ylim)/md.dy) + 1));
  [X Y] = meshgrid(x, y);
  slip_c = dc3dm.mConstInterp(rid, slip, X, Y);
  traction_c = 1e-6*dc3dm.mConstInterp(rid, traction, X, Y);
  slip_1 = dc3dm.mLinterp(rid, slip, bdy_vals, X, Y);
  traction_1 = 1e-6*dc3dm.mLinterpWExtrap(rid, traction, X, Y);
  slip_3 = dc3dm.mCinterp(rid, slip, bdy_vals, X, Y);
  traction_3 = 1e-6*dc3dm.mCinterpWExtrap(rid, traction, X, Y);

  dc3dm.mClear(rid);

  clf;
  img = @(im) imagesc(x, y, im);
  h(1) = subplot(321); img(slip_c); title('slip, mesh resolution');
  h(2) = subplot(322); img(traction_c); title('traction, mesh resolution');
  h(3) = subplot(323); img(slip_1); title('slip, linear interp');
  h(4) = subplot(324); img(traction_1); title('traction, linear interp');
  h(5) = subplot(325); img(slip_3); title('slip, cubic interp');
  h(6) = subplot(326); img(traction_3); title('traction, cubic interp');
  linkaxes(h); zoom on;
  for (i = 1:numel(h))
    subplot(h(i));
    if (mod(i, 2) == 0) caxis(135*[-1 1]); else caxis([-1 1]); end
    colorbar;
    axis equal; axis xy; axis tight;
    % Draw the mesh.
    draw_rects_r(rs, 0, 'k');
  end
  saveas(gcf, 'figures/test1_fig4.png')

end


% -------------- PRIVATE ---------------------

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

% make_base_fn()
% generates a file name based on the inputted options
function bfn = make_base_fn (o)
  bfn = sprintf('%sdc3t_rf%1.2flf%1.2fvf%1.2fnbr%du%d', o.dir, o.rfac, ...
                o.len_fac, o.vary_fac, o.neighborhood, o.do_uniform);
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
  if (o.want_free_surface) o.depth_min = 0; end
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

function addpaths()
  addpath /auto/home/jmwick/EarthquakeSimulation/dc3dm-main/matlab;
end

function c = CC (v)
% Cell-centered from node-centered
  c = 0.5*(v(1:end-1) + v(2:end));
end

function h = draw_rect (xlo, ylo, xhi, yhi, filled, clr, varargin)
  if (filled)
    h = fill([xlo xhi xhi xlo], [ylo ylo yhi yhi], clr);
  else
    h = line([xlo xhi; xhi xhi; xhi xlo; xlo xlo]',...
             [ylo ylo; ylo yhi; yhi yhi; yhi ylo]');
    set(h, 'color', clr, varargin{:});
  end
end

function h = draw_rects_r (r, varargin)
  h = [];
  for (i = 1:size(r, 2))
    h1 = draw_rect(r(1,i), r(2,i), r(1,i) + r(3,i), r(2,i) + r(4,i),...
                   varargin{:});
    h = [h; h1(:)];
  end
end
