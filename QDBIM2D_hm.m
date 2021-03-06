% QDBIM2D_hm
% basically QDBIM2D but with hmmvp stress kernel

function varargout = exmb (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% --------------------- Public -----------------------------

% sets up hmatrix
function r = build()
  addpaths();

  % properties
  c.lambdaZ=40e3;
  c.n = 400; n = c.n;
  c.dz = c.lambdaZ/c.n;
  c.tol = 1.0e-8;

  y3 = (0:c.n-1)'*c.dz;
  W = ones(c.n,1)*c.dz;

  % reference friction coefficient
  ss.fo = 0.6*ones(size(y3));

  % frictional parameters
  % a is 'direct effect' - strength parameter
  % b is 'evolution effect' - weakening parameter

  % all vs
  %ss.a = 1e-3*ones(size(y3));
  %ss.b = ss.a - 2.1e-4*ones(size(y3)); % make + for with vw region

  % with vw
  ss.a=1e-2+Ramp((y3-15e3)/3e3)*(0.025-0.01);
  ss.b=0.015*ones(size(y3));

  % effective normal stress (MPa)
  ss.sigma=50.0*ones(size(y3));

  % characteristic weakening distance (m)
  ss.L=8e-3*ones(size(y3));

  % plate rate (m/s)
  ss.Vpl=1e-9*ones(size(y3));

  % reference slip rate (m/s)
  ss.Vo=1e-6*ones(size(y3));

  rho = 2670;
  Vs = 3464;
  G = rho*Vs^2/1e6;
  c.G = G;

  % Radiation damping coefficient
  ss.eta = G./(2*Vs);

  % housekeeping
  c.command = 'compress';
  c.write_hmat_filename = getFname(c);
  c.write_hd_filename = [c.write_hmat_filename '-hm'];
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.err_method  = 'mrem-fro';
  c.allow_overwrite = 1;

  c.greens_fn = 'okadaS12';

  % edges
  x = zeros(1,n);
  y = zeros(1,n);
  z = linspace(0,c.lambdaZ-c.dz,n);
  c.Y = [x; y; z];

  % centers
  z_c = z+ (c.dz/2);
  c.X = [x; y; z_c];

  kvf('Write', c.kvf, c, 4);

  disp('run this in a shell:')
  cmd = ['     ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)

  r.c = c;
  r.ss = ss;

end

% runs numerical solution
function solve(r)

  r.ss.dgf = 4;

  Y0 = zeros(r.c.n*r.ss.dgf,1);
  Y0(1:r.ss.dgf:end) = zeros(r.c.n,1);
  Y0(2:r.ss.dgf:end) = max(r.ss.a).*r.ss.sigma.*asinh(r.ss.Vpl./r.ss.Vo/2.*exp((r.ss.fo+r.ss.b.*log(r.ss.Vo./r.ss.Vpl))./max(r.ss.a))) + r.ss.eta.*r.ss.Vpl;
  Y0(3:r.ss.dgf:end)=r.ss.a./r.ss.b.*log(2*r.ss.Vo./r.ss.Vpl.*sinh((Y0(2:r.ss.dgf:end)-r.ss.eta.*r.ss.Vpl)./r.ss.a./r.ss.sigma))-r.ss.fo./r.ss.b;
  Y0(4:r.ss.dgf:end)=log(r.ss.Vpl./r.ss.Vo);

  % load kernel as hm
  hm = hmmvp('init', r.c.write_hmat_filename);

  yp=@(t,y) DieterichRuinaRegAging(t,y,r,hm);

  tic % val had err at 1e-8
  options=odeset('Refine',1,'RelTol',1e-6,'InitialStep',1e-5);
  [t,Y]=ode45(yp,[0 1e10],Y0,options); % 500*3.15e7
  disp('Done solving');
  toc

  % figures
  V = Y(:,4:r.ss.dgf:end)';
  disp(size(Y));
  disp(size(V));
  disp(std(V(:)));

  clf;
  imagesc(V); colorbar;
  title('QDBIM2D slip velocity')
  xlabel('time steps')
  ylabel('depth')
  saveas(gcf, 'figures/QDBIM2D_v.png')

end

function plot_hm(r)
  addpaths();

  hm_fname = r.c.write_hmat_filename;
  hm = hmmvp('init', hm_fname, 4);
  m = hmmvp('getm', hm);
  n = hmmvp('getn', hm);
  rs = (1:1:m); cs = (1:1:n);
  fullM = hmmvp('extract', hm, rs, cs);

  clf;
  imagesc(fullM); title('h-matrix'); colorbar;
  saveas(gcf, 'figures/qdbim2dhm_fullHM.png')

end

function benchmark(r)
  addpaths();

  % do hmmvp
  hm_fname = r.c.write_hmat_filename;
  hm = hmmvp('init', hm_fname, 32);
  m = hmmvp('getm', hm);
  n = hmmvp('getn', hm);

  x = zeros(m,1);
  start = int8((m/2)-(m/4));
  fin = int8((m/2)+(m/4));
  x(start:fin)=1;

  disp('hm mvp:');
  tic
  y_hm = hmmvp('mvp', hm, x);
  toc

  % check against sparse matrix
  s12h=@(x2,x3,y2,y3,W, G) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-W)./((x2-y2).^2+(x3-y3-W).^2)-(x3+y3+W)./((x2-y2).^2+(x3+y3+W).^2) ...
    )/2/pi;

  rho = 2670;
  Vs = 3464;
  G = rho*Vs^2/1e6;
  y3=(0:m-1)'*r.c.dz;

  K=zeros(m,m);
  for k=1:m
    K(:,k)=s12h(0,y3+r.c.dz/2,0,y3(k),r.c.dz, G);
  end

  disp('dense mvp:');
  tic
  y_d = K*x;
  toc

  %err metric.
  rs = (1:1:m); cs = (1:1:n);
  fullhm = hmmvp('extract', hm, rs, cs);
  diff = K-fullhm;
  frobD = norm(diff, 'fro');
  frobK = norm(K, 'fro');

  r.c.tol
  r.c.n
  err = frobD / frobK ;
  disp('err:')
  disp(err)

  % figure
  clf;
  subplot(211); plot(y_hm); title('hm mvp');
  subplot(212); plot(y_d); title('dense mvp');
  saveas(gcf, 'figures/qdbim2dhm_benchmark.png')

end

% -------------------- Private -----------------------------

function addpaths()
  addpath('../hmmvp-okada/matlab');
end

function f = getFname(p)
  f = sprintf('./tmp/QDBIM_tol%f_lz%d_n%d', p.tol, p.lambdaZ, p.n);
end

% boxcar function
function bc = BC(x)
  boxc=@(x) (x+0.5>=0)-(x-0.5>=0);
  bc = boxc(x);
end

% Heaviside function
function hs = HS(x)
  hs = 0+x>=0;
end

% ramp funciton
function r = Ramp(x)
  r = x.*BC(x-1/2)+HS(x-1);
end
