% ViscoelasticCycles.m

function varargout = vc (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% --------------------- Public ---------------------------

% set up h-matrices
% there will be 11 kernels
function r = build()
  addpaths();

  % we'll have 9 kernels:
  % - 1 for fault-fault interaction (s12)
  % - 4 for shear-shear interaction (1212, 1213, 1312, 1313)
  % - 2 for fault-shear interaction (1212, 1213)
  % - 2 for shear-fault interaction (1212, 1213)

  % fault mesh
  lambdaZ = 40e3; %(m)

  % visco mesh
  transition = 35e3;
  vL = 25e3;
  vW = 15e3;

  dz = 100;
  % fault has 400 elems
  % visco is 250x150

  % set up kvf general propts
  c.lambdaZ = lambdaZ;
  c.dz = dz;
  c.tol = 1e-3;
  c.G = 30e3;
  c.command = 'compress';
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  % - s12 kernel for fault-fault interaction; 1D -
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/VC_ff-s12';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.W = lambdaZ;
  c.L = 0;

  n = lambdaZ/dz;
  x = zeros(1,n);
  y = zeros(1,n) + 15000;
  z = linspace(c.dz,lambdaZ,n);
  c.X = [x;y;z];

  kvf('Write', c.kvf, c, 4);
  disp('run these in a shell:')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.s12 = c.write_hmat_filename;

  % - shear 1212 kernel for shear-shear interaction -
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_ss-shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.transition=transition;
  c.L = vL;
  c.W = vW;
  n = vL/dz;
  m = vW/dz;
  x = zeros(1,n*m);
  y = linspace(c.dz,vL,n);
  z = linspace(c.dz,vW,m);
  [Y Z] = ndgrid(y,z);
  c.X = [x;Y(:)';Z(:)'];

  kvf('Write', c.kvf, c, 4);
  cmd = ['    ...hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1212 = c.write_hmat_filename;

  % - shear 1213 kernel for shear-shear interaction -
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/VC_ss-shear1213';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is the same for all shear-shear interactions
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1213 = c.write_hmat_filename;

  % - shear 1312 kernel for shear-shear interaction -
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/VC_ss-shear1312';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp' c.kvf];
  disp(cmd)
  r.ss1312 = c.write_hmat_filename;

  % - shear 1313 kernel for shear-shear interaction -
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/VC_ss-shear1313';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp' c.kvf];
  disp(cmd)
  r.ss1313 = c.write_hmat_filename;

  % - shear 1212 kernel for fault-shear interaction -
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_fs-shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % update geometry
  fn = lambdaZ/dz;                %fault
  fx = zeros(1,fn);
  fy = zeros(1,fn) + 15000;
  fz = linspace(c.dz,lambdaZ,fn);
  sn = vL/dz;                     %shear
  sm = vW/dz;
  sx = zeros(1,sn*sm);
  sy = linspace(c.dz,vL,sn);
  sz = linspace(c.dz,vW,sm);
  [Ys Zs] = ndgrid(sy,sz);
  bx = [fx sx];                   %together
  by = [fy Ys(:)'];
  bz = [fz Zx(:)'];
  c.X = [bx;by;bz];
  % write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp' c.kvf];
  disp(cmd)
  r.fs1212 = c.write_hmat_filename;

  % - shear 1213 kernel for fault-shear interaction -
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/VC_fs-shear1213';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is same as fs1212
  % %write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp' c.kvf];
  disp(cmd)
  r.fs1213 = c.write_hmat_filename;

  % - shear 1212 kernel for shear-fault interaction -
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_sf-shear1212';
  c.kvf = [write_hmat_filename '.kvf'];
  % update geometry - need to adjust fault locs to be receiver
  fn = lambdaZ/dz;                %fault
  fx = zeros(1,fn);
  fy = zeros(1,fn) + 15000 + 0.5*dz;
  fz = linspace(c.dz,lambdaZ,fn);
  bx = [fx sx];                   %together
  by = [fy Ys];
  bz = [fz Zx];
  c.X = [bx;by;bz];
  % write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp' c.kvf];
  disp(cmd)
  r.sf1212 = c.write_hmat_filename;

  % - shear 1313 kernel for shear-fault interaction -
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/Vc_sf-shear1313';
  c.kvf = [write_hmat_filename '.kvf'];
  % geometry is same as sf1212
  % write
  kvf('Write', c.kvf, c, 4);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp' c.kvf];
  disp(cmd)
  r.sf1313 = c.write_hmat_filename;

end

% -------------------- Private ---------------------------

function addpaths()
  addpath('../hmmvp-okada/matlab')
end

function f = getFname(p)
  f = sprintf('./tmp/VC_%s_tol%f_lz%d_n%d', p.greens_fn, p.tol, p.lambdaZ, p.n);
end

% boxcar function
function bc = BC(x)
  boxc=@(x) (x+0.5>=0)-(x-0.5>=0);
  bc = boxc(x);
end
