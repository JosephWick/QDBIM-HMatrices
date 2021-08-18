function varargout = testKernels(varargin)

% - testing the cpp kernels I wrote for use in hmmvp against those
%   valere wrote in ViscoelasticCycles in matlab
% - less concerned with the setup and more if the math matches

  [varargout{1:nargout}] = feval(varargin{:});

end

% kernels to test:
% - GreensFnOkadaS12 (already tested but may was well here too)
% - GreensFnOkadaS13
% - GreensFnShear1212
% - GreensFnShear1213
% - GreensFnShear1312
% - GreensFnShear1313

% ----------------- Public -------------------------

% builds kvf's for all of the kernels
function b = build()

  addpaths();

  % general fields
  c.command = 'compress';
  c.err_method = 'mrem-fro';
  c.tol = 1e-3;
  c.allow_overwrite = 1;
  c.halfspace = 1;

  rho = 2670;
  Vs = 3464;
  G = rho*Vs^2/1e6;
  c.G = G;

  c.mu = 30e9;
  c.nu = 0.25;
  c.W = 10000;
  c.L = 10000;
  c.dz = 100;
  n = c.L/c.dz;

  x = zeros(1,n);
  y = linspace(c.dz, c.L, n);
  z = linspace(c.dz, c.W, n);
  [X Y] = ndgrid(x,y);
  c.X = [X(:)';Y(:)';z];

  % per kernel now

  %s12
  c.write_hmat_filename = './tmp/tk_s12';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'okadaS12';
  kvf('Write', c.kvf, c, 4);

  %s13
  c.write_hmat_filename = './tmp/tk_s13';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'okadaS13';
  kvf('Write', c.kvf, c, 4);

  %shear1212
  c.write_hmat_filename = './tmp/tk_shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'shear1212';
  kvf('Write', c.kvf, c, 4);

  %shear1213
  c.write_hmat_filename = './tmp/tk_shear1213';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'shear1213';
  kvf('Write', c.kvf, c, 4);

  %shear1312
  c.write_hmat_filename = './tmp/tk_shear1312';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'shear1312';
  kvf('Write', c.kvf, c, 4);

  %shear1313
  c.write_hmat_filename = './tmp/tk_shear1313';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'shear1313';
  kvf('Write', c.kvf, c, 4);

  disp('run these shell commands:')
  fprintf('   ../hmmvp-okada/bin/hmmvpbuild_omp ./tmp/tk_s12.kvf \n')
  fprintf('   ../hmmvp-okada/bin/hmmvpbuild_omp ./tmp/tk_s13.kvf \n')
  fprintf('   ../hmmvp-okada/bin/hmmvpbuild_omp ./tmp/tk_shear1212.kvf \n')
  fprintf('   ../hmmvp-okada/bin/hmmvpbuild_omp ./tmp/tk_shear1213.kvf \n')
  fprintf('   ../hmmvp-okada/bin/hmmvpbuild_omp ./tmp/tk_shear1312.kvf \n')
  fprintf('   ../hmmvp-okada/bin/hmmvpbuild_omp ./tmp/tk_shear1313.kvf \n')

  b=c;

end

function test(b)
  addpaths();

  %s12
  hm_fname = './tmp/tk_s12';
  hm = hmmvp('init', hm_fname, 4);
  m = hmmvp('getm', hm);
  n = hmmvp('getn', hm);
  rs = (1:1:m); cs = (1:1:n);
  fullM = hmmvp('extract', hm, rs, cs);

  disp('S12:')
  for idx = 1:5
    i = randi([1,n], 1,1);
    j = randi([1,n], 1,1);

    hmVal = fullM(i,j);



    fprintf('hm: %f, d: %f\n', hmVal, dVal)
  end


end

% ---------------- Private -------------------------

function addpaths()
  addpath('../hmmvp-okada/matlab');
end

function s = s12h(x2, x3, y2, y3, Wf)
  rho = 2670;
  Vs = 3464;
  G = rho*Vs^2/1e6;

  s12=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

  s = s12(x2, x3, y2, y3, Wf);

end
