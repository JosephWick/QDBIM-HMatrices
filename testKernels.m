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
  c.X = [x;y;z];

  % per kernel now

  %s12
  c.write_hmat_filename = './tmp/tk_s12';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'OkadaS12';
  kvf('Write', c.kvf, c, 4);

  %s13
  c.write_hmat_filename = './tmp/tk_s13';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'OkadaS13';
  kvf('Write', c.kvf, c, 4);

  %shear1212
  c.write_hmat_filename = './tmp/tk_shear1212';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'shear1212';
  kvf('Write', c.kvf, c, 4);

  %shear1213
  c.write_hmat_filename = './tmp/tk_shear1213';
  c.kvf = [c.write_hmat_filename 'kvf'];
  c.greens_fn = 'shear1213';
  kvf('Write', c.kvf, c, 4);

  %shear1312
  c.write_hmat_filename = './tmp/tk_shear1312';
  c.kvf = [c.write_hmat_filename 'kvf'];
  c.greens_fn = 'shear1312';
  kvf('Write', c.kvf, c, 4);

  %shear1313
  c.write_hmat_filename = './tmp/tk_shear1313';
  c.kvf = [c.write_hmat_filename '.kvf'];
  c.greens_fn = 'shear1313';
  kvf('Write', c.kvf, c, 4);


end


% ---------------- Private -------------------------

function addpaths()
  addpath('../hmmvp-okada/matlab');
end
