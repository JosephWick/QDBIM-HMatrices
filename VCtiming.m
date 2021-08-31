function varargout = vct (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ------------------------- Public ----------------------------------

% time single MVP action using hmmvp and with dense matrices
function time(b)
  addpaths();

  % Slip velocities
  V = (2.*b.ss.Vs.*b.ss.a.*b.ss.sigmab./G).*...
       Lambert_W(G*b.ss.Vo./(2*b.ss.Vs.*b.ss.a.*b.ss.sigmab).* ...
       exp((tauF-b.ss.mu0.*b.ss.sigmab-b.ss.sigmab.*b.ss.b.*th)./ ...
       (b.ss.sigmab.*b.ss.a)));

  % slices for kernels
  % greater than ss.M
  gM = (b.ss.M+1:1: (b.ss.Nx*ss.Nz)+b.ss.M);
  % less than or equal to ss.M
  lM = (1:1:b.ss.M);

  fs1212 = hmmvp('init', b.fs1212, 32);
  m = hmmvp('getm', fs1212);)
  X = rand(m,1);

  % hmmvp mvp
  disp('hmmvp:')
  tic
  p = hmmvp('mvp', fs1212, X, lM, gM);
  toc

  % dense mvp
  fs1212_d = hmmvp('extract', b.fs1212, (1:1:m), (1:1:m));
  tic
  pd = fs1212_d(lm,gm)*X(lm);


end

% ----------------------- Private ----------------------------------
function addpaths()
  addpath('../hmmvp-okada/matlab')
end
