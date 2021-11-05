% BP1_visco.m

function varargout = bp1v (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ---------------------------- Public ----------------------------

% set up h-matrices/kernels
% there will be 9 Kernels

function r = build()
  addpaths();

  % 9 kernels:
  % - 1 for fault-fault interaction (S12)
  % - 4 for shear-shear interaction (1212, 1213, 1312, 1313)
  % - 4 for fault-shear interaction (1212, 1213, 1312, 1313))
  % - 2 for shear-fault interaction (s12, s13)

  % problem specifications
  % - fault goes 40 km deep
  %  - shear TBD 


end

% ---------------------------- Private ---------------------------

% add paths to hmmvp
function addpaths()
  addpath('../hmmvp-okada/matlab')
end

% boxcar function
function bc = BC(x)
  boxc=@(x) (x+0.5>=0)-(x-0.5>=0);
  bc = boxc(x);
end
