% compareKerns.m

function varargout = ck (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% ------------------------- Public ------------------------------

% load h-matrices
function hm = load(r)
  addpaths();

  hm.s12    = hmmvp('init', r.s12,     4);
  hm.ss1212 = hmmvp('init', r.ss1212, 32);
  hm.ss1213 = hmmvp('init', r.ss1213, 32);
  hm.ss1312 = hmmvp('init', r.ss1312, 32);
  hm.ss1313 = hmmvp('init', r.ss1313, 32);
  hm.fs1212 = hmmvp('init', r.fs1212, 32);
  hm.fs1312 = hmmvp('init', r.fs1312, 32);
  hm.sf12   = hmmvp('init', r.sf12,   32);
  hm.sf13   = hmmvp('init', r.sf13,   32);

end

% compShearKerns()
% - builds dense kernels and compares them (subtraction) to hm
function compShearKerns(hm, r)
  addpaths();
  G = 3e30;

  r.ss.k1212 = zeros(length(r.ss.shearY_chat)*length(r.ss.shearZ_chat));
  r.ss.k1213 = zeros(length(r.ss.shearY_chat)*length(r.ss.shearZ_chat));
  r.ss.k1312 = zeros(length(r.ss.shearY_chat)*length(r.ss.shearZ_chat));
  r.ss.k1313 = zeros(length(r.ss.shearY_chat)*length(r.ss.shearZ_chat));

  disp('begining kernels')

  L = r.ss.shearYhat(2:end) - r.ss.shearYhat(1:end-1);
  W = r.ss.shearZhat(2:end) - r.ss.shearZhat(1:end-1);
  L(length(L)+1) = L(1);
  W(length(W)+1) = W(1);

  % fields from shear zones
  for ky=1:length(r.ss.shearY_chat)
    for kz=1:length(r.ss.shearZ_chat)
      r.ss.k1212(:,(kz-1)*r.ss.Ny+ky) = s1212(G, r.ss.shearZ_c(kz)+r.ss.transition, L(ky), W(kz), r.ss.shearY_c'-r.ss.shearYhat(ky)', r.ss.shearZ_c');
      r.ss.k1213(:,(kz-1)*r.ss.Ny+ky) = s1213(G, r.ss.shearZ_c(kz)+r.ss.transition, L(ky), W(kz), r.ss.shearY_c'-r.ss.shearYhat(ky)', r.ss.shearZ_c');
      r.ss.k1312(:,(kz-1)*r.ss.Ny+ky) = s1312(G, r.ss.shearZ_c(kz)+r.ss.transition, L(ky), W(kz), r.ss.shearY_c'-r.ss.shearYhat(ky)', r.ss.shearZ_c');
      r.ss.k1313(:,(kz-1)*r.ss.Ny+ky) = s1313(G, r.ss.shearZ_c(kz)+r.ss.transition, L(ky), W(kz), r.ss.shearY_c'-r.ss.shearYhat(ky)', r.ss.shearZ_c');
    end
  end
  disp('kernels created')

  m = hmmvp('getm', hm.s1212);
  n = hmmvp('getn', hm.s1212);

  disp('beginning extraction')

  rs = (1:1:m); cs = (1:1:n);
  s1212HM = hmmvp('extract', hm.s1212, rs, cs);
  s1213HM = hmmvp('extract', hm.s1213, rs, cs);
  s1312HM = hmmvp('extract', hm.s1312, rs, cs);
  s1313HM = hmmvp('extract', hm.s1313, rs, cs);

  disp('hms extracted')
  disp('begin subtraction')

  s1212Diff = s1212HM - r.ss.k1212;
  s1213Diff = s1213HM - r.ss.k1213;
  s1312Diff = s1312HM - r.ss.k1312;
  s1313Diff = s1313HM - r.ss.k1313;

  disp('subtraction complete')
  disp('begin figurse')

  clf;
  imagesc(s1212Diff); colorbar;
  saveas(gcf, 'figures/diff_s1212')

  clf;
  imagesc(s1213Diff); colorbar;
  saveas(gcf, 'figures/diff_s1213')

  clf;
  imagesc(s1312Diff); colorbar;
  saveas(gcf, 'figures/diff_s1312');

  clf;
  imagesc(s1313Diff); colorbar;
  saveas(gcf, 'figures/diff_s1313')

  disp('figures created')

end

% ------------------------ Private -------------------------------
function addpaths()
  addpath('../hmmvp-okada/matlab')
end

% - Kernels
function k = s1312(G, D,L,W,x2,x3)
    k = G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    +log((x2-L/2).^2+(x3+D+W).^2) - log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    -log((x2-L/2).^2+(x3+D).^2) + log((x2+L/2).^2+(x3+D).^2));
end

function k = s1212(G, D,L,W,x2,x3)
    k = G/pi*( ...
    atan((x3-D)./(x2+L/2))-atan((x3-D)./(x2-L/2)) ...
    +atan((x3-D-W)./(x2-L/2))-atan((x3-D-W)./(x2+L/2)) ...
    -atan((x3+D+W)./(x2-L/2))-atan((x3+D)./(x2+L/2)) ...
    +atan((x3+D)./(x2-L/2))+atan((x3+D+W)./(x2+L/2)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);
end

function k = s1213(G, D,L,W,x2,x3)
    k = G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    -log((x2-L/2).^2+(x3+D+W).^2) + log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    +log((x2-L/2).^2+(x3+D).^2) - log((x2+L/2).^2+(x3+D).^2));
end

function k = s1313(G, D,L,W,x2,x3)
    k = G/pi*( ...
    atan((x2+L/2)./(x3-D))  -atan((x2-L/2)./(x3-D)) ...
    -atan((x2+L/2)./(x3-D-W))+atan((x2-L/2)./(x3-D-W)) ...
    +atan((x2+L/2)./(x3+D))  -atan((x2-L/2)./(x3+D)) ...
    -atan((x2+L/2)./(x3+D+W))+atan((x2-L/2)./(x3+D+W)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);
end

function bc = boxc(x) 
  bc = (x+0.5>=0)-(x-0.5>=0);
end
