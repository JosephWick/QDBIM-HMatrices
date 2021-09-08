% ViscoelasticCycles.m

function varargout = vc (varargin)

  [varargout{1:nargout}] = feval(varargin{:});

end

% --------------------- Public ---------------------------

% set up h-matrices
% there will be 9 kernels
function r = build()
  addpaths();

  % we'll have 9 kernels:
  % - 1 for fault-fault interaction (s12)
  % - 4 for shear-shear interaction (1212, 1213, 1312, 1313)
  % - 4 for fault-shear interaction (1212, 1213, 1312, 1313))
  % - 2 for shear-fault interaction (s12, s13)

  % problem specifications
  % 200km long (L), 200km deep (W)
  % transition depth is 40 km
  % fault goes 40km deep

  % x1 (x) is in/out of page
  % x2 (y) is left/right of page
  % x3 (z) is up/down

  % ---     general params      ---
  probL = 200e3;
  probW = 200e3;
  lambdaZ = 35e3; % fault depth
  transition = 40e3; %where shear zone starts

  ss.Ny = 50; %num elems for shear mesh
  ss.Nz = 50;

  tag = string(probL/1000);

  % fault mesh
  ss.M = 400;
  ss.dz = lambdaZ/ss.M;
  % fault patch edges (top left)
  faultX = zeros(1,ss.M);
  faultY = zeros(1,ss.M);
  faultZ = linspace(0,lambdaZ-ss.dz, ss.M);
  %tops of fault patches
  ss.y3f = faultZ;
  % fault patch centers
  faultX_c = faultX;
  faultY_c = faultY;
  faultZ_c = faultZ+ss.dz/2;

  % shear mesh
  eps = 1e-12;

  % shear patch edges
  nc = (-ss.Nz/2:ss.Nz/2);
  shearZhat = transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*transition;
  shearYhat = tan(nc*pi/(2.5*max(nc)))*32e3;
  shearX = zeros(1,ss.Ny*ss.Nz);

  % shear patch centers
  shearX_c = shearX;
  shearY_chat = zeros(1,ss.Ny);
  shearZ_chat = zeros(1,ss.Nz);
  for idx=(1:length(shearZhat)-1)
    shearZ_chat(idx) = shearZhat(idx) + (shearZhat(idx+1) - shearZhat(idx))/2;
    shearY_chat(idx) = shearYhat(idx) + (shearYhat(idx+1) - shearYhat(idx))/2;
  end

  % grid and flatten
  shearZhat(end)=[]; shearYhat(end)=[];
  [shearZ shearY] = ndgrid(shearZhat, shearYhat);
  shearY = shearY(:)';
  shearZ = shearZ(:)';

  [shearZ_c shearY_c] = ndgrid(shearZ_chat, shearY_chat);
  shearZ_c = shearZ_c(:)';
  shearY_c = shearY_c(:)';

  c.X = [shearX_c; shearY_c; shearZ_c];
  c.Y = [shearX; shearY; shearZ];

  y2 = probL/2;
  ss.Nx = ss.Nz;

  ss.polesz = transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*transition;
  % center of shear zone (x3)
  ss.x3c=(ss.polesz(2:end)+ss.polesz(1:end-1))/2;
  W=ss.polesz(2:end)-ss.polesz(1:end-1);
  ss.polesxc=(2*y2/1e3:(2*y2)/(26e3):2*y2/1e3)'*1e3;
  edges= floor(ss.Nx-length(ss.polesxc)+1)/2;
  ss.polesxl=flipud(min(ss.polesxc)-tan((0:edges)'*pi/(2.2*(edges)+eps))*transition);
  ss.polesxr=max(ss.polesxc)+tan((0:edges)'*pi/(2.2*(edges)+eps))*transition;
  ss.polesx=[ss.polesxl(1:end-1);ss.polesxc;ss.polesxr(2:end)];
  % center of shear zone (x2)
  ss.x2c=(ss.polesx(2:end)+ss.polesx(1:end-1))/2;

  % combo mesh
  comboX = [faultX shearX];
  comboY = [faultY shearY];
  comboZ = [faultZ shearZ];

  comboX_c = [faultX_c shearX_c];
  comboY_c = [faultY_c shearY_c];
  comboZ_c = [faultZ_c shearZ_c];

  c.eta=3;
  % ---       kvf params      ---
  c.command = 'compress';
  c.lambdaZ = lambdaZ; ss.lambdaZ = lambdaZ;
  c.dz = ss.dz;
  c.tol = 1e-10;
  c.G = 30e3;
  c.command = 'compress';
  c.allow_overwrite = 1;
  c.err_method = 'mrem-fro';

  % ---       s12 kernel for fault-fault interaction      ---
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/VC_ff-s12_200';
  c.write_hd_filename = './tmp/VC_ff-s12-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.W = lambdaZ;
  c.L = 0;
  ss.Wf = 0;
  c.Y = [faultX; faultY; faultZ];
  c.X = [faultX_c; faultY_c; faultZ_c];

  kvf('Write', c.kvf, c, 4);
  disp('run these in a shell:')
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.s12 = c.write_hmat_filename;

  % ---       shear1212 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_ss-shear1212_200';
  c.write_hd_filename = './tmp/VC_ss-shear1212-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.transition = transition;
  c.L = probL;
  c.W = probW;
  c.Y = [shearX; shearY; shearZ];
  c.X = [shearX_c; shearY_c; shearZ_c];

  clf;
  imagesc(c.X); colorbar;
  saveas(gcf, 'figures/VCX.png')

  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1212 = c.write_hmat_filename;
  r.c=c;
  % ---       Shear 1213 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1213';
  c.write_hmat_filename = './tmp/VC_ss-shear1213_200';
  c.write_hd_filename = './tmp/VC_ss-shear1213-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  % geometry is the same for all shear-shear interaction
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1213 = c.write_hmat_filename;

  % ---     Shear 1312 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/VC_ss-shear1312_200';
  c.write_hd_filename = './tmp/VC_ss-shear1312-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1312 = c.write_hmat_filename;

  % ---     shear 1313 kernel for shear-shear interaction ---
  c.greens_fn = 'shear1313';
  c.write_hmat_filename = './tmp/VC_ss-shear1313_200';
  c.write_hd_filename = './tmp/VC_ss-shear1313-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.ss1313 = c.write_hmat_filename;

  % ---     shear 1212 kernels for fault-shear interaction ---
  c.greens_fn = 'shear1212';
  c.write_hmat_filename = './tmp/VC_fs-shear1212_200';
  c.write_hd_filename = './tmp/VC_fs-shear1212-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];

  c.X = [comboX_c; comboY_c; comboZ_c];
  c.Y = [comboX; comboY; comboZ];

  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1212 = c.write_hmat_filename;

  % ---       shear 1312 kernel for fault-shear interaction ---
  c.greens_fn = 'shear1312';
  c.write_hmat_filename = './tmp/VC_fs-shear1312_200';
  c.write_hd_filename = './tmp/VC_fs-shear1312-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.fs1312 = c.write_hmat_filename;

  % ---       s12 kernel for shear-fault interaction ---
  c.greens_fn = 'okadaS12';
  c.write_hmat_filename = './tmp/VC_sf-12_200';
  c.write_hd_filename = './tmp/VC_sf-s12-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.sf12 = c.write_hmat_filename;

  % ---     s13 kernel for shear-fault interaction ---
  c.greens_fn = 'okadaS13';
  c.write_hmat_filename = './tmp/VC_sf-s13_200';
  c.write_hd_filename = './tmp/VC_sf-s13-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 32);
  cmd = ['    ../hmmvp-okada/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd);
  r.sf13 = c.write_hmat_filename;

  r.ss = ss;

end

% runs numerical solution
function out = solve(r)
  addpaths();

  % rigidity (MPa)
  G=30e3;

  % - frictional parameters -

  % Confining pressure (MPa) and Temperature (K)
  k = 3.138; % thermal conductivity (W/m/K)
  Cp = 1171 ; % specific heat (J/kg/K)
  Rm = 3330 ; % mantle density (kg/m^3)

  Pconf       = Rm*9.8*r.ss.x3c/1e6;  % Shear zones
  Pconf_fault = Rm*9.8*(r.ss.y3f+r.ss.dz); % Faults

  Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
  Age_plate = 2e15; % seconds
  r.ss.Tprof  = 300+1380*erf(r.ss.x3c/(sqrt(4* Kappa * Age_plate)));  % Kelvin

  % default friction properties (velocity-weakening friction)
  % effective confining pressure (MPa)
  r.ss.sigmab = 1000;

  % frictional parameters
  r.ss.a = 1e-3*ones(size(r.ss.y3f));
  r.ss.b = r.ss.a+2.1e-4*ones(size(r.ss.y3f));

  r.ss.mu0 = 0.2*ones(size(r.ss.y3f));
  % characteristic weakening distance (m)
  r.ss.L = 0.012*ones(size(r.ss.y3f));
  % plate velocity (m/s)
  r.ss.V_plate = 1e-9*ones(size(r.ss.y3f));
  % reference slip rate (m/s)
  r.ss.Vo = 1e-6*ones(size(r.ss.y3f));
  % shear wave speed (m/s)
  r.ss.Vs = 3e3*ones(size(r.ss.y3f));

  % Velocity-strengthening at edges
  top    = floor(5e3/(r.ss.lambdaZ/r.ss.M));
  bottom = ceil(15e3/(r.ss.lambdaZ/r.ss.M));
  r.ss.b(1:top)      = r.ss.a(1:top)-2.1e-4;
  r.ss.b(bottom:end) = r.ss.a(bottom:end)-2.1e-4;

  % - Rheology -
  r.ss.e12p_plate = 1e-14; %1e-14*ones(length(ss.x2c)*length(ss.x3c),1);
  r.ss.e13p_plate = 0.0; %zeros(length(ss.x2c)*length(ss.x3c),1);
  % Rheological Parameters
  % Reference Strain Rate (for stress in MPa)
  r.ss.Adif = 1e6*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.Adis = 90 *ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Power-Law Exponent
  r.ss.n = 3.5*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Activation Energy Wet Oliving (J/mol)
  r.ss.Qdif = 335e3*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.Qdis = 480e3*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Activation Volume (m^3/mol)
  r.ss.Voldif = 4e-6*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.Voldis = 11e-6*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Grain size (m)
  r.ss.d    = 1e-2*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.pexp = 3*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Water fugacity (H/10^6 Si)
  r.ss.COH = 1000*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  r.ss.R   = 1.2*ones(length(r.ss.x3c)*length(r.ss.x2c),1);
  % Pressure (Pa)
  r.ss.P = repmat(1e6*Pconf',length(r.ss.x2c),1);
  r.ss.P = reshape(r.ss.P,[length(r.ss.x2c)*length(r.ss.x3c),1]);
  % Temperature (K)
  Te0 = repmat(r.ss.Tprof',length(r.ss.x2c),1);
  Te0 = reshape(Te0,[length(r.ss.x2c)*length(r.ss.x3c),1]);
  % Coefficients for dislocation and diffusion creep
  r.ss.Const_dis = r.ss.Adis.*exp(-(r.ss.Qdis+r.ss.P.*r.ss.Voldis)./(8.314.*Te0)).*r.ss.COH.^(r.ss.R);
  r.ss.Const_diff=r.ss.Adif.*exp(-(r.ss.Qdif+r.ss.P.*r.ss.Voldif)./(8.314.*Te0)).*r.ss.COH.^(r.ss.R).*r.ss.d.^(-r.ss.pexp);
  % Strengh profile
  s120 = (r.ss.e12p_plate./r.ss.Const_dis).^(1./r.ss.n);
  s130 = zeros(size(s120));
  e120 = zeros(size(s120));
  e130 = zeros(size(s120));
  % Fault Strength
  r.ss.strength = r.ss.sigmab.*(r.ss.mu0+(r.ss.a-r.ss.b).*log(r.ss.V_plate./r.ss.Vo))+G*r.ss.V_plate./(2*r.ss.Vs);

  % - Numerical Solution -

  % state parameters
  r.ss.dgfF=3;
  r.ss.dgfS=4;
  %% Initialize State Vector
  Y0=zeros(r.ss.M*r.ss.dgfF+length(r.ss.x2c)*length(r.ss.x3c)*r.ss.dgfS,1);

  % Fault patches
  Y0(1:r.ss.dgfF:r.ss.M*r.ss.dgfF)=zeros(size(r.ss.y3f));
  Y0(2:r.ss.dgfF:r.ss.M*r.ss.dgfF)=r.ss.strength;
  Y0(3:r.ss.dgfF:r.ss.M*r.ss.dgfF)=log(r.ss.Vo./r.ss.V_plate);

  % Shear zones
  Y0(r.ss.M*r.ss.dgfF+1:r.ss.dgfS:end)=s120;
  Y0(r.ss.M*r.ss.dgfF+2:r.ss.dgfS:end)=s130;
  Y0(r.ss.M*r.ss.dgfF+3:r.ss.dgfS:end)=e120;
  Y0(r.ss.M*r.ss.dgfF+4:r.ss.dgfS:end)=e130;

  % get the h-matrices
  hm.s12    = hmmvp('init', r.s12,     4);
  hm.ss1212 = hmmvp('init', r.ss1212, 32);
  hm.ss1213 = hmmvp('init', r.ss1213, 32);
  hm.ss1312 = hmmvp('init', r.ss1312, 32);
  hm.ss1313 = hmmvp('init', r.ss1313, 32);
  hm.fs1212 = hmmvp('init', r.fs1212, 32);
  hm.fs1312 = hmmvp('init', r.fs1312, 32);
  hm.sf12   = hmmvp('init', r.sf12,   32);
  hm.sf13   = hmmvp('init', r.sf13,   32);

  % initialize the function handle with
  % set constitutive parameters
  yp=@(t,y) odeViscoelastic(t,y,r.ss, hm);
  tic
  % Solve the system
  options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6);
  [t,Y]=ode45(yp,[0 1e7],Y0,options); %1e10
  toc
  % Compute the instantaneous derivative
  Yp=zeros(length(t)-1,size(Y,2));
  for k=1:length(t)-1
    Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
  end

  % - figures -
  % Strain rate at center
  Ep=sqrt(Yp(:,r.ss.M*r.ss.dgfF+floor(length(r.ss.x2c)/2)*r.ss.dgfS+3:r.ss.dgfS*length(r.ss.x2c):end)'.^2 +...
        Yp(:,r.ss.M*r.ss.dgfF+floor(length(r.ss.x2c)/2)*r.ss.dgfS+4:r.ss.dgfS*length(r.ss.x2c):end)'.^2);

  Epall = sqrt( Yp(:,r.ss.M*r.ss.dgfF+3:r.ss.dgfS:end)'.^2 +...
                Yp(:,r.ss.M*r.ss.dgfF+4:r.ss.dgfS:end)'.^2);

  % Velocity
  V=Yp(:,1:r.ss.dgfF:r.ss.M*r.ss.dgfF);

  % make movie
  movie=true;
  if movie
    clf;
    fig = figure;
    fname = 'figures/strain.gif';
    for idx = 1:size(Epall, 2)
      oneE = Epall(:,idx);
      oneEsq = reshape(oneE, [r.ss.Ny, r.ss.Nz]);
      imagesc(oneEsq); colorbar;
      title(idx)
      drawnow
      frame = getframe(fig);
      im{idx} = frame2im(frame);

      [A,map] = rgb2ind(im{idx},256);
      if idx==1
        imwrite(A,map,fname,'gif','LoopCount',Inf,'DelayTime',0.1);
      else
        imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',0.1);
      end

    end

    % velocity movie
    clf;
    fig = figure;
    fname='figures/faultV.gif';
    for idx = 1:size(V,1)
      oneV = V(idx,:)';
      imagesc(oneV); colorbar;
      title(idx)
      drawnow
      frame = getframe(fig);
      im{idx} = frame2im(frame);

      [A, map] = rgb2ind(im{idx}, 256);
      if idx==1
        imwrite(A,map,fname,'gif','LoopCount',Inf,'DelayTime',0.1);
      else
        imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',0.1);
      end
    end
  end

  % Maximum Velocity
  Vmax=zeros(length(t)-1,1);
  for ts=1:length(t)-1
      Vmax(ts)=max(V(ts,:));
  end

  %%
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %                    Function of Time                   %
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  figure(1);clf;set(gcf,'name','Time Evolution')
  f1=subplot(5,1,1);cla;
  pcolor(t(1:end-1)/3.15e7,r.ss.y3f/1e3,log10(V')), shading flat
  set(gca,'YDir','reverse');

  h=colorbar('Location','NorthOutside');
  caxis([min(min(log10(V))) max(max(log10(V)))+0.1]);
  colormap(f1,parula);
  title(h,'Slip Rate West (m/s)')
  xlabel('Time (yr)')
  ylabel('Depth (km)');

  f2=subplot(5,1,2);cla;
  pcolor(t(1:end-1)/3.15e7,r.ss.y3f/1e3,log10(V')), shading flat
  set(gca,'YDir','reverse');

  h=colorbar('Location','NorthOutside');
  caxis([min(min(log10(V))) max(max(log10(V)))+0.1]);
  colormap(f2,parula);
  title(h,'Slip Rate East (m/s)')
  xlabel('Time (yr)')
  ylabel('Depth (km)');

  f3=subplot(5,1,3);cla;
  %pcolor(t(1:end-1)/3.15e7,r.ss.x3c/1e3,log10(Ep)), shading flat
  set(gca,'YDir','reverse');

  %caxis([log10(min(min(Ep))) log10(max(max(Ep)))]);
  h1=colorbar('Location','NorthOutside');
  colormap(f3,hot);
  title(h1,'Strain Rate (1/s)')
  xlabel('Time (Yr)')
  ylabel('Depth (km)');

  subplot(5,1,4);cla;
  plot(t(1:end-1)/3.15e7,log10(Vmax))
  xlabel('Time (Yr)')
  ylabel('Velocity (m/s) log10')
  title('Maximum slip rates on faults')

  subplot(5,1,5);cla;
  plot(t(1:end-1)/3.15e7,log10(V(:,floor((top+bottom)/2))))
  xlabel('Time (Yr)')
  ylabel('Velocity (m/s) log10')
  title('Time series at center of seismogenic zones')

  % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %                Function of Time Steps                 %
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  f1=figure(2);clf;set(gcf,'name','Time Step Evolution')

  subplot(5,1,1);cla;
  pcolor(1:length(t)-1,r.ss.y3f/1e3,log10(V')), shading flat
  set(gca,'YDir','reverse');

  h=colorbar('Location','NorthOutside');
  caxis([min(min(log10(V))) max(max(log10(V)))+0.1]);
  colormap(f1,parula);
  title(h,'Slip Rate West (m/s)')
  xlabel('Time Steps')
  ylabel('Depth (km)');

  out.t = t;
  out.V = V;
  out.E = Ep;
  out.Eall = Epall;

  f3=subplot(5,1,3);cla;
  pcolor(1:length(t)-1, r.ss.x3c(1:end)/1e3, log10(Ep)), shading flat
  set(gca,'YDir','reverse');

  %caxis([log10(min(min(Ep))) log10(max(max(Ep)))]);
  h1=colorbar('Location','NorthOutside');
  colormap(f3,hot);
  title(h1,'Strain Rate (1/s)')
  xlabel('Time Steps')
  ylabel('Depth (km)');

  subplot(5,1,4);cla;
  plot(1:length(t)-1,log10(Vmax))
  xlabel('Time Steps')
  ylabel('Velocity (m/s) log10')
  title('Maximum slip rates on faults')

  subplot(5,1,5);cla;
  plot(1:length(t)-1,log10(V(:,floor((top+bottom)/2))))
  xlabel('Time Steps')
  ylabel('Velocity (m/s) log10')
  title('Time series at center of seismogenic zones')

  saveas(f1, 'figures/VC_f1.png')
  saveas(f2, 'figures/VC_f2.png')
  saveas(f3, 'figures/VC_f3.png')

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
