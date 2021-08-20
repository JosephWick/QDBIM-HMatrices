function [Yp]= odeViscoelastic(~,Y,ss)
% function OdefunViscoelastic describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state vector y is
%
%        /        s          \
%    y = |       tau         |
%        | log(theta Vo / L) |
%        |       ...         |
%        |       s12         |
%        |       s13         |
%        |       e12         |
%        \       e13         /
%
% Instead of integrating numerically the aging law
%
%    d theta / dt = 1 - V theta / L
%
% as is, we operate the following change of variable
%
%    phi = ln (theta Vo / L)
%
% Considering the following relationships
%
%    d phi = d theta / theta
%
%    d theta = d phi theta = d phi exp(phi) L / Vo
%
%    d theta / dt = d phi exp(phi) / dt L / Vo
%
%    d phi exp(phi) / dt L / Vo = 1 - V theta / L = 1 - V exp(phi) / Vo
%
% we obtain the evolution law for the new variable
%
%    d phi / dt = ( Vo exp(-phi) - V ) / L
%
%  The state vector is split with initial cells considering slip on a fault
%  and the latter portion considering strain in finite volumes with (s12,s13)
%  and (e12,e13) being the 2D antiplane stress and strain components.

G=30e3; % MPa

% Shear stress on faults
tauF = Y(2:ss.dgfF:ss.M*ss.dgfF);

% State variables
th   = Y(3:ss.dgfF:ss.M*ss.dgfF);

% Slip velocities
V = (2.*ss.Vs.*ss.a.*ss.sigmab./G).*...
     Lambert_W(G*ss.Vo./(2*ss.Vs.*ss.a.*ss.sigmab).*...
     exp((tauF-ss.mu0.*ss.sigmab-ss.sigmab.*ss.b.*th_W)./(ss.sigmab.*ss.a)));

% Shear stress in zones of distributed deformation
tau12=Y(2*ss.M*ss.dgfF+1:ss.dgfS:end);
tau13=Y(2*ss.M*ss.dgfF+2:ss.dgfS:end);
tau=sqrt(tau12.^2+tau13.^2);

% Dislocation strain rate
Aeff = ss.Const_dis.* (tau).^(ss.n-1);
e12p = tau12 .* Aeff;
e13p = tau13 .* Aeff;

% Initiate state derivative
Yp=zeros(size(Y));

% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                     West Fault                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Slip velocity
Yp(1:ss.dgfF:ss.M*ss.dgfF)=V;

% Shear stress rate on fault due to fault and shear zones
Yp(2:ss.dgfF:ss.M*ss.dgfF)=  hmmvp('mvp', hm.s12, (V-ss.V_plate)) + ...
                             hmmvp('mvp', hm.fs1212, (e12p-ss.e12p_plate)) + ...
                             hmmvp('mvp', hm.fs1312, (e13p-ss.e13p_plate));

% Rate of state
Yp(3:ss.dgfF:ss.M*ss.dgfF)=(ss.Vo.*exp(-th)-V)./ss.L;


% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                     Shear Zones                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Stress rate due to shear zones and fault slip velocity
Yp(2*ss.M*ss.dgfF+1 : ss.dgfS : end) =  hmmvp('mvp', hm.ss1212, (e12p-ss.e12p_plate)) + ...
                                        hmmvp('mvp', hm.ss1312, (e13p-ss.e13p_plate)) + ...
                                        hmmvp('mvp', hm.sf12, (V-ss.V_plate));

Yp(2*ss.M*ss.dgfF+2 : ss.dgfS : end) = hmmvp('mvp', hm.ss1213, (e12p-ss.e12p_plate)) + ...
                                       hmmvp('mvp', hm.ss1313, (e13p-ss.e13p_plate)) + ...
                                       hmmvp('mvp', hm.sf13, (V_W-ss.V_plate));

% Strain rate
Yp(2*ss.M*ss.dgfF+3 : ss.dgfS : end) = e12p;
Yp(2*ss.M*ss.dgfF+4 : ss.dgfS : end) = e13p;


end