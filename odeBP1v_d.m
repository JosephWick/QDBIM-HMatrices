function [Yp] = odeBp1v(~,Y,ss,)
% This file describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state
% vector y is
%
%        /        s          \
%    y = |       tau         |
%        | log(theta Vo / D_rs) |
%        \    log(V / Vo)    /

% based on the regularized form of Dieterich-Ruina rate-and-state friction
% and using the aging law
%
% Velocity is determined by the balance of shear resistance and shear stress
% Resistance : tau = a sigma asinh( V/2Vo exp( (fo + b psi ) / a))
% Stress     : tau = tauo + f(z,t) - eta*(V-Vpl)
% Where we use radiation damping to approximate the inertial terms
% as the quasi-dynamic approximation with f(z,t) = K(delta - Vpl*t)
%
% This is done by taking the time-derivative of both equations and equating
% the two.
%
%   a sigma  alpha/ sqrt(1 + alpha^2) ( 1/V dV/dt + b/a dPhi/dt)
%
%                   = K( V - Vpl) - eta dV/dt

% where
%      alpha = V / 2Vo exp( ( fo + b phi) / a)
%
%
%
% Note this form assumes no time variation in sigma, the frictional
% parameters, or the loading plate rate. Loading is done purely through
% back slip at plate rate Vpl.
%
% Instead of directly integrating numerically the aging law
%
%    d theta / dt = 1 - V theta / D_rs
%
% as is, we operate the following change of variable
%
%    phi = ln (theta Vo / D_rs)
%
% where we obtain the evolution law for the new variable
%
%    d phi / dt = ( Vo exp(-phi) - V ) / D_rs
%
% Given the regularized form of Dieterich-Ruina R+S we can express,
%    1 dV      K (V - Vpl) - b sigma dphi / dt Q
%    - --  =  ----------------------------------
%    V dt           a sigma Q  + eta V
%
% where
%
%                            1
%    Q = -----------------------------------------------
%        /                                             \(1/2)
%        |  1 + [2 Vo / V exp(-(fo + b phi) / a )]^2   |
%        \                                             /
%
% Note that d/dt log(V/Vo) = 1/V dV/dt, and is much more efficient to
% integrate than dV/dt alone

G = 30e3;

% slices for kernels
% greater than ss.M
em = hmmvp('getm', hm.ss1212);
gM = (ss.M+1:1:em);
% less than or equal to ss.M
lM = (1:1:ss.M);

% State variable
th = Y(3:ss.dgfF:ss.M*ss.dgfF);

% Slip rate
V = ss.Vo.* exp(-Y(4:ss.dgfF:ss.M*ss.dgfF));

% Initialize Time Derivative
Yp=zeros(size(Y));

% Slip
Yp(1:ss.dgfF:ss.M*ss.dgfF)=V;

% Shear stress in zones of distributed deformation
tau12 = Y(ss.M*ss.dgfF+1:ss.dgfS:end);
tau13 = Y(ss.M*ss.dgfF+2:ss.dgfS:end);
tau = sqrt(tau12.^2 + tau13.^2);

% Dislocation strain rate
Aeff = ss.Const_dis.* (tau).^(ss.n-1);
e12p = tau12 .* Aeff;
e13p = tau13 .* Aeff;

% ---       FAULT       ---
% Shear stress rate on fault due to fault and shear zones
%TODO
Yp(2:ss.dgfF:ss.M*ss.dgfF) = ss. + t2(1:ss.M) + t3(1:ss.M);

% Rate of state
%TODO
Yp(3:ss.dgfF:ss.M*ss.dgfF) = (ss.Vo.*exp(-th)-V)./ss.L;

% ---       SHEAR         ---
% Stress rate due to shear zones and fault slip velocity
Yp(ss.M*ss.dgfF+1:ss.dgfS:end) = ss.k1212*(e12p-ss.e12p_plate) + ...
  ss.k1312*(e13p-ss.e13p_plate);

Yp(ss.M*ss.dgfF+2:ss.dgfS:end) = ss.k1213*(e12p-ss.e12p_plate) + ...
  ss.k1313*(e13p-ss.e13p_plate);

% Strain rate
Yp(ss.M*ss.dgfF+3 : ss.dgfS : end) = e12p;
Yp(ss.M*ss.dgfF+4 : ss.dgfS : end) = e13p;


end
