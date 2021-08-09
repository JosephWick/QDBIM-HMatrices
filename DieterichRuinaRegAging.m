function [yp] = DieterichRuinaRegAging(~,y,r)
% This file describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state
% vector y is
%
%        /        s          \
%    y = |       tau         |
%        | log(theta Vo / L) |
%        \    log(V / Vo)    /

% based on the regularized form of Dieterich-Ruina rate-and-state friction
% and using the aging law
%
% Velocity is determined by the balance of shear resistance and shear stress
% Resistance : tau = a sigma asinh( V/2Vo exp( (fo + b psi ) / a))
% Stress     : tau = tauo + f(z,t) - eta*(V-Vpl)
% Where we use radiation damping for approximating the inertial terms
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
%    d theta / dt = 1 - V theta / L
%
% as is, we operate the following change of variable
%
%    phi = ln (theta Vo / L)
%
% where we obtain the evolution law for the new variable
%
%    d phi / dt = ( Vo exp(-phi) - V ) / L
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

% State variable
th=y(3:r.ss.dgf:end);

% Slip rate
V = r.ss.Vo.* exp(y(4:r.ss.dgf:end));

% Initialize Time Derivative
yp=zeros(size(y));

% Slip
yp(1:r.ss.dgf:end)=V;

% State Variable
dth = (r.ss.Vo.*exp(-th)-V)./r.ss.L;
yp(3:r.ss.dgf:end)=dth;

% Slip Velocity
K = hmmvp('init', r.kvf, 4);

func = hmmvp('mvp', K, (V-r.ss.Vpl));
f1=2*r.ss.Vo./V.*exp(-(r.ss.fo+r.ss.b.*th)./r.ss.a);
f2=1./sqrt(1+f1.^2);

yp(4:r.ss.dgf:end) = (func - r.ss.b.*r.ss.sigma.*dth.*f2)./ ...
                    (r.ss.a.*r.ss.sigma.*f2 + r.ss.eta.*V);


% Evolution of shear stress
yp(2:r.ss.dgf:end)=func - r.ss.eta.*V.*yp(4:r.ss.dgf:end);



end
