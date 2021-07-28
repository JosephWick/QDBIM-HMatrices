%

addpath('../okada_wrapper')

source_depth = 3.0;
obs_depth = 3.0;
dislocation = [1.0, 0.0, 0.0];

mu = 30e9;
nu = 0.25;
alpha = 0.67;

dip = 90;

n = 100;

x = linspace(-0.5*n, 0.5*n, n);
y = linspace(-0.5*n, 0.5*n, n);
S_xy = zeros(n, n);
U_xy = zeros(n,n);
for i = 1:100
  for j = 1:100
    [sucess, u, grad_u] = DC3Dwrapper(alpha, [x(i), y(j), -obs_depth], ...
                                      source_depth, dip, ...
                                      [-0.7, 0.7], ...
                                      [20, 0.7], ...
                                      [1.0, 0.0, 0.0]);

    S_xy(i,j) = ( grad_u(1,2) + grad_u(2,1) ) * mu;
    %U_xy(i,j) = grad

  end
end

clf;
imagesc(S_xy); title('okada wrapper'); colorbar;
saveas(gcf, 'figures/test_ow.png');
