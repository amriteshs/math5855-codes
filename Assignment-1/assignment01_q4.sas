proc iml;

mu = {172, 104, 105, 94};
cov = {1000 -80 1100 275, -80 500 90 -90, 1100 90 1500 200, 275 -90 200 720};
cor = {1 -0.11317 0.89815 0.3241, -0.11317 1 0.1039 -0.15, 0.89815 0.1039 1 0.19245, 0.3241 -0.15 0.19245 1};

print 'Q4(i)';

mu_y = mu[1, ];
mu_x = mu[{2 3 4}, ];
x_x = {100, 100, 100};
sigma_yy = cov[1, 1];
sigma_yx = cov[1, {2 3 4}];
sigma_xy = sigma_yx`;
sigma_xx = cov[{2 3 4}, {2 3 4}];
sigma_xx_inv = inv(sigma_xx);

mu_est = mu_y + (sigma_yx * sigma_xx_inv * (x_x - mu_x));
cov_est = sigma_yy - (sigma_yx * sigma_xx_inv * sigma_xy);

beta = sigma_xx_inv * sigma_xy;
b = mu_y - (beta` * mu_x);

lin_apprx = sum((beta` || b) # {100 100 100 1});

print mu_est[label='Estimated Mean Vector for (X1) | (X2, X3, X4)'];
print cov_est[label='Estimated Covariance Matrix for (X1) | (X2, X3, X4)'];
print lin_apprx[label='Best Linear Approximation of Sales for (X1) | (X2, X3, X4)'];

print 'Q4(ii)';

r_yx = cor[1, {2 3 4}];
r_xy = r_yx`;
r_xx = cor[{2 3 4}, {2 3 4}];
r_xx_inv = inv(r_xx);

rho_sq_1234 = r_yx * r_xx_inv * r_xy;
rho_sq_12 = cor[1, 2] * cor[1, 2];

if rho_sq_1234 > 0.81 then
    corr_1234 = 'Strong';
else if rho_sq_1234 < 0.0625 then
    corr_1234 = 'Weak';
else
    corr_1234 = 'Moderate';
	
if rho_sq_12 > 0.81 then
    corr_12 = 'Strong';
else if rho_sq_12 < 0.0625 then
    corr_12 = 'Weak';
else
    corr_12 = 'Moderate';

print rho_sq_1234[label='Multiple Correlation Coefficient for (X1) | (X2, X3, X4)'];
print corr_1234[label='Nature of Correlation Between (X1) and (X2, X3, X4)'];
print corr_12[label='Nature of Correlation Between (X1) and (X2)'];

print 'Q4(iii)';

mu_y1 = mu[{1 2}, ];
mu_x1 = mu[{3 4}, ];
x_x1 = {100, 100};
sigma_yy1 = cov[{1 2}, {1 2}];
sigma_yx1 = cov[{1 2}, {3 4}];
sigma_xy1 = sigma_yx1`;
sigma_xx1 = cov[{3 4}, {3 4}];
sigma_xx_inv1 = inv(sigma_xx1);

mu_est1 = mu_y1 + (sigma_yx1 * sigma_xx_inv1 * (x_x1 - mu_x1));
cov_est1 = sigma_yy1 - (sigma_yx1 * sigma_xx_inv1 * sigma_xy1);

print mu_est1[label='Estimated Mean Vector for (X1, X2) | (X3, X4)'];
print cov_est1[label='Estimated Covariance Matrix for (X1, X2) | (X3, X4)'];

quit;
