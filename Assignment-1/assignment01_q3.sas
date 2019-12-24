data iris; 
infile '/folders/myfolders/iris.dat'; 
input x1 x2 x3 x4;
run;

proc iml; 
use iris var {x1 x2 x3 x4};
read all var _num_ into x;

print 'Q3(i)';

n = nrow(x);
p = ncol(x);
dfchi = p * (p + 1) * (p + 2) / 6;

q = i(n) - (1 / n) * j(n, n, 1);
s = (1 / n) * x` * q * x; 
s_inv = inv(s);
g_matrix = q * x * s_inv * x` * q;
beta1hat = sum(g_matrix # g_matrix # g_matrix) / (n * n);
beta2hat = trace(g_matrix # g_matrix) / n;
kappa1 = n * beta1hat / 6;
kappa2 = (beta2hat - p * (p + 2)) / sqrt(8 * p * (p + 2) / n);
pvalskew = 1 - probchi(kappa1, dfchi);
pvalkurt = 2 * (1 - probnorm(abs(kappa2)));

if pvalskew > 0.05 & pvalkurt > 0.05 then
    hypothesis = 'Accepted';
else
    hypothesis = 'Rejected';

print pvalskew[label='P-value of Multivariate Skewness'];
print pvalkurt[label='P-value of Multivariate Kurtosis'];
print hypothesis[label='Hypothesis for Multivariate Normality'];

print 'Q3(ii)';

mu = x[:, ]`;
cov = ((x` - mu) * (x` - mu)`) / (n - 1);

print mu[label='Estimate of Mean Vector'];
print cov[label='Sample Covariance Matrix'];

print 'Q3(iii)';

mu_y = mu[{3 4}, ];
mu_x = mu[{1 2}, ];
x_x = {1.3, 1.5};
sigma_yy = cov[{3 4}, {3 4}];
sigma_yx = cov[{3 4}, {1 2}];
sigma_xy = sigma_yx`;
sigma_xx = cov[{1 2}, {1 2}];
sigma_xx_inv = inv(sigma_xx);

mu_est = mu_y + (sigma_yx * sigma_xx_inv * (x_x - mu_x));
cov_est = sigma_yy - (sigma_yx * sigma_xx_inv * sigma_xy);

print mu_est[label='Estimate of Mean Vector for (X3, X4) | (X1, X2)'];
print cov_est[label='Estimate of Covariance Matrix for (X3, X4) | (X1, X2)'];

print 'Q3(iv)';

cov_inv = inv(cov);
mu0 = {1.2, 0.6, 4.0, 1.6};

t_squared = n * (mu - mu0)` * cov_inv * (mu - mu0);
f_critical_val = p * (n - 1) * quantile('F', 0.95, p, n - p) / (n - p);

if t_squared <= f_critical_val then
    null_hypothesis = 'Accepted';
else
    null_hypothesis = 'Rejected';
	
print cov_inv[label='Inverse of Sample Covariance Matrix'];
print t_squared[label='Hotelling T^2'];
print f_critical_val[label='F Critical Value'];
print null_hypothesis[label='Null Hypothesis H0'];

print 'Q3(v)';

y_conv = {1 0 -1 0, 0 1 0 -1};
p_new = nrow(y_conv);
mu_test = {0, 0};

mu1 = y_conv * mu;
cov1 = y_conv * cov * y_conv`;
cov_inv1 = inv(cov1);

t_squared1 = n * (mu1 - mu_test)` * cov_inv1 * (mu1 - mu_test);
f_critical_val1 = p_new * (n - 1) * quantile('F', 0.95, p_new, n - p_new) / (n - p_new);

if t_squared1 <= f_critical_val1 then
    null_hypothesis1 = 'Accepted';
else
    null_hypothesis1 = 'Rejected';

print t_squared1[label='Hotelling T^2'];
print f_critical_val1[label='F Critical Value'];
print null_hypothesis1[label='Null Hypothesis H0'];

quit;
