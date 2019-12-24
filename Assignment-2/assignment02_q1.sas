data brothers;
infile '/folders/myfolders/brothers.dat';
input x1 x2 x3 x4;
run;

proc iml;
use brothers var {x1 x2 x3 x4};
read all var _num_ into x;

print 'Q1(i)';

start partial_correlation(cov, j, i);
	cov_est = cov[j, j] - (cov[j, i] * inv(cov[i, i]) * cov[i, j]);
	return cov_est[1, 2] / sqrt(cov_est[1, 1] * cov_est[2, 2]);
finish;

n = nrow(x);
p = ncol(x);

mu = x[:, ]`;
cov = ((x` - mu) * (x` - mu)`) / (n - 1);

r_34_12 = partial_correlation(cov, {3 4}, {1 2});
r_21_34 = partial_correlation(cov, {1 2}, {3 4});

print r_34_12[label='R(34.12)'];
print r_21_34[label='R(21.34)'];

print 'Q1(ii)';

r_12 = cov[1, 2] / sqrt(cov[1, 1] * cov[2, 2]);
print r_12[label='R(12)'];

print 'Q1(iii)';

fisher_z = 0.5 * log((1 + r_21_34) / (1 - r_21_34));
std_error = 1 / sqrt(n - (p - 2) - 3);

z_L = fisher_z + quantile('Normal', (0.05 / 2), 0, 1) * std_error;
z_U = fisher_z + quantile('Normal', 1 - (0.05 / 2), 0, 1) * std_error;

ci = ((exp(2 * z_L) - 1) / (exp(2 * z_L) + 1)) || ((exp(2 * z_U) - 1) / (exp(2 * z_U) + 1));

print fisher_z[label="Fisher's Z"];
print ci[label='Confidence Interval' colname={'Lower limit' 'Upper limit'}];

print 'Q1(iv)';

start r(cov, i, j);
	return cov[i, j] / sqrt(cov[i, i] * cov[j, j]);
finish;

r = sqrt(((r(cov, 1, 3) ** 2) + (r(cov, 2, 3) ** 2) - (2 * r(cov, 1, 3) * r(cov, 2, 3) * r(cov, 1, 2))) / (1 - (r(cov, 1, 2) ** 2)));

p_new = 3;
f = ((r ** 2) * (n - p_new)) / ((1 - r ** 2) * (p_new - 1));
f_critical = quantile('F', 0.95, p_new - 1, n - p_new);

if f <= f_critical then
	null_hypothesis = 'Not Rejected';
else
	null_hypothesis = 'Rejected';

print r[label='Multiple Correlation Coefficient'];
print f[label='Test Statistic F'];
print f_critical[label='Critical Value of F'];
print null_hypothesis[label='Null Hypothesis H0'];

print 'Q1(v)';

r_34 = cov[3, 4] / sqrt(cov[3, 3] * cov[4, 4]);

t = r_34 * sqrt((n - 2) / (1 - r_34 ** 2));
t_critical = quantile('T', 1 - (0.05 / 2), n - 2);

if t <= t_critical then
	null_hypothesis1 = 'Not Rejected';
else
	null_hypothesis1 = 'Rejected';

print t[label='Test Statistic T'];
print t_critical[label='Critical Value of T'];
print null_hypothesis1[label='Null Hypothesis H0'];

print 'Q1(i) [TESTING]';

proc corr data = brothers nosimple noprob;
var x3 x4;
partial x1 x2;
run;

proc corr data = brothers nosimple noprob;
var x1 x2;
partial x3 x4;
run;

quit;