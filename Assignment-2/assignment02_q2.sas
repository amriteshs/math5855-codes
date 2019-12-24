data soil;
infile '/folders/myfolders/soil.dat';
input x1 x2 x3 x4;
run;

proc iml;
use soil var {x1 x2 x3 x4};
read all var _num_ into x;

print 'Q2(i)';

y = {1 -1 0 0, 0 1 -1 0, 0 0 1 -1} * x`;

n = ncol(y);
p = nrow(y);

mu = y[, :];
cov = ((y - mu) * (y - mu)`) / (n - 1);
mu0 = {0, 0, 0};

t_squared = n * (mu - mu0)` * inv(cov) * (mu - mu0);
f_critical = p * (n - 1) * quantile('F', 0.95, p, n - p) / (n - p);

if t_squared <= f_critical then
	null_hypothesis = 'Not Rejected';
else
	null_hypothesis = 'Rejected';

print t_squared[label="Hotelling's T^2"];
print f_critical[label='Critical Value of F'];
print null_hypothesis[label='Null Hypothesis H0'];

print 'Q2(ii)';

y1 = {1 -2 0 0, 0 1 -2 0, 0 0 1 -2} * x`;

n1 = ncol(y1);
p1 = nrow(y1);

mu1 = y1[, :];
cov1 = ((y1 - mu1) * (y1 - mu1)`) / (n1 - 1);

t_squared1 = n1 * (mu1 - mu0)` * inv(cov1) * (mu1 - mu0);
f_critical1 = p1 * (n1 - 1) * quantile('F', 0.95, p1, n1 - p1) / (n1 - p1);

if t_squared1 <= f_critical1 then
	null_hypothesis1 = 'Not Rejected';
else
	null_hypothesis1 = 'Rejected';

print t_squared1[label="Hotelling's T^2"];
print f_critical1[label='Critical Value of F'];
print null_hypothesis1[label='Null Hypothesis H0'];

quit;