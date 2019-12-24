proc iml;

x = {2.9 1.2 0.75, 4.6 1.6 1.2, 0.6 0.15 0.25};
cov = {2.30 0.25 0.47, 0.25 0.60 0.03, 0.47 0.03 0.59};

print 'Q3(i)';

start d_hat(x, cov, i);
	return (x[i, ] * inv(cov)) || (-0.5 * x[i, ] * inv(cov) * x[i, ]`);
finish;

d1_hat = d_hat(x, cov, 1);
d2_hat = d_hat(x, cov, 2);
d3_hat = d_hat(x, cov, 3);

print d1_hat[label='Linear Discriminant Score for Anxiety State Group' colname={'Coefficient of A' 'Coefficient of B' 'Coefficient of C' 'Constant term'}];
print d2_hat[label='Linear Discriminant Score for Obsession Group' colname={'Coefficient of A' 'Coefficient of B' 'Coefficient of C' 'Constant term'}];
print d3_hat[label='Linear Discriminant Score for Normal Group' colname={'Coefficient of A' 'Coefficient of B' 'Coefficient of C' 'Constant term'}];

print 'Q3(ii)';

individuals = {3.00 1.20 1.00, 4.00 1.40 1.32, 1.00 0.50 0.33};
d_hat_scores = d1_hat` || d2_hat` || d3_hat`;

start individual_group(individuals, scores, i);
	max_d_index = ((individuals[i, ] || {1}) * scores)[, <:>];
	
	if max_d_index = 1 then
		group = 'Anxiety State';
	else if max_d_index = 2 then
		group = 'Obsession';
	else if max_d_index = 3 then
		group = 'Normal';
	
	return group;
finish;

group_mary = individual_group(individuals, d_hat_scores, 1);
group_fred = individual_group(individuals, d_hat_scores, 2);
group_gise = individual_group(individuals, d_hat_scores, 3);

print group_mary[label="Mary's Group"];
print group_fred[label="Fred's Group"];
print group_gise[label="Giselda's Group"];

print 'Q3(iii)';

x1 = (x[1, ]` || x[2, ]`)`;
fisher_d_hat = ((x1[1, ] - x1[2, ]) * inv(cov)) || (-0.5 * (x1[1, ] - x1[2, ]) * inv(cov) * (x1[1, ] + x1[2, ])`);

delta_hat = sqrt((x1[1, ] - x1[2, ]) * inv(cov) * (x1[1, ] - x1[2, ])`);
p21 = cdf('Normal', (-0.5 * delta_hat));
p12 = p21;

print fisher_d_hat[label="Linear Discriminant Function ('Anxiety State' if >= 0, else 'Obsession')" colname={'Coefficient of A' 'Coefficient of B' 'Coefficient of C' 'Constant term'}];
print p21[label='P(2|1)'];
print p12[label='P(1|2)'];

quit;