data bank;
infile '/folders/myfolders/bank.dat';
input x1-x6;
forge = 2 - (_n_ < 101);
run;

proc iml;
print 'Q1(i)';

proc princomp cov data=bank;
var x1-x6;
ods select eigenvalues eigenvectors;
run;

proc iml;
print 'Q1(ii)';

proc princomp data=bank;
var x1-x6;
ods select eigenvalues eigenvectors;
run;

proc iml;
print 'Q1(iii)';

proc princomp cov data=bank plots(only)=(score(ncomp=2));
var x1-x6;
id forge;
ods select scoreplot;
run;

proc iml;
print 'Q1(iv)';

proc discrim data=bank method=normal pool=yes crosslisterr;
var x1-x6;
class forge;
ods select lineardiscfunc postcrossval;
run;

quit;