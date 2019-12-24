data pricecons;
infile '/folders/myfolders/price-cons.dat';
input y x1-x4;
run;

proc iml;
print 'Q3';

proc cancorr data=pricecons;
var x1-x2;
with x3-x4;
run;

quit;