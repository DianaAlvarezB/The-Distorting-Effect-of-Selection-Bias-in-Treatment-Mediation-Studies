*Collider in Prevention code draft for restricting the range of the collider (U) variable in the 4-variable model;

/*TITLE 'SIMULATION 4-VARIABLE MEDIATOR MODEL WITH A COLLIDER (U) CAUSED BY M AND Y';
OPTIONS PS=59 LS=80 REPLACE NONOTES;
FILENAME NULLOG DUMMY 'C:\NULL';
PROC PRINTTO LOG=NULLOG;*/

data list;
run;
data correlations;
run;

%MACRO SIMULATE(NREP,NOBS,B1,B2,B3,B4,B5,FILE,TYPE,ERROR,combo);
DATA SUMMARY; 
SET _NULL_;  
%DO rep=1 %TO &NREP;
*GENERATE*;
DATA SIM;
DO I=1 TO &NOBS;
X = RAND('BERNOULLI', .5);
M = &B1*X + (&error)*RANNOR(0);
Y = &B2*X + &B3*M + (&error)*RANNOR(0);
U = &B4*M + &B5*Y + (&error)*RANNOR(0);
combo=&combo;
OUTPUT;
END;

ods graphics off;
ods exclude all;ods noresults; options nonotes;

ods output Quantiles=perc;       
proc univariate data = SIM;
var U; 
quit;

*this saves the individual values for each percentile as a macro variable;
data percentiles;
set perc;
if Quantile="10%" then call symputx("P10", Estimate);
if Quantile="25% Q1" then call symputx("P25", Estimate);
if Quantile="50% Median" then call symputx("P50", Estimate);
if Quantile="75% Q3" then call symputx("P75", Estimate);
if Quantile="90%" then call symputx("P90", Estimate);
if Quantile="95%" then call symputx("P95", Estimate);
run;
%put P10 = &P10;
%put P25 = &P25;
%put P50 = &P50;
%put P75 = &P75;
%put P90 = &P90;
%put P95 = &P95;

data CUT_75;
set sim;
IF U < &P75 THEN DELETE;
run;

proc surveyselect data=CUT_75 noprint method=srs n=200 out=s_75th;
run;

data CUT_90;
set sim;
IF U < &P90 THEN DELETE;
run;

proc surveyselect data=CUT_90 noprint method=srs n=200 out=s_90th;
run;

data CUT_95;
set sim;
IF U < &P95 THEN DELETE;
run;

proc surveyselect data=CUT_95 noprint method=srs n=200 out=s_95th;
run;

*pulling a random sample with no selection on U;
proc surveyselect data=sim noprint method=srs n=200 out=s_none;
run;

/*other cutoff values;
data CUT_10;
set sim;
IF U < &P10 THEN DELETE;
run;

proc surveyselect data=CUT_10 noprint method=srs n=200 out=s_10th;
run;

data CUT_25;
set sim;
IF U < &P25 THEN DELETE;
run;

proc surveyselect data=CUT_25 noprint method=srs n=200 out=s_25th;
run;

data CUT_50;
set sim;
IF U < &P50 THEN DELETE;
run;

proc surveyselect data=CUT_50 noprint method=srs n=200 out=s_50th;
run;*/


*mediation analysis for each sample and saving the coefficients;
*75th perc;
proc reg data=s_75th;
model y = x m;
ODS output ParameterEstimates=reg_75th;
run;

data reg_75th;
set reg_75th;
if Variable='X' then call symputx("cp", Estimate);
if Variable='X' then call symputx("cp_se", StdErr);
if Variable='M' then call symputx("b", Estimate);
if Variable='M' then call symputx("b_se", StdErr);
run;

proc reg data=s_75th;
model m = x;
ODS output ParameterEstimates=reg2_75th;
run;

data reg2_75th;
set reg2_75th;
if Variable='X' then call symputx("a", Estimate);
if Variable='X' then call symputx("a_se", StdErr);
run;

proc reg data=s_75th;
model y = x;
ODS output ParameterEstimates=reg3_75th;
run;

data reg3_75th;
set reg3_75th;
if Variable='X' then call symputx("c", Estimate);
if Variable='X' then call symputx("c_se", StdErr);
run;

data coef_75;
combo=&combo;
cutoff='75th';
a_true=&B1;
b_true=&B3;
cp_true=&B2;
d_true=&B4;
e_true=&B5;
ab_true=a_true*b_true;
c_true=ab_true+cp_true;
rep=&rep;
a=&a;
b=&b;
c=&c;
cp=&cp;
a_stderr=&a_se;
b_stderr=&b_se;
c_stderr=&c_se;
cp_stderr=&cp_se;
ab=&a*&b;
a_z = a/a_stderr;
if abs(a_z) > abs(1.96) then a_sig=1;
else a_sig=0;
b_z=b/b_stderr;
if abs(b_z) > abs(1.96) then b_sig=1;
else b_sig=0;
c_z=c/c_stderr;
if abs(c_z) > abs(1.96) then c_sig=1;
else c_sig=0;
cp_z=cp/cp_stderr;
if abs(cp_z) > abs(1.96) then cp_sig=1;
else cp_sig=0;
ab_stderr=SQRT(&A*&A*&b_se*&b_se+&B*&B*&a_se*&a_se);
ab_z=ab/ab_stderr;
if abs(ab_z) > abs(1.96) then ab_sig=1;
else ab_sig=0;
bias_a = a_true - a;
rb_a = bias_a/a_true;
bias_b = b_true - b;
rb_b = bias_b/b_true;
bias_cp = cp_true - cp;
rb_cp = bias_cp/cp_true;
bias_ab = ab_true - ab;
rb_ab = bias_ab/ab_true;
bias_c = c_true - c;
rb_c = bias_c/c_true;
run;

*90th perc;
proc reg data=s_90th;
model y = x m;
ODS output ParameterEstimates=reg_90th;
run;

data reg_90th;
set reg_90th;
if Variable='X' then call symputx("cp", Estimate);
if Variable='X' then call symputx("cp_se", StdErr);
if Variable='M' then call symputx("b", Estimate);
if Variable='M' then call symputx("b_se", StdErr);
run;

proc reg data=s_90th;
model m = x;
ODS output ParameterEstimates=reg2_90th;
run;

data reg2_90th;
set reg2_90th;
if Variable='X' then call symputx("a", Estimate);
if Variable='X' then call symputx("a_se", StdErr);
run;

proc reg data=s_90th;
model y = x;
ODS output ParameterEstimates=reg3_90th;
run;

data reg3_90th;
set reg3_90th;
if Variable='X' then call symputx("c", Estimate);
if Variable='X' then call symputx("c_se", StdErr);
run;

data coef_90;
combo=&combo;
cutoff='90th';
a_true=&B1;
b_true=&B3;
cp_true=&B2;
d_true=&B4;
e_true=&B5;
ab_true=a_true*b_true;
c_true=ab_true+cp_true;
rep=&rep;
a=&a;
b=&b;
c=&c;
cp=&cp;
a_stderr=&a_se;
b_stderr=&b_se;
c_stderr=&c_se;
cp_stderr=&cp_se;
ab=&a*&b;
a_z = a/a_stderr;
if abs(a_z) > abs(1.96) then a_sig=1;
else a_sig=0;
b_z=b/b_stderr;
if abs(b_z) > abs(1.96) then b_sig=1;
else b_sig=0;
c_z=c/c_stderr;
if abs(c_z) > abs(1.96) then c_sig=1;
else c_sig=0;
cp_z=cp/cp_stderr;
if abs(cp_z) > abs(1.96) then cp_sig=1;
else cp_sig=0;
ab_stderr=SQRT(&A*&A*&b_se*&b_se+&B*&B*&a_se*&a_se);
ab_z=ab/ab_stderr;
if abs(ab_z) > abs(1.96) then ab_sig=1;
else ab_sig=0;
bias_a = a_true - a;
rb_a = bias_a/a_true;
bias_b = b_true - b;
rb_b = bias_b/b_true;
bias_cp = cp_true - cp;
rb_cp = bias_cp/cp_true;
bias_ab = ab_true - ab;
rb_ab = bias_ab/ab_true;
bias_c = c_true - c;
rb_c = bias_c/c_true;
run;

*95th perc;
proc reg data=s_95th;
model y = x m;
ODS output ParameterEstimates=reg_95th;
run;

data reg_95th;
set reg_95th;
if Variable='X' then call symputx("cp", Estimate);
if Variable='X' then call symputx("cp_se", StdErr);
if Variable='M' then call symputx("b", Estimate);
if Variable='M' then call symputx("b_se", StdErr);
run;

proc reg data=s_95th;
model m = x;
ODS output ParameterEstimates=reg2_95th;
run;

data reg2_95th;
set reg2_95th;
if Variable='X' then call symputx("a", Estimate);
if Variable='X' then call symputx("a_se", StdErr);
run;

proc reg data=s_95th;
model y = x;
ODS output ParameterEstimates=reg3_95th;
run;

data reg3_95th;
set reg3_95th;
if Variable='X' then call symputx("c", Estimate);
if Variable='X' then call symputx("c_se", StdErr);
run;

data coef_95;
combo=&combo;
cutoff='95th';
a_true=&B1;
b_true=&B3;
cp_true=&B2;
d_true=&B4;
e_true=&B5;
ab_true=a_true*b_true;
c_true=ab_true+cp_true;
rep=&rep;
a=&a;
b=&b;
c=&c;
cp=&cp;
a_stderr=&a_se;
b_stderr=&b_se;
c_stderr=&c_se;
cp_stderr=&cp_se;
ab=&a*&b;
a_z = a/a_stderr;
if abs(a_z) > abs(1.96) then a_sig=1;
else a_sig=0;
b_z=b/b_stderr;
if abs(b_z) > abs(1.96) then b_sig=1;
else b_sig=0;
c_z=c/c_stderr;
if abs(c_z) > abs(1.96) then c_sig=1;
else c_sig=0;
cp_z=cp/cp_stderr;
if abs(cp_z) > abs(1.96) then cp_sig=1;
else cp_sig=0;
ab_stderr=SQRT(&A*&A*&b_se*&b_se+&B*&B*&a_se*&a_se);
ab_z=ab/ab_stderr;
if abs(ab_z) > abs(1.96) then ab_sig=1;
else ab_sig=0;
bias_a = a_true - a;
rb_a = bias_a/a_true;
bias_b = b_true - b;
rb_b = bias_b/b_true;
bias_cp = cp_true - cp;
rb_cp = bias_cp/cp_true;
bias_ab = ab_true - ab;
rb_ab = bias_ab/ab_true;
bias_c = c_true - c;
rb_c = bias_c/c_true;
run;

*Population;
proc reg data=s_None;
model y = x m;
ODS output ParameterEstimates=reg_None;
run;

data reg_None;
set reg_None;
if Variable='X' then call symputx("cp", Estimate);
if Variable='X' then call symputx("cp_se", StdErr);
if Variable='M' then call symputx("b", Estimate);
if Variable='M' then call symputx("b_se", StdErr);
run;

proc reg data=s_None;
model m = x;
ODS output ParameterEstimates=reg2_None;
run;

data reg2_None;
set reg2_None;
if Variable='X' then call symputx("a", Estimate);
if Variable='X' then call symputx("a_se", StdErr);
run;

proc reg data=s_None;
model y = x;
ODS output ParameterEstimates=reg3_None;
run;

data reg3_None;
set reg3_None;
if Variable='X' then call symputx("c", Estimate);
if Variable='X' then call symputx("c_se", StdErr);
run;

data coef_none;
combo=&combo;
cutoff='None';
a_true=&B1;
b_true=&B3;
cp_true=&B2;
d_true=&B4;
e_true=&B5;
ab_true=a_true*b_true;
c_true=ab_true+cp_true;
rep=&rep;
a=&a;
b=&b;
c=&c;
cp=&cp;
a_stderr=&a_se;
b_stderr=&b_se;
c_stderr=&c_se;
cp_stderr=&cp_se;
ab=&a*&b;
a_z = a/a_stderr;
if abs(a_z) > abs(1.96) then a_sig=1;
else a_sig=0;
b_z=b/b_stderr;
if abs(b_z) > abs(1.96) then b_sig=1;
else b_sig=0;
c_z=c/c_stderr;
if abs(c_z) > abs(1.96) then c_sig=1;
else c_sig=0;
cp_z=cp/cp_stderr;
if abs(cp_z) > abs(1.96) then cp_sig=1;
else cp_sig=0;
ab_stderr=SQRT(&A*&A*&b_se*&b_se+&B*&B*&a_se*&a_se);
ab_z=ab/ab_stderr;
if abs(ab_z) > abs(1.96) then ab_sig=1;
else ab_sig=0;
bias_a = a_true - a;
rb_a = bias_a/a_true;
bias_b = b_true - b;
rb_b = bias_b/b_true;
bias_cp = cp_true - cp;
rb_cp = bias_cp/cp_true;
bias_ab = ab_true - ab;
rb_ab = bias_ab/ab_true;
bias_c = c_true - c;
rb_c = bias_c/c_true;
run;

/*OTHER PERCENTILE CUTOFFS;
*10th perc;
proc reg data=s_10th;
model y = x m;
ODS output ParameterEstimates=reg_10th;
run;

data reg_10th;
set reg_10th;
if Variable='X' then call symputx("cp", Estimate);
if Variable='X' then call symputx("cp_se", StdErr);
if Variable='M' then call symputx("b", Estimate);
if Variable='M' then call symputx("b_se", StdErr);
run;

proc reg data=s_10th;
model m = x;
ODS output ParameterEstimates=reg2_10th;
run;

data reg2_10th;
set reg2_10th;
if Variable='X' then call symputx("a", Estimate);
if Variable='X' then call symputx("a_se", StdErr);
run;

data coef_10;
combo=&combo;
cutoff='10th';
a_true=&B1;
b_true=&B3;
cp_true=&B2;
d_true=&B4;
e_true=&B5;
ab_true=a_true*b_true;
c_true=ab_true+cp_true;
rep=&rep;
a=&a;
b=&b;
c=&c;
cp=&cp;
a_stderr=&a_se;
b_stderr=&b_se;
c_stderr=&c_se;
cp_stderr=&cp_se;
ab=&a*&b;
a_z = a/a_stderr;
if abs(a_z) > abs(1.96) then a_sig=1;
else a_sig=0;
b_z=b/b_stderr;
if abs(b_z) > abs(1.96) then b_sig=1;
else b_sig=0;
c_z=c/c_stderr;
if abs(c_z) > abs(1.96) then c_sig=1;
else c_sig=0;
cp_z=cp/cp_stderr;
if abs(cp_z) > abs(1.96) then cp_sig=1;
else cp_sig=0;
ab_stderr=SQRT(&A*&A*&b_se*&b_se+&B*&B*&a_se*&a_se);
ab_z=ab/ab_stderr;
if abs(ab_z) > abs(1.96) then ab_sig=1;
else ab_sig=0;
bias_a = a_true - a;
rb_a = bias_a/a_true;
bias_b = b_true - b;
rb_b = bias_b/b_true;
bias_cp = cp_true - cp;
rb_cp = bias_cp/cp_true;
bias_ab = ab_true - ab;
rb_ab = bias_ab/ab_true;
bias_c = c_true - c;
rb_c = bias_c/c_true;
run;

*25th perc;
proc reg data=s_25th;
model y = x m;
ODS output ParameterEstimates=reg_25th;
run;

data reg_25th;
set reg_25th;
if Variable='X' then call symputx("cp", Estimate);
if Variable='X' then call symputx("cp_se", StdErr);
if Variable='M' then call symputx("b", Estimate);
if Variable='M' then call symputx("b_se", StdErr);
run;

proc reg data=s_25th;
model m = x;
ODS output ParameterEstimates=reg2_25th;
run;

data reg2_25th;
set reg2_25th;
if Variable='X' then call symputx("a", Estimate);
if Variable='X' then call symputx("a_se", StdErr);
run;

data coef_25;
combo=&combo;
cutoff='25th';
a_true=&B1;
b_true=&B3;
cp_true=&B2;
d_true=&B4;
e_true=&B5;
ab_true=a_true*b_true;
c_true=ab_true+cp_true;
rep=&rep;
a=&a;
b=&b;
c=&c;
cp=&cp;
a_stderr=&a_se;
b_stderr=&b_se;
c_stderr=&c_se;
cp_stderr=&cp_se;
ab=&a*&b;
a_z = a/a_stderr;
if abs(a_z) > abs(1.96) then a_sig=1;
else a_sig=0;
b_z=b/b_stderr;
if abs(b_z) > abs(1.96) then b_sig=1;
else b_sig=0;
c_z=c/c_stderr;
if abs(c_z) > abs(1.96) then c_sig=1;
else c_sig=0;
cp_z=cp/cp_stderr;
if abs(cp_z) > abs(1.96) then cp_sig=1;
else cp_sig=0;
ab_stderr=SQRT(&A*&A*&b_se*&b_se+&B*&B*&a_se*&a_se);
ab_z=ab/ab_stderr;
if abs(ab_z) > abs(1.96) then ab_sig=1;
else ab_sig=0;
bias_a = a_true - a;
rb_a = bias_a/a_true;
bias_b = b_true - b;
rb_b = bias_b/b_true;
bias_cp = cp_true - cp;
rb_cp = bias_cp/cp_true;
bias_ab = ab_true - ab;
rb_ab = bias_ab/ab_true;
bias_c = c_true - c;
rb_c = bias_c/c_true;
run;

*50th perc;
proc reg data=s_50th;
model y = x m;
ODS output ParameterEstimates=reg_50th;
run;

data reg_50th;
set reg_50th;
if Variable='X' then call symputx("cp", Estimate);
if Variable='X' then call symputx("cp_se", StdErr);
if Variable='M' then call symputx("b", Estimate);
if Variable='M' then call symputx("b_se", StdErr);
run;

proc reg data=s_50th;
model m = x;
ODS output ParameterEstimates=reg2_50th;
run;

data reg2_50th;
set reg2_50th;
if Variable='X' then call symputx("a", Estimate);
if Variable='X' then call symputx("a_se", StdErr);
run;

data coef_50;
combo=&combo;
cutoff='50th';
a_true=&B1;
b_true=&B3;
cp_true=&B2;
d_true=&B4;
e_true=&B5;
ab_true=a_true*b_true;
c_true=ab_true+cp_true;
rep=&rep;
a=&a;
b=&b;
c=&c;
cp=&cp;
a_stderr=&a_se;
b_stderr=&b_se;
c_stderr=&c_se;
cp_stderr=&cp_se;
ab=&a*&b;
a_z = a/a_stderr;
if abs(a_z) > abs(1.96) then a_sig=1;
else a_sig=0;
b_z=b/b_stderr;
if abs(b_z) > abs(1.96) then b_sig=1;
else b_sig=0;
c_z=c/c_stderr;
if abs(c_z) > abs(1.96) then c_sig=1;
else c_sig=0;
cp_z=cp/cp_stderr;
if abs(cp_z) > abs(1.96) then cp_sig=1;
else cp_sig=0;
ab_stderr=SQRT(&A*&A*&b_se*&b_se+&B*&B*&a_se*&a_se);
ab_z=ab/ab_stderr;
if abs(ab_z) > abs(1.96) then ab_sig=1;
else ab_sig=0;
bias_a = a_true - a;
rb_a = bias_a/a_true;
bias_b = b_true - b;
rb_b = bias_b/b_true;
bias_cp = cp_true - cp;
rb_cp = bias_cp/cp_true;
bias_ab = ab_true - ab;
rb_ab = bias_ab/ab_true;
bias_c = c_true - c;
rb_c = bias_c/c_true;
run;*/

data list;
set list coef_75 coef_90 coef_95 coef_none;
run;

/*for if you want to save all the cutoffs;
data list;
set list coef_10 coef_25 coef_50 coef_75 coef_90 coef_95 coef_none;
run;*/


%END;


%MEND;

*9 nonzero parameter combinations;
%SIMULATE(NREP=1000,NOBS=500000,B1=.28,  B2=.14,     B3=.14, B4=.59,     B5= .59,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=1);
%SIMULATE(NREP=1000,NOBS=500000,B1=.78,  B2=.14,     B3=.39, B4=.59,     B5= .59,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=2);
%SIMULATE(NREP=1000,NOBS=500000,B1=1.18,  B2=.14,     B3=.59, B4=.59,     B5= .59,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=3);
%SIMULATE(NREP=1000,NOBS=500000,B1=.28,  B2=.14,     B3=.14, B4=.39,     B5= .39,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=4);
%SIMULATE(NREP=1000,NOBS=500000,B1=.78,  B2=.14,     B3=.39, B4=.39,     B5= .39,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=5);
%SIMULATE(NREP=1000,NOBS=500000,B1=1.18,  B2=.14,     B3=.59, B4=.39,     B5= .39,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=6);
%SIMULATE(NREP=1000,NOBS=500000,B1=.28,  B2=.14,     B3=.14, B4=.14,     B5=.14,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=7);
%SIMULATE(NREP=1000,NOBS=500000,B1=.78,  B2=.14,     B3=.39, B4=.14,     B5=.14,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=8);
%SIMULATE(NREP=1000,NOBS=500000,B1=1.18,  B2=.14,     B3=.59, B4=.14,     B5=.14,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=9);

*3 zero parameter combinations;
/*%SIMULATE(NREP=1000,NOBS=500000,B1=.78,  B2=.14,     B3=0, B4=.59,     B5= .59,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=37);
%SIMULATE(NREP=1000,NOBS=500000,B1=0,  B2=.14,     B3=.39, B4=.59,     B5= .59,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=38);
%SIMULATE(NREP=1000,NOBS=500000,B1=0,  B2=.14,     B3=0, B4=.59,     B5= .59,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=39);*/

*alchol treatment example;
*%SIMULATE(NREP=1000,NOBS=500000,B1=-.28,  B2=-.14,     B3=.36, B4=.59,     B5= .59,     FILE = TEMP, TYPE = 'CCCC', ERROR = 1, combo=0);

*final dataset with all estimates saved;
data full_list;
set list;
if _n_=1 then delete;
run;
