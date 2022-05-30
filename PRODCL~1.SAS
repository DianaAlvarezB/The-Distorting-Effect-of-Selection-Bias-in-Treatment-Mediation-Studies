/*Prodclin.sas - MacKinnon, Fritz, Williams, and Lockwood, 2005       */
/*Dept. of Psychology, Arizona State University                       */
/*Necessary files:prodclin.sas, prodclin.exe                          */

/*The prodclin program calculates asymmetric confidence intervals for */
/*the mediated effect using the distribution of the product of two    */
/*random variables and the Fortran program Prodclin.exe.  To run, open*/
/* in SAS or copy and paste the code into a SAS editor window.  Then, */
/*the location of the Fortran program needs to be entered in four     */
/*places(both files should be in the same location).  Next, values for*/
/*the a and b paths, as well as the standarderrors of a (sea) and b   */
/*(seb) must be input in the macro call line at the end of the        */
/*program.  Lastly, values for rho (correlation between a and b, in   */
/*most cases equal to zero) and alpha (the type 1 error rate) must be */
/*specified in the beginning of the program. Output consists of the   */
/*user specified values for a, sea, b, seb, rho, alpha, ab, the first */
/*order standard error (sobelse), the lower and upper standardized    */
/*critical values from Prodclin [low,high], the asymmetric confidence */
/*limits [prodlow,produp] and the corresponding normal theory         */
/*confidence limits [normlow,normup]. If either prodlow or produp is  */
/*equal to 99999, there was a problem with convergence in the Fortran */
/*program.  First, check the values you input and then, if the problem*/
/*still occurs, try using fewer decimal places.                       */


/*All places that need user input have been labeled with comments.    */
*Change file address to match the location of the file prodclin.exe; 

*data from the simulation results;


data allperc;
set "C:\Users\sjlamp\Documents\collider in mediation\negative coefficients\results_zeros_allperc.sas7bdat";
run;


ods graphics off;
ods exclude all;ods noresults; options nonotes;

options noxwait;

proc sort data=allperc;
by combo rep cutoff;
run;

data results;
run;

%macro sim;
data trial;
call streaminit(52712);
%do row =1 %to 21000;

*saving values of a, b, ab, and std errors;
data allperc;
set allperc;
if _n_ = &row then call symputx("a", a);
if _n_ = &row then call symputx("b", b);
if _n_ = &row then call symputx("sea", a_stderr);
if _n_ = &row then call symputx("seb", b_stderr);
run;

*Input values for rho and alpha;
%let rho=0;
%let alpha=.05;
%let a = &a;
%let b = &b;
%let sea = &sea;
%let seb = &seb;

data data1; 
*Change file address to match the location of the file prodclin.exe;
file "C:\Users\sjlamp\Documents\PRODCLIN\raw.txt";
a=&a; sea=&sea; b=&b; seb=&seb; rho=&rho; alpha=&alpha;
put a @; put sea @; put b @; put seb @; put rho @; put alpha @;data data1;

*Change file address to match the location of the file prodclin.exe; 
X cd C:\Users\sjlamp\Documents\PRODCLIN\;
*Change file address to match the location of the file prodclin.exe;
X call "C:\Users\sjlamp\Documents\PRODCLIN\ProdClin2_Sas.exe";

data data2; 
*Change file address to match the location of the file prodclin.exe;
infile "C:\Users\sjlamp\Documents\PRODCLIN\critval.txt";
input lcrit ucrit;
a=&a; sea=&sea; b=&b; seb=&seb; rho=&rho; alpha=&alpha;
da=a/sea;
db=b/seb;
sedadb=sqrt(da*da+db*db+1);
dadb=da*db;
low=(lcrit-dadb)/(sedadb);
high=(ucrit-dadb)/(sedadb);
ab=a*b;
sobelse=sqrt(a*a*seb*seb+b*b*sea*sea);
prodlow=ab+low*sobelSE;
produp=ab+high*sobelSE;
nl=probit(alpha/2);
normlow=ab+nl*sobelse;
normup=ab-nl*sobelse;
if lcrit=99999 then prodlow=99999;
if ucrit=99999 then produp=99999;


data results;
set results data2;
run;

run;
%end;
%mend;
%sim;
run;

data results1;
set results;
if _n_ = 1 then delete;
run;

data info;
set allperc;
keep combo cutoff rep b_sig ab_sig;
run;

data blah;
merge results1 info;
run;

data blah;
set blah;
if lcrit*ucrit > 0 then d_power = 1;
else d_power=0;
run;

ods graphics on;
ods exclude none;ods results;
options notes;

proc sort data=blah;
by combo cutoff;
run;

proc means data=blah mean;
var b_sig ab_sig d_power;
by combo cutoff;
output out=sim_results mean=b_sig sobel distprod;
run;

data sim_results;
set sim_results;
drop _type_ _freq_;
run;
