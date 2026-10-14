data muscle;
    input obs Y X;
    datalines;
1 106 43
2 106 41
3 97 47
4 113 46
5 96 45
6 119 41
7 92 47
8 112 41
9 92 48
10 102 48
11 107 42
12 107 47
13 102 43
14 115 44
15 101 42
16 87 55
17 91 57
18 97 56
19 82 59
20 78 57
21 95 54
22 98 53
23 94 52
24 96 53
25 100 54
26 84 60
27 70 59
28 104 51
29 76 59
30 93 57
31 73 68
32 73 63
33 76 60
34 80 63
35 84 63
36 71 64
37 64 66
38 88 65
39 79 60
40 88 65
41 73 65
42 74 69
43 76 61
44 87 70
45 70 68
46 69 78
47 54 78
48 62 78
49 78 72
50 65 70
51 64 73
52 74 76
53 87 78
54 63 78
55 82 71
56 80 75
57 52 77
58 56 76
59 70 72
60 74 76
;
run;

proc print data=muscle;
run;

/*(a)*/

ods graphics off;

proc univariate data=muscle plot;
    var X;
    title "Stem-and-Leaf Plot of Ages";
run;

ods graphics on;

 

/* Part (b): Obtain residuals and make a dot plot */

proc reg data=muscle;
    model Y = X;
    output out=muscle_resid
        predicted=Yhat
        residual=e;
run;
quit;

/* Print residuals */
proc print data=muscle_resid;
    var Y X Yhat e;
    title "Residuals from the Regression of Muscle Mass on Age";
run;

/* Dot plot of residuals */
proc sgplot data=muscle_resid;
    dot e;
    xaxis label="Residual";
    title "Dot Plot of Residuals";
run;

/* Part (c): Residual plots */

proc reg data=muscle;
    model Y = X;
    output out=muscle_resid
        predicted=Yhat
        residual=e;
run;
quit;

/* Residuals versus fitted values */
proc sgplot data=muscle_resid;
    scatter x=Yhat y=e;
    refline 0 / axis=y lineattrs=(pattern=shortdash);
    xaxis label="Fitted Values";
    yaxis label="Residuals";
    title "Residuals versus Fitted Values";
run;

/* Residuals versus age */
proc sgplot data=muscle_resid;
    scatter x=X y=e;
    refline 0 / axis=y lineattrs=(pattern=shortdash);
    xaxis label="Age";
    yaxis label="Residuals";
    title "Residuals versus Age";
run;

/* Part (d): Normal probability plot and Shapiro-Wilk test */

proc reg data=muscle;
    model Y = X;
    output out=muscle_resid
        predicted=Yhat
        residual=e;
run;
quit;

/* Normal probability plot and Shapiro-Wilk test */
proc univariate data=muscle_resid normal;
    var e;
    qqplot e / normal(mu=est sigma=est);
    title "Normal Probability Plot and Shapiro-Wilk Test for Residuals";
run;

/*\(H_0\): residuals are normally distributed.
\(H_a\): residuals are not normally distributed.
If Pr < W is small, for example less than 0.05, reject normality.
If Pr < W is not small, we do not have strong evidence against normality.*/

 /* Part (e): Breusch-Pagan test manually */

/* Step 1: Fit original regression and save residuals */
proc reg data=muscle;
    model Y = X;
    output out=muscle_resid
        predicted=Yhat
        residual=e;
run;
quit;

/* Step 2: Compute squared residuals */
data muscle_bp;
    set muscle_resid;
    e2 = e**2;
run;

/* Step 3: Auxiliary regression: e_i^2 on X */
ods output FitStatistics=fit_bp;

proc reg data=muscle_bp;
    model e2 = X;
    title "Auxiliary Regression for Breusch-Pagan Test";
run;
quit;

ods output close;

/* Step 4: Extract R-square */
data _null_;
    set fit_bp;
    if Label2 = "R-Square" then call symputx("Rsq", nValue2);
run;

/* Step 5: Compute BP = n R^2 and p-value */
proc sql noprint;
    select count(*) into :n
    from muscle_bp
    where e2 is not missing and X is not missing;
quit;

data bp_result;
    n = &n;
    Rsq = &Rsq;
    BP = n * Rsq;
    df = 1;
    p_value = 1 - probchi(BP, df);
run;

/* Step 6: Print chi-square statistic and p-value */
proc print data=bp_result noobs;
    var n Rsq BP df p_value;
    title "Breusch-Pagan Test Statistic and P-value";
run;

/*\(H_0\): constant error variance.
\(H_a\): error variance depends on \(X\).
If p_value < alpha, reject \(H_0\).
If p_value >= alpha, do not reject \(H_0\).*/
