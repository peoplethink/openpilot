/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7783346681631977471);
void inv_err_fun(double *nom_x, double *true_x, double *out_3370115840177616313);
void H_mod_fun(double *state, double *out_212879511712864462);
void f_fun(double *state, double dt, double *out_4638309030321542643);
void F_fun(double *state, double dt, double *out_7325456146553586026);
void h_3(double *state, double *unused, double *out_802485685206233605);
void H_3(double *state, double *unused, double *out_5539980980987807745);
void h_4(double *state, double *unused, double *out_7027348924825988923);
void H_4(double *state, double *unused, double *out_2172946816918086839);
void h_9(double *state, double *unused, double *out_2930814106272619331);
void H_9(double *state, double *unused, double *out_6673049544758445784);
void h_10(double *state, double *unused, double *out_6364753259239910793);
void H_10(double *state, double *unused, double *out_5267622721369518965);
void h_12(double *state, double *unused, double *out_3131420646843713785);
void H_12(double *state, double *unused, double *out_6263659756908986840);
void h_31(double *state, double *unused, double *out_5547007004235775228);
void H_31(double *state, double *unused, double *out_64332741002999008);
void h_32(double *state, double *unused, double *out_3800963807657689828);
void H_32(double *state, double *unused, double *out_1074388609318319310);
void h_13(double *state, double *unused, double *out_4818214050445037326);
void H_13(double *state, double *unused, double *out_3569746897862172482);
void h_14(double *state, double *unused, double *out_2930814106272619331);
void H_14(double *state, double *unused, double *out_6673049544758445784);
void h_19(double *state, double *unused, double *out_963171697598985061);
void H_19(double *state, double *unused, double *out_4668531848327094932);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);