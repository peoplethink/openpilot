/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8096049958724316807);
void inv_err_fun(double *nom_x, double *true_x, double *out_8279088223888228505);
void H_mod_fun(double *state, double *out_2457549075545957044);
void f_fun(double *state, double dt, double *out_5188444040622817444);
void F_fun(double *state, double dt, double *out_6880204027575935328);
void h_25(double *state, double *unused, double *out_3624333546785423039);
void H_25(double *state, double *unused, double *out_8668137839260624228);
void h_24(double *state, double *unused, double *out_1640118378038141378);
void H_24(double *state, double *unused, double *out_1409310185087476488);
void h_30(double *state, double *unused, double *out_2429398138357197661);
void H_30(double *state, double *unused, double *out_924817957936797882);
void h_26(double *state, double *unused, double *out_3599589183117062072);
void H_26(double *state, double *unused, double *out_5430464877248331228);
void h_27(double *state, double *unused, double *out_5324631439891927468);
void H_27(double *state, double *unused, double *out_7574479069461655864);
void h_29(double *state, double *unused, double *out_7901847793471347442);
void H_29(double *state, double *unused, double *out_1364363635470041804);
void h_28(double *state, double *unused, double *out_751004511462432978);
void H_28(double *state, double *unused, double *out_8869473968987790654);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
