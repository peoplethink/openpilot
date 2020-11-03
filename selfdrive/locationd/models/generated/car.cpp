
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8096049958724316807) {
   out_8096049958724316807[0] = delta_x[0] + nom_x[0];
   out_8096049958724316807[1] = delta_x[1] + nom_x[1];
   out_8096049958724316807[2] = delta_x[2] + nom_x[2];
   out_8096049958724316807[3] = delta_x[3] + nom_x[3];
   out_8096049958724316807[4] = delta_x[4] + nom_x[4];
   out_8096049958724316807[5] = delta_x[5] + nom_x[5];
   out_8096049958724316807[6] = delta_x[6] + nom_x[6];
   out_8096049958724316807[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8279088223888228505) {
   out_8279088223888228505[0] = -nom_x[0] + true_x[0];
   out_8279088223888228505[1] = -nom_x[1] + true_x[1];
   out_8279088223888228505[2] = -nom_x[2] + true_x[2];
   out_8279088223888228505[3] = -nom_x[3] + true_x[3];
   out_8279088223888228505[4] = -nom_x[4] + true_x[4];
   out_8279088223888228505[5] = -nom_x[5] + true_x[5];
   out_8279088223888228505[6] = -nom_x[6] + true_x[6];
   out_8279088223888228505[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_2457549075545957044) {
   out_2457549075545957044[0] = 1.0;
   out_2457549075545957044[1] = 0.0;
   out_2457549075545957044[2] = 0.0;
   out_2457549075545957044[3] = 0.0;
   out_2457549075545957044[4] = 0.0;
   out_2457549075545957044[5] = 0.0;
   out_2457549075545957044[6] = 0.0;
   out_2457549075545957044[7] = 0.0;
   out_2457549075545957044[8] = 0.0;
   out_2457549075545957044[9] = 1.0;
   out_2457549075545957044[10] = 0.0;
   out_2457549075545957044[11] = 0.0;
   out_2457549075545957044[12] = 0.0;
   out_2457549075545957044[13] = 0.0;
   out_2457549075545957044[14] = 0.0;
   out_2457549075545957044[15] = 0.0;
   out_2457549075545957044[16] = 0.0;
   out_2457549075545957044[17] = 0.0;
   out_2457549075545957044[18] = 1.0;
   out_2457549075545957044[19] = 0.0;
   out_2457549075545957044[20] = 0.0;
   out_2457549075545957044[21] = 0.0;
   out_2457549075545957044[22] = 0.0;
   out_2457549075545957044[23] = 0.0;
   out_2457549075545957044[24] = 0.0;
   out_2457549075545957044[25] = 0.0;
   out_2457549075545957044[26] = 0.0;
   out_2457549075545957044[27] = 1.0;
   out_2457549075545957044[28] = 0.0;
   out_2457549075545957044[29] = 0.0;
   out_2457549075545957044[30] = 0.0;
   out_2457549075545957044[31] = 0.0;
   out_2457549075545957044[32] = 0.0;
   out_2457549075545957044[33] = 0.0;
   out_2457549075545957044[34] = 0.0;
   out_2457549075545957044[35] = 0.0;
   out_2457549075545957044[36] = 1.0;
   out_2457549075545957044[37] = 0.0;
   out_2457549075545957044[38] = 0.0;
   out_2457549075545957044[39] = 0.0;
   out_2457549075545957044[40] = 0.0;
   out_2457549075545957044[41] = 0.0;
   out_2457549075545957044[42] = 0.0;
   out_2457549075545957044[43] = 0.0;
   out_2457549075545957044[44] = 0.0;
   out_2457549075545957044[45] = 1.0;
   out_2457549075545957044[46] = 0.0;
   out_2457549075545957044[47] = 0.0;
   out_2457549075545957044[48] = 0.0;
   out_2457549075545957044[49] = 0.0;
   out_2457549075545957044[50] = 0.0;
   out_2457549075545957044[51] = 0.0;
   out_2457549075545957044[52] = 0.0;
   out_2457549075545957044[53] = 0.0;
   out_2457549075545957044[54] = 1.0;
   out_2457549075545957044[55] = 0.0;
   out_2457549075545957044[56] = 0.0;
   out_2457549075545957044[57] = 0.0;
   out_2457549075545957044[58] = 0.0;
   out_2457549075545957044[59] = 0.0;
   out_2457549075545957044[60] = 0.0;
   out_2457549075545957044[61] = 0.0;
   out_2457549075545957044[62] = 0.0;
   out_2457549075545957044[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_5188444040622817444) {
   out_5188444040622817444[0] = state[0];
   out_5188444040622817444[1] = state[1];
   out_5188444040622817444[2] = state[2];
   out_5188444040622817444[3] = state[3];
   out_5188444040622817444[4] = state[4];
   out_5188444040622817444[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5188444040622817444[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5188444040622817444[7] = state[7];
}
void F_fun(double *state, double dt, double *out_6880204027575935328) {
   out_6880204027575935328[0] = 1;
   out_6880204027575935328[1] = 0;
   out_6880204027575935328[2] = 0;
   out_6880204027575935328[3] = 0;
   out_6880204027575935328[4] = 0;
   out_6880204027575935328[5] = 0;
   out_6880204027575935328[6] = 0;
   out_6880204027575935328[7] = 0;
   out_6880204027575935328[8] = 0;
   out_6880204027575935328[9] = 1;
   out_6880204027575935328[10] = 0;
   out_6880204027575935328[11] = 0;
   out_6880204027575935328[12] = 0;
   out_6880204027575935328[13] = 0;
   out_6880204027575935328[14] = 0;
   out_6880204027575935328[15] = 0;
   out_6880204027575935328[16] = 0;
   out_6880204027575935328[17] = 0;
   out_6880204027575935328[18] = 1;
   out_6880204027575935328[19] = 0;
   out_6880204027575935328[20] = 0;
   out_6880204027575935328[21] = 0;
   out_6880204027575935328[22] = 0;
   out_6880204027575935328[23] = 0;
   out_6880204027575935328[24] = 0;
   out_6880204027575935328[25] = 0;
   out_6880204027575935328[26] = 0;
   out_6880204027575935328[27] = 1;
   out_6880204027575935328[28] = 0;
   out_6880204027575935328[29] = 0;
   out_6880204027575935328[30] = 0;
   out_6880204027575935328[31] = 0;
   out_6880204027575935328[32] = 0;
   out_6880204027575935328[33] = 0;
   out_6880204027575935328[34] = 0;
   out_6880204027575935328[35] = 0;
   out_6880204027575935328[36] = 1;
   out_6880204027575935328[37] = 0;
   out_6880204027575935328[38] = 0;
   out_6880204027575935328[39] = 0;
   out_6880204027575935328[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6880204027575935328[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6880204027575935328[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6880204027575935328[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6880204027575935328[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6880204027575935328[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6880204027575935328[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6880204027575935328[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6880204027575935328[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6880204027575935328[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6880204027575935328[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6880204027575935328[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6880204027575935328[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6880204027575935328[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6880204027575935328[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6880204027575935328[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6880204027575935328[56] = 0;
   out_6880204027575935328[57] = 0;
   out_6880204027575935328[58] = 0;
   out_6880204027575935328[59] = 0;
   out_6880204027575935328[60] = 0;
   out_6880204027575935328[61] = 0;
   out_6880204027575935328[62] = 0;
   out_6880204027575935328[63] = 1;
}
void h_25(double *state, double *unused, double *out_3624333546785423039) {
   out_3624333546785423039[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8668137839260624228) {
   out_8668137839260624228[0] = 0;
   out_8668137839260624228[1] = 0;
   out_8668137839260624228[2] = 0;
   out_8668137839260624228[3] = 0;
   out_8668137839260624228[4] = 0;
   out_8668137839260624228[5] = 0;
   out_8668137839260624228[6] = 1;
   out_8668137839260624228[7] = 0;
}
void h_24(double *state, double *unused, double *out_1640118378038141378) {
   out_1640118378038141378[0] = state[4];
   out_1640118378038141378[1] = state[5];
}
void H_24(double *state, double *unused, double *out_1409310185087476488) {
   out_1409310185087476488[0] = 0;
   out_1409310185087476488[1] = 0;
   out_1409310185087476488[2] = 0;
   out_1409310185087476488[3] = 0;
   out_1409310185087476488[4] = 1;
   out_1409310185087476488[5] = 0;
   out_1409310185087476488[6] = 0;
   out_1409310185087476488[7] = 0;
   out_1409310185087476488[8] = 0;
   out_1409310185087476488[9] = 0;
   out_1409310185087476488[10] = 0;
   out_1409310185087476488[11] = 0;
   out_1409310185087476488[12] = 0;
   out_1409310185087476488[13] = 1;
   out_1409310185087476488[14] = 0;
   out_1409310185087476488[15] = 0;
}
void h_30(double *state, double *unused, double *out_2429398138357197661) {
   out_2429398138357197661[0] = state[4];
}
void H_30(double *state, double *unused, double *out_924817957936797882) {
   out_924817957936797882[0] = 0;
   out_924817957936797882[1] = 0;
   out_924817957936797882[2] = 0;
   out_924817957936797882[3] = 0;
   out_924817957936797882[4] = 1;
   out_924817957936797882[5] = 0;
   out_924817957936797882[6] = 0;
   out_924817957936797882[7] = 0;
}
void h_26(double *state, double *unused, double *out_3599589183117062072) {
   out_3599589183117062072[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5430464877248331228) {
   out_5430464877248331228[0] = 0;
   out_5430464877248331228[1] = 0;
   out_5430464877248331228[2] = 0;
   out_5430464877248331228[3] = 0;
   out_5430464877248331228[4] = 0;
   out_5430464877248331228[5] = 0;
   out_5430464877248331228[6] = 0;
   out_5430464877248331228[7] = 1;
}
void h_27(double *state, double *unused, double *out_5324631439891927468) {
   out_5324631439891927468[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7574479069461655864) {
   out_7574479069461655864[0] = 0;
   out_7574479069461655864[1] = 0;
   out_7574479069461655864[2] = 0;
   out_7574479069461655864[3] = 1;
   out_7574479069461655864[4] = 0;
   out_7574479069461655864[5] = 0;
   out_7574479069461655864[6] = 0;
   out_7574479069461655864[7] = 0;
}
void h_29(double *state, double *unused, double *out_7901847793471347442) {
   out_7901847793471347442[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1364363635470041804) {
   out_1364363635470041804[0] = 0;
   out_1364363635470041804[1] = 1;
   out_1364363635470041804[2] = 0;
   out_1364363635470041804[3] = 0;
   out_1364363635470041804[4] = 0;
   out_1364363635470041804[5] = 0;
   out_1364363635470041804[6] = 0;
   out_1364363635470041804[7] = 0;
}
void h_28(double *state, double *unused, double *out_751004511462432978) {
   out_751004511462432978[0] = state[5];
   out_751004511462432978[1] = state[6];
}
void H_28(double *state, double *unused, double *out_8869473968987790654) {
   out_8869473968987790654[0] = 0;
   out_8869473968987790654[1] = 0;
   out_8869473968987790654[2] = 0;
   out_8869473968987790654[3] = 0;
   out_8869473968987790654[4] = 0;
   out_8869473968987790654[5] = 1;
   out_8869473968987790654[6] = 0;
   out_8869473968987790654[7] = 0;
   out_8869473968987790654[8] = 0;
   out_8869473968987790654[9] = 0;
   out_8869473968987790654[10] = 0;
   out_8869473968987790654[11] = 0;
   out_8869473968987790654[12] = 0;
   out_8869473968987790654[13] = 0;
   out_8869473968987790654[14] = 1;
   out_8869473968987790654[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
