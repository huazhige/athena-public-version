#ifndef INTERPOLATION_HPP_
#define INTERPOLATION_HPP_

#define SQR(x) ( (x)*(x) )

/* 2-nd order plm for non-uniform grid
template<typename T>
inline T interp_plm(T phim1, T phi, T phip1, Real dxl, Real dx) {
  Real dwl = phi - phim1;
  Real dwr = phip1 - phi;
  Real dw2 = dwl*dwr;
  Real dwm = 2.0*dw2/(dwl + dwr);
  if (dw2 <= 0.0) dwm = 0.0;
  return phi - dxl/dx*dwm;
}*/

// limiter
template<typename T>
inline T minmod(T a, T b) {
  return a*std::max(0., std::min(1., b/a));
}

template<typename T>
inline T superbee(T a, T b) {
  Real r = b/a;
  return a*std::max(0., std::min(1., 2.*r), std::min(2., r));
}

template<typename T>
inline T vanleer(T a, T b) {
  Real r = b/a;
  return a*(r + fabs(r))/(1. + fabs(r));
}

template<typename T>
inline T mclimiter(T a, T b) {
  Real r = b/a;
  Real c = (1. + r)/2.;
  return std::max(0., std::min(std::min(c, 2.), 2.*r));
}

// 2-rd polynomial interpolation
template<typename T>
inline T interp_bp2(T phi, T phip1) {
  return 1.5*phi - 0.5*phip1;
}

template<typename T>
inline T interp_cp2(T phim1, T phi) {
  return 0.5*phim1 + 0.5*phi;
}

// 3-rd polynomial interpolation
template<typename T>
inline T interp_cp3(T phim1, T phi, T phip1) {
  return 1./3.*phim1 + 5./6.*phi - 1./6.*phip1;
};

template<typename T>
inline T interp_bp3(T phi, T phip1, T phip2) {
  return 11./6.*phi - 7./6.*phip1 + 1./3.*phip2;
};

// WENO 3 interpolation
template<typename T>
inline T interp_weno3(T phim1, T phi, T phip1) {
    T p0 = (-1.0/2.0) * phim1 + (3.0/2.0) * phi;
    T p1 = (1.0/2.0) * phi + (1.0/2.0) * phip1;

    T beta1 = (phip1 - phi) * (phip1 - phi);
    T beta0 = (phi - phim1) * (phi - phim1);

    T alpha0 = (1.0/3.0) /((beta0 + 1e-10) * (beta0 + 1.0e-10));
    T alpha1 = (2.0/3.0)/((beta1 + 1e-10) * (beta1 + 1.0e-10));

    T alpha_sum_inv = 1.0/(alpha0 + alpha1);

    T w0 = alpha0*alpha_sum_inv;
    T w1 = alpha1*alpha_sum_inv;

    return w0*p0 + w1*p1;
};

// 4-th polynomial interpolation
template<typename T>
inline T interp_cp4(T phim2, T phim1, T phi, T phip1) {
  return -1./12*phim2 + 7./12.*phim1 + 7./12.*phi - 1./12.*phip1;
};

// 5-th polynomial interpolation
template<typename T>
inline T interp_cp5(T phim2, T phim1, T phi, T phip1, T phip2) {
  return -1./20*phim2 + 9./20.*phim1 + 47./60.*phi - 13./60.*phip1 + 1./30.*phip2;
};

// WENO 5 interpolation
template<typename T>
inline T interp_weno5(T phim2, T phim1, T phi, T phip1, T phip2) {
  T p0 = (1./3.)*phim2 - (7./6.)*phim1 + (11./6.)*phi;
  T p1 = (-1./6.)*phim1 + (5./6.)*phi + (1./3.)*phip1;
  T p2 = (1./3.)*phi + (5./6.)*phip1 - (1./6.)*phip2;

  T beta0 = 13./12.*SQR(phim2 - 2.*phim1 + phi) + .25*SQR(phim2 - 4.*phim1 + 3.*phi);
  T beta1 = 13./12.*SQR(phim1 - 2.*phi + phip1) + .25*SQR(phim1 - phip1);
  T beta2 = 13./12.*SQR(phi - 2.*phip1 + phip2) + .25*SQR(3.*phi - 4.*phip1 + phip2);

  T alpha0 = .1/SQR(beta0 + 1e-6);
  T alpha1 = .6/SQR(beta1 + 1e-6);
  T alpha2 = .3/SQR(beta2 + 1e-6);

  return (alpha0*p0 + alpha1*p1 + alpha2*p2)/(alpha0 + alpha1 + alpha2);
};

// 3rd order polynomial with inflection point
template<typename T>
inline T inflection3_cell1(T f1, T f2, T f3) {
  return 7./3.*f1 - 5./3.*f2 + 1./3.*f3;
}

template<typename T>
inline T inflection3_cell2(T f1, T f2, T f3) {
  return 10./3.*f1 - 8./3.*f2 + 1./3.*f3;
}

template<typename T>
inline T inflection3_cell3(T f1, T f2, T f3) {
  return 10./3.*f1 - 5./3.*f2 - 2./3.*f3;
}

#undef SQR

#endif
