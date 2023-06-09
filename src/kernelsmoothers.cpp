#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Uniform kernel, unscaled (with roughness = int x^2 f(x) dx = 1/3)
arma::vec kuniform2(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i = 0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 1.0) k = 0.5;
    ax[i] = k;
  }
  return ax;
}

// Uniform kernel convolution
arma::vec kuniform2conv(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 2.0) k = 0.5 - 0.25*ax[i];
    ax[i] = k;
  }
  return ax;
}

// 4th-order uniform
arma::vec kuniform4(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i = 0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] <= 0.5)
      k =  1.16666666666666667;
    else if (ax[i] < 1.0)
      k = -1.66666666666666667e-1;
    ax[i] = k;
  }
  return ax;
}

// 4th-order uniform convolution
// x|-->-65/36*x + 25/18 on (0, 1/2], x|-->-49/36*x + 7/6 on (1/2, 1], x|-->5/12*x - 11/18 on (1, 3/2], x|-->-1/36*x + 1/18 on (3/2, 2]; x)
arma::vec kuniform4conv(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i = 0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] <= 0.5)
      k = -1.80555555555555556*ax[i] + 1.38888888888888888;
    else if (ax[i] <= 1.0)
      k = -1.36111111111111111*ax[i] + 1.16666666666666667;
    else if (ax[i] <= 1.5)
      k =  4.16666666666666667e-1*ax[i] - 6.11111111111111111e-1;
    else if (ax[i] < 2.0)
      k = -2.77777777777777778e-2*ax[i] + 5.55555555555555556e-2;
    ax[i] = k;
  }
  return ax;
}

// 6th-order uniform
// arma::vec kuniform6(arma::vec x) {
//   arma::vec ax = arma::abs(x);
//   for (arma::uword i = 0; i < ax.n_elem; i++) {
//     double k = 0.0;
//     if (ax[i] <= 0.33333333333333333)
//       k = 1.85;
//     else if (ax[i] <= 0.66666666666666667)
//       k = -0.4;
//     else if (ax[i] < 1.0)
//       k = 0.05;
//     ax[i] = k;
//   }
//   return ax;
// }

// 6th-order uniform convolution
// arma::vec kuniform6conv(arma::vec x) {
//   arma::vec ax = arma::abs(x);
//   for (arma::uword i = 0; i < ax.n_elem; i++) {
//     double k = 0.0;
//     if (ax[i] <= 0.33333333333333333)
//       k = -5.2675*ax[i] + 2.39;
//     else if (ax[i] <= 0.66666666666666667)
//       k = -3.19745*ax[i] + 1.7;
//     else if (ax[i] <= 1.0)
//       k =  1.64*ax[i] - 1.525;
//     else if (ax[i] < 1.33333333333333333)
//       k = -0.385*ax[i] + 0.5;
//     else if (ax[i] < 1.66666666666666667)
//       k = 0.0425*ax[i] - 0.07;
//     else if (ax[i] < 2.0)
//       k = -0.0025*ax[i] + 0.005;
//     ax[i] = k;
//   }
//   return ax;
// }

// Triangular kernel, unscaled (with roughness = int x^2 f(x) dx = 1/6)
arma::vec ktriangular2(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 1.0) k = 1.0 - ax[i];
    ax[i] = k;
  }
  return ax;
}

// Triangular kernel convolution (2nd order)
arma::vec ktriangular2conv(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 2.0) {
      if (ax[i] < 1.0)
        k = ax[i]*ax[i]*(0.5*ax[i] - 1.0) + 0.666666666666666667;
      else {
        double xm2 = ax[i] - 2.0;
        k = -0.1666666666666666667*xm2*xm2*xm2;
      }
    }
    ax[i] = k;
  }
  return ax;
}

// Triangular 4th-order, piecewise linear
// y = a-b|x|) if |x|<(a+b)/2b and -b+b*|x| otherwise (intersection at (a-b)/2a)
//  b=sqrt(2)+1
arma::vec ktriangular4(arma::vec x) {
  double b = 2.41421356237309492343; // sqrt(2)+1
  double a = 1.64599349813977635648; // -b + sqrt(2b*(b+1))
  double c = 0.840896415253714613058; // (a+b)/2/b = sqrt((b+1)/2/b)
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < c)
      k = a-b*ax[i];
    else if (ax[i] < 1)
      k = -b + b*ax[i];
    ax[i] = k;
  }
  return ax;
}

// Triangular kernel convolution (4th order)
// Computed in Sage:
// f4 = piecewise([[(-1, 0), (1+x)*(12/7 - 30/7*x^2)], [(0, 1), (1-x)*(12/7 - 30/7*x^2)]])
// g4 = f4.convolution(f4); g4
arma::vec ktriangular4conv(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 2.0) {
      double y = ax[i];
      double y2 = y*y;
      if (y < 0.1591035847462854570)
        k = 1.2627507209715991835 + 0.97140452079103168293*y2*(-6.0+7.0*y);
      else if (y < 0.84089641525371454303)
        k = 1.2784002043808324125 - 0.29508103354532229150*y - 3.9737798267869814275*y2 + 2.9142135623730950488*y2*y;
      else if (y < 1.0)
        k = 5.8992048750628466607 - 16.780362407783892584*y + 15.630634076279361623*y2 - 4.8570226039551584147*y2*y;
      else if (y < 1.6817928305074290861)
        k = 2.0135867918987199289 -5.1235081582915123891*y + 3.9737798267869814275*y2 - 0.97140452079103168293*y2*y;
      else if (y < 1.8408964152537145430)
        k = -16.469631890829337064 + 27.847054590185628197*y -15.630634076279361623*y2 +2.9142135623730950488*y2*y;
      else {
        double ym2 = y - 2.0;
        k = -0.97140452079103168293 * ym2*ym2*ym2;
      }
    }
    ax[i] = k;
  }
  return ax;
}

// arma::vec ktriangular6(arma::vec x) {
//   arma::vec ax = arma::abs(x);
//   for (arma::uword i=0; i < ax.n_elem; i++) {
//     double k = 0.0;
//     if (ax[i] < 1.0) {
//       double x2 = ax[i]*ax[i];
//       k = (1.0 - ax[i])*(2.39385065885797932 -15.3733528550512446*x2 + 17.5256222547584173*x2*x2);
//     }
//     ax[i] = k;
//   }
//   return ax;
// }

// Triangular kernel convolution (6th order)
// f6 = piecewise([[(-1, 0), (1+x)*(1635/683 - 10500/683*x^2 + 11970/683*x^4)], [(0, 1), (1-x)*(1635/683 - 10500/683*x^2 + 11970/683*x^4)]])
// g6 = f6.convolution(f6); g6
// arma::vec ktriangular6conv(arma::vec x) {
//   arma::vec ax = arma::abs(x);
//   for (arma::uword i=0; i < ax.n_elem; i++) {
//     double k = 0.0;
//     if (ax[i] < 2.0) {
//       double y = ax[i];
//       double y2 = y*y;
//       if (ax[i] < 1.0)
//         k = (((((0.332410644390133736*y -0.48753561177219612)*y2 -3.69500674185243394*y  +5.13195380812838042)*y2 +16.1897708198907146*y -24.6594882194435456)*y2 -18.1743192229613122*y  +33.7205378904968853)*y2 +5.354708256786334*y  -15.688312050230552)*y2 + 1.97780167865207379;
//       else
//         k = (((((0.110803548130044574*y +0.48753561177219612)*y2 +0.581621431587883153*y -5.13195380812838042)*y2 +1.44601480420760198*y +24.6594882194435456)*y2 -38.5049593881099028*y  +1.33070661901995546)*y2 +32.2810666489456324*y -15.0583936598719372)*y2 -4.78770131771595864*y + 2.80989963906388507;
//     }
//     ax[i] = k;
//   }
//   return ax;
// }


// Epanechnikov kernel, unscaled (with roughness = int x^2 f(x) dx = 1/5)
arma::vec kepanechnikov2(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 1.0) k = 0.75*(1.0 - ax[i]*ax[i]);
    ax[i] = k;
  }
  return ax;
}

// Epanechnikov kernel convolution
arma::vec kepanechnikov2conv(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 2.0) {
      double y = ax[i];
      double ym2 = 2.0 - y;
      double ym2_3 = ym2*ym2*ym2;
      k = 0.01875 * ym2_3 * (y*(y + 6) + 4);
    }
    ax[i] = k;
  }
  return ax;
}

arma::vec kepanechnikov4(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 0.402102389929217485243) {
      k = 1.51301275841878699284 - 5.30812951240326835745*ax[i]*ax[i];
    } else if (ax[i] < 1) {
      double xm = ax[i] - 0.804204779858434970485;
      k = 5.30812951240326835745*xm*xm - 0.203491222723821640894;
    }
    ax[i] = k;
  }
  return ax;
}

// Epanechnikov kernel convolution (4th order)
arma::vec kepanechnikov4conv(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 2.0) {
      double y = ax[i];
      double y2 = y*y;
      double y3 = y2*y;
      if (y < 0.59789761007078252799)
        k = 0.93920796401488528437*(-216.36559038442943064-40.0*y2+20.0*y3-5.0*y2*y3 + 5.2011568408621214332*(1.0+3.0*y2)+12.934906558943076212*(4.0+(-6.0+y)*y2)+8.0420477985843494401*(20.0+y2*(12.0+(-4.0+y)*y)));
      else if (y < 0.80420477985843494401)
        k = 0.23480199100372132109*(-2.0913975961073576591*(8.0+7.0*y)-25.869813117886152424*(-14.0+15.0*y+y3)+10.402313681724242866*(-5.0+y*(12.0+y))+64.336382388674795521*(3.0+y*(8.0+y*(-3.0+2.0*y)))-4.0*(120.0+y*(60.0+y*(-40.0+20.0*y+y3))));
      else if (y < 1.4021023899292174720)
        k = 0.23480199100372132109*(544.0-10.402313681724242866*(-1.0+y)*(-5.0+7.0*y)+2.0913975961073576591*(-8.0+9.0*y)+25.869813117886152424*(14.0+3.0*y*(-5.0+y2))+4.0*y*(-60.0+y*(40.0-20.0*y+3.0*y3))-64.336382388674795521*(13.0+y*(-8.0+y*(3.0+(-2.0+y)*y))));
      else {
        double ym2 = y - 2.0;
        k = 0.93920796401488528437*ym2*ym2*ym2*(-16.934906558943076212-6.0*y-y2+8.0420477985843494401*(2.0+y));
      }
    }
    ax[i] = k;
  }
  return ax;
}


// arma::vec kepanechnikov6(arma::vec x) {
//   arma::vec ax = arma::abs(x);
//   for (arma::uword i=0; i < ax.n_elem; i++) {
//     double k = 0.0;
//     if (ax[i] < 1.0) {
//       double x2 = ax[i]*ax[i];
//       k = 0.75*(1.0 - x2)*(2.734375 -16.406250*x2 + 18.046875*x2*x2);
//     }
//     ax[i] = k;
//   }
//   return ax;
// }

// Epanechnikov kernel convolution (6th order)
// f6 = piecewise([[(-1, 1), 0.75*(1.0 - x^2)*(2.734375 -16.406250*x^2 + 18.046875*x^4)]])
// g6 = f6.convolution(f6); g6
// arma::vec kepanechnikov6conv(arma::vec x) {
//   arma::vec ax = arma::abs(x);
//   for (arma::uword i=0; i < ax.n_elem; i++) {
//     double k = 0.0;
//     if (ax[i] < 2.0) {
//       double y = ax[i];
//       double y2 = y*y;
//       k = y2*(y*(y2*(y2*(y2*(y2*(-0.0152514531062199511*y2 + 0.3028106689453125) -2.6019287109375) +14.996337890625) -13.53515625*y -19.73876953125) +25.83984375*y +7.177734375) -14.35546875) +1.89302884615384626;
//     }
//     ax[i] = k;
//   }
//   return ax;
// }

// Quartic kernel, unscaled (with roughness = int x^2 f(x) dx = 1/7)
arma::vec kquartic2(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 1.0) {
      double x2 = ax[i]*ax[i];
      double y = 1.0 - x2;
      k = 0.9375 * y*y;
    }
    ax[i] = k;
  }
  return ax;
}

// Quartic kernel convolution
arma::vec kquartic2conv(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 2.0) {
      double y = ax[i];
      double y2 = y*y;
      double ym2 = 2.0 - y;
      double ym2_2 = ym2*ym2;
      double ym2_5 = ym2_2*ym2_2*ym2;
      k = 0.0013950892857142857 * ym2_5 * (16.0 + y*(2.0 + y)*(20.0 + 8.0*y + y2));
    }
    ax[i] = k;
  }
  return ax;
}

arma::vec kquartic4(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 1.0) {
      double x2 = ax[i]*ax[i];
      double y = 1.0 - x2;
      k = 0.9375 * y*y * (1.75 - 5.25*x2);
    }
    ax[i] = k;
  }
  return ax;
}

// Quartic kernel convolution (4th order)
arma::vec kquartic4conv(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 2.0) {
      double y = ax[i];
      double y2 = y*y;
      k = y2*(y2*(y2*(y*(y2*(y2*(-0.00201672107189685302*y2 + 0.0489390980113636395) + -0.52978515625) + 4.1015625) + -4.921875) + -5.7421875*y + 11.484375) + -5.96590909090909083) + 1.40734265734265729;
    }
    ax[i] = k;
  }
  return ax;
}

// arma::vec kquartic6(arma::vec x) {
//   arma::vec ax = arma::abs(x);
//   for (arma::uword i=0; i < ax.n_elem; i++) {
//     double k = 0.0;
//     if (ax[i] < 1.0) {
//       double x2 = ax[i]*ax[i];
//       double y = 1.0 - x2;
//       k = 0.9375 * y*y * (2.4609375 -18.0468750*x2 + 23.4609375*x2*x2);
//     }
//     ax[i] = k;
//   }
//   return ax;
// }

// Quartic kernel convolution (6th order)
// arma::vec kquartic6conv(arma::vec x) {
//   arma::vec ax = arma::abs(x);
//   for (arma::uword i=0; i < ax.n_elem; i++) {
//     double k = 0.0;
//     if (ax[i] < 2.0) {
//       double y = ax[i];
//       double y2 = y*y;
//       double y4 = y2*y2;
//       k = y2*(y2*(y2*(y*(y2*(y2*(y2*(y2*(-0.00221108689027674048*y2 + 0.0594806671142578125) + -0.706280928391676666) + 4.97955322265625) + -27.05108642578125) + 21.99462890625*y + 47.757568359375) + -60.908203125) + -29.06982421875*y + 58.1396484375) + -17.2265625) + 2.07119626696832571;
//     }
//     ax[i] = k;
//   }
//   return ax;
// }


// Gaussian kernel unscaled (with roughness 1)
arma::vec kgaussian2(arma::vec x) {
  return 0.39894228040143268 * arma::exp(-0.5 * (x % x));
}

// Gaussian kernel convolution
arma::vec kgaussian2conv(arma::vec x) {
  return 0.28209479177387814 * arma::exp(-0.25 * (x % x));
}

arma::vec kgaussian4(arma::vec x) {
  arma::vec x2 = x % x;
  return 0.39894228040143268 * arma::exp(-0.5 * (x2)) * (1.5 - 0.5*x2);
}

arma::vec kgaussian4conv(arma::vec x) {
  arma::vec x2 = x % x;
  arma::vec x4 = x2 % x2;
  return 0.28209479177387814 * arma::exp(-0.25 * x2) * (1.6875 -0.4375*x2 + 0.015625*x4);
}

// arma::vec kgaussian6(arma::vec x) {
//   arma::vec x2 = x % x;
//   arma::vec x4 = x2 % x2;
//   return 0.39894228040143268 * arma::exp(-0.5 * (x2)) * (1.875 -1.25*x2 + 0.125*x4);
// }

// arma::vec kgaussian6conv(arma::vec x) {
//   arma::vec x2 = x % x;
//   arma::vec x4 = x2 % x2;
//   arma::vec x6 = x4 % x2;
//   arma::vec x8 = x4 % x4;
//   return 0.28209479177387814 * arma::exp(-0.25 * x2) * (2.2119140625 -1.181640625*x2 + 0.14111328125*x4 -0.00537109375*x6 + 0.00006103515625*x8);
// }

arma::vec kernelOrder(arma::vec x, std::string kernel, int order, bool convolution = false) {
  if (!convolution) {
    if (kernel == "gaussian") {
      if (order == 2)
        return kgaussian2(x);
      else if (order == 4)
        return kgaussian4(x);
      // else if (order == 6)
      //   return kgaussian6(x);
    } else if (kernel == "triangular") {
      if (order == 2)
        return ktriangular2(x);
      else if (order == 4)
        return ktriangular4(x);
      // else if (order == 6)
      //   return ktriangular6(x);
    } else if (kernel == "epanechnikov") {
      if (order == 2)
        return kepanechnikov2(x);
      else if (order == 4)
        return kepanechnikov4(x);
      // else if (order == 6)
      //   return kepanechnikov6(x);
    } else if (kernel == "quartic") {
      if (order == 2)
        return kquartic2(x);
      else if (order == 4)
        return kquartic4(x);
      // else if (order == 6)
      //   return kquartic6(x);
    } else if (kernel == "uniform") {
      if (order == 2)
        return kuniform2(x);
      else if (order == 4)
        return kuniform4(x);
      // else if (order == 6)
      //   return kuniform6(x);
    }
  } else { // Convolution
    if (kernel == "gaussian") {
      if (order == 2)
        return kgaussian2conv(x);
      else if (order == 4)
        return kgaussian4conv(x);
      // else if (order == 6)
      //   return kgaussian6conv(x);
    } else if (kernel == "triangular") {
      if (order == 2)
        return ktriangular2conv(x);
      else if (order == 4)
        return ktriangular4conv(x);
      // else if (order == 6)
      //   return ktriangular6conv(x);
    } else if (kernel == "epanechnikov") {
      if (order == 2)
        return kepanechnikov2conv(x);
      else if (order == 4)
        return kepanechnikov4conv(x);
      // else if (order == 6)
      //   return kepanechnikov6conv(x);
    } else if (kernel == "quartic") {
      if (order == 2)
        return kquartic2conv(x);
      else if (order == 4)
        return kquartic4conv(x);
      // else if (order == 6)
      //   return kquartic6conv(x);
    } else if (kernel == "uniform") {
      if (order == 2)
        return kuniform2conv(x);
      else if (order == 4)
        return kuniform4conv(x);
      // else if (order == 6)
      //   return kuniform6conv(x);
    }
  }
  throw std::runtime_error("kernelOrder: Invalid kernel type (should be gaussian, uniform, triangular, epanechnikov, or quartic) or order (should be 2 or 4).");
  return x;
}

// In the innermost loop, using ksum += dnorm((xgrid[i] - x[j]) / bw, 0, 1, 0) gives terrible performance
// [[Rcpp::export]]
arma::mat kernelWeightsOneCPP(arma::vec x, arma::vec xgrid, double bw, std::string kernel = "gaussian", int order = 2, bool convolution = false) {
  arma::uword n = xgrid.n_elem;
  arma::uword nx = x.n_elem;
  arma::mat kw(n, nx);
  arma::vec xs = x / bw; // Scaling by the bandwidth
  arma::vec gs = xgrid / bw;

  for (arma::uword i=0; i < n; i++) {
    kw.row(i) = kernelOrder(gs[i] - xs, kernel, order, convolution).t();
  }

  return kw;
}

// [[Rcpp::export]]
arma::mat kernelWeightsCPP(arma::mat x, arma::mat xgrid, arma::vec bw, std::string kernel = "gaussian", int order = 2, bool convolution = false) {
  arma::uword d = x.n_cols;
  // The product kernel matrix starts with the first dimension (there is at least one column or row)
  arma::mat pk = kernelWeightsOneCPP(x.col(0), xgrid.col(0), bw[0], kernel);
  if (d > 1) { // We need to compute the product kernel starting from the 2nd till the last dimension
    for (arma::uword k=1; k < d; k++) {
      pk %= kernelWeightsOneCPP(x.col(k), xgrid.col(k), bw[k], kernel, order, convolution);
    }
  }
  return pk;
}

// [[Rcpp::export]]
Rcpp::NumericVector kernelDensityCPP(arma::mat x, arma::mat xgrid, arma::vec bw, std::string kernel = "gaussian", int order = 2, bool convolution = false) {
  arma::uword d = x.n_cols;
  double nb = (double)x.n_rows; // n*prod(b) in the denominator
  for (arma::uword k=0; k < d; k++) {
    nb *= bw[k];
  }

  arma::mat kw = kernelWeightsCPP(x, xgrid, bw, kernel, order);
  arma::vec out = sum(kw, 1) / nb;
  Rcpp::NumericVector rout = Rcpp::NumericVector(out.begin(), out.end());
  return rout;
}

// [[Rcpp::export]]
Rcpp::NumericVector kernelSmoothCPP(arma::mat x, arma::vec y, arma::mat xgrid, arma::vec bw, std::string kernel = "gaussian", int order = 2, bool LOO = false, bool convolution = false) {
  arma::mat kw = kernelWeightsCPP(x, xgrid, bw, kernel, order, convolution);

  // LOO: setting diagonal elements to zero, assuming x = xgrid! (R makes sure it happens, though.)
  if (LOO) {
    for (arma::uword i=0; i < x.n_rows; i++) {
      kw(i, i) = 0;
    }
  }

  arma::vec ksum = sum(kw, 1); // Nadaraya--Watson denominator
  kw.each_row() %= y.t(); // Nadaraya--Watson numerator: y_i * w_ij (in place to save memory)
  arma::vec ysum = sum(kw, 1);
  arma::vec out = ysum / ksum;
  Rcpp::NumericVector rout = Rcpp::NumericVector(out.begin(), out.end());
  return rout;
}

