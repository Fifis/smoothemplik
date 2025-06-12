#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "kernelsm.h"
using namespace Rcpp;
using namespace RcppParallel;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom RcppParallel RcppParallelLibs

// All auxiliary computations were done in Mathematica 14 and Sage

// Uniform kernel, unscaled (with roughness = int x^2 f(x) dx = 1/3)
arma::vec kuniform2(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) { val = (val <= 1.0) ? 0.5 : 0.0; });
  return ax;
}

// Uniform kernel convolution
// u2 = piecewise([[(-1, 1), 1/2]])
// u2c = u2.convolution(u2); u2c
// piecewise(x|-->1/4*x + 1/2 on (-2, 0], x|-->-1/4*x + 1/2 on (0, 2]; x)
arma::vec kuniform2conv(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    val = (val < 2.0) ? 0.5 - 0.25 * val : 0.0;
  });
  return ax;
}

// 4th-order uniform
arma::vec kuniform4(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    if (val < 0.666666666666666667)
      val = 0.95;
    else if (val <= 1.0)
      val = -0.4;
    else
      val = 0.0;
  });
  return ax;
}

// 4th-order uniform convolution
// f1[x_] := 19/20*Boole[RealAbs[x] <= 2/3] - 2/5*Boole[RealAbs[x] > 2/3 && RealAbs[x] <= 1]
// fc1[x_] := Simplify[Integrate[f1[t]*f1[x - t], {t, -1, 1}]]
arma::vec kuniform4conv(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    if (val <= 0.333333333333333333)
      val = 1.31 - 1.9825*val;
    else if (val <= 1.33333333333333333)
      val = 0.95 - 0.9025*val;
    else if (val < 1.66666666666666667)
      val = -1.48 + 0.92*val;
    else if (val <= 2.0)
      val = 0.32 - 0.16*val;
    else
      val = 0.0;
  });
  return ax;
}

// Triangular kernel, unscaled (with roughness = int x^2 f(x) dx = 1/6)
arma::vec ktriangular2(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    val = (val < 1.0) ? 1.0 - val : 0.0;
  });
  return ax;
}

// Triangular kernel convolution (2nd order)
// t2 = piecewise([[(-1, 0), (1+x)], [(0, 1), (1-x)]])
// t2c = t2.convolution(t2); t2c
// x|-->1/2*x^3 - x^2 + 2/3 on (0, 1], x|-->-1/6*x^3 + x^2 - 2*x + 4/3 on (1, 2]; x)
arma::vec ktriangular2conv(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    double k = 0.0;
    if (val < 2.0) {
      if (val < 1.0)
        k = val * val * (0.5 * val - 1.0) + 0.66666666666666666667;
      else {
        double xm2 = val - 2.0;
        k = -0.166666666666666666667 * xm2 * xm2 * xm2;
      }
    }
    val = k;
  });
  return ax;
}

// Triangular 4th-order, piecewise linear
// y = a-b|x| if |x|<(a+b)/2b and -b+b*|x| otherwise (intersection at (a-b)/2a)
//  b=sqrt(2)+1
arma::vec ktriangular4(const arma::vec& x) {
  const double a = 1.64599349813977617069556; // -b + sqrt(2b*(b+1))
  const double b = 2.41421356237309504880169; // sqrt(2)+1
  const double t = 0.84089641525371454303113; // (a+b)/2/b = sqrt((b+1)/2/b)

  arma::vec ax = arma::abs(x);
  ax.for_each([a, b, t](arma::vec::elem_type& val) {
    if (val < t)
      val = a - b * val;
    else if (val < 1)
      val = -b + b * val;
    else
      val = 0.0;
  });
  return ax;
}

// Triangular kernel convolution (4th order)
// t4 = piecewise([[(-1, 0), (1+x)*(12/7 - 30/7*x^2)], [(0, 1), (1-x)*(12/7 - 30/7*x^2)]])
// t4c = t4.convolution(t4); t4c
// However, this code was obtained in Mathematica in a convenient Horner form
// Avoiding mutually exclusive checks yields a 3x speed-up, so this is a MORE efficient solution
// because there are more correct CPU guesses and a mispredicted branch flushes the pipeline
arma::vec ktriangular4conv(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& v) {
    if (v < 2.0) {

      double seg1 =  1.2627507209715991834795    + v*v*(-5.828427124746190097603 + 6.799831645537221780537*v);
      double seg2 =  1.2784002043808324124710    + v*(-0.29508103354532229149832 + v*(-3.9737798267869814274534 + 2.9142135623730950488017*v));
      double seg3 =  5.89920487506284666066      + v*(-16.78036240778389258431 + (15.63063407627936162266 - 4.85702260395515841467*v)*v);
      double seg4 =  2.0135867918987199289201    + v*(-5.123508158291512389102 + (3.9737798267869814274534 - 0.9714045207910316829339*v)*v);
      double seg5 = -16.46963189082933706382     + v*(27.84705459018562819652 + v*(-15.63063407627936162266 + 2.914213562373095048802*v));
      double seg6 =  7.7712361663282534634711699 + v*(-11.6568542494923801952067549 + (5.82842712474619009760337745 - 0.971404520791031682933896241*v)*v);

      v = seg1 * (v <  0.159103584746285456968875) +
        seg2 * (v >= 0.159103584746285456968875) * (v < 0.84089641525371454303113) +
        seg3 * (v >= 0.84089641525371454303113)  * (v < 1.0) +
        seg4 * (v >= 1.00000000000000000000000)  * (v < 1.68179283050742908606225) +
        seg5 * (v >= 1.68179283050742908606225)  * (v < 1.84089641525371454303113) +
        seg6 * (v >= 1.84089641525371454303113);
    } else {
      v = 0.0; // Set to zero if v is >= 2.0
    }
  });
  return ax;
}


// Epanechnikov kernel, unscaled (with roughness = int x^2 f(x) dx = 1/5)
arma::vec kepanechnikov2(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    val = (val < 1.0) ? 0.75 * (1.0 - val * val) : 0.0;
  });
  return ax;
}

// Epanechnikov kernel convolution
// e2 = piecewise([[(-1, 1), 3/4*(1-x^2)]])
// e2c = e2.convolution(e2); e2c
arma::vec kepanechnikov2conv(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    if (val < 2.0) {
      double ym2 = 2.0 - val;
      val = 0.01875 * ym2 * ym2 * ym2 * (val * (val + 6) + 4);
    } else {
      val = 0.0;
    }
  });
  return ax;
}

arma::vec kepanechnikov4(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    if (val < 0.402102389929217485243) {
      val = 1.5130127584187869928425 - 5.30812951240326826778399 * val * val;
    } else if (val < 1) {
      val = 3.2295167395613956861352 - 8.5376462519646639539192 * val + 5.30812951240326826778399 * val * val;
    } else {
      val = 0.0;
    }
  });
  return ax;
}

// Epanechnikov kernel convolution (4th order)
// e4 = piecewise([[(-1, -0.402102389929217472006029),                         3.2295167395613956861352 + 8.5376462519646639539192*x + 5.30812951240326826778399*x^2],
//                [(-0.402102389929217472006029, 0.402102389929217472006029), 1.5130127584187869928425 - 5.30812951240326826778399*x^2],
//                [(0.402102389929217472006029, 1),                           3.2295167395613956861352 - 8.5376462519646639539192*x + 5.30812951240326826778399*x^2]])
//  e4c = e4.convolution(e4); e4c
arma::vec kepanechnikov4conv(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 2.0) {
      double y = ax[i];
      double y2 = y*y;
      double y3 = y2*y;
      double y4 = y2*y2;
      double y5 = y2*y3;
      if (y < 0.59789761007078252799)
        k = 1.3300581033563908047756 -  5.1669542301715031721552*y2 + 0.7201051765702296976554*y3 +  7.553155339418797088324*y4 - 4.6960398200744264218529*y5;
      else if (y < 0.80420477985843494401)
        k = 1.5130127584187869928425 - 0.7438894894558010581981*y -  5.3081295124032682677840*y2 + 5.3541784504036264841143*y3 -  0.9392079640148852843706*y5;
      else if (y < 1.4021023899292174720)
        k = 0.24928068482791296534562 + 7.113139623286665983254*y -  24.848001220516091366474*y2 + 29.651312958299051211193*y3 -  15.106310678837594176648*y4 + 2.8176238920446558531117*y5;
      else {
        k = 6.3927074493573063443 - 24.932321616671849800*y +  35.323084963090895541*y2 - 23.577029331325203325*y3 + 7.5531553394187983377*y4 - 0.93920796401488549078*y5;
      }
    }
    ax[i] = k;
  }
  return ax;
}


// Quartic kernel, unscaled (with roughness = int x^2 f(x) dx = 1/7)
arma::vec kquartic2(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    if (val < 1.0) {
      double y = 1.0 - val * val;
      val = 0.9375 * y * y;
    } else {
      val = 0.0;
    }
  });
  return ax;
}

// Quartic kernel convolution
arma::vec kquartic2conv(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    if (val < 2.0) {
      double y = val;
      double y2 = y * y;
      double ym2 = 2.0 - y;
      double ym2_2 = ym2 * ym2;
      double ym2_5 = ym2_2 * ym2_2 * ym2;
      val = 0.0013950892857142857 * ym2_5 * (16.0 + y * (2.0 + y) * (20.0 + 8.0 * y + y2));
    } else {
      val = 0.0;
    }
  });
  return ax;
}

arma::vec kquartic4(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (arma::uword i=0; i < ax.n_elem; i++) {
    double k = 0.0;
    if (ax[i] < 1.0) {
      double y = ax[i];
      if (y < 0.78473987365599073452) {
        k = 4.6563049882118414382631 * (-0.9245053574195424742574 + y) * (-0.6139406996720486380031 + y) * (0.6139406996720486380031 + y) * (0.9245053574195424742574 + y);
      } else {
        k = -123.7643807892723135558 * (-1 + y) * (-1 + y) * (-0.56947974731198146905 + y) * (-0.56947974731198146905 + y);
      }
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
      if (y < 0.21526012634400926547470)
        k = -50.49126109160526797810 * (-0.64935947085522228899 + y) * (0.46974783056122619477 - 1.26085718948377406347*y + y2) * (0.48792231357963262276 - 0.59385649192870404089*y + y2) * (0.25611084691786996550 + 0.88802200926364339091*y + y2) * (0.68907513661844702964 + 1.36014627339129782358*y + y2);
      else if (y < 1.5694797473119814690506)
        k = -0.034414565306740439421181 * (-1.84111179474874436 + y) * (-1.4377963552304881514 + y) * (-0.75994894568687385283 + y) * (0.4788491554147262815930 + y) * (0.949424897648337898003 + y) * (3.45344786332473319 - 3.70808010848103044 * y + y2) * (12.18722190215235760461 + 6.31866315108407263441*y + y2);
      else if (y < 1.7847398736559907345253)
        k = 26.14316088796202325336 * (-1.43779712926478619915 +  y) * (3.692596063915834557912 - 3.83497472010219620819*y +  y2) * (3.64175067444915645841 - 3.816151698746538734088*y +  y2) * (2.21858428493804009 - 2.6397857343145376826*y +  y2) * (0.9394117034107813899 - 1.902369830015286824*y +  y2);
      else {
        double ym2 = -2 + y;
        double ym25 = ym2*ym2*ym2*ym2*ym2;
        k = -24.3136856383365042853 * ym25 * (1.9854151897194796069 - 2.8134555537564699353*y + y2)*(0.6360558906432458046 - 1.311862172051363286*y + y2);
      }
    }
    ax[i] = k;
  }
  return ax;
}


// Gaussian kernel unscaled (with roughness 1)
// NOTE: since pnorm(8.29) = 2^-53 but pnorm(8.3) = 0,
// returning zero outside the reasonable range is a good idea.
// 1. The ratio phi(x) / phi(0) < 2^-53 when abs(x) > 2*sqrt(106*log(2)) ~ 8.572
// Solve[PDF[NormalDistribution[], x]/PDF[NormalDistribution[], 0] == 2^{-53}, x, Reals]
// 2. Tail area: 2*(1 - \int_{-x}^x Phi(x)) < 2^-53 when abs(x) > ~ 8.2924
// Solve[CDF[NormalDistribution[], x] == 2^{-54}, x, Reals]
arma::vec kgaussian2(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    val = (val > 8.2924) ? 0 : 0.39894228040143268 * std::exp(-0.5 * val * val); });
  return ax;
}

// Gaussian kernel convolution
arma::vec kgaussian2conv(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    val = (val > 11.7272) ? 0 : 0.28209479177387814 * std::exp(-0.25 * val * val); });
  return ax;
}

// The cut-off level is such that the area of two tails is < MachEps:
// NIntegrate[PDF[NormalDistribution[], t]*(3/2 - t^2/2), {t, -Infinity, -8.713}] ~ -5.54079*10^-17
arma::vec kgaussian4(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    val = (val > 8.713) ? 0 : 0.39894228040143268 * std::exp(-0.5 * val*val) * (1.5 - 0.5 * val*val); });
  return ax;
}

// NIntegrate[PDF[NormalDistribution[0, Sqrt[2]],  t]*(1.6875 - 0.4375*t^2 + 0.015625*t^4), {t, -Infinity, -12.68}] ~ 5.43501*10^-17
// N[Factor[(108 - 28*x^2 + x^4)/64/2/Sqrt[Pi], Extension -> {Sqrt[22]}], 20]
arma::vec kgaussian4conv(const arma::vec& x) {
  arma::vec ax = arma::abs(x);
  ax.for_each([](arma::vec::elem_type& val) {
    val = (val > 12.68) ? 0 : -0.0044077311214668460 * std::exp(-0.25 * val*val) * (23.380831519646859 - val*val) * (-4.6191684803531409 + val*val);
  });
  return ax;
}


// [[Rcpp::export]]
arma::vec kernelFunCPP(arma::vec x, std::string kernel, int order, bool convolution = false) {
  if (!convolution) {
    if (kernel == "gaussian") {
      switch (order) {
      case 2: return kgaussian2(x);
      case 4: return kgaussian4(x);
      }
    } else if (kernel == "triangular") {
      switch (order) {
      case 2: return ktriangular2(x);
      case 4: return ktriangular4(x);
      }
    } else if (kernel == "epanechnikov") {
      switch (order) {
      case 2: return kepanechnikov2(x);
      case 4: return kepanechnikov4(x);
      }
    } else if (kernel == "quartic") {
      switch (order) {
      case 2: return kquartic2(x);
      case 4: return kquartic4(x);
      }
    } else if (kernel == "uniform") {
      switch (order) {
      case 2: return kuniform2(x);
      case 4: return kuniform4(x);
      }
    }
  } else { // Convolution
    if (kernel == "gaussian") {
      switch (order) {
      case 2: return kgaussian2conv(x);
      case 4: return kgaussian4conv(x);
      }
    } else if (kernel == "triangular") {
      switch (order) {
      case 2: return ktriangular2conv(x);
      case 4: return ktriangular4conv(x);
      }
    } else if (kernel == "epanechnikov") {
      switch (order) {
      case 2: return kepanechnikov2conv(x);
      case 4: return kepanechnikov4conv(x);
      }
    } else if (kernel == "quartic") {
      switch (order) {
      case 2: return kquartic2conv(x);
      case 4: return kquartic4conv(x);
      }
    } else if (kernel == "uniform") {
      switch (order) {
      case 2: return kuniform2conv(x);
      case 4: return kuniform4conv(x);
      }
    }
  }
  throw std::runtime_error("kernelFunCPP: Invalid kernel type (should be gaussian, uniform, triangular, epanechnikov, or quartic) or order (should be 2 or 4).");
  return x;
}

// [[Rcpp::export]]
arma::mat kernelWeightsOneCPP(arma::vec x, arma::vec xout, arma::vec bw, std::string kernel = "gaussian", int order = 2, bool convolution = false) {
  arma::uword ng = xout.n_elem;
  arma::uword nb = bw.n_elem;
  arma::uword nx = x.n_elem;
  // Rcout << "ngrid: " << ng << ", nbw: " << nb << "\n";
  if (nb != ng)
    Rcpp::stop("kernelWeightsOneCPP: bw and xout must have the same length.");
  // arma::vec xs = x / bw; // Scaling by the bandwidth
  // arma::vec gs = xout / bw;
  const arma::vec inv_bw = 1.0 / bw;
  arma::mat kw(ng, nx);

  for (arma::uword i = 0; i < nx; i++)
    kw.col(i) = kernelFunCPP((xout - x[i]) % inv_bw, kernel, order, convolution);

  return kw;
}

// [[Rcpp::export]]
arma::sp_mat sparseKernelWeightsOneCPP(arma::vec x, arma::vec xout, arma::vec bw, std::string kernel = "gaussian", int order = 2, bool convolution = false) {
  arma::uword ng = xout.n_elem;
  arma::uword nb = bw.n_elem;
  arma::uword nx = x.n_elem;
  if (nb != ng)
    Rcpp::stop("sparseKernelWeightsOneCPP: bw and xout must have the same length.");

  std::vector<arma::uword> rows;
  std::vector<arma::uword> cols;
  std::vector<double>      vals;
  rows.reserve(ng); cols.reserve(ng); vals.reserve(ng);   // coarse guess

  const arma::vec inv_bw = 1.0 / bw;
  double maxX = 1.0;
  if (kernel == "gaussian") {
    if ((order == 2) & (!convolution)) maxX =  8.2924;
    if ((order == 2) & ( convolution)) maxX =  11.7272;
    if ((order == 4) & (!convolution)) maxX =  8.713;
    if ((order == 4) & ( convolution)) maxX =  12.68;
  }

  for (arma::uword i = 0; i < nx; ++i) {
    arma::vec  u   = (xout - x[i]) % inv_bw;
    arma::uvec nz  = arma::find(arma::abs(u) < maxX); // cheap range test
    if (nz.is_empty()) continue;

    arma::vec k = kernelFunCPP(u.elem(nz), kernel, order, convolution);
    rows.insert(rows.end(), nz.begin(), nz.end());
    cols.insert(cols.end(), nz.n_elem, i);                  // repeat i
    vals.insert(vals.end(), k.begin(), k.end());
  }

  arma::umat loc(2, rows.size());
  for (size_t j = 0; j < rows.size(); ++j) {
    loc(0, j) = rows[j];
    loc(1, j) = cols[j];
  }

  return arma::sp_mat(loc, arma::vec(vals), ng, nx);
}

// [[Rcpp::export]]
arma::mat kernelWeightsCPP(arma::mat x, arma::mat xout, arma::mat bw, std::string kernel = "gaussian", int order = 2, bool convolution = false) {
  arma::uword d = x.n_cols;
  // The product kernel matrix starts with the first dimension (there is at least one column or row)
  // We need to compute the product kernel starting from the 2nd till the last dimension
  // The loop from k=1, k<d ensures that
  arma::mat pk = kernelWeightsOneCPP(x.col(0), xout.col(0), bw.col(0), kernel, order, convolution);
  for (arma::uword k = 1; k < d; ++k)
    pk %= kernelWeightsOneCPP(x.col(k), xout.col(k), bw.col(k), kernel, order, convolution);
  return pk;
}

// [[Rcpp::export]]
arma::sp_mat sparseKernelWeightsCPP(arma::mat x, arma::mat xout, arma::mat bw, std::string kernel = "gaussian", int order = 2, bool convolution = false) {
  arma::uword d = x.n_cols;
  // The product kernel matrix starts with the first dimension (there is at least one column or row)
  // We need to compute the product kernel starting from the 2nd till the last dimension
  // The loop from k=1, k<d ensures that
  arma::sp_mat pk = sparseKernelWeightsOneCPP(x.col(0), xout.col(0), bw.col(0), kernel, order, convolution);
  for (arma::uword k = 1; k < d; ++k)
    pk %= sparseKernelWeightsOneCPP(x.col(k), xout.col(k), bw.col(k), kernel, order, convolution);
  return pk;
}

// [[Rcpp::export]]
NumericVector kernelDensityCPP(arma::mat x, arma::mat xout, arma::vec weights, arma::mat bw, std::string kernel = "gaussian", int order = 2, bool convolution = false, int chunks = 0) {
  arma::uword ng = xout.n_rows;
  arma::uword nbw = bw.n_rows;
  if (nbw != ng)
    Rcpp::stop("kernelDensityCPP: bw and xout must have the same number of rows.");
  arma::uword nx = x.n_rows;
  // arma::uword d = x.n_cols;
  if (chunks == 0) { // Auto-selection of matrix slices
    long long memuse = static_cast<long long>(ng) * static_cast<long long>(nx) * 8;
    int neededblks = std::max(1, static_cast<int>(memuse / 536870912)); // 2^29 = 512 MB RAM target use; not allowing more that double that amount
    chunks = neededblks;
  }

  const double n = arma::sum(weights);
  arma::vec prod_bw = arma::prod(bw, 1);
  arma::vec nb   = n * prod_bw; // n*prod(b) in the denominator

  arma::vec out(ng);
  if (chunks == 1) {
    arma::mat kw = kernelWeightsCPP(x, xout, bw, kernel, order, convolution);
    kw.each_row() %= weights.t();
    out = arma::sum(kw, 1) / nb;
  } else { // Splitting the job into chunks
    arma::uvec inds = round(arma::linspace<arma::uvec>(1, ng+1, chunks+1));
    arma::uvec starts = inds.subvec(0, inds.n_elem-2) - 1;
    arma::uvec ends = inds.subvec(1, inds.n_elem-1) - 2;

    KernelDensityWorker worker(x, xout, weights, bw, kernel, order, convolution, starts, ends, out, nb);
    parallelFor(0, chunks, worker);
  }

  NumericVector rout = Rcpp::NumericVector(out.begin(), out.end());
  return rout;
}

// [[Rcpp::export]]
NumericVector kernelSmoothCPP(arma::mat x, arma::vec y, arma::mat xout, arma::vec weights, arma::mat bw, std::string kernel = "gaussian", int order = 2, bool LOO = false, bool convolution = false, int chunks = 0) {
  arma::uword ng = xout.n_rows;
  arma::uword nb = bw.n_rows;
  if (nb != ng)
    Rcpp::stop("kernelSmoothCPP: bw and xout must have the same number of rows.");
  arma::uword nx = x.n_rows;
  if (chunks == 0) { // Auto-selection of matrix slices
    long long memuse = static_cast<long long>(ng) * static_cast<long long>(nx) * 8;
    int neededblks = std::max(1, static_cast<int>(memuse / 536870912)); // 2^29 = 512 MB RAM target use; not allowing more that double that amount
    chunks = neededblks;
  }

  arma::vec ysum(ng); // Nadaraya--Watson numerator
  arma::vec ksum(ng); // Nadaraya--Watson denominator
  arma::vec out(ng);
  if (chunks == 1) {
    arma::mat kw = kernelWeightsCPP(x, xout, bw, kernel, order, convolution);
    kw.each_row() %= weights.t();

    // LOO: setting diagonal elements to zero, assuming x = xout (R makes sure it happens, though.)
    if (LOO) kw.diag().zeros();
    ksum = sum(kw, 1);
    kw.each_row() %= y.t();      // Nadaraya--Watson numerator: y_i * w_ij (in place to save memory)
    ysum = sum(kw, 1);
  } else {
    arma::uvec inds = round(arma::linspace<arma::uvec>(1, ng+1, chunks+1));
    arma::uvec starts = inds.subvec(0, inds.n_elem-2) - 1;
    arma::uvec ends = inds.subvec(1, inds.n_elem-1) - 2;
    KernelSmoothWorker worker(x, y, xout, weights, bw, kernel, order, convolution, LOO, starts, ends, ksum, ysum);
    parallelFor(0, chunks, worker);
  }

  out = ysum / ksum;
  NumericVector rout = Rcpp::NumericVector(out.begin(), out.end());
  return rout;
}

