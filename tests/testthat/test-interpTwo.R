test_that("interpolation works and is monotone", {
  xx <- -4:4  # Global data for EL evaluation
  xt <- c(-3, 2.5)  # Transition at these points
  w <- c(9, rep(1, 8))
  w <- w / sum(w)
  f <- Vectorize(function(m) -2*EL0(xx, mu = m, ct = w, chull.fail = "none", SEL = TRUE)$logelr)
  xl <- seq(xt[1]-1, xt[1]+1, 0.1)
  xr <- seq(xt[2]-1, xt[2]+1, 0.1)
  fl <- f(xl)
  fr <- f(xr)
  wm <- weighted.mean(xx, w)
  wv <- weighted.mean((xx-wm)^2, w) / sum(w)
  # xseq <- seq(min(xx), max(xx), 0.1)
  # plot(xseq, f(xseq))
  # lines(xseq, (xseq - wm)^2 / wv, col = 2, lty = 2)
  fil <- interpTwo(xl, f, mean = wm, var = wv, at = xt[1], gap = 0.5)
  fir <- interpTwo(xr, f, mean = wm, var = wv, at = xt[2], gap = 0.5)
  # plot(xl, fl); lines(xl, fil)
  # plot(xr, fr); lines(xr, fir)
  expect_equal(tail(fl, 6), tail(fil, 6))  # 6 because gap = 0.5 and step = 0.1
  expect_equal(head(fr, 6), head(fr, 6))
  expect_true(all(sign(diff(fil)) == -1))
  expect_true(all(sign(diff(fir)) == 1))
  expect_true(all((fil <= fl)[is.finite(fl)]))  # Extrapolation to a lower function
  expect_true(all((fir >= fr)[1:16]))  # Extrapolation to a higher function; ignore the end

  # Miror case with branches in the opposite direction
  xt <- c(-2.5, 3)  # Transition at these points
  xl <- seq(xt[1]-1, xt[1]+1, 0.1)
  xr <- seq(xt[2]-1, xt[2]+1, 0.1)
  w <- c(rep(1, 8), 9)
  w <- w / sum(w)
  fl <- f(xl)
  fr <- f(xr)
  wm <- weighted.mean(xx, w)
  wv <- weighted.mean((xx-wm)^2, w) / sum(w)
  # plot(xseq, f(xseq))
  # lines(xseq, (xseq - wm)^2 / wv, col = 2, lty = 2)
  fil <- interpTwo(xl, f, mean = wm, var = wv, at = xt[1], gap = 0.5)
  fir <- interpTwo(xr, f, mean = wm, var = wv, at = xt[2], gap = 0.5)
  # plot(xl, fl); lines(xl, fil)
  # plot(xr, fr); lines(xr, fir)
  expect_equal(tail(fl, 6), tail(fil, 6))  # 6 because gap = 0.5 and step = 0.1
  expect_equal(head(fr, 6), head(fr, 6))
  expect_true(all(sign(diff(fil)) == -1))
  expect_true(all(sign(diff(fir)) == 1))
  expect_true(all(tail(fil >= fl, 16)))  # Extrapolation to a higher function
  expect_true(all((fir <= fr)[is.finite(fr)]))  # Extrapolation to a lower function
})



