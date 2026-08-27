## Description: basic function for the paper Spherically Embedded Time Series with Unknown Trend and Periodic Components

#############################################################################################
#############################################################################################
## Basic functions for optimal transport
#############################################################################################
#############################################################################################

# - usage: SpheGeoDist(a,b)
# - Description: geodesic distance between two vectors after projecting to the sphere
# - Inputs:
#    a, b: vectors
# - Output: numeric scalar.
# - Assumptions: 'a' and 'b' have the same length.

SpheGeoDist <- function(y1,y2) {
  if (abs(length(y1) - length(y2)) > 0) {
    stop("y1 and y2 should be of the same length.")
  }
  if ( !isTRUE( all.equal(l2norm(y1),1) ) ) {
    stop("y1 is not a unit vector.")
  }
  if ( !isTRUE( all.equal(l2norm(y2),1) ) ) {
    stop("y2 is not a unit vector.")
  }
  y1 = y1 / l2norm(y1)
  y2 = y2 / l2norm(y2)
  if (sum(y1 * y2) > 1){
    return(0)
  } else if (sum(y1*y2) < -1){
    return(pi)
  } else return(acos(sum(y1 * y2)))
}

# - usage: unit_norm(a)
# - Description: project the vector to the sphere
# - Inputs:
#    a: numeric vector.
# - Output: numeric spherical vector.
# - Assumptions: NULL.

unit_norm <- function(x) x / sqrt(sum(x^2))

# - usage: geo_dist(a,b)
# - Description: geodesic distance between two spherical vectors
# - Inputs:
#    a, b: spherical vectors
# - Output: numeric scalar.
# - Assumptions: 'a' and 'b' have the same length.

geo_dist <- function(u, v) {
  val <- sum(u * v)
  val <- min(1, max(-1, val))
  acos(val)
}

# - usage: log_rotation(a,b,c)
# - Description: rotation operator given in (2.1) of the paper
# - Inputs:
#    a, b: spherical vectors
#    c: tolerance
# - Output: list of rotation operator, angle for the rotation, and basis of span. 
# - Assumptions: 'a' and 'b' have the same length.


log_rotation <- function(x, y, tol = 1e-12) { ## move x to y, mimic y-x
  x <- as.numeric(x); y <- as.numeric(y)
  stopifnot(length(x) == length(y))
  
  # ensure unit length
  x <- x / sqrt(sum(x^2))
  y <- y / sqrt(sum(y^2))
  
  cxy <- sum(x * y)
  cxy <- max(min(cxy, 1), -1)
  
  # identical or antipodal cases
  if (1 - abs(cxy) < tol) {
    if (cxy > 0) {
      return(list(L = matrix(0, length(x), length(x)), theta = 0, u1 = x, u2 = NA))
    } else {
      # antipodal: choose any unit u2 perpendicular to x
      j <- which.min(abs(x))
      e <- rep(0, length(x)); e[j] <- 1
      v <- e - sum(e * x) * x
      vn <- sqrt(sum(v^2))
      if (vn < tol) {
        if (length(x) < 2) stop("Need dimension >= 2 for antipodal case")
        v <- c(-x[2], x[1], if (length(x) > 2) rep(0, length(x) - 2) else NULL)
        vn <- sqrt(sum(v^2))
      }
      u2 <- v / vn
      theta <- pi
      Q <- tcrossprod(u2, x) - tcrossprod(x, u2)
      return(list(L = theta * Q, theta = theta, u1 = x, u2 = u2))
    }
  }
  
  # general case
  u1 <- x
  v <- y - cxy * u1
  u2 <- v / sqrt(sum(v^2))
  theta <- acos(cxy)
  Q <- tcrossprod(u2, u1) - tcrossprod(u1, u2)
  
  list(L = theta * Q, theta = theta, u1 = u1, u2 = u2)
}

# - usage: exp_from_log_rotation_predict(a,b)
# - Description: exponential operator given in (2.2) of the paper
# - Inputs:
#    a: rotation operator
#    b: angle for the rotation
# - Output: exponential map 
# - Assumptions: 'a' and 'b' have the same length.

exp_from_log_rotation_predict <- function(q,theta) {
  th <- theta
  L <- q
  if (is.null(th)) {
    eig <- eigen(L, symmetric = FALSE)
    V <- eig$vectors; D <- diag(exp(eig$values))
    Re(V %*% D %*% solve(V))
  } else {
    I <- diag(nrow(L))
    K <- L / max(th, 1e-16)
    I + sin(th) * K + (1 - cos(th)) * (K %*% K)
  }
}


#############################################################################################
#############################################################################################
## Basic functions for simulated data generation
#############################################################################################
#############################################################################################

# - usage: rSkewSym(a,b)
# - Description: random rotation operator generation
# - Inputs:
#    a: dimension of sphere
#    b: concertration parameter for von-Mises Fisher distribution
# - Output: rotation operator
# - Assumptions: NULL.

rSkewSym <- function(d,kappa) {
  m2 <- c(1:d) / sqrt(sum((c(1:d))^2))
  m1 <- rmovMF(1, kappa*m2)
  opt_trans <- log_rotation(m2,m1) 
  M <- opt_trans$L
  return(M)
}


# - usage: AR_generate_function_plot(a,b,c,d,e)
# - Description: spherical AR generation for the data generation plot based on (3.8) in the paper 
# - Inputs:
#    a: seed
#    b: sample size
#    c: dimension of sphere
#    d: AR coefficient
#    e: burnin number
# - Output: list of skew-symmetric operators 
# - Assumptions: NULL.

AR_generate_function_plot <- function(seed_set, sample_size, d, coef,burnin=500) {
  set.seed(seed_set)
  matrix_mean <- 0 
  matrix_list <- lapply(1:3, function(i) rSkewSym(d, d^3))
  for (t in 4:(sample_size + burnin + 3)) {
    matrix_list[[t]] <- matrix_mean + coef[1] * (matrix_list[[t-1]] - matrix_mean) +
      coef[2] * (matrix_list[[t-2]] - matrix_mean) + coef[3] * (matrix_list[[t - 3]] - matrix_mean) + rSkewSym(d, d^4)
  }
  return(matrix_list[(4 + burnin):(sample_size + burnin + 3)])
}

# - usage: AR_generate_function(a,b,c,d,e)
# - Description: spherical AR generation for the simulation based on (3.8) in the paper 
# - Inputs:
#    a: seed
#    b: sample size
#    c: dimension of sphere
#    d: AR coefficient
#    e: burnin number
# - Output: list of skew-symmetric operators 
# - Assumptions: NULL.

AR_generate_function <- function(seed_set, sample_size, d, coef,burnin=500) {
  set.seed(seed_set)
  matrix_mean <- 0 
  matrix_list <- lapply(1:length(coef), function(i) rSkewSym(d, 5*d))
  for (t in (length(coef)+1):(sample_size + burnin + length(coef))) {
    matrix_save <- matrix_mean
    for (coef_num in 1:length(coef)) {
      matrix_save <- matrix_save + coef[coef_num] * (matrix_list[[t-coef_num]] - matrix_mean)
    }
    matrix_list[[t]] <- matrix_save + rSkewSym(d, 20*d)
  }
  return(matrix_list[(length(coef) + 1 + burnin):(sample_size + burnin + length(coef))])
}


# - usage: intrinsic_mean_sphere(a,b,c)
# - Description: Frechet mean calculation
# - Inputs:
#    a: matrix of spherical data
#    b: tolerence
#    c: maximum number of iterations
# - Output: numeric spherical vector
# - Assumptions: NULL. 

intrinsic_mean_sphere <- function(y, tol = 1e-8, max_iter = 1000) {
  n <- nrow(y)
  d <- ncol(y)
  
  # Initialize with extrinsic mean
  mu <- colSums(y)
  mu <- mu / sqrt(sum(mu^2))
  
  for (iter in 1:max_iter) {
    tangent_sum <- rep(0, d)
    for (i in 1:n) {
      theta <- acos(sum(mu * y[i,]))
      if (theta > 1e-12) {
        tangent <- (theta / sin(theta)) * (y[i,] - cos(theta) * mu)
        tangent_sum <- tangent_sum + tangent
      }
    }
    tangent_mean <- tangent_sum / n
    norm_t <- sqrt(sum(tangent_mean^2))
    
    # Check convergence
    if (norm_t < tol) break
    
    # Update mu using exponential map
    mu <- cos(norm_t) * mu + sin(norm_t) * (tangent_mean / norm_t)
  }
  
  return(mu)
}

# - usage: norm_vec(a)
# - Description: spherical norm calculation
# - Inputs:
#    a: numeric spherical vector.
# - Output: numeric scalar.
# - Assumptions: NULL.

norm_vec <- function(x) {
  sqrt(sum(x^2))
}


# - usage: get_rotation_matrix(a,b)
# - Description: rotation operator given in (2.1) of the paper
# - Inputs:
#    a, b: spherical vectors
# - Output: rotation operator
# - Assumptions: 'a' and 'b' have the same length.

get_rotation_matrix <- function(g1, g3) { ## g3-g1
  p <- length(g1)
  theta <- geo_dist(g1, g3)
  
  # If g1 and g3 are effectively the same, return Identity
  if (theta < 1e-10) {
    return(diag(p))
  }
  
  u1 <- g1
  proj <- g3 - sum(g3 * u1) * u1
  u2 <- proj / norm_vec(proj)
  Q <- (u2 %*% t(u1)) - (u1 %*% t(u2))
  Q2 <- Q %*% Q
  I <- diag(p)
  R <- I + sin(theta) * Q + (1 - cos(theta)) * Q2
  
  return(R)
}

# - usage: periodic_component_function(a,b,c)
# - Description: periodic component generation for simulation
# - Inputs:
#    a: sample size
#    b: dimension of sphere
#    c: period
# - Output: numeric matrix 
# - Assumptions: NULL.

periodic_component_function <- function(sample_size,d,period){
  ## ensure the intrinsice mean of g(t) over a period is exactly mu
  mu <- c(1:d) / sqrt(sum((c(1:d))^2))
  
  raw_basis <- pracma::nullspace(t(mu))
  u <- raw_basis[, 1]
  v <- raw_basis[, 2]
  theta <- 0.3 
  
  # 4. Construct the points g(1)...g(12)
  g_period <- matrix(0, nrow = period, ncol = d)
  for (t in 1:period) {
    phi <- 2 * pi * (t / period)
    g_period[t, ] <- cos(theta) * mu + sin(theta) * (cos(phi) * u + sin(phi) * v)
  }
  return(g_period)
}

# - usage: period_generate_funtion(a)
# - Description: periodic spherical time series generation
# - Inputs:
#    a: list of skew-symmetric operators 
# - Output: numeric matrix 
# - Assumptions: NULL.

period_generate_funtion <- function(AR_list){
  sample_size <- length(AR_list)
  d <- dim(AR_list[[1]])[1]
  sphere_matrix <- matrix(NA,nrow = sample_size,ncol = d)
  mu <- c(1:d) / sqrt(sum((c(1:d))^2))
  
  ## get R_t^{(2)}
  for (t in 1:sample_size) {
    angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
    sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]] ,angle_n) %*% mu
  }
  
  ## get R_t^{(1)}
  periodic_function_matrix <- periodic_component_function(sample_size,d,period=12)
  
  for (t in 1:sample_size) {
    current_ind <- t + 12 - floor((t + 12 - 1)/12) * 12
    sphere_matrix[t,] <- get_rotation_matrix(mu,periodic_function_matrix[current_ind,]) %*% sphere_matrix[t,]
  }
  
  return(sphere_matrix)
}

# - usage: trend_component_function(a,b,c,d)
# - Description: trend component generation
# - Inputs:
#    a: scaled time point
#    b: sphere dimension
#    c: spherical vector indices for trend generation
#    d: scale of trend
# - Output: spherical vector
# - Assumptions: NULL.

trend_component_function <- function(u, d, inc_idx = c(1,2,3), A = 0.6) {
  # intrinsic mean
  mu <- c(1:d) / sqrt(sum((c(1:d))^2))
  
  # project onto tangent space at mu
  proj_tangent <- function(v) {
    v - sum(v * mu) * mu
  }
  
  # trend direction: increasing in inc_idx
  trend_raw <- rep(0, d)
  trend_raw[inc_idx] <- 1
  v_dir <- proj_tangent(trend_raw)
  v_dir <- v_dir / sqrt(sum(v_dir^2))
  
  v <- A * (u - 0.5) * v_dir
  
  # exponential map on sphere
  nv <- sqrt(sum(v^2))
  if (nv < 1e-10) return(mu)
  
  cos(nv) * mu + sin(nv) * v / nv
}

# - usage: trend_generate_funtion(a)
# - Description: spherical time series with trend generation
# - Inputs:
#    a: numeric matrix of spherical data
# - Output: numeric matrix 
# - Assumptions: NULL.

trend_generate_funtion <- function(period_matrix){
  sample_size <- nrow(period_matrix)
  d <- ncol(period_matrix)
  
  mu <- c(1:d) / sqrt(sum((c(1:d))^2))
  
  u_grid <- (1:sample_size) / sample_size
  
  trend_function_matrix <- t(sapply(u_grid, trend_component_function,
                                    d,
                                    inc_idx = c(1,2,3),
                                    A = 0.6))
  
  
  y_matrix <- matrix(NA,nrow = sample_size,ncol = d)
  for (t in 1:sample_size) {
    y_matrix[t,] <- get_rotation_matrix(mu,trend_function_matrix[t,]) %*% period_matrix[t,]
  }
  
  return(y_matrix)
}



#############################################################################################
#############################################################################################
## Basic functions for global and local Frechet regression, 
## mostly downloaded from <https://github.com/functionaldata/tFrechet>
#############################################################################################
#############################################################################################


# - Description: same with the function in <https://github.com/functionaldata/tFrechet>

l2norm <- function(x){
  as.numeric(sqrt(crossprod(x)))
}

# - Description: same with the function in <https://github.com/functionaldata/tFrechet>

SpheGeoGrad <- function(x, y, tol=1e-8) {
  dot <- sum(x * y)
  dot <- max(min(dot, 1), -1)  # clamp
  tmp <- 1 - dot^2
  
  # if nearly coincident or antipodal, return 0 gradient
  if (tmp < tol) {
    return(rep(0, length(x)))
  }
  
  return(-(1/sqrt(tmp)) * x)
}

# - Description: same with the function in <https://github.com/functionaldata/tFrechet>

SpheGeoHess <- function(x, y, tol=1e-8) {
  dot <- sum(x * y)
  dot <- max(min(dot, 1), -1)
  tmp <- 1 - dot^2
  p <- length(x)
  
  # near coincident points: limit exists
  if (tmp < tol) {
    return(2 * (diag(p) - outer(y, y)))
  }
  
  return(- dot * (tmp)^(-1.5) * (x %*% t(x)))
}


## local Frechet regression using moving window cross validation

# - Description: same with the function in <https://github.com/functionaldata/tFrechet>

SetBwRange <- function(xin, xout, kernel_type) {
  xinSt <- unique(sort(xin))
  bw.min <- max(diff(xinSt), xinSt[2] - min(xout), max(xout) -
                  xinSt[length(xin)-1])*1.1 / (ifelse(kernel_type == "gauss", 3, 1) *
                                                 ifelse(kernel_type == "gausvar", 2.5, 1))
  bw.max <- diff(range(xin))/3
  if (bw.max < bw.min) {
    if (bw.min > bw.max*3/2) {
      bw.max <- bw.min*1.01
    } else bw.max <- bw.max*3/2
  }
  return(list(min=bw.min, max = bw.max))
}

# - Description: same with the function in <https://github.com/functionaldata/tFrechet>

kerFctn <- function(kernel_type){
  if (kernel_type=='gauss'){
    ker <- function(x){
      dnorm(x)
    }
  } else if(kernel_type=='rect'){
    ker <- function(x){
      as.numeric((x<=1) & (x>=-1))
    }
  } else if(kernel_type=='epan'){
    ker <- function(x){
      n <- 1
      (2*n+1) / (4*n) * (1-x^(2*n)) * (abs(x)<=1)
    }
  } else if(kernel_type=='gausvar'){
    ker <- function(x) {
      dnorm(x)*(1.25-0.25*x^2)
    }
  } else if(kernel_type=='quar'){
    ker <- function(x) {
      (15/16)*(1-x^2)^2 * (abs(x)<=1)
    }
  } else {
    stop('Unavailable kernel')
  }
  return(ker)
}

# - Description: same with the function in <https://github.com/functionaldata/tFrechet>

LocSpheReg <- function(xin=NULL, yin=NULL, xout=NULL, optns=list()){
  
  if (is.null(xin))
    stop ("xin has no default and must be input by users.")
  if (is.null(yin))
    stop ("yin has no default and must be input by users.")
  if (is.null(xout))
    xout <- xin
  if (!is.vector(xin) | !is.numeric(xin))
    stop("xin should be a numerical vector.")
  if (!is.matrix(yin) | !is.numeric(yin))
    stop("yin should be a numerical matrix.")
  if (!is.vector(xout) | !is.numeric(xout))
    stop("xout should be a numerical vector.")
  if (length(xin)!=nrow(yin))
    stop("The length of xin should be the same as the number of rows in yin.")
  if (sum(abs(rowSums(yin^2) - rep(1,nrow(yin))) > 1e-6)){
    yin = yin / sqrt(rowSums(yin^2))
    warning("Each row of yin has been standardized to enforce sum of squares equal to 1.")
  }
  
  if (is.null(optns$bw)){
    optns$bw <- "CV" #max(sort(xin)[-1] - sort(xin)[-length(xin)]) * 1.2
  }
  if (is.character(optns$bw)) {
    if (optns$bw != "CV") {
      warning("Incorrect input for optns$bw.")
    }
  } else if (!is.numeric(optns$bw)) {
    stop("Mis-specified optns$bw.")
  }
  if (length(optns$bw) > 1)
    stop("bw should be of length 1.")
  
  if (is.null(optns$kernel))
    optns$kernel <- "gauss"
  
  if (is.numeric(optns$bw)) {
    bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$kernel)
    if (optns$bw < bwRange$min | optns$bw > bwRange$max) {
      optns$bw <- "CV"
      warning("optns$bw is too small or too large; reset to be chosen by CV.")
    }
  } 
  if (optns$bw == "CV") {
    optns$bw <- bwCV_sphe_moving_window(xin = xin, yin = yin, xout = xout, optns = optns)
  }
  yout <- LocSpheGeoReg(xin = xin, yin = yin, xout = xout, optns = optns)
  res <- list(xout = xout, yout = yout, xin = xin, yin = yin, optns = optns)
  class(res) <- "spheReg"
  return(res)
}

# - Description: same with the function in <https://github.com/functionaldata/tFrechet>

LocSpheGeoReg <- function(xin, yin, xout, optns = list()) {
  k = length(xout)
  n = length(xin)
  m = ncol(yin)
  
  bw <- optns$bw
  ker <- kerFctn(optns$kernel)
  
  yout = sapply(1:k, function(j){
    mu0 = mean(ker((xout[j] - xin) / bw))
    mu1 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j]))
    mu2 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j])^2)
    s = ker((xout[j] - xin) / bw) * (mu2 - mu1 * (xin - xout[j])) /
      (mu0 * mu2 - mu1^2)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if (sum(sapply(1:n, function(i) sum(yin[i,]*y0))[ker((xout[j] - xin) / bw)>0] > 1-1e-8)){
      y0 = y0 + rnorm(3) * 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
        return(list(value = Inf))
      }
      f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
      res = sapply(1:n, function(i){
        grad_i = SpheGeoGrad(yin[i,], y)
        return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
      }, simplify = "array")
      h = 2 * apply(res, 1:2, mean)
      return(list(value=f, gradient=g, hessian=h))
    }
    res = trust::trust(objFctn, y0, 0.1, 1e5)
    return(res$argument / l2norm(res$argument))
  })
  return(t(yout))
}

# - Description: same with the function in <https://github.com/functionaldata/tFrechet>

bwCV_sphe <- function(xin, yin, xout, optns) {
  yin <- yin[order(xin),]
  xin <- sort(xin)
  compareRange <- (xin > min(xin) + diff(range(xin))/5) & (xin < max(xin) - diff(range(xin))/5)
  
  objFctn <- function(bw) {
    optns1 <- optns
    optns1$bw <- bw
    folds <- numeric(length(xin))
    n <- sum(compareRange)
    numFolds <- ifelse(n > 30, 10, sum(compareRange))
    
    tmp <- c(sapply(1:ceiling(n/numFolds), function(i)
      sample(x = seq_len(numFolds), size = numFolds, replace = FALSE)))
    tmp <- tmp[1:n]
    repIdx <- which(diff(tmp) == 0)
    for (i in which(diff(tmp) == 0)) {
      s <- tmp[i]
      tmp[i] <- tmp[i-1]
      tmp[i-1] <- s
    }

    
    folds[compareRange] <- tmp
    
    yout <- lapply(seq_len(numFolds), function(foldidx) {
      testidx <- which(folds == foldidx)
      res <- LocSpheGeoReg(xin = xin[-testidx], yin = yin[-testidx,], xout = xin[testidx], optns = optns1)
      res # each row is a spherical vector
    })
    yout <- do.call(rbind, yout)
    yinMatch <- yin[which(compareRange)[order(tmp)],]
    mean(sapply(1:nrow(yout), function(i) SpheGeoDist(yout[i,], yinMatch[i,])^2))
  }
  bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$ker)

  res <- optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
  res$minimum
}

# - usage: bwCV_sphe_moving_window(a,b,c,d)
# - Description: rolling-window cross validation for local Frechet regression
# - Inputs:
#    a: numeric matrix for scaled time indices
#    b: numeric matrix for spherical data
#    c: numeric matrix for scaled time indices for fitted model
#    d: option for kernel choice
# - Output: numeric scalar for the chosen bandwidth 
# - Assumptions: length of a should be equal to the row number of b

bwCV_sphe_moving_window <- function(xin, yin, xout, optns) {
  ord <- order(xin)
  yin <- yin[ord, ]
  xin <- xin[ord]
  
  n <- length(xin)

  window_size <- ceiling(n * 0.7) 
  
  objFctn <- function(bw) {
    optns1 <- optns
    optns1$bw <- bw
    
    test_indices <- (window_size + 1):n
    
    errors <- sapply(test_indices, function(i) {
      train_idx <- (i - window_size):(i - 1)
      test_idx <- i
      res <- LocSpheGeoReg(
        xin = xin[train_idx], 
        yin = yin[train_idx, , drop = FALSE], 
        xout = xin[test_idx], 
        optns = optns1
      )
      geo_dist(res, yin[test_idx, ])^2
    })

    mean(errors)
  }
  
  # Set bandwidth search range
  bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$ker)
  
  # Optimize bandwidth
  res <- optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
  
  return(res$minimum)
}



# - Description: same with the function in <https://github.com/functionaldata/tFrechet>

GloSpheGeoReg <- function(xin, yin, xout = NULL) {
  if (is.null(xout))
    xout <- xin
  
  k = nrow(xout)
  n = nrow(xin)
  m = ncol(yin)
  
  xbar <- colMeans(xin)
  invSigma <- diag(1/xbar) + rep(1,length(xbar)) %*% t(rep(1,length(xbar))) / length(xbar)
  
  
  yout = sapply(1:k, function(j){
    s <- 1 + t(t(xin) - xbar) %*% invSigma %*% (xout[j,] - xbar)
    s <- as.vector(s)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if ( any( sapply( 1:n, function(i) isTRUE( all.equal( sum(yin[i,]*y0), 1 ) ) ) ) ){
      y0[1] = y0[1] + 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
        return(list(value = Inf))
      }
      f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
      res = sapply(1:n, function(i){
        grad_i = SpheGeoGrad(yin[i,], y)
        return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
      }, simplify = "array")
      h = 2 * apply(res, 1:2, mean)
      return(list(value=f, gradient=g, hessian=h))
    }
    res = trust::trust(objFctn, y0, 0.1, 1)
    return(res$argument)
  })
  return(t(yout))
}

#############################################################################################
#############################################################################################
## period estimation
#############################################################################################
#############################################################################################

# - usage: covariate_function(a,b)
# - Description: covariate matrix for period estimation, given by (3.4) in the paper
# - Inputs:
#    a: period
#    b: sample size
# - Output: numeric matrix 
# - Assumptions: NULL

covariate_function <- function(period,sample_size){
  x_matrix <- rep(0,period)
  x_matrix[1] <- 1
  for (i in 2:sample_size) {
    x_vec <- rep(0,period)
    current_ind <- i + period - floor((i + period -1)/period) * period
    x_vec[current_ind] <- 1
    x_matrix <- rbind(x_matrix,x_vec)
  }
  return(x_matrix)
}

# - usage: period_est_sphere_function(a,b,c,d)
# - Description: period estimation for spherical time series
# - Inputs:
#    a: numeric matrix of spherical data
#    b: upper bound of period candidates
#    c: sample size
#    d: tuning parameter for the penalty term
# - Output: Vector of penalized RSS in (3.2) of the paper
# - Assumptions: NULL

period_est_sphere_function <- function(data,period_max,sample_size,lambda){
  y_matrix <- data
  
  loss_vec <- rep(0,period_max)
  for (p in 1:period_max) {
    x <- covariate_function(p,sample_size)
    
    res <- GloSpheGeoReg(xin = x, yin = y_matrix)
    
    loss_value <- 0
    for (i in 1:sample_size) {
      single_diff <- acos(t(y_matrix[i,]) %*% res[i,])
      
      loss_value <- loss_value + single_diff
    }
    loss_vec[p] <- loss_value + lambda * p
    print(paste('Period',p,': loss',loss_vec[p]))
  }
  return(loss_vec)
}

# - usage: BIC_function_new(a,b,c,d,e)
# - Description: BIC function for tuning parameter selection
# - Inputs:
#    a: RSS vector
#    b: sequence of tuning parameter
#    c: sample size
#    d: selected bandwidth
#    e: constant in penalty 
# - Output: Vector of BIC values
# - Assumptions: NULL

BIC_function_new <- function(loss_seq,lambda_seq,sample_size,bw,cst_set){
  
  BIC_vec <- NULL
  for (lambda in lambda_seq) {
    loss_vec_now <- loss_seq
    loss_vec_now <- loss_vec_now + seq(1,length(loss_seq)) * lambda
    p <- which(loss_vec_now==min(loss_vec_now))
    
    loss_mean <- loss_seq[p] / sample_size
    
    select_rate <- max(bw^2*log(1/bw),sqrt((log(1/bw))^2/(sample_size*bw)))
    
    current_BIC <- log(loss_mean) + p * select_rate / cst_set 
    BIC_vec <- c(BIC_vec,current_BIC)
    
  }
  
  return(BIC_vec)
}

#############################################################################################
#############################################################################################
## AR model fitting
#############################################################################################
#############################################################################################

# - usage: get_spherical_residuals(a,b,c)
# - Description: spherical residual based on spherical removal operation
# - Inputs:
#    a,b,c: numric spherical vector
# - Output: numric spherical vector
# - Assumptions: lengths of a and b are the same

get_spherical_residuals <- function(Y, M1, mu = NULL) {
  
  if (is.null(mu)) {
    mu <- intrinsic_mean_sphere(Y)
  }
  
  n <- nrow(Y)
  p <- ncol(Y)
  residuals <- matrix(NA, nrow = n, ncol = p)
  
  for (t in 1:n) {
    y_t <- Y[t, ]
    m_t <- M1[t, ]
    
    R <- get_rotation_matrix(g1 = m_t, g3 = mu)

    res_t <- R %*% y_t
    
    residuals[t, ] <- as.vector(res_t)
  }
  
  return(list(residuals = residuals, mu = mu))
}




# - usage: AR_fit_function(a,b)
# - Description: spherical AR coefficient estimation
# - Inputs:
#    a: numeric matric of spherical data
#    b: spherical AR order
# - Output: list of AR coefficient estimation and the intercept estimation
# - Assumptions: NULL


AR_fit_function <- function(residual_sphere_matrix,AR_order){
  sample_size <- nrow(residual_sphere_matrix)
  d <- ncol(residual_sphere_matrix)
  
  mu <- intrinsic_mean_sphere(residual_sphere_matrix)
  T_list <- list()
  for (t in 1:sample_size) {
    opt_trans <- log_rotation(mu,residual_sphere_matrix[t,]) ## m1-m2
    T_list[[t]] <- opt_trans$L
  }
  # 1. Estimate the sample mean operator
  mu_T_hat <- Reduce("+", T_list) / sample_size
  T_centered <- lapply(T_list, function(x) x - mu_T_hat)
  # 2. Compute sample autocovariances lambda_k
  get_lambda_k <- function(k) {
    if (k == 0) {
      sum_val <- sum(sapply(T_centered, function(x) sum(x * x)))
      return(sum_val / sample_size)
    }
    
    sum_val <- 0
    for (t in 1:(sample_size - k)) {
      sum_val <- sum_val + sum(T_centered[[t]] * T_centered[[t + k]])
    }
    return(sum_val / (sample_size - k))
  }
  
  # Compute lambdas from 0 to p
  lambdas <- sapply(0:AR_order, get_lambda_k)
  
  # 3. Construct the Yule-Walker system
  Gamma_matrix <- matrix(0, nrow = AR_order, ncol = AR_order)
  for (i in 1:AR_order) {
    for (j in 1:AR_order) {
      Gamma_matrix[i, j] <- lambdas[abs(i - j) + 1]
    }
  }
  
  rhs_vector <- lambdas[2:(AR_order + 1)]
  
  # 4. Solve for alpha coefficients
  alpha_hat <- solve(Gamma_matrix) %*% rhs_vector
  return(list(alpha_hat,mu_T_hat))
}

# - usage: difference_AR_fit_function(a,b)
# - Description: differencing-based spherical AR coefficient estimation
# - Inputs:
#    a: numeric matric of spherical data
#    b: spherical AR order
# - Output: list of AR coefficient estimation and the intercept estimation
# - Assumptions: NULL

difference_AR_fit_function <- function(residual_sphere_matrix,AR_order){
  sample_size <- nrow(residual_sphere_matrix)-1
  d <- ncol(residual_sphere_matrix)
  
  T_list <- list()
  for (t in 1:sample_size) {
    opt_trans <- log_rotation(residual_sphere_matrix[t,],residual_sphere_matrix[(t+1),]) ## m1-m2
    T_list[[t]] <- opt_trans$L
  }
  # 1. Estimate the sample mean operator
  mu_T_hat <- Reduce("+", T_list) / sample_size
  T_centered <- lapply(T_list, function(x) x - mu_T_hat)
  # 2. Compute sample autocovariances lambda_k
  # Using the Frobenius inner product for matrices (corresponds to H \otimes H)
  # cite: 137, 209
  get_lambda_k <- function(k) {
    if (k == 0) {
      sum_val <- sum(sapply(T_centered, function(x) sum(x * x)))
      return(sum_val / sample_size)
    }
    
    sum_val <- 0
    for (t in 1:(sample_size - k)) {
      sum_val <- sum_val + sum(T_centered[[t]] * T_centered[[t + k]])
    }
    return(sum_val / (sample_size - k))
  }
  
  # Compute lambdas from 0 to p
  lambdas <- sapply(0:AR_order, get_lambda_k)
  
  # 3. Construct the Yule-Walker system
  Gamma_matrix <- matrix(0, nrow = AR_order, ncol = AR_order)
  for (i in 1:AR_order) {
    for (j in 1:AR_order) {
      Gamma_matrix[i, j] <- lambdas[abs(i - j) + 1]
    }
  }
  
  rhs_vector <- lambdas[2:(AR_order + 1)]
  
  # 4. Solve for alpha coefficients
  alpha_hat <- solve(Gamma_matrix) %*% rhs_vector
  return(list(alpha_hat,mu_T_hat))
}


# - usage: order_select_function(a,b,c)
# - Description: spherical AR order selection
# - Inputs:
#    a: numeric matrix of spherical data
#    b: spherical AR order candidates
#    c: number of steps ahead
# - Output: optimal order
# - Assumptions: NULL

order_select_function <- function(residual_sphere_matrix, candidate.order=c(1:10), predict.step.num){
  sample_size <- nrow(residual_sphere_matrix)
  train_num <- round(sample_size*0.7) ## 0.6,0.7,0.8
  test_num <- sample_size - train_num
  test_round_num <- floor(test_num/predict.step.num)
  test_num <- test_round_num*predict.step.num
  train_num <- sample_size - test_num
  
  test_results <- NULL
  for (p_candidate in candidate.order) {
    rolling_predict_angle_save <- NULL
    
    for (i in 1:test_round_num) {
      
      train_select <- ((i-1)*predict.step.num + 1):((i-1)*predict.step.num + train_num)
      train_sphere <- residual_sphere_matrix[train_select,]
      mean_est <- intrinsic_mean_sphere(train_sphere)
      
      AR_fit_results <- AR_fit_function(train_sphere, p_candidate)
      AR_coef_est <- as.vector(AR_fit_results[[1]])
      intercept <- AR_fit_results[[2]] 

      prediction_buffer <- train_sphere[(nrow(train_sphere) - length(AR_coef_est) + 1):nrow(train_sphere), , drop = FALSE]
      
      for (l in 1:predict.step.num) {
        predict_proj <- intercept
        for (k in 1:length(AR_coef_est)) {
          prev_point <- prediction_buffer[(nrow(prediction_buffer) - k + 1), ]
          opt_trans <- log_rotation(mean_est, prev_point)
          
          predict_proj <- predict_proj + AR_coef_est[k] * (opt_trans$L - intercept)
        }
        angle_n <- sqrt(sum(predict_proj^2))
        sphere_predict <- exp_from_log_rotation_predict(predict_proj, angle_n) %*% mean_est
        
        current_test_idx <- train_num + (i-1)*predict.step.num + l
        
        actual_val <- residual_sphere_matrix[current_test_idx, ]
        rolling_predict_angle_save <- c(rolling_predict_angle_save, acos(pmin(pmax(sum(actual_val * sphere_predict), -1), 1)))
        
        prediction_buffer <- rbind(prediction_buffer[-1, ], as.vector(sphere_predict))
      }
    }
    
    test_results <- c(test_results,mean(rolling_predict_angle_save))
  }
  min_ind <- which(test_results==min(test_results))
  return(candidate.order[min_ind])
}

# - usage: difference_order_select_function(a,b,c)
# - Description: differencing-based spherical AR order selection
# - Inputs:
#    a: numeric matrix of spherical data
#    b: spherical AR order candidates
#    c: number of steps ahead
# - Output: optimal order
# - Assumptions: NULL

difference_order_select_function <- function(residual_sphere_matrix, candidate.order=c(1:10), predict.step.num){
  sample_size <- nrow(residual_sphere_matrix)
  train_num <- round(sample_size*0.7) ## 0.6,0.7,0.8
  test_num <- sample_size - train_num
  test_round_num <- floor(test_num/predict.step.num)
  test_num <- test_round_num*predict.step.num
  train_num <- sample_size - test_num
  
  test_results <- NULL
  for (p_candidate in candidate.order) {
    rolling_predict_angle_save <- NULL
    
    
    for (i in 1:test_round_num) {
      train_select <- ((i-1)*predict.step.num + 1):((i-1)*predict.step.num + train_num)
      train_sphere <- residual_sphere_matrix[train_select,]
      
      AR_fit_results <- difference_AR_fit_function(train_sphere, p_candidate)
      AR_coef_est <- as.vector(AR_fit_results[[1]])
      intercept <- AR_fit_results[[2]]
      
      prediction_buffer <- train_sphere[(nrow(train_sphere) - length(AR_coef_est)):nrow(train_sphere), , drop = FALSE]
      
      for (l in 1:predict.step.num) {
        predict_proj <- intercept
        for (k in 1:length(AR_coef_est)) {
          p_prev <- prediction_buffer[nrow(prediction_buffer) - k, ]
          p_curr <- prediction_buffer[nrow(prediction_buffer) - k + 1, ]
          opt_trans <- log_rotation(p_prev, p_curr)
          predict_proj <- predict_proj + AR_coef_est[k] * (opt_trans$L - intercept)
        }
        
        base_point <- prediction_buffer[nrow(prediction_buffer), ]
        angle_n <- sqrt(sum((predict_proj %*% base_point)^2))
        sphere_predict <- exp_from_log_rotation_predict(predict_proj, angle_n) %*% base_point
        
        idx <- ((i-1)*predict.step.num + l)
        rolling_predict_angle_save <- c(rolling_predict_angle_save, acos(pmin(pmax(sum(residual_sphere_matrix[train_num + idx, ] * sphere_predict), -1), 1)))
        
        prediction_buffer <- rbind(prediction_buffer[-1, ], as.vector(sphere_predict))
      }
    }
    
    test_results <- c(test_results,mean(rolling_predict_angle_save))
  }
  min_ind <- which(test_results==min(test_results))
  return(candidate.order[min_ind])
}

#############################################################################################
#############################################################################################
## estimation simulation workflow
#############################################################################################
#############################################################################################

# - usage: full_workflow_simulation_function(a,b,c)
# - Description: estimation simulation workflow
# - Inputs:
#    a: numeric matrix of spherical data
#    b: upper bound of spherical AR order candidates
#    c: upper bound of period candidates
# - Output: list of estimated trend component, periodic component, final residuals, AR coefficients, 
#           penalized RSS for de-trended time series, penalized RSS for final residuals, selected bandwidth
# - Assumptions: NULL

full_workflow_simulation_function <- function(y_matrix,AR_order,period_upper){
  sample_size <- nrow(y_matrix)
  d <- ncol(y_matrix) 
  
  ## local fit
  x_input <- c(1:sample_size)/sample_size
  trend_component <- LocSpheReg(xin = x_input, yin = y_matrix)
  bw_current <- trend_component$optns$bw
  trend_component <- trend_component$yout
  
  mu <- intrinsic_mean_sphere(y_matrix)
  remove_trend_result <- get_spherical_residuals(y_matrix, trend_component,mu = mu)
  remove_trend_data <- remove_trend_result$residuals
  
  ## period est
  loss_vec_save_np <- period_est_sphere_function(remove_trend_data,period_upper,sample_size,0)
  
  ## remove periodic component
  x <- covariate_function(12,sample_size)
  periodic_component <- GloSpheGeoReg(xin = x, yin = remove_trend_data)
  
  final_result <- get_spherical_residuals(remove_trend_data, periodic_component)
  final_resid_data <- final_result$residuals
  
  ## check whether period=1
  loss_vec_save_np_1 <- period_est_sphere_function(final_resid_data,period_upper,sample_size,0)
  
  ## AR modeling
  AR_coef <- AR_fit_function(final_resid_data,AR_order)
  AR_coef <- as.vector(AR_coef[[1]])
  
  return(list(trend_component,periodic_component,final_resid_data,AR_coef,loss_vec_save_np,loss_vec_save_np_1,bw_current))
}

# - usage: extrinsic_full_workflow_simulation_function(a,b)
# - Description: extrinsic estimation simulation workflow
# - Inputs:
#    a: numeric matrix of spherical data
#    b: upper bound of spherical AR order candidates
# - Output: list of estimated periodic component, final residuals, AR coefficient
# - Assumptions: NULL

extrinsic_full_workflow_simulation_function <- function(y_matrix,AR_order){
  sample_size <- nrow(y_matrix)
  d <- ncol(y_matrix) 
  
  ## remove periodic component
  x <- covariate_function(12,sample_size)
  periodic_component <- GloSpheGeoReg(xin = x, yin = y_matrix)
  
  final_result <- y_matrix-periodic_component
  row_norms <- sqrt(rowSums(final_result^2))
  final_resid_data <- final_result/row_norms

  ## AR modeling
  AR_coef <- AR_fit_function(final_resid_data,AR_order)
  AR_coef <- as.vector(AR_coef[[1]])
  return(list(periodic_component,final_resid_data,AR_coef))
}

# - usage: full_workflow_MSE_function(a,b,c,d,e)
# - Description: estimation simulation workflow to get MSE
# - Inputs:
#    a: sample size
#    b: seed
#    c: true spherical AR coefficients
#    d: upper bound of period candidates
#    e: sphere dimension 
# - Output: list of MSEs of estimated trend component, periodic component, final residuals, AR coefficients, 
#           penalized RSS for de-trended time series, penalized RSS for final residuals, selected bandwidth
# - Assumptions: NULL

full_workflow_MSE_function <- function(sample.size,rep,true_AR_coef,period_upper,d_set = 7){

  AR_list <- AR_generate_function(seed_set=rep,sample_size=sample.size,d=d_set,coef=true_AR_coef)
  
  sphere_matrix <- matrix(NA,nrow = sample.size,ncol = d_set)
  for (t in 1:sample.size) {
    mu <- c(1:d_set) / sqrt(sum((c(1:d_set))^2))
    angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
    sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]] ,angle_n) %*% mu
  }
  
  period_matrix <- period_generate_funtion(AR_list)
  trend_matrix <- trend_generate_funtion(period_matrix)
  
  output_list <- full_workflow_simulation_function(trend_matrix,AR_order=length(true_AR_coef),period_upper)
  
  
  ## true trend and periodic component
  
  u_grid <- (1:sample.size) / sample.size
  
  trend_function_matrix <- t(sapply(u_grid, trend_component_function,
                                    d=d_set,
                                    inc_idx = c(1,2,3),
                                    A = 0.6))
  
  periodic_function_matrix <- periodic_component_function(sample.size,d=d_set,period=12)
  
  trend_MSE_save <- NULL
  periodic_MSE_save <- NULL
  residual_MSE_save <- NULL
  
  for (k in 1:sample.size) {
    trend_MSE_save <- c(trend_MSE_save,(acos(t(trend_function_matrix[k,]) %*% output_list[[1]][k,]))^2)
    current_ind <- k + 12 - floor((k + 12 - 1)/12) * 12
    periodic_MSE_save <- c(periodic_MSE_save,(acos(t(periodic_function_matrix[current_ind,]) %*% output_list[[2]][k,]))^2)
    residual_MSE_save <- c(residual_MSE_save,(acos(t(sphere_matrix[k,]) %*% output_list[[3]][k,]))^2)
  }
  
  results_vec <- c(mean(trend_MSE_save),mean(periodic_MSE_save),mean(residual_MSE_save))
  
  return(list(results_vec,output_list[[4]],output_list[[5]],output_list[[6]],output_list[[7]]))
}

# - usage: extrinsic_full_workflow_MSE_function(a,b,c,d)
# - Description: extrinsic estimation simulation workflow to get MSE
# - Inputs:
#    a: sample size
#    b: seed
#    c: true spherical AR coefficients
#    e: sphere dimension 
# - Output: list of MSEs of estimated periodic component, final residuals, and AR coefficients.
# - Assumptions: NULL

extrinsic_full_workflow_MSE_function <- function(sample.size,rep,true_AR_coef,d_set = 7){
  
  AR_list <- AR_generate_function(seed_set=rep,sample_size=sample.size,d=d_set,coef=true_AR_coef)
  
  sphere_matrix <- matrix(NA,nrow = sample.size,ncol = d_set)
  for (t in 1:sample.size) {
    mu <- c(1:d_set) / sqrt(sum((c(1:d_set))^2))
    angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
    sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]] ,angle_n) %*% mu
  }
  
  period_matrix <- period_generate_funtion(AR_list)
  trend_matrix <- trend_generate_funtion(period_matrix)
  
  ## true trend and periodic component
  
  u_grid <- (1:sample.size) / sample.size
  
  trend_function_matrix <- t(sapply(u_grid, trend_component_function,
                                    d=d_set,
                                    inc_idx = c(1,2,3),
                                    A = 0.6))
  
  periodic_function_matrix <- periodic_component_function(sample.size,d=d_set,period=12)
  
  final_result <- trend_matrix-trend_function_matrix
  row_norms <- sqrt(rowSums(final_result^2))
  period_matrix <- final_result/row_norms
  
  output_list <- extrinsic_full_workflow_simulation_function(period_matrix,AR_order=length(true_AR_coef))
  
  
  

  periodic_MSE_save <- NULL
  residual_MSE_save <- NULL
  
  for (k in 1:sample.size) {
    current_ind <- k + 12 - floor((k + 12 - 1)/12) * 12
    periodic_MSE_save <- c(periodic_MSE_save,(acos(t(periodic_function_matrix[current_ind,]) %*% output_list[[1]][k,]))^2)
    residual_MSE_save <- c(residual_MSE_save,(acos(t(sphere_matrix[k,]) %*% output_list[[2]][k,]))^2)
  }
  
  results_vec <- c(mean(periodic_MSE_save),mean(residual_MSE_save))
  
  return(list(results_vec,output_list[[3]]))
}




#############################################################################################
#############################################################################################
## basic functions for prediction
#############################################################################################
#############################################################################################

# - usage: KL(a,b)
# - Description: KL divergence
# - Inputs:
#    a, b: two distributional vector
# - Output: numeric scalar
# - Assumptions: lengths of a and b are the same

KL <- function(p, q) {
  # Only indices where p > 0
  idx <- which(p > 0)
  sum(p[idx] * log(p[idx] / q[idx]))
}

# - usage: JSD(a,b)
# - Description: JSD divergence
# - Inputs:
#    a, b: two distributional vector
# - Output: numeric scalar
# - Assumptions: lengths of a and b are the same

JSD <- function(p, q) {
  m <- 0.5 * (p + q)
  0.5 * KL(p, m) + 0.5 * KL(q, m)
}

# - usage: stationary_AR_predict_function(a,b,c,d,e,f)
# - Description: prediction based on spherical AR model
# - Inputs:
#    a: numeric matrix of spherical data
#    b: spherical AR order candidates
#    c: test number
#    d: number of steps ahead
#    e: choice for moving window or expanding window validation
#    f: distance measure for prediction error
# - Output: list of prediction errors, predictions, and selected AR order
# - Assumptions: NULL

stationary_AR_predict_function <- function(data,
                                           candidate.order = c(1:20),
                                           test.num,
                                           predict.step.num,
                                           eval.method,
                                           error.measure = "sphere") {
  
  sample_size <- nrow(data)
  test_round_num <- test.num/predict.step.num
  test_num <- test.num
  train_num <- sample_size - test_num
  
  test_results <- NULL
  sphere_predict_matrix <- matrix(NA, nrow = test_num, ncol = ncol(data))
  order_select_vec <- NULL
  
  for (i in 1:test_round_num) {
    print(paste("Round:", i))
    
    # Define training set for this window
    if (eval.method=='moving'){
      ## moving window
      train_select <- ((i-1)*predict.step.num + 1):((i-1)*predict.step.num + train_num)
    } else{
      ## expanding window
      train_select <- 1:((i-1)*predict.step.num + train_num)
    }
    
    train_sphere <- data[train_select,]
    mean_est <- intrinsic_mean_sphere(train_sphere)
    
    # Model Estimation
    select_order <- order_select_function(train_sphere, candidate.order, predict.step.num)
    order_select_vec <- c(order_select_vec,select_order)
    AR_fit_results <- AR_fit_function(train_sphere, select_order)
    AR_coef_est <- as.vector(AR_fit_results[[1]])
    intercept <- AR_fit_results[[2]] 
    
    # k-step recursive prediction
    prediction_buffer <- train_sphere[(nrow(train_sphere) - length(AR_coef_est) + 1):nrow(train_sphere), , drop = FALSE]
    
    for (l in 1:predict.step.num) {
      predict_proj <- intercept
      for (k in 1:length(AR_coef_est)) {
        prev_point <- prediction_buffer[(nrow(prediction_buffer) - k + 1), ]
        opt_trans <- log_rotation(mean_est, prev_point)

        predict_proj <- predict_proj + AR_coef_est[k] * (opt_trans$L - intercept)
      }
      
      angle_n <- sqrt(sum(predict_proj^2))
      sphere_predict <- exp_from_log_rotation_predict(predict_proj, angle_n) %*% mean_est
      
      current_test_idx <- train_num + (i-1)*predict.step.num + l
      sphere_predict_matrix[((i-1)*predict.step.num + l), ] <- sphere_predict
      
      actual_val <- data[current_test_idx, ]
      
      if (error.measure=='JCD'){
        predict_error_current <- JSD(sphere_predict^2,actual_val^2)
      } else{
        predict_error_current <- acos(pmin(pmax(sum(actual_val * sphere_predict), -1), 1))
      }
      
      test_results <- c(test_results, predict_error_current)
      prediction_buffer <- rbind(prediction_buffer[-1, ], as.vector(sphere_predict))
    }
  }
  
  return(list(errors = test_results, predictions = sphere_predict_matrix,order = order_select_vec))
}

# - usage: difference_AR_predict_function(a,b,c,d,e,f)
# - Description: prediction based on differenncing-based spherical AR model
# - Inputs:
#    a: numeric matrix of spherical data
#    b: spherical AR order candidates
#    c: test number
#    d: number of steps ahead
#    e: choice for moving window or expanding window validation
#    f: distance measure for prediction error
# - Output: list of prediction errors, predictions, and selected AR order
# - Assumptions: NULL

difference_AR_predict_function <- function(data,
                                           candidate.order = c(1:20),
                                           test.num,
                                           predict.step.num,
                                           eval.method,
                                           error.measure = "sphere") {
  
  sample_size <- nrow(data)
  test_round_num <- test.num/predict.step.num
  test_num <- test.num
  train_num <- sample_size - test_num
  
  test_results <- NULL
  sphere_predict_matrix <- matrix(NA, nrow = test_num, ncol = ncol(data))
  order_select_vec <- NULL
  
  for (i in 1:test_round_num) {
    print(paste("Round:", i))
    
    # Define training set for this window
    if (eval.method=='moving'){
      ## moving window
      train_select <- ((i-1)*predict.step.num + 1):((i-1)*predict.step.num + train_num)
    } else{
      ## expanding window
      train_select <- 1:((i-1)*predict.step.num + train_num)
    }
    train_sphere <- data[train_select,]
    
    select_order <- difference_order_select_function(train_sphere, candidate.order, predict.step.num)
    order_select_vec <- c(order_select_vec,select_order)
    AR_fit_results <- difference_AR_fit_function(train_sphere, select_order)
    AR_coef_est <- as.vector(AR_fit_results[[1]])
    intercept <- AR_fit_results[[2]]
    
    prediction_buffer <- train_sphere[(nrow(train_sphere) - length(AR_coef_est)):nrow(train_sphere), , drop = FALSE]
    
    for (l in 1:predict.step.num) {
      predict_proj <- intercept
      for (k in 1:length(AR_coef_est)) {
        p_prev <- prediction_buffer[nrow(prediction_buffer) - k, ]
        p_curr <- prediction_buffer[nrow(prediction_buffer) - k + 1, ]
        opt_trans <- log_rotation(p_prev, p_curr)
        predict_proj <- predict_proj + AR_coef_est[k] * (opt_trans$L - intercept)
      }
      
      base_point <- prediction_buffer[nrow(prediction_buffer), ]
      angle_n <- sqrt(sum((predict_proj %*% base_point)^2))
      sphere_predict <- exp_from_log_rotation_predict(predict_proj, angle_n) %*% base_point
      
      # Storage and recursive update
      idx <- ((i-1)*predict.step.num + l)
      sphere_predict_matrix[idx, ] <- sphere_predict
      
      if (error.measure=='JCD'){
        predict_error_current <- JSD(sphere_predict^2,(data[train_num + idx, ])^2)
      } else{
        predict_error_current <- acos(pmin(pmax(sum(data[train_num + idx, ] * sphere_predict), -1), 1))
      }
      
      test_results <- c(test_results, predict_error_current)
      ##test_results <- c(test_results, (acos(pmin(pmax(sum(data[train_num + idx, ] * sphere_predict), -1), 1)))^2)
      
      prediction_buffer <- rbind(prediction_buffer[-1, ], as.vector(sphere_predict))
    }
  }
  return(list(errors = test_results, predictions = sphere_predict_matrix,order = order_select_vec))
}

# - usage: PTAR_predict_function(a,b,c,d,e,f,g)
# - Description: prediction based on peridic-trend-spherical AR model
# - Inputs:
#    a: numeric matrix of spherical data
#    b: period
#    c: spherical AR order candidates
#    d: test number
#    e: number of steps ahead
#    f: choice for moving window or expanding window validation
#    g: distance measure for prediction error
# - Output: list of prediction errors, predictions, and selected AR order
# - Assumptions: NULL

PTAR_predict_function <- function(data,
                                  period_est,
                                  candidate.order = c(1:20),
                                  test.num,
                                  predict.step.num,
                                  eval.method,
                                  error.measure = "sphere") {
  
  sample_size <- nrow(data)
  test_round_num <- test.num/predict.step.num
  test_num <- test.num
  train_num <- sample_size - test_num
  
  test_results <- NULL
  sphere_predict_matrix <- matrix(NA, nrow = test_num, ncol = ncol(data))
  order_select_vec <- NULL
  
  for (i in 1:test_round_num) {
    print(paste("Round:", i))
    
    if (eval.method=='moving'){
      train_select <- ((i-1)*predict.step.num + 1):((i-1)*predict.step.num + train_num)
    } else{
      train_select <- 1:((i-1)*predict.step.num + train_num)
    }
    train_sphere <- data[train_select,]
    current_train_num <- length(train_select)
    
    x_input_full <- c(1:(current_train_num + predict.step.num)) / (current_train_num)
    trend_comp <- LocSpheReg(xin = x_input_full[1:current_train_num], yin = train_sphere, xout = x_input_full)$yout
    
    mu_Y <- intrinsic_mean_sphere(train_sphere)
    rem_trend <- get_spherical_residuals(train_sphere, trend_comp[1:current_train_num,], mu = mu_Y)
    mean_est1 <- intrinsic_mean_sphere(rem_trend$residuals)
    
    x_per <- covariate_function(period_est, (current_train_num + predict.step.num))
    periodic_comp <- GloSpheGeoReg(xin = x_per[1:current_train_num,], yin = rem_trend$residuals, xout = x_per)
    
    final_resid_train <- get_spherical_residuals(rem_trend$residuals, periodic_comp[1:current_train_num,])$residuals
    mean_est2 <- intrinsic_mean_sphere(final_resid_train)
    
    select_order <- order_select_function(final_resid_train, candidate.order, predict.step.num)
    order_select_vec <- c(order_select_vec,select_order)
    AR_fit <- AR_fit_function(final_resid_train, select_order)
    AR_coef <- as.vector(AR_fit[[1]])
    intercept <- AR_fit[[2]]
    
    resid_buffer <- final_resid_train[(nrow(final_resid_train) - length(AR_coef) + 1):nrow(final_resid_train), ,  drop = FALSE]
    
    
    for (l in 1:predict.step.num) {
      predict_proj <- intercept
      for (k in 1:length(AR_coef)) {
        opt_trans <- log_rotation(mean_est2, resid_buffer[(nrow(resid_buffer) - k + 1), ])
        predict_proj <- predict_proj + AR_coef[k] * (opt_trans$L - intercept)
      }
      
      angle_n <- sqrt(sum((predict_proj %*% mean_est2)^2))
      resid_pred <- exp_from_log_rotation_predict(predict_proj, angle_n) %*% mean_est2
      
      plus_per <- get_rotation_matrix(mean_est1, periodic_comp[current_train_num + l, ]) %*% resid_pred
      sphere_pred <- get_rotation_matrix(mu_Y, trend_comp[current_train_num + l, ]) %*% plus_per
      
      # Storage and Update
      idx <- ((i-1)*predict.step.num + l)
      sphere_predict_matrix[idx, ] <- sphere_pred
      
      if (error.measure=='JCD'){
        predict_error_current <- JSD(sphere_pred^2,(data[train_num + idx, ])^2)
      } else{
        predict_error_current <- acos(pmin(pmax(sum(data[train_num + idx, ] * sphere_pred), -1), 1))
      }
      
      test_results <- c(test_results, predict_error_current)
      
      resid_buffer <- rbind(resid_buffer[-1, ], as.vector(resid_pred))
    }
  }
  return(list(errors = test_results, predictions = sphere_predict_matrix,order = order_select_vec))
}



#############################################################################################
#############################################################################################
## prediction simulation workflow
#############################################################################################
#############################################################################################

# - usage: full_workflow_prediction_function(a,b,c,d,e)
# - Description: prediction simulation workflow
# - Inputs:
#    a: sample size
#    b: seed 
#    c: true spherical AR coefficients
#    d: upper bound of spherical AR order candidates
#    e: sphere dimension
# - Output: numeric matrix of prediction errors
# - Assumptions: NULL

full_workflow_prediction_function <- function(sample.size,rep,true_AR_coef,order_upper_set=c(1:20),d_set=7){
  
  ## data generation
  AR_list <- AR_generate_function(seed_set=rep,sample_size=sample.size,d=d_set,coef=true_AR_coef)
  
  sphere_matrix <- matrix(NA,nrow = sample.size,ncol = d_set)
  for (t in 1:sample.size) {
    mu <- c(1:d_set) / sqrt(sum((c(1:d_set))^2))
    angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
    sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]] ,angle_n) %*% mu
  }
  
  period_matrix <- period_generate_funtion(AR_list)
  y_matrix <- trend_generate_funtion(period_matrix)
  
  ## prediction procedure
  test_num <- round(sample.size*0.1)
  prediction_error_matrix <- matrix(NA,nrow = 3,ncol = test_num)
  
  for (i in 1:test_num) {
    y_matrix_current <- y_matrix[i:nrow(y_matrix),]
    zhu_AR_predict <- stationary_AR_predict_function(
      data = y_matrix_current,
      candidate.order = order_upper_set,
      test.num = (test_num+1-i),
      predict.step.num = (test_num+1-i),
      eval.method = 'moving'
    )
    zhu_DAR_predict <- difference_AR_predict_function(
      data = y_matrix_current,
      candidate.order = order_upper_set,
      test.num = (test_num+1-i),
      predict.step.num = (test_num+1-i),
      eval.method = 'moving'
    )
    PTAR_predict <- PTAR_predict_function(
      data = y_matrix_current,
      period_est = 12,
      candidate.order = order_upper_set,
      test.num = (test_num+1-i),
      predict.step.num = (test_num+1-i),
      eval.method = 'moving'
    )
    prediction_error_matrix[1,i] <- mean(zhu_AR_predict[[1]]) 
    prediction_error_matrix[2,i] <- mean(zhu_DAR_predict[[1]]) 
    prediction_error_matrix[3,i] <- mean(PTAR_predict[[1]]) 
  }
  
  
  return(prediction_error_matrix)
}




