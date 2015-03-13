binom_diff_conf <- function (n1, x1, n2, x2, step=0.1, level=3.841459) {
#   3.841459 == qchisq(.95, df=1)

  if ((x1 >= n1) | (x2 >=n2)) {
    stop("n (number of successes) should not be larger than x (number of tries).")
  }
  
  mle <- (x1 / n1) - (x2 / n2)   
  eps <- .Machine$double.eps
  
  llr_test <- function(theta, n1=n1, x1=x1, n2=n2, x2=x2)
  {
    llr <- function(p, n1, x1, n2, x2, theta) {
      npllik <-  function(x,n,p) -2*(x*log((p*n)/x)+(n-x)*log(((1-p)*n)/(n-x)))
      return (npllik(x1,n1,p) + npllik(x2,n2,p-theta)) 
    }
    
#     we must ensure that 0 < const < 1 and 0 < const-theta < 1

    u_bound <- min(1-eps, 1+theta-eps)
    l_bound <- max(eps, theta+eps)
    
    optimized_llr <- optimize(f = llr, 
                              lower = l_bound, 
                              upper = u_bound, 
                              n1 = n1, 
                              x1 = x1, 
                              n2 = n2, 
                              x2 = x2,
                              theta = theta)
    
    c_star <- optimized_llr$minimum
    optimized_p <- optimized_llr$objective
    p_value <- 1 - pchisq(optimized_p, df=1)
    list(neg_2_llr = optimized_p,
         c_star = c_star,
         p_value = p_value)
  }
  
  l_bound <- u_bound <- mle
  l_step <- u_step <- step
  l_llr_test_val <- u_llr_test_val <- 0
  
  for (i in 1:8) {
    l_bound <- l_bound - l_step
    while ((l_llr_test_val < level) & (l_bound > -1)) {
      l_llr_test_val <- llr_test(theta = l_bound, 
                                 n1 = n1, 
                                 x1 = x1, 
                                 n2 = n2, 
                                 x2 = x2)$neg_2_llr
      l_bound <- l_bound - l_step
    }
    
    l_bound <- l_bound + l_step
    l_step <- l_step / 10
    
    if (l_bound != mle) {
      l_llr_test_val <- llr_test(theta = l_bound, 
                                 n1 = n1, 
                                 x1 = x1, 
                                 n2 = n2, 
                                 x2 = x2)$neg_2_llr
    }
  }
  
  for (i in 1:8) {
    u_bound = u_bound + u_step
    while ((u_llr_test_val < level) & (u_bound < 1)) {
      u_llr_test_val <- llr_test(theta = u_bound, 
                                 n1 = n1, 
                                 x1 = x1, 
                                 n2 = n2, 
                                 x2 = x2)$neg_2_llr
      u_bound <- u_bound + u_step
    }
    u_bound <- u_bound - u_step
    u_step <- u_step / 10
    if (l_bound != mle) {
      u_llr_test_val <- llr_test(theta = u_bound, 
                                n1 = n1, 
                                x1 = x1, 
                                n2 = n2, 
                                x2 = x2)$neg_2_llr
    }
  }
  
  return(list(l_bound=l_bound, 
              u_bound=u_bound, 
              l_last_step_size=l_step, 
              u_last_step_size=u_step,
              u_llr_test_val=u_llr_test_val, 
              l_llr_test_val=l_llr_test_val) )
}

