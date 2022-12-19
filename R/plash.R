# likelihood of full model
plash_lik <- function(Y, LL, FF, cc, size) {

  Lambda <- exp(crossprod(LL, FF)) - outer(cc, size)
  lik <- dpois(drop(Y), drop(Lambda), log = TRUE)
  return(sum(lik))

}

#' Title
#'
#' @param Y
#' @param K
#' @param tol
#' @param update_L
#' @param update_F
#' @param update_c
#'
#' @return
#' @export
#'
plash <- function(
  Y, K, tol = 1e-8, update_L = T, update_F = F, update_c = F,
  true_LL = NULL, true_FF = NULL, true_cc = NULL
  ) {

  n <- nrow(Y)
  p <- ncol(Y)

  if (update_c) {

    cc <- rep(0, n)

  } else {

    cc <- true_cc

  }

  if (update_L) {

    LL <- matrix(
      data = runif(K * n), nrow = K, ncol = n
    )

  } else {

    LL <- true_LL

  }

  if (update_F) {

    FF <- matrix(
      data = runif(K * n), nrow = K, ncol = p
    )

  } else {

    FF <- true_FF

  }

  current_lik <- plash_lik(Y, LL, FF, cc)
  converged <- FALSE

  t <- 1

  while (!converged) {

    print(t)
    print(current_lik)

    FF_T <- t(FF)

    print("Updating L...")
    if (update_L) {

      # Update L
      for (i in 1:n) {

        LL[, i] <- solve_pois_reg_offset_b(
          X = FF_T, y = Y[i, ], c = cc[i], b_init = LL[, i]
        )

      }

    }

    LL_T <- t(LL)

    print("Updating F...")
    if (update_F) {

      # Update F
      for (j in 1:p) {

        FF[, j] <- solve_pois_reg_offset_b(
          X = LL_T, y = Y[, j], c = cc, b_init = FF[, j]
        )

      }

    }

    # get updated transpose of F
    FF_T <- t(FF)

    print("updating c...")
    if (update_c) {

      # Update c
      for (i in 1:n) {

        cc[i] <- solve_pois_reg_offset_c(
          X = FF_T, y = Y[i, ], b = LL[, i], cc[i]
        )

      }

    }

    new_lik <- plash_lik(Y, LL, FF, cc)
    if (new_lik < current_lik) {

      converged <- TRUE

    } else {

      rel_improvement <- abs((new_lik - current_lik) / current_lik)
      if (rel_improvement < tol) {

        converged <- TRUE

      } else {

        current_lik <- new_lik

      }

    }

    t <- t + 1

  }

  return(
    list(
      LL = LL, FF = FF, cc = cc
    )
  )

}

#' Title
#'
#' @param Y
#' @param K
#' @param tol
#' @param update_L
#' @param update_F
#' @param update_c
#'
#' @return
#' @export
#'
plash_parallel <- function(
  Y, K, cluster, tol = 1e-8, update_L = T, update_F = T, update_c = T,
  true_LL = NULL, true_FF = NULL, true_cc = NULL
) {

  n <- nrow(Y)
  p <- ncol(Y)

  library(foreach)

  doParallel::registerDoParallel(cl = cluster)

  if (update_c) {

    cc <- rep(0, n)

  } else {

    cc <- true_cc

  }

  if (update_L) {

    LL <- matrix(
      data = runif(K * n), nrow = K, ncol = n
    )

  } else {

    LL <- true_LL

  }

  if (update_F) {

    FF <- matrix(
      data = runif(K * p), nrow = K, ncol = p
    )

  } else {

    FF <- true_FF

  }

  current_lik <- plash_lik(Y, LL, FF, cc)
  converged <- FALSE

  t <- 1

  while (!converged) {

    print(t)
    print(current_lik)

    FF_T <- t(FF)

    print("Updating L...")
    if (update_L) {

      LL <- foreach::foreach(
        i = 1:n,
        .combine = 'cbind'
      ) %dopar% {

        solve_pois_reg_offset_b(
          X_T = FF, X = FF_T, y = Y[i, ], c = cc[i], b_init = LL[, i]
        )

      }

    }

    LL_T <- t(LL)

    print("Updating F...")
    if (update_F) {

      FF <- foreach::foreach(
        j = 1:p,
        .combine = 'cbind'
      ) %dopar% {

        solve_pois_reg_offset_b(
          X_T = LL, X = LL_T, y = Y[, j], c = cc, b_init = FF[, j]
        )

      }

    }

    print("updating c...")
    if (update_c) {

      cc <- foreach::foreach(
        i = 1:n,
        .combine = 'c'
      ) %dopar% {

        solve_pois_reg_offset_c(
          X_T = FF, y = Y[i, ], b = LL[, i], cc[i]
        )

      }

    }
    
    print(max(cc))

    new_lik <- plash_lik(Y, LL, FF, cc)
    if (new_lik < current_lik) {

      converged <- TRUE

    } else {

      rel_improvement <- abs((new_lik - current_lik) / current_lik)
      if (rel_improvement < tol) {

        converged <- TRUE

      } else {

        current_lik <- new_lik

      }

    }

    t <- t + 1

  }

  return(
    list(
      LL = LL, FF = FF, cc = cc, lik = new_lik
    )
  )

}

#' Title
#'
#' @param Y
#' @param K
#' @param tol
#' @param update_L
#' @param update_F
#' @param update_c
#'
#' @return
#' @export
#'
plash_parallel_optim_v2 <- function(
  Y, K, cluster, tol = 1e-8, update_L = T, update_F = T, update_c = T,
  true_LL = NULL, true_FF = NULL, true_cc = NULL
) {
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  library(foreach)
  
  doParallel::registerDoParallel(cl = cluster)
  
  if (update_L) {
    
    LL <- matrix(
      data = runif(K * n), nrow = K, ncol = n
    )
    
    cc <- rep(0, n)
    
  } else {
    
    LL <- true_LL
    cc <- true_cc
    
  }
  
  if (update_F) {
    
    FF <- matrix(
      data = runif(K * p), nrow = K, ncol = p
    )
    
  } else {
    
    FF <- true_FF
    
  }
  
  current_lik <- plash_lik(Y, LL, FF, cc)
  converged <- FALSE
  
  t <- 1
  
  while (!converged) {
    
    print(t)
    print(current_lik)
    
    FF_T <- t(FF)
    
    print("Updating L and c...")
    if (update_L) {
      
      LL_cc <- foreach::foreach(
        i = 1:n,
        .combine = 'cbind'
      ) %dopar% {
        
        solve_pois_reg_unknown_offset(X = FF_T, y = Y[i, ], init = c(cc[i], LL[, i]))
        
      }
      
    }
    
    LL <- LL_cc[-1, ]
    cc <- pmax(LL_cc[1, ], 0)
    LL_T <- t(LL)
    print(glue::glue("max c = {max(cc)}"))
    
    print("Updating F...")
    if (update_F) {
      
      FF <- foreach::foreach(
        j = 1:p,
        .combine = 'cbind'
      ) %dopar% {
        
        solve_pois_reg_offset_b(
          X_T = LL, X = LL_T, y = Y[, j], c = cc, b_init = FF[, j]
        )
        
      }
      
    }
    
    new_lik <- plash_lik(Y, LL, FF, cc)
    if (new_lik < current_lik) {
      
      converged <- TRUE
      
    } else {
      
      rel_improvement <- abs((new_lik - current_lik) / current_lik)
      if (rel_improvement < tol) {
        
        converged <- TRUE
        
      } else {
        
        current_lik <- new_lik
        
      }
      
    }
    
    t <- t + 1
    
  }
  
  return(
    list(
      LL = LL, FF = FF, cc = cc, lik = new_lik
    )
  )
  
}

# Write function to fit L and F with c 0 and a potential offset and intercept
# Then, compare this to the offset of glmpca to ensure that things are consistent. 
#' Title
#'
#' @param Y
#' @param K
#' @param tol
#' @param update_L
#' @param update_F
#' @param update_c
#'
#' @return
#' @export
#'
plash_parallel_glmpca <- function(
  Y, K, cluster, offset = FALSE, offset_vec = NULL, intercept = FALSE,
  tol = 1e-8, update_L = T, update_F = T, update_c = F,
  true_LL = NULL, true_FF = NULL, true_cc = NULL
) {
  
  K_total <- K + offset + intercept
  K_fixed_LL <- as.numeric(offset)
  K_fixed_FF <- K_total - K
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  library(foreach)
  
  doParallel::registerDoParallel(cl = cluster)
  
  if (update_L) {
    
    LL <- matrix(
      data = runif(K_total * n), nrow = K_total, ncol = n
    )
    
    cc <- rep(0, n)
    
  } else {
    
    LL <- true_LL
    cc <- true_cc
    
  }
  
  if (update_F) {
    
    FF <- matrix(
      data = runif(K_total * p), nrow = K_total, ncol = p
    )
    
  } else {
    
    FF <- true_FF
    
  }
  
  if (offset && intercept) {
    
    # Offset
    LL[1, ] <- 1
    FF[1, ] <- offset_vec
    
    # Intercept
    FF[2, ] <- 1
    
  } else if (offset) {
    
    # Offset
    LL[1, ] <- 1
    FF[1, ] <- offset_vec
    
  } else if (intercept) {
    
    # Intercept
    FF[1, ] <- 1
    
  }
  
  current_lik <- plash_lik(Y, LL, FF, cc)
  converged <- FALSE
  
  t <- 1
  
  while (!converged) {
    
    print(t)
    print(current_lik)
    
    FF_T <- t(FF)
    
    print("Updating L...")
    if (update_L) {
      
      LL <- foreach::foreach(
        i = 1:n,
        .combine = 'cbind'
      ) %do% {
        
        if (K_fixed_LL > 0) {
          
          LL_fixed_b <- LL[1:K_fixed_LL, i]
          
        } else {
          
          LL_fixed_b <- NULL
          
        }
        
        solve_pois_reg_offset_fixed_b(
          X_T = FF, X = FF_T, y = Y[i, ], c = cc[i],
          fixed_b = LL_fixed_b, b_init = LL[(K_fixed_LL + 1):K_total, i]
        )
        
      }
      
    }
    
    LL_T <- t(LL)
    
    print("Updating F...")
    if (update_F) {
      
      FF <- foreach::foreach(
        j = 1:p,
        .combine = 'cbind'
      ) %do% {
        
        if (K_fixed_FF > 0) {
          
          FF_fixed_b <- FF[1:K_fixed_FF, j]
          
        } else {
          
          FF_fixed_b <- NULL
          
        }
        
        solve_pois_reg_offset_fixed_b(
          X_T = LL, X = LL_T, y = Y[, j], c = cc,
          fixed_b = FF_fixed_b, b_init = FF[(K_fixed_FF + 1):K_total, j]
        )
        
      }
      
    }
    
    new_lik <- plash_lik(Y, LL, FF, cc)
    if (new_lik < current_lik) {
      
      converged <- TRUE
      
    } else {
      
      rel_improvement <- abs((new_lik - current_lik) / current_lik)
      if (rel_improvement < tol) {
        
        converged <- TRUE
        
      } else {
        
        current_lik <- new_lik
        
      }
      
    }
    
    t <- t + 1
    
  }
  
  return(
    list(
      LL = LL, FF = FF, cc = cc, lik = new_lik
    )
  )
  
}

# Write function to fit L and F with c 0 and a potential offset and intercept
# Then, compare this to the offset of glmpca to ensure that things are consistent. 
#' Title
#'
#' @param Y
#' @param K
#' @param tol
#' @param update_L
#' @param update_F
#' @param update_c
#'
#' @return
#' @export
#'
plash_omni <- function(
  Y, K, offset = FALSE, intercept = FALSE,
  constant_offset = TRUE, parallel = TRUE,
  tol = 1e-8, update_L = T, update_F = T, update_c = T, update_L_c_simul = FALSE,
  init_LL = NULL, init_FF = NULL, init_cc = NULL, update_c_first = TRUE,
  min_iter = 3, max_iter = 30
) {
  
  if (parallel) {
    
    `%loopdo%` <- foreach::`%dopar%`
    
  } else {
    
    `%loopdo%` <- foreach::`%do%`
    
  }
  
  K_total <- K + offset + intercept
  K_fixed_LL <- as.numeric(offset)
  K_fixed_FF <- K_total - K
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  cm <- Matrix::colMeans(Y)
  
  if (offset) {
    
    offset_vec <- log(cm)
    
  }
  
  if (constant_offset && offset) {
    
    constant_offset_vec <- cm
    
  } else {
    
    constant_offset_vec <- rep(1, p)
    
  }
  
  if (is.null(init_LL)) {
    
    LL <- matrix(
      data = runif(K_total * n), nrow = K_total, ncol = n
    )
    
  } else {
    
    LL <- init_LL
    
  }
  
  if (is.null(init_FF)) {
    
    FF <- matrix(
      data = runif(K_total * p), nrow = K_total, ncol = p
    )
    
  } else {
    
    FF <- init_FF
    
  }
  
  if (is.null(init_cc)) {
    
    cc <- rep(0, n)
    
  } else {
    
    cc <- init_cc
    
  }
  
  if (offset && intercept) {
    
    # Offset
    LL[1, ] <- 1
    FF[1, ] <- offset_vec
    
    # Intercept
    FF[2, ] <- 1
    
  } else if (offset) {
    
    # Offset
    LL[1, ] <- 1
    FF[1, ] <- offset_vec
    
  } else if (intercept) {
    
    # Intercept
    FF[1, ] <- 1
    
  }
  
  current_lik <- plash_lik(Y, LL, FF, cc, constant_offset_vec)
  converged <- FALSE
  
  t <- 1
  
  while ((!converged && t <= max_iter) || t <= min_iter) {
    
    print(t)
    print(current_lik)
    
    FF_T <- t(FF)
    
    if (!update_L_c_simul) {
      
      if (update_c && update_c_first) {
        
        print("Updating c...")
        cc <- foreach::foreach(
          i = 1:n,
          .combine = 'c'
        ) %loopdo% {
          
          solve_pois_reg_offset_c(
            X_T = FF, y = Y[i, ], b = LL[, i], c_init = cc[i], size = constant_offset_vec
          )
          
        }
        
      }
      
      print("Updating L...")
      if (update_L) {
        
        LL <- foreach::foreach(
          i = 1:n,
          .combine = 'cbind'
        ) %loopdo% {
          
          if (K_fixed_LL > 0) {
            
            LL_fixed_b <- LL[1:K_fixed_LL, i]
            
          } else {
            
            LL_fixed_b <- NULL
            
          }
          
          solve_pois_reg_offset_fixed_b(
            X_T = FF, X = FF_T, y = Y[i, ], c = cc[i] * constant_offset_vec,
            fixed_b = LL_fixed_b, b_init = LL[(K_fixed_LL + 1):K_total, i]
          )
          
        }
        
      }
      
    } else{
      
      if (update_L && update_c) {
        
        print("Updating L and c")
        
        LL_cc <- foreach::foreach(
          i = 1:n,
          .combine = 'cbind'
        ) %loopdo% {
          
          if (K_fixed_LL > 0) {
            
            LL_fixed_b <- LL[1:K_fixed_LL, i]
            
          } else {
            
            LL_fixed_b <- NULL
            
          }
          
          solve_pois_reg_unknown_offset_fixed_b(
            X = FF_T, y = Y[i, ], init = c(cc[i], LL[(K_fixed_LL + 1):K_total, i]), 
            fixed_b = LL_fixed_b, size = constant_offset_vec
          )
          
        }
        
      }
      
      LL <- LL_cc[-1, ]
      cc <- pmax(LL_cc[1, ], 0)
        
    }
    
    LL_T <- t(LL)
    
    print("Updating F...")
    if (update_F) {
      
      FF <- foreach::foreach(
        j = 1:p,
        .combine = 'cbind'
      ) %loopdo% {
        
        if (K_fixed_FF > 0) {
          
          FF_fixed_b <- FF[1:K_fixed_FF, j]
          
        } else {
          
          FF_fixed_b <- NULL
          
        }
        
        solve_pois_reg_offset_fixed_b(
          X_T = LL, X = LL_T, y = Y[, j], c = cc * constant_offset_vec[j],
          fixed_b = FF_fixed_b, b_init = FF[(K_fixed_FF + 1):K_total, j]
        )
        
      }
      
    }
    
    if (!update_L_c_simul && !update_c_first && update_c) {
      
      print("updating c...")
        
      cc <- foreach::foreach(
        i = 1:n,
        .combine = 'c'
      ) %loopdo% {
          
        solve_pois_reg_offset_c(
          X_T = FF, y = Y[i, ], b = LL[, i], c_init = cc[i], size = constant_offset_vec
        )
          
        }
      
    }
    
    new_lik <- plash_lik(Y, LL, FF, cc, constant_offset_vec)
    if (new_lik == -Inf) {
      
      return(
        list(
          LL = LL, FF = FF, cc = cc, lik = new_lik
        )
      )
      
    }
    if (new_lik < current_lik && t >= min_iter) {
      
      converged <- TRUE
      
    } else {
      
      rel_improvement <- abs((new_lik - current_lik) / current_lik)
      if (rel_improvement < tol && t >= min_iter) {
        
        converged <- TRUE
        
      } else {
        
        current_lik <- new_lik
        
      }
      
    }
    
    t <- t + 1
    
  }
  
  return(
    list(
      LL = LL, FF = FF, cc = cc, lik = new_lik, size = constant_offset_vec
    )
  )
  
}

get_feasible_init <- function(
    Y, K, cc, offset = FALSE, intercept = FALSE, constant_offset = TRUE  
  ) {
  
  K_total <- K + offset + intercept
  K_fixed_LL <- as.numeric(offset)
  K_fixed_FF <- K_total - K
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  cm <- Matrix::colMeans(Y)
  
  if (offset) {
    
    offset_vec <- log(cm)
    
  } else {
    
    offset_vec <- 0
    
  }
  
  if (constant_offset && offset) {
    
    constant_offset_vec <- cm
    
  } else {
    
    constant_offset_vec <- rep(1, p)
    
  }
  
  LL <- matrix(
    data = 1, nrow = K_total, ncol = n
  )
  
  FF <- matrix(
    data = 1, nrow = K_total, ncol = p
  )
  
  if (offset && intercept) {
    
    # Offset
    LL[1, ] <- 1
    FF[1, ] <- offset_vec
    
    min_offset <- min(offset_vec)
    
    # Intercept + rest of rows set
    FF[2:K_total, ] <- 1
    
    for (j in 1:n) {
      
      LL[2:K_total, j] <- ((log(cc[j]) - min_offset) / K) + 1e-10
      
    }
    
  } else if (offset) {
    
    min_offset <- min(offset_vec)
    # Offset
    LL[1, ] <- 1
    FF[1, ] <- offset_vec
    
    for (j in 1:n) {
      
      LL[2:K_total, j] <- ((log(cc[j]) - min_offset) / K) + 1e-10
      
    }
    
  } else if (intercept) {
    
    # Intercept
    FF[1, ] <- 1
    
    for (j in 1:n) {
      
      LL[, j] <- (log(cc[j]) / K) + 1e-10
      
    }
    
  }
  
  # Now, check for feasibility
  Lambda <- exp(crossprod(LL, FF)) - outer(cc, constant_offset_vec)
  if (all(drop(Lambda) >= 0)) {
    
    print("SUCCESS!!!!!!!!!!!!")
    
  } else {
    
    print("You messed up :(")
    
  }
  
  return(list(
    LL = LL, FF = FF
  ))
  
}
