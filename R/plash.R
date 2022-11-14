# likelihood of full model
plash_lik <- function(Y, LL, FF, cc) {

  n <- nrow(Y)
  p <- ncol(Y)
  Lambda <- exp(crossprod(LL, FF)) - matrix(data = rep(cc, p), ncol = p)
  lik <- dpois(x = drop(Y), lambda = drop(Lambda), log = TRUE)
  return(mean(lik))

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
  Y, K, cluster, tol = 1e-8, update_L = T, update_F = F, update_c = F,
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
