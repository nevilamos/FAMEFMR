#' Growth Stage Optimisation (GSO)
#'
#' @details Internal GSO function.
#' @details Calcuates the GSO using a non-linear optimisation funciton. The input is the
#' species lists from DataGen() function
#'
#' @param spp the species lists from DataGen() function
#'
#' @return returns a data.frame with the result of the optimisation with a 'geom' column
#' @export
#'

gso <- function(spp) {
  x0 <-
    const_vector(ncol(spp), 1 / ncol(spp)) ## *** Corresponds to the number of columns in the input ***
  run <-
    nloptr::nloptr(
      x0 = x0,
      eval_f = geomean.obj,
      eval_grad_f = geomean.grad,
      lb = const_vector(length(x0), 0),
      ub = const_vector(length(x0), 1),
      eval_g_ineq = eval_cs,
      eval_jac_g_ineq = eval_cs_jac,
      opts = list(
        ## nloptr Algorithm to use.
        "algorithm" = "NLOPT_LD_MMA",
        ## Termination conditions
        "maxeval" = 1000,
        "ftol_rel" = 1.0e-15,
        "xtol_rel" = 1.0e-8,
        ## Suppress output during optimization
        "print_level" = 0
      ),
      spp = spp
    )
  ## Find the objective value corresponding to the optimal solution
  res <- rbind(c(run$solution, geomean.fun(run$solution, spp)))
  colnames(res) <- c(names(spp), "geom")
  res <- data.frame(res)
  return (res)
}
