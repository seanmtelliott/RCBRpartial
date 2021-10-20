#' getbound
#'
#' calls Rmosek to solve for two linear programming problem for the identified set
#'
#' @param obj defines the linear coefficients of the linear objective function
#' @param Acons defines the linear constraints in the form of Acons%*%x = bcons
#' @param bcons defines the linear constraints in the form of Acons%*%x = bcons
#' @param verb defines how much information of the optimization routine to be printed. Default is set at 1 for no print, change to 5 for full output print
#' @return bound is the lower and upper bound of the value function, phat is the optimized solution for lower and upper bound, status records the status of the two optimization problem
#' @examples
#'  ## min/max_x obj'x s.t. Acons %*% x = bcons
#' @export
getbound <- function(obj, Acons, bcons, verb = 1){
  # solve for upper and lower bound of the parameter of interest
  # use Rmosek
  lb = list()
  lb$sense = "min"
  lb$c = obj
  blx = c(rep(0, ncol(Acons)))
  bux = c(rep(1, ncol(Acons)))
  lb$bx = rbind(blx, bux)
  lb$A = as(Acons, 'dgCMatrix')
  buc = c(bcons)
  blc = c(bcons)
  lb$bc = rbind(blc, buc)
  rlb = mosek(lb, list(verbose = verb))
  lbstatus = rlb$sol$bas$solsta
  if (lbstatus == "OPTIMAL"){
    lbphat = rlb$sol$bas$xx
    lowerbound = obj%*% rlb$sol$bas$xx
  }else{
    lbphat = rep(NA, length(obj))
    lowerbound = NA
  }
  ub = list()
  ub$sense = "max"
  ub$c = obj
  blx = c(rep(0, ncol(Acons)))
  bux = c(rep(1, ncol(Acons)))
  ub$bx = rbind(blx, bux)
  ub$A = as(Acons, 'dgCMatrix')
  buc = c(bcons)
  blc = c(bcons)
  ub$bc = rbind(blc, buc)
  rub = mosek(ub, list(verbose = verb))
  ubstatus = rub$sol$bas$solsta
  if (ubstatus == "OPTIMAL"){
    ubphat = rub$sol$bas$xx
    upperbound = obj%*%rub$sol$bas$xx
  }else{
    ubphat = rep(NA, length(obj))
    upperbound = NA
  }
  list(bound = c(lowerbound, upperbound), phat = cbind(lbphat, ubphat), status = c(lbstatus, ubstatus))
}

#' getbound_gurobi
#'
#' calls gurobi to solve for two linear programming problem for the identified set
#'
#' @param obj defines the linear coefficients of the linear objective function
#' @param Acons defines the linear constraints in the form of Acons%*%x = bcons
#' @param bcons defines the linear constraints in the form of Acons%*%x = bcons
#' @param verb defines how much information of the optimization routine to be printed. Default is set at 1 for no print, change to 5 for full output print
#' @return bound is the lower and upper bound of the value function, phat is the optimized solution for lower and upper bound, status records the status of the two optimization problem
#' @examples
#'  ## min/max_x obj'x s.t. Acons %*% x = bcons
#' @export
getbound_gurobi <- function(obj, Acons, bcons, verb = 1){
  # solve for upper and lower bound of the parameter of interest
  # use Gurobi
  LB = list()
  LB$A = Acons
  LB$obj = obj
  LB$sense = rep("=", nrow(Acons))
  LB$rhs = bcons
  LB$modelsense = "min"
  LB$lb = rep(0, ncol(Acons))
  LB$ub = rep(1, ncol(Acons))
  if (verb == 1){
    par <- list(OutputFlag = 0)
  }else{
    par = NULL
  }
  rLB = gurobi(LB, params = par)
  if (rLB$status =="OPTIMAL"){
    lowerbound = rLB$objval
    lbphat = rLB$x
  }else{
    lbphat = rep(NA, length(obj))
    lowerbound = NA
  }
  UB = list()
  UB$A = Acons
  UB$obj = obj
  UB$sense = rep("=", nrow(Acons))
  UB$rhs = bcons
  UB$modelsense = "max"
  UB$lb = rep(0, ncol(Acons))
  UB$ub = rep(1, ncol(Acons))
  if (verb == 1){
    par <- list(OutputFlag = 0)
  }else{
    par = NULL
  }
  rUB = gurobi(UB, params = par)
  if (rUB$status =="OPTIMAL"){
    upperbound = rUB$objval
    ubphat = rUB$x
  }else{
    ubphat = rep(NA, length(obj))
    upperbound = NA
  }
  list(bound = c(lowerbound, upperbound), phat = cbind(lbphat, ubphat))
}

#' getbound_consistent
#'
#' calls Rmosek to solve for two linear programming problem for the identified set with a small slack introduced on the constraints for consistent plug-in estimator of the identified set
#'
#' @param obj defines the linear coefficients of the linear objective function
#' @param constraints output from make_constraints() - defines the linear constraints in the form of A%*%x - b <= b_n
#' @param verb defines how much information of the optimization routine to be printed. Default is set at 1 for no print, change to 5 for full output print
#' @param slack defines the amount to slack b_n in the equality constraints to get consistent plug-in estimator of the identified set.
#' @return bound is the lower and upper bound of the value function, phat is the optimized solution for lower and upper bound, status records the status of the two optimization problem
#' @examples
#'  ## min/max_x obj'x
#'  s.t. Acons %*% x - bcons <= slack
#'       Acons %*% x - bcons >= -slack
#' @export
getbound_consistent <- function(obj, constraints, slack, verb = 1){
  # solve for bound
  # use Rmosek

  Acons = constraints$Aconsbig
  bcons = constraints$bconsbig
  Acons = Matrix(Acons, sparse = TRUE)

  lb = list()
  lb$sense = "min"
  lb$c = obj
  blx = c(rep(0, ncol(Acons)))
  bux = c(rep(1, ncol(Acons)))
  lb$bx = rbind(blx, bux)
  lb$A = as(rbind(Acons,-Acons), 'dgCMatrix')
  buc = c(bcons+slack, -bcons+slack)
  blc = rep(-Inf, 2*length(bcons))
  lb$bc = rbind(blc, buc)
  rlb = mosek(lb, list(verbose = verb))
  lbstatus = rlb$sol$bas$solsta
  if (lbstatus == "OPTIMAL"){
    lbphat = rlb$sol$bas$xx
    lowerbound = obj%*% rlb$sol$bas$xx
  }else{
    lbphat = rep(NA, length(obj))
    lowerbound = NA
  }
  ub = list()
  ub$sense = "max"
  ub$c = obj
  blx = c(rep(0, ncol(Acons)))
  bux = c(rep(1, ncol(Acons)))
  ub$bx = rbind(blx, bux)
  ub$A = as(rbind(Acons,-Acons), 'dgCMatrix')
  buc = c(bcons + slack, -bcons+slack)
  blc = rep(-Inf, 2*length(bcons))
  ub$bc = rbind(blc, buc)
  rub = mosek(ub, list(verbose = verb))
  ubstatus = rub$sol$bas$solsta
  if (ubstatus == "OPTIMAL"){
    ubphat = rub$sol$bas$xx
    upperbound = obj%*%rub$sol$bas$xx
  }else{
    ubphat = rep(NA, length(obj))
    upperbound = NA
  }
  list(bound = c(lowerbound, upperbound), phat = cbind(lbphat, ubphat), status = c(lbstatus, ubstatus))
}


#' gen_perturb
#'
#' generate perturbation to apply Cho and Russell (2020)
#'
#' @param obj defines the linear coefficients of the linear objective function
#' @param Acons defines the linear constraints in the form of Acons%*%x - bcons <= b_n
#' @return per_b1 and per_b2 is the perturbation of two inequality (implied from the equality constraint) constraints, each of length equals to the number of constraints (i.e. nrows(Acons)). per_psi is the perturbation of the objective, of length equals to the number of variables (i.e. ncol(Acons)). per_r1 and per_r2 is the perturbation on the support of the variables, each of length the number of variables.
#' @export
gen_perturb <- function(Acons, perturb_support = c(-1e-06, 1e-06, 0, 1e-06)){
  # perturb support contains four values: LB and UB for perturb_psi;  LB and UB for perturb_b, perturb_psi is the perturbation in the objective, perturb_b is the perturbation in the constraints.
  # need to convert all equality constraints of the form Acons %*% x = bcons, into Acons %*% x - bcons <= perturb_b1, -Acons %*% x + bcons <= perturb_b2
  dtheta = ncol(Acons)
  k <- nrow(Acons)
  perturb_b <- matrix(runif(2*k, perturb_support[3], perturb_support[4]), ncol = 2)
  perturb_b1 <- perturb_b[,1]
  perturb_b2 <- perturb_b[,2]
  perturb_range <- matrix(runif(2*dtheta, perturb_support[3], perturb_support[4]), ncol = 2)
  perturb_r1 <- perturb_range[,1]
  perturb_r2 <- perturb_range[,2]
  perturb_psi <- runif(dtheta, perturb_support[1], perturb_support[2])
  return(list(per_b1 = perturb_b1, per_b2 = perturb_b2, per_psi = perturb_psi, per_r1 = perturb_r1, per_r2 = perturb_r2))
}

#' getbound_per
#'
#' Use Rmosek to solve four perturbed LP problems as proposed in Cho and Russell (2020) to get bias-corrected estimator for the identified set.
#'
#' @param obj defines the linear coefficients of the linear objective function
#' @param Acons defines the linear constraints in the form of Acons%*%x = bcons
#' @param pervar defines the perturbation as in Cho and Russell (2020).
#' @return optimal value of the objective functions for the lower_minus, lower_plus, upper_minus, upper_plus LP problems. Definition see Cho and Russell (2020). phat are the associated optimal solutions and status is the solution status of the optimization problem.
#' #' @examples
#'  Suppose the original LP is min/max_x c'x s.t. A%*%x = b, Cho and Russell (2020) solves four perturbed LP. Details see Cho and Russell (2020).
#' @export
getbound_per <- function(obj, Acons, bcons, pervar, verb = 1){
  # solves for the four perturbed LPs.
  # use Rmosek
  lb_m = list()
  lb_m$sense = "min"
  lb_m$c = obj - pervar$per_psi
  blx = c(rep(0, ncol(Acons))) -pervar$per_r2
  bux = c(rep(1, ncol(Acons))) + pervar$per_r1
  lb_m$bx = rbind(blx, bux)
  lb_m$A = as(rbind(Acons, -Acons), 'dgCMatrix')
  buc = c(bcons+pervar$per_b1, -bcons+pervar$per_b2)
  blc = c(rep(-Inf, 2*length(bcons)))
  lb_m$bc = rbind(blc, buc)
  lb_m$iparam <- list(INFEAS_REPORT_AUTO="ON")
  rlb_m = mosek(lb_m, list(verbose = verb))
  lbstatus_m = rlb_m$sol$bas$solsta
  if (lbstatus_m == "OPTIMAL"){
    lbphat_m = rlb_m$sol$bas$xx
    lowerbound_m = (obj - pervar$per_psi)%*% rlb_m$sol$bas$xx
  }else{
    lbphat_m = rep(NA, length(obj))
    lowerbound_m = NA
  }
  ub_m = list()
  ub_m$sense = "max"
  ub_m$c = obj- pervar$per_psi
  blx = c(rep(0, ncol(Acons))) -pervar$per_r2
  bux = c(rep(1, ncol(Acons))) + pervar$per_r1
  ub_m$bx = rbind(blx, bux)
  ub_m$A = as(rbind(Acons, -Acons), 'dgCMatrix')
  buc = c(bcons+pervar$per_b1, -bcons+pervar$per_b2)
  blc = c(rep(-Inf, 2*length(bcons)))
  ub_m$bc = rbind(blc, buc)
  ub_m$iparam <- list(INFEAS_REPORT_AUTO="ON")
  rub_m = mosek(ub_m, list(verbose = verb))
  ubstatus_m = rub_m$sol$bas$solsta
  if (ubstatus_m == "OPTIMAL"){
    ubphat_m = rub_m$sol$bas$xx
    upperbound_m = (obj- pervar$per_psi)%*%rub_m$sol$bas$xx
  }else{
    ubphat_m = rep(NA, length(obj))
    upperbound_m = NA
  }
  lb_p = list()
  lb_p$sense = "min"
  lb_p$c = obj + pervar$per_psi
  blx = c(rep(0, ncol(Acons))) -pervar$per_r2
  bux = c(rep(1, ncol(Acons))) + pervar$per_r1
  lb_p$bx = rbind(blx, bux)
  lb_p$A = as(rbind(Acons, -Acons), 'dgCMatrix')
  buc = c(bcons+pervar$per_b1, -bcons+pervar$per_b2)
  blc = c(rep(-Inf, 2*length(bcons)))
  lb_p$bc = rbind(blc, buc)
  lb_p$iparam <- list(INFEAS_REPORT_AUTO="ON")
  rlb_p = mosek(lb_p, list(verbose = verb))
  lbstatus_p = rlb_p$sol$bas$solsta
  if (lbstatus_p == "OPTIMAL"){
    lbphat_p = rlb_p$sol$bas$xx
    lowerbound_p = (obj + pervar$per_psi)%*% rlb_p$sol$bas$xx
  }else{
    lbphat_p = rep(NA, length(obj))
    lowerbound_p = NA
  }
  ub_p = list()
  ub_p$sense = "max"
  ub_p$c = obj + pervar$per_psi
  blx = c(rep(0, ncol(Acons))) -pervar$per_r2
  bux = c(rep(1, ncol(Acons))) + pervar$per_r1
  ub_p$bx = rbind(blx, bux)
  ub_p$A = as(rbind(Acons, -Acons), 'dgCMatrix')
  buc = c(bcons+pervar$per_b1, -bcons+pervar$per_b2)
  blc = c(rep(-Inf, 2*length(bcons)))
  ub_p$bc = rbind(blc, buc)
  ub_p$iparam <- list(INFEAS_REPORT_AUTO="ON")
  rub_p = mosek(ub_p, list(verbose = verb))
  ubstatus_p = rub_p$sol$bas$solsta
  if (ubstatus_p == "OPTIMAL"){
    ubphat_p = rub_p$sol$bas$xx
    upperbound_p = (obj + pervar$per_psi)%*%rub_p$sol$bas$xx
  }else{
    ubphat_p = rep(NA, length(obj))
    upperbound_p = NA
  }
  list(bound = c(lowerbound_m, lowerbound_p, upperbound_m, upperbound_p), phat = cbind(lbphat_m, lbphat_p, ubphat_m, ubphat_p), status = c(lbstatus_m, lbstatus_p, ubstatus_m, ubstatus_p))
}

#' getbound_per_gurobi
#'
#' Use gurobi to solve four perturbed LP problems as proposed in Cho and Russell (2020) to get bias-corrected estimator for the identified set.
#'
#' @param obj defines the linear coefficients of the linear objective function
#' @param Acons defines the linear constraints in the form of Acons%*%x = bcons
#' @param pervar defines the perturbation as in Cho and Russell (2020).
#' @return optimal value of the objective functions for the lower_minus, lower_plus, upper_minus, upper_plus LP problems. Definition see Cho and Russell (2020). phat are the associated optimal solutions and status is the solution status of the optimization problem.
#' #' @examples
#'  Suppose the original LP is min/max_x c'x s.t. A%*%x = b, Cho and Russell (2020) solves four perturbed LP. Details see Cho and Russell (2020).
#' @export
getbound_per_gurobi <- function(obj, Acons, bcons, pervar, verb = 1){
  # solves for the four perturbed LPs.
  # use Gurobi
  if (verb == 1){
    par <- list(OutputFlag = 0)
  }else{
    par = NULL
  }
  constraintform = rep("<=", 2*length(bcons))
  lb_m = list()
  lb_m$modelsense = "min"
  lb_m$obj = obj - pervar$per_psi
  lb_m$lb = c(rep(0, ncol(Acons))) -pervar$per_r2
  lb_m$ub = c(rep(1, ncol(Acons))) + pervar$per_r1
  lb_m$A = as(rbind(Acons, -Acons), 'dgCMatrix')
  lb_m$rhs = c(bcons + pervar$per_b1, -bcons+pervar$per_b2)
  lb_m$sense = constraintform
  rlb_m = gurobi(lb_m, params = par)
  if (rlb_m$status == "OPTIMAL"){
    lbphat_m = rlb_m$x
    lowerbound_m = rlb_m$objval
  }else{
    lbphat_m = rep(NA, length(obj))
    lowerbound_m = NA
  }
  ub_m = list()
  ub_m$modelsense = "max"
  ub_m$obj = obj- pervar$per_psi
  ub_m$lb = c(rep(0, ncol(Acons))) -pervar$per_r2
  ub_m$ub = c(rep(1, ncol(Acons))) + pervar$per_r1
  ub_m$A = as(rbind(Acons, -Acons), 'dgCMatrix')
  ub_m$rhs = c(bcons+pervar$per_b1, -bcons+pervar$per_b2)
  ub_m$sense = constraintform
  rub_m = gurobi(ub_m, params = par)
  if (rub_m$status == "OPTIMAL"){
    ubphat_m = rub_m$x
    upperbound_m = rub_m$objval
  }else{
    ubphat_m = rep(NA, length(obj))
    upperbound_m = NA
  }
  lb_p = list()
  lb_p$modelsense = "min"
  lb_p$obj = obj + pervar$per_psi
  lb_p$lb = c(rep(0, ncol(Acons))) -pervar$per_r2
  lb_p$ub = c(rep(1, ncol(Acons))) + pervar$per_r1
  lb_p$A = as(rbind(Acons, -Acons), 'dgCMatrix')
  lb_p$rhs = c(bcons+pervar$per_b1, -bcons+pervar$per_b2)
  lb_p$sense = constraintform
  rlb_p = gurobi(lb_p, params = par)
  if (rlb_p$status == "OPTIMAL"){
    lbphat_p = rlb_p$x
    lowerbound_p = rlb_p$objval
  }else{
    lbphat_p = rep(NA, length(obj))
    lowerbound_p = NA
  }
  ub_p = list()
  ub_p$modelsense = "max"
  ub_p$obj = obj + pervar$per_psi
  ub_p$lb = c(rep(0, ncol(Acons))) -pervar$per_r2
  ub_p$ub = c(rep(1, ncol(Acons))) + pervar$per_r1
  ub_p$A = as(rbind(Acons, -Acons), 'dgCMatrix')
  ub_p$rhs = c(bcons+pervar$per_b1, -bcons+pervar$per_b2)
  ub_p$sense = constraintform
  rub_p = gurobi(ub_p, params = par)
  if (rub_p$status == "OPTIMAL"){
    ubphat_p = rub_p$x
    upperbound_p = rub_p$objval
  }else{
    ubphat_p = rep(NA, length(obj))
    upperbound_p = NA
  }
  list(bound = c(lowerbound_m, lowerbound_p, upperbound_m, upperbound_p), phat = cbind(lbphat_m, lbphat_p, ubphat_m, ubphat_p))
}

#' getbound_per_relax
#'
#' In the original LP, we allow relaxation of the IV assumption: Aextra %*% x = 0 to be
#' Aextra %*% x <= bextra and
#' Aextra %*% x >= -bextra
#' while maintaining the equality constraints Acons %*% x = bcons.
#' And then we solve the four perturbed LP problems as in Cho and Russell (2020) using Rmosek.
#' @param obj defines the linear coefficients of the linear objective function
#' @param Acons defines the linear constraints in the form of Acons%*%x = bcons
#' @param bcons defines the linear constraints in the form of Acons%*%x = bcons
#' @param Aextra defines the relaxed linear equality constraints in the form Aextra %*% x <= bextra and Aextra %% x >= -bextra
#' @param bextra defines the amount of relaxation of the equality constraints.
#' @param pervar defines the perturbation as in Cho and Russell (2020).
#' @return optimal value of the objective functions for the lower_minus, lower_plus, upper_minus, upper_plus LP problems. Definition see Cho and Russell (2020). phat are the associated optimal solutions and status is the solution status of the optimization problem.
#' #' @examples
#' @examples
#' Original LP is defined as
#'  ## min/max_x obj'x
#'  s.t. Acons %*% x - bcons = 0
#'       Aextra %*% x <= bextra
#'       Aextra %*% x >= -bextra
#' We then define the four perturbed LP as in Cho and Russell (2020).
getbound_per_relax <- function(obj, Acons, bcons, Aextra, bextra, pervar, verb = 1){
  # solves for the four relaxed perturbed LPs.
  # Acons and bcons contains model and add up constraint [these constraints are not relaxed]
  # Aextra and bextra contains indep or/and mono constraint of the form Aextra %*% x <= bextra, -Aextra %*% x <= bextra [these are relaxed from equality to inequality with a slack
  lb_m = list()
  lb_m$sense = "min"
  lb_m$c = obj - pervar$per_psi
  blx = c(rep(0, ncol(Acons))) -pervar$per_r2
  bux = c(rep(1, ncol(Acons))) + pervar$per_r1
  lb_m$bx = rbind(blx, bux)
  lb_m$A = as(rbind(Acons, -Acons, Aextra, -Aextra), 'dgCMatrix')
  buc = c(bcons+pervar$per_b1[1:nrow(Acons)], -bcons+pervar$per_b2[1:nrow(Acons)], bextra + pervar$per_b1[-c(1:nrow(Acons))], bextra + pervar$per_b2[-c(1:nrow(Acons))])
  blc = c(rep(-Inf, 2*length(bcons)), rep(-Inf, 2*length(bextra)))
  lb_m$bc = rbind(blc, buc)
  lb_m$iparam <- list(INFEAS_REPORT_AUTO="ON")
  rlb_m = mosek(lb_m, list(verbose = verb))
  lbstatus_m = rlb_m$sol$bas$solsta
  if (lbstatus_m == "OPTIMAL"){
    lbphat_m = rlb_m$sol$bas$xx
    lowerbound_m = (obj - pervar$per_psi)%*% rlb_m$sol$bas$xx
  }else{
    lbphat_m = rep(NA, length(obj))
    lowerbound_m = NA
  }
  ub_m = list()
  ub_m$sense = "max"
  ub_m$c = obj- pervar$per_psi
  blx = c(rep(0, ncol(Acons))) -pervar$per_r2
  bux = c(rep(1, ncol(Acons))) + pervar$per_r1
  ub_m$bx = rbind(blx, bux)
  ub_m$A = as(rbind(Acons, -Acons, Aextra, -Aextra), 'dgCMatrix')
  buc = c(bcons+pervar$per_b1[1:nrow(Acons)], -bcons+pervar$per_b2[1:nrow(Acons)], bextra + pervar$per_b1[-c(1:nrow(Acons))], bextra + pervar$per_b2[-c(1:nrow(Acons))])
  blc = c(rep(-Inf, 2*length(bcons)), rep(-Inf, 2*length(bextra)))
  ub_m$bc = rbind(blc, buc)
  ub_m$iparam <- list(INFEAS_REPORT_AUTO="ON")
  rub_m = mosek(ub_m, list(verbose = verb))
  ubstatus_m = rub_m$sol$bas$solsta
  if (ubstatus_m == "OPTIMAL"){
    ubphat_m = rub_m$sol$bas$xx
    upperbound_m = (obj- pervar$per_psi)%*%rub_m$sol$bas$xx
  }else{
    ubphat_m = rep(NA, length(obj))
    upperbound_m = NA
  }
  lb_p = list()
  lb_p$sense = "min"
  lb_p$c = obj + pervar$per_psi
  blx = c(rep(0, ncol(Acons))) -pervar$per_r2
  bux = c(rep(1, ncol(Acons))) + pervar$per_r1
  lb_p$bx = rbind(blx, bux)
  lb_p$A = as(rbind(Acons, -Acons, Aextra, -Aextra), 'dgCMatrix')
  buc = c(bcons+pervar$per_b1[1:nrow(Acons)], -bcons+pervar$per_b2[1:nrow(Acons)], bextra + pervar$per_b1[-c(1:nrow(Acons))], bextra + pervar$per_b2[-c(1:nrow(Acons))])
  blc = c(rep(-Inf, 2*length(bcons)), rep(-Inf, 2*length(bextra)))
  lb_p$bc = rbind(blc, buc)
  lb_p$iparam <- list(INFEAS_REPORT_AUTO="ON")
  rlb_p = mosek(lb_p, list(verbose = verb))
  lbstatus_p = rlb_p$sol$bas$solsta
  if (lbstatus_p == "OPTIMAL"){
    lbphat_p = rlb_p$sol$bas$xx
    lowerbound_p = (obj + pervar$per_psi)%*% rlb_p$sol$bas$xx
  }else{
    lbphat_p = rep(NA, length(obj))
    lowerbound_p = NA
  }
  ub_p = list()
  ub_p$sense = "max"
  ub_p$c = obj + pervar$per_psi
  blx = c(rep(0, ncol(Acons))) -pervar$per_r2
  bux = c(rep(1, ncol(Acons))) + pervar$per_r1
  ub_p$bx = rbind(blx, bux)
  ub_p$A = as(rbind(Acons, -Acons, Aextra, -Aextra), 'dgCMatrix')
  buc = c(bcons+pervar$per_b1[1:nrow(Acons)], -bcons+pervar$per_b2[1:nrow(Acons)], bextra + pervar$per_b1[-c(1:nrow(Acons))], bextra + pervar$per_b2[-c(1:nrow(Acons))])
  blc = c(rep(-Inf, 2*length(bcons)), rep(-Inf, 2*length(bextra)))
  ub_p$bc = rbind(blc, buc)
  ub_p$iparam <- list(INFEAS_REPORT_AUTO="ON")
  rub_p = mosek(ub_p, list(verbose = verb))
  ubstatus_p = rub_p$sol$bas$solsta
  if (ubstatus_p == "OPTIMAL"){
    ubphat_p = rub_p$sol$bas$xx
    upperbound_p = (obj + pervar$per_psi)%*%rub_p$sol$bas$xx
  }else{
    ubphat_p = rep(NA, length(obj))
    upperbound_p = NA
  }
  list(bound = c(lowerbound_m, lowerbound_p, upperbound_m, upperbound_p), phat = cbind(lbphat_m, lbphat_p, ubphat_m, ubphat_p), status = c(lbstatus_m, lbstatus_p, ubstatus_m, ubstatus_p))
}


#' ATEobj
#'
#' Constructs the objective function to bound the ATE
#'
#' @param processed_data A bcDAT object which is the output of the process_data() function.
#' @param constraints A list of constraints which is the output of the make_constraints() function. Only the hyperplane arrangement is used here.
#' @return A numerical vector containing values of the objective function used to bound the ATE.
#' @examples
#' # process_data1 is output from process_data()
#' # cons0 is output from make_constraints()
#' objlam_ATE = ATEobj(processed_data = process_data1,constraints = cons0)
#' @export
ATEobj <- function(processed_data,constraints){
  # zcol should be the columns in the XZsupport that belongs to the exogenous variables
  # if P_yxzw has NA, set it to zero so that they drop out of the objective

  jset1 <- processed_data$counterfactual$jset1
  jset0 <- processed_data$counterfactual$jset0
  fsub <- constraints$f
  YXZWsupport <- processed_data$support %>% select(-c(P_all,P_exo,P_y1))
  if(processed_data$case_list$insvars==F){
  XZsupport <- YXZWsupport %>% select(-c(y)) %>% unique()
  }else if(processed_data$case_list$insvars==T){
  XZsupport <- YXZWsupport %>% select(-c(y)) %>% select(-all_of(ivlist)) %>% unique()
  }
  P_yxzw <- processed_data$support$P_all
  zcol <- which(names(XZsupport)%in%processed_data$exolist)

  P_yxzw[is.na(P_yxzw)] = 0
  ATEcoef1 = NULL
  for (j in 1:length(jset1)){
    counf.index = jset1[j]
    sign = 0.5 * (fsub$SignVector[counf.index,]+1)
    signe = NULL
    subset = which(apply(YXZWsupport, 1, function(x) identical(as.numeric(x[zcol+1]), as.numeric(XZsupport[counf.index,zcol]))))
    for (q in 1:(length(subset)/2)){
      signe = c(signe, kronecker(sign,c(P_yxzw[subset[(q-1)*2 + c(1:2)]])))
    }
    ATEcoef1 = c(ATEcoef1,signe)
  }

  ATEcoef0 = NULL
  for (j in 1:length(jset0)){
    counf.index = jset0[j]
    sign = 0.5 * (fsub$SignVector[counf.index,]+1)
    signe = NULL
    subset = which(apply(YXZWsupport, 1, function(x) identical(as.numeric(x[zcol+1]), as.numeric(XZsupport[counf.index,zcol]))))
    for (q in 1:(length(subset)/2)){
      signe = c(signe, kronecker(sign,c(P_yxzw[subset[(q-1)*2 + c(1:2)]])))
    }
    ATEcoef0 = c(ATEcoef0,signe)
  }
  return(ATEcoef1 - ATEcoef0)
}

#' CCPy0obj
#'
#' Constructs the objective function to bound the CCP for Y=0
#'
#' @param processed_data A bcDAT object which is the output of the process_data() function.
#' @param constraints A list of constraints which is the output of the make_constraints() function. Only the hyperplane arrangement is used here.
#' @return A numerical vector containing values of the objective function used to bound the ATE.
#' @examples
#' objlam_CCPy0 = CCPy0obj(processed_data = process_data1,constraints = cons0)
#' @export
CCPy0obj = function(processed_data,constraints){
  # zcol should be the columns in the XZsupport that belongs to the exogenous variables
  # CCP_y0: P(varphi(1, z, theta,beta) >=0| Y = 0, X = 0) = \sum_{z, w} P(varphi(1, z,theta,beta)>=0|Y = 0, X = 0, Z = z, W = w) P(Z = z, W = w|Y = 0, X = 0)
  # ATEcoef1 takes the form P(varphi(1, z, theta,beta> =0)) = \sum_{y,x,z,w} P(varphi(1,z,theta,beta)>=0|Y = y, X = x, Z = z, W = w) P(Y = y, X = x, Z = z, W = w)
  # hence one strategy to get the CCP_y0 objective is to set P(Y = 0, X = 0, Z = z, W = w) to be its value and the rest to zero, and divide all vector by P(Y = 0, X = 0), which is simply to compute
  # If P_yxzw has NA, meaning the values (y,x,z,w) is not realized in the data, then set it to zero so that it drops out of the objective

  jset1 <- processed_data$counterfactual$jset1
  jset0 <- processed_data$counterfactual$jset0
  fsub <- constraints$f
  datsub <- processed_data$data
  YXZWsupport <- processed_data$support %>% select(-c(P_all,P_exo,P_y1))
  if(processed_data$case_list$insvars==F){
    XZsupport <- YXZWsupport %>% select(-c(y)) %>% unique()
  }else if(processed_data$case_list$insvars==T){
    XZsupport <- YXZWsupport %>% select(-c(y)) %>% select(-all_of(ivlist)) %>% unique()
  }
  P_yxzw <- processed_data$support$P_all
  zcol <- which(names(XZsupport)%in%processed_data$exolist)
  names(datsub) <- names(YXZWsupport)

  CCPcoefy0 = NULL
  P_y0x0 = sum(datsub$y==0 & datsub$ins==0)/nrow(datsub)
  P_yxzw[is.na(P_yxzw)] = 0
  for (j in 1:length(jset1)){
    counf.index = jset1[j]
    sign = 0.5 * (fsub$SignVector[counf.index,]+1)
    signe = NULL
    subset = which(apply(YXZWsupport, 1, function(x) identical(as.numeric(x[zcol +1]), as.numeric(XZsupport[counf.index,zcol]))))
    for (q in 1:(length(subset)/2)){
      tempq = subset[(q-1)*2 + c(1:2)]
      tempp = P_yxzw[tempq]
      tempd = YXZWsupport[tempq,]
      tempi = which(tempd[,1]!=0|tempd[,2]!=0)
      tempp[tempi]=0  # set P(Y = y, X = x, Z = z, W = w) = 0 for y!= 1 or x!=0
      tempp = tempp/P_y0x0  # normalize by P(Y = 1, X = 0)
      signe = c(signe, kronecker(sign,c(tempp)))
    }
    CCPcoefy0 = c(CCPcoefy0,signe)
  }
  return(CCPcoefy0)
}
