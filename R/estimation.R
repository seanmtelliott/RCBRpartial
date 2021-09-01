
#' getbound
#'
#' description
#'
#' @param obj
#' @param Acons description
#' @param bcons description
#' @param verb description
#' @return something returned
#' @examples
#'  ## Example here
#' example code
#' @export
getbound <- function(obj, Acons, bcons, verb = 1){
  # solve for bound
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


getbound_gurobi <- function(obj, Acons, bcons, verb = 1){
  # solve for bound
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

gen_perturb <- function(Acons, perturb_support = c(-1e-06, 1e-06, 0, 1e-06)){
  # perturb support contains four values: LB and UB for perturb_psi;  LB and UB for perturb_b, perturb_psi is the perturbation in the pobjective, perturb_b is the perturbation in the constraints.
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
