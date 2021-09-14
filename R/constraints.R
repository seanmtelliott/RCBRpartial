#' Make Constrains
#'
#'
#' This function constructs constraints and bounds. The behaviour of this function is determined by the formula provided by the user
#' in the process_data() function.
#'
#' @param data A bcDAT object - the output from process_data()
#' @param f The arrangement sign vector.
#' @param eliminate A list containing the monotonicity constraints.
#' @return The A matrix and b vector for the system of linear equations Ax=b
#' @examples
#' ## The general case where no functional form is imposed.
#' fsub <- list(SignVector = t(expand.grid(rep(list(c(-1,1)), nrow(XZsupport)))))
#' eliset <- list( c(1, 5), c(1,7),c(3, 7),c(3,5))
#' constraints <- make_constraints(data = processed_data, f = fsub, eliminate = eliset)
#' @export
make_constraints = function(data,f,eliminate){

  # Rename objects from data processing module
  OP <- data$support$P_y1 %>% na.omit() %>% unlist() %>% as.numeric()
  pyxz <- data$support$P_exo %>% unlist() %>% as.numeric()
  ind <- data$exolist

  if(data$case_list$insvars==F){
    Hsupport <- data$support %>% select(-c(y,P_all,P_exo,P_y1)) %>% unique()
    Arrangesupport <- data$support %>% select(-c(y,P_all,P_exo,P_y1)) %>% unique()
    excl <- NULL
    ivexo = TRUE
  }else if(data$case_list$insvars==T){
    Hsupport <- data$support  %>% select(-c(y,P_all,P_exo,P_y1)) %>% unique()
    Arrangesupport <- data$support %>% select(-data$ivlist) %>% select(-c(y,P_all,P_exo,P_y1)) %>% unique()
    excl <- data$ivlist
    ivexo <- ifelse(data$ivlist %in% data$exolist,T,F)
  }


  # eliminate provides restrictions on the signVector that needs to be impose to have zero probability, through monotonicity assumptions
  # eliminate encodes in the form of a list as: (i, j) , which means sign in row i > sign in row j of SignVector
  # if OP constrains NA, set to zero.
  OP[is.na(OP)] = 0
  AA1 = NULL
  AA2 = NULL
  AA3 = NULL
  AA4 = NULL
  bb1 = NULL
  bb2 = NULL
  bb3 = NULL
  bb4 = NULL
  if (!length(eliminate) & !is.list(eliminate)) stop("eliminate needs to be a list")
  Mstar = ncol(f$SignVector)
  AA1 = makecons.obs_Y2(Hsupport, Arrangesupport, f = f, OP = OP, excl = excl)
  bb1 = rep(1, length(OP)*2)  # changed
  AA3 =kronecker( diag(1, nrow(Hsupport)), rbind(c(rep(c(1,0), Mstar)), c(rep(c(0,1), Mstar)))) # sums to one constraint on P(cell | Y, X, Z)
  bb3 = rep(1, 2*nrow(Hsupport))
  AA3 = Matrix(AA3, sparse=TRUE)
  if (!length(ind) & !length(eliminate)){
    Aconsbig = rbind(AA1, AA3)
    bconsbig = c(bb1,bb3)
  }
  if (length(ind) > 0 & !length(eliminate)){
    AA2 = makecons.ind.test(Hsupport, ind = ind, f= f, pyxz = pyxz, ivexo)
    bb2 = rep(0, nrow(AA2))
    Aconsbig = rbind(AA1, AA2, AA3)
    bconsbig = c(bb1, bb2, bb3)
  }
  if (length(ind) > 0 & length(eliminate)>0){
    AA2 = makecons.ind.test(Hsupport, ind = ind,f = f, pyxz = pyxz,ivexo)
    bb2 = rep(0, nrow(AA2))
    eli = unlist(lapply(eliminate, function(x) which(f$SignVector[x[1],] > f$SignVector[x[2],])))
    celleli = rep(0, Mstar)
    celleli[unique(eli)] = 1
    celleli = kronecker(celleli, c(1,1))
    AA4 = t(kronecker(diag(ncol(AA1)/(2*Mstar)), celleli))
    bb4 = rep(0, nrow(AA4))
    Aconsbig = rbind(AA1, AA2, AA3, AA4)
    bconsbig = c(bb1, bb2, bb3, bb4)
  }
  if (!length(ind)  & length(eliminate) >0){
    eli = unlist(lapply(eliminate, function(x) which(f$SignVector[x[1],] > f$SignVector[x[2],])))
    celleli = rep(0, Mstar)
    celleli[unique(eli)] = 1
    celleli = kronecker(celleli, c(1,1))
    AA4 = t(kronecker(diag(ncol(AA1)/(2*Mstar)), celleli))
    bb4 = rep(0, nrow(AA4))
    Aconsbig = rbind(AA1, AA3, AA4)
    bconsbig = c(bb1, bb3, bb4)
  }
  list(Aconsbig = Aconsbig, bconsbig = bconsbig, AA.model = AA1, AA.indep = AA2, AA.addup = AA3, AA.mono = AA4, bb.model = bb1, bb.indep = bb2, bb.addup = bb3, bb.mono = bb4)
}


makecons.ind1 = function(Hsupport, ind, f, pyxz, ivexo=TRUE){
  # connstraints due to independence of variables Z
  # allows for both exogeneous variable or exclusion restriction [this influences the relevant arrangement since exclusion restriction do not enter the choice equation
  # Hsupport is the full support of (X, Z), as a data frame with variable names
  # f is the revalent arrangement: conditional on some beta value
  # PYXZ is the conditional probablity of (Y, X) | Z [i.e. P_exo]
  # ind is the list of variable names that are assumed to be independent
  indcol = which(names(Hsupport)%in% ind)
  if (!length(indcol)){
    Acons = NULL
  }
  if (length(indcol)>0){
    indset = Hsupport[,indcol,drop=FALSE]
    uniqueind = unique(indset)  # unique values of the independent variables
    if (dim(uniqueind)[1] < 2) stop("check support of ind variable")
    #Acons = matrix(0, (nrow(uniqueind)-1)*ncol(f$SignVector), ncol(f$SignVector) * length(pyxz))
    pmat = matrix(0, nrow(indset), 2* nrow(indset))
    mark1 = apply(indset, 1, function(x) (sum(x==uniqueind[1,])==length(x)))  # find the entry where the iv variable Z takes the first value Z = z1
    pmat[1,which(kronecker(mark1, c(1,1))==1)] = pyxz[which(kronecker(mark1, c(1,1))==1)]
    markj = list()
    for (j in 2:nrow(uniqueind)){
      markj[[j]] = apply(indset, 1, function(x) (sum(x == uniqueind[j,])== length(x)))
      pmat[j,which(kronecker(markj[[j]], c(1,1))==1)] = -pyxz[which(kronecker(markj[[j]], c(1,1))==1)]
    }
    Mcell = ncol(f$SignVector)
    Acons= Matrix(nrow = (nrow(uniqueind)-1)*Mcell, ncol = Mcell * length(pyxz), data = 0, sparse=TRUE)
    if (ivexo){
      # if iv is assumed to be exogenous
      for (j in 1:(nrow(uniqueind)-1)){
        pos1 = c(1,0.5*nrow(Hsupport)+1) # here we use the fact the endo variable is binary, take 2 values. If not, we would have to read out all odd entries of pmat[1,] not zero. we have nrow(Hsupport) number of support point, first half for X = 0, last half for X = 1
        pos2 = j + pos1
        tp1 = list(matrix(pmat[1,c(1:2)+(pos1[1]-1)*2],nrow=1))
        tp11 = rep(tp1, Mcell)
        temp1a = bdiag(tp11)
        tp2 = list(matrix(pmat[1,c(1:2)+(pos1[2]-1)*2],nrow=1))
        tp22 = rep(tp2, Mcell)
        temp1b = bdiag(tp22)
        tp3 = list(matrix(pmat[j+1, c(1:2)+(pos2[1]-1)*2], nrow=1))
        tp33 = rep(tp3, Mcell)
        temp1c = bdiag(tp33)
        tp4 = list(matrix(pmat[j+1, c(1:2)+(pos2[2]-1)*2], nrow=1))
        tp44 = rep(tp4, Mcell)
        temp1d = bdiag(tp44)
        Acons[(j-1)*Mcell + (1:Mcell),(pos1[1]-1)*2*Mcell + (1:(2*Mcell))] = temp1a
        Acons[(j-1)*Mcell + (1:Mcell),(pos1[2]-1)*2*Mcell + (1:(2*Mcell))] = temp1b
        Acons[(j-1)*Mcell + (1:Mcell),(pos2[1]-1)*2*Mcell + (1:(2*Mcell))] = temp1c
        Acons[(j-1)*Mcell + (1:Mcell),(pos2[2]-1)*2*Mcell + (1:(2*Mcell))] = temp1d
      }
    }else{
      # if not imposing IV exo, e.g., IV is conditional indep
      for (j in 2:nrow(uniqueind)){
        temp1 = matrix(unlist(lapply(asplit(matrix(pmat[1,],nrow = 2), 2), function(x) t(kronecker(diag(ncol(f$SignVector)), x)) )), nrow = ncol(f$SignVector))
        temp2 = matrix(unlist(lapply(asplit(matrix(pmat[j,],nrow = 2), 2), function(x) t(kronecker(diag(ncol(f$SignVector)), x)) )), nrow = ncol(f$SignVector))
        #Acons = rbind(Acons, temp1 + temp2)
        Acons[((j-2)*ncol(f$SignVector) + 1:ncol(f$SignVector)),] = temp1 + temp2
      }
    }
  }
  if (sum(is.na(pyxz))>0){
    Narows = which(rowSums(is.na(Acons))>0)
    Acons[Narows,] = 0
  }
  Acons
}

makecons.ind.test = function(Hsupport, ind, f, pyxz, ivexo=TRUE){
  # connstraints due to independence of variables Z
  # allows for both exogeneous variable or exclusion restriction [this influences the relevant arrangement since exclusion restriction do not enter the choice equation
  # Hsupport is the full support of (X, Z), as a data frame with variable names
  # f is the revalent arrangement: conditional on some beta value
  # PYXZ is the conditional probablity of (Y, X) | Z [i.e. P_exo]
  # ind is the list of variable names that are assumed to be independent
  indcol = which(names(Hsupport)%in% ind)
  if (!length(indcol)){
    Acons = NULL
  }
  if (length(indcol)>0){
    indset = Hsupport[,indcol,drop=FALSE]
    uniqueind = unique(indset)  # unique values of the independent variables
    if (dim(uniqueind)[1] < 2) stop("check support of ind variable")
    #Acons = matrix(0, (nrow(uniqueind)-1)*ncol(f$SignVector), ncol(f$SignVector) * length(pyxz))
    pmat = matrix(0, nrow(indset), 2* nrow(indset))
    mark1 = apply(indset, 1, function(x) (sum(x==uniqueind[1,])==length(x)))  # find the entry where the iv variable Z takes the first value Z = z1
    pmat[1,which(kronecker(mark1, c(1,1))==1)] = pyxz[which(kronecker(mark1, c(1,1))==1)]
    markj = list()
    for (j in 2:nrow(uniqueind)){
      markj[[j]] = apply(indset, 1, function(x) (sum(x == uniqueind[j,])== length(x)))
      pmat[j,which(kronecker(markj[[j]], c(1,1))==1)] = -pyxz[which(kronecker(markj[[j]], c(1,1))==1)]
    }
    Mcell = ncol(f$SignVector)
    Acons= Matrix(nrow = (nrow(uniqueind)-1)*Mcell, ncol = Mcell * length(pyxz), data = 0, sparse=TRUE)
    if (ivexo){
      # if iv is assumed to be exogenous
      # A matrix construction for makecons.ind1()

      mat_rows <- ncol(f$SignVector)

      #endo_var_count <- nrow(data$support %>% select(all.vars(formula(data$formula,rhs=1))[[2]]) %>% unique()) # Take number of unique points of X

      # By defining the zero matrix ahead of time it saves us having to calculate it on each iteration
      zero_block <- bdiag(replicate(mat_rows,t(c(0,0)),simplify = FALSE))

      # Select probabilities in first row and generate first matrix
      unique_z <- which(rep(c(1:(nrow(uniqueind) * 2)), nrow(unique(Hsupport %>% select(-ind)))) %in% c(1,2))
      initial_row <- matrix(pmat[1,unique_z],nrow=2)
      initial_mat <- apply(initial_row,2,function(x) bdiag(replicate(mat_rows, t(x), simplify = FALSE)))

      mat_list <- list()
      for(i in 1:(nrow(uniqueind)-1)){
        unique_z <- unique_z + 2
        matlistsub <- vector(mode = "list", length = ncol(pmat)/2)
        psub <- matrix(pmat[(i+1),unique_z],nrow=2) # is possible that probablilites are zero
        psub_mat <- apply(psub,2,function(x) bdiag(replicate(mat_rows, t(x), simplify = FALSE)))

        # Where to insert initial rows

        initial_values <- c(1:length(matlistsub))[seq(1,length(matlistsub),nrow(uniqueind))]
        h=1
        for(j in initial_values){
          matlistsub[[j]] <- initial_mat[[h]]
          h=h+1
        }

        # Inserting remaining non-zero matricies

        remaining_values <- initial_values + i

        h=1
        for(j in remaining_values){
          matlistsub[[j]] <- psub_mat[[h]]
          h=h+1
        }

        # Insert zero block everywhere else

        for(j in 1:length(matlistsub)){
          if(is.null(matlistsub[[j]])==T){
            matlistsub[[j]] <- zero_block
          }
        }

        mat_list[[i]] <- do.call(cbind,matlistsub)

      }

      Acons <- do.call(rbind,mat_list)


    }else{
      # if not imposing IV exo, e.g., IV is conditional indep
      for (j in 2:nrow(uniqueind)){
        temp1 = matrix(unlist(lapply(asplit(matrix(pmat[1,],nrow = 2), 2), function(x) t(kronecker(diag(ncol(f$SignVector)), x)) )), nrow = ncol(f$SignVector))
        temp2 = matrix(unlist(lapply(asplit(matrix(pmat[j,],nrow = 2), 2), function(x) t(kronecker(diag(ncol(f$SignVector)), x)) )), nrow = ncol(f$SignVector))
        #Acons = rbind(Acons, temp1 + temp2)
        Acons[((j-2)*ncol(f$SignVector) + 1:ncol(f$SignVector)),] = temp1 + temp2
      }
    }
  }
  if (sum(is.na(pyxz))>0){
    Narows = which(rowSums(is.na(Acons))>0)
    Acons[Narows,] = 0
  }
  Acons
}

makecons.cind <- function(XZWsupport, uniquecat, uniquecatindex, zvalues, f, P_all){
  # uniquecat is a list defining grouping of W values
  # uniquecatindex is the list of the correspondng index value of W values
  # P_all = P(Y = y, X = x, Z = z, W = w)
  Mstar = ncol(f$SignVector)
  #uniquecat = list()
  #uniquecat[[1]] = sort(unique(XZWsupport[which(XZWsupport[,4] <= 10),4]))
  #uniquecat[[2]] = sort(unique(XZWsupport[which(XZWsupport[,4]>10 & XZWsupport[,4]<=49),4]))
  #uniquecat[[3]] = sort(unique(XZWsupport[which(XZWsupport[,4]>49 & XZWsupport[,4]<=100),4]))
  #uniquecat[[4]] = sort(unique(XZWsupport[which(XZWsupport[,4] > 100),4]))
  #uniquecatindex = list()
  #uniquecatindex[[1]] = c(1,2)
  #uniquecatindex[[2]] = c(3,4,5)
  #uniquecatindex[[3]] = c(6,7,8,9)
  #uniquecatindex[[4]] = c(10,11)

  Acindep = NULL
  for (cc in 1:length(uniquecat)){
    Ablock = NULL
    for (k in 1:nrow(zvalues)){
      p = matrix(0, nrow(XZWsupport), 2*nrow(XZWsupport))
      mark1 = rep(FALSE, nrow(XZWsupport))
      mark1[which(XZWsupport[,2]==zvalues[k,1] & XZWsupport[,3]==zvalues[k,2] & XZWsupport[,4]==uniquecat[[cc]][1])] = TRUE
      p[1,which(kronecker(mark1, c(1,1))==1)] = P_all[which(kronecker(mark1, c(1,1))==1)]/(sum(datsub$married==zvalues[k,1]&datsub$healthy==zvalues[k,2]&datsub$emp ==uniquecat[[cc]][1])/nrow(datsub))
      markj = list()
      for (s in 2:length(uniquecat[[cc]])){
        markjt = rep(FALSE, nrow(XZWsupport))
        markjt[which(XZWsupport[,2]==zvalues[k,1] & XZWsupport[,3]==zvalues[k,2] & XZWsupport[,4]==uniquecat[[cc]][s])] = TRUE
        markj[[s]] = markjt
        p[s,which(kronecker(markj[[s]], c(1,1))==1)] = - P_all[which(kronecker(markj[[s]], c(1,1))==1)]/(sum(datsub$married==zvalues[k,1]&datsub$healthy==zvalues[k,2]&datsub$emp ==uniquecat[[cc]][s])/nrow(datsub))
      }
      Ablockz = Matrix(nrow = (length(uniquecat[[cc]])-1) * Mstar, ncol = Mstar * 2* nrow(XZWsupport), data = 0, sparse = TRUE)
      for (s in 2:length(uniquecat[[cc]])){
        pos1 = c(2*uniquecatindex[[cc]][1]-1 + (k-1)*22, 2*uniquecatindex[[cc]][1]-1+ (k-1)*22+nrow(XZWsupport))
        pos2 = pos1 + (s-1)*2
        tp1 = list(matrix(p[1,c(pos1[1],pos1[1]+1)],nrow=1))
        tp11 = rep(tp1, Mstar)
        temp1a = bdiag(tp11)
        tp2 = list(matrix(p[1,c(pos1[2],pos1[2]+1)],nrow=1))
        tp22 = rep(tp2, Mstar)
        temp1b = bdiag(tp22)
        tp3 = list(matrix(p[s, c(pos2[1], pos2[1]+1)], nrow=1))
        tp33 = rep(tp3, Mstar)
        temp1c = bdiag(tp33)
        tp4 = list(matrix(p[s, c(pos2[2],pos2[2]+1)], nrow=1))
        tp44 = rep(tp4, Mstar)
        temp1d = bdiag(tp44)
        Ablockz[((s-2)*Mstar+1:Mstar),(pos1[1]-1)*Mstar + (1:(2*Mstar))] = temp1a
        Ablockz[((s-2)*Mstar+1:Mstar), (pos1[2]-1)*Mstar + (1:(2*Mstar))] = temp1b
        Ablockz[((s-2)*Mstar+1:Mstar),(pos2[1]-1)*Mstar + (1:(2*Mstar))] = temp1c
        Ablockz[((s-2)*Mstar+1:Mstar),(pos2[2]-1)*Mstar + (1:(2*Mstar))] = temp1d
      }
      Ablock = rbind(Ablock, Ablockz)
    }
    Acindep = rbind(Acindep, Ablock)
  }
  if (sum(is.na(P_all))>0){
    Narows = which(rowSums(is.na(Acindep))>0)
    Acindep[Narows,] = 0
  }
  return(Acindep)
}

makecons.obs_Y2 = function(Hsupport, Arrangesupport, f, OP, excl){
  # version stated in the paper
  # Hsupport is the full support of all observable except for Y: (X, Z, W where W is instruments)
  # Arrangesupport is the support of all variables involved in the arrangement. (e.g. If Hsupport is on X,Z,W, and W is excl, then Arrangesupport will be support of X, Z, arranged in the same way as how Hyperplane arrangement is done.)
  # f is the relevant arrangement
  # obsp is the observed probability
  # cond.X is the value of X we condition on
  # OP: observed choice probability (with encoding P(Y = 1| Hsupport), mainly need to know which is the NA entry
  # output: constraint matrix on the vector of variable P(cell | X = x, Z = z)
  # if excl is of positive length, implies some of the Hsupport is used as instrument and does not enter into the choice equation, (i.e. existence of W in Hsupport)
  # if the input OP has NAs (meaning some values of Hsupport is not realized in the data), then Acons will have NAs corresponding to the constraints on the P(cells | Y = y, X = x, Z = z), set the NA to zero will mean that we anihilate the constraints on P(cells | Y = y, X = x, Z = z) for (y,x,z) not realized in the data.
  if (!length(excl)){
    matAl = list()
    for (j in 1:nrow(Hsupport)){
      if (!is.na(OP[j])){
        matAl[[j]] = rbind(kronecker(0.5*(f$SignVector[j,]+1), c(0, 1)), kronecker(-0.5*(f$SignVector[j,]-1), c(1,0)))  # both needed!
      }else{
        matAl[[j]] = rbind(kronecker(0.5*(f$SignVector[j,]+1), c(0, 0)), kronecker(-0.5*(f$SignVector[j,]-1), c(0,0)))
      }
    }
    Acons = bdiag(matAl)
  }else{
    Hsupportsub = Arrangesupport
    excl.index = which(names(Hsupport)%in%excl)
    matAl = list()
    for (j in 1:nrow(Hsupportsub)){
      matAindex  = apply(Hsupport, 1, function(x) (sum(x[-excl.index]==Hsupportsub[j,])==length(x[-excl.index])))
      for (k in 1:sum(matAindex)){
        tempinx = which(matAindex)[k]
        if (!is.na(OP[tempinx])){
          matAl[[tempinx]] = rbind(kronecker(0.5*(f$SignVector[j,]+1), c(0, 1)), kronecker(-0.5*(f$SignVector[j,]-1), c(1,0)))
        }else{
          matAl[[tempinx]] = rbind(kronecker(0.5*(f$SignVector[j,]+1), c(0, 0)), kronecker(-0.5*(f$SignVector[j,]-1), c(0,0)))
        }
      }
    }
    Acons = bdiag(matAl)
  }
  Acons[is.na(Acons)] = 0
  Acons
}

asplit <- function (x, MARGIN)
{
  dl <- length(dim(x))
  if (!dl)
    stop("dim(x) must have a positive length")
  if (is.object(x))
    x <- if (dl == 2L)
      as.matrix(x)
  else as.array(x)
  d <- dim(x)
  dn <- dimnames(x)
  ds <- seq_len(dl)
  if (is.character(MARGIN)) {
    if (is.null(dnn <- names(dn)))
      stop("'x' must have named dimnames")
    MARGIN <- match(MARGIN, dnn)
    if (anyNA(MARGIN))
      stop("not all elements of 'MARGIN' are names of dimensions")
  }
  s.call <- ds[-MARGIN]
  s.ans <- ds[MARGIN]
  d.call <- d[-MARGIN]
  d.ans <- d[MARGIN]
  dn.call <- dn[-MARGIN]
  dn.ans <- dn[MARGIN]
  d2 <- prod(d.ans)
  newx <- aperm(x, c(s.call, s.ans))
  dim(newx) <- c(prod(d.call), d2)
  ans <- vector("list", d2)
  for (i in seq_len(d2)) {
    ans[[i]] <- array(newx[, i], d.call, dn.call)
  }
  array(ans, d.ans, dn.ans)
}
