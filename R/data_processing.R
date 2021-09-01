#' Process data
#'
#'
#' Description here
#'
#' @param formula Path to the input file
#' @param data Path to the input file
#' @param support Path to the input file
#' @return An object
#' @export
process_data <- function(formula,data,support){
 # This function returns an object containing all relevant components necessary to compute ATE and CCP
 # data: data frame
 # formula: y ~ endogenous | exogenous | instrument
 # assumption: a list of assumptions - could refer to those in the paper i.e. A1, A2, etc...

  case_list <- parse.formula(formula)

    YXZW <- as.data.frame(data[all.vars(formula)])
    XZW <- YXZW[,-1]

    if(case_list$exovars==F){
      exovars <- c()
    }else if(case_list$exovars==T){
      formula <- as.Formula(formula)
      exolist <- all.vars(formula(formula,rhs=2)[[3]])
      exocol <- which(names(XZW)%in%exolist)
    }

    formula <- as.Formula(formula)

    exolist <- all.vars(formula(formula,rhs=2)[[3]])
    if(case_list$insvars==T){
    ivlist <- all.vars(formula(formula,rhs=3)[[3]])}
    else if(case_list$insvars==F){
      ivlist <- NULL
    }

    pvec <- pvec.process(data = list(YXZW,XZW), exolist = exolist, ivlist = ivlist, support = support)

    output <- list()
    output[["formula"]] <- formula
    output[["case_list"]] <- case_list
    output[["data"]] <- YXZW
    output[["support"]] <- pvec$p_support
    output[["exolist"]] <- exolist
    output[["ivlist"]] <- ivlist

  class(output) <- "bcDAT"

  return(output)

}

parse.formula <- function(formula){

  # Formula parsing
  # Four separate cases
  # (1) y ~ x1 + x2 : no exogenous variables, no instruments
  # (2) y ~ x1 + x2 | z1 + z2 : exogenous variables, no instruments
  # (3) y ~ x1 + x2 | -1 | w1 + w2 : no exogenous variables, instruments
  # (4) y ~ x1 + x2 | z1 + z2 | w1 + w2 : exogenous variables, instruments

  case_list <- list()

  formula <- as.Formula(formula)

  if(length(formula)[2]==1){
    exovars <- F
    insvars <- F
  }else if(length(formula)[2]==2){
    exovars <- T
    insvars <- F
  }else if(length(formula)[2]==3 & paste(as.character(unlist(formula(formula,rhs=2)[[3]])),collapse='') == "-1"){
    exovars <- F
    insvars <- T
  }else if(length(formula)[2]==3 & paste(as.character(unlist(formula(formula,rhs=2)[[3]])),collapse='') != "-1"){
    exovars <- T
    insvars <- T
  }

  case_list$exovars <- exovars
  case_list$insvars <- insvars

  return(case_list)

}

pvec.process = function(data, exolist, ivlist, support){

  fakedat = cbind(rep(c(0,1), nrow(support)),support[rep(seq_len(nrow(support)), each = 2),])
  names(fakedat) = c("y", names(support))
  XZ = data[[2]]
  YXZ = data[[1]]
  datsub = as.data.frame(YXZ)
  names(datsub) = names(fakedat)
  fakedatsub = rbind(datsub, fakedat)
  indepcol = which(names(support)%in%exolist)

  # Compute frequency table of exovars - to calculate what needs to be removed at final step

  exo_fake <- fakedat[,exolist] %>% as.data.frame()
  names(exo_fake) <- exolist
  exo_fake_sub <- fakedatsub[,exolist] %>% as.data.frame()
  names(exo_fake_sub) <- exolist
  fake_sub <- fakedatsub[,names(support)] %>% as.data.frame()
  names(fake_sub) <- names(support)

  exovars_sym <- syms(names(exo_fake))
  freq_remove <- as_tibble(exo_fake) %>%
    group_by(!!!exovars_sym) %>%
    summarize(total_remove = n() , .groups = 'drop')

  P_yxzw = rep(NA, nrow(support)*2) #P(Y = y, X = x, Z = z)  make sure even missing support point is included
  for (j in 1:nrow(support)){
    row_matches <- which(is.na(row.match(datsub[-1],as.data.frame(support[j,])))==F)
    temp = datsub[row_matches,]
    if (nrow(temp)>0){
      P_yxzw[(j-1)*2 + 1] = sum(temp$y==0)/nrow(datsub)
      P_yxzw[(j-1)*2 + 2] = sum(temp$y==1)/nrow(datsub)
    }else{
      P_yxzw[(j-1)*2 + 1] = NA
      P_yxzw[(j-1)*2 + 2] = NA
    }
  }
  # Assign P_yxzw based on order of the support
  support_P_yxzw <- support %>% slice(rep(1:n(), each = 2)) %>% mutate(y = ifelse(row_number() %% 2 ==0,1,0))
  support_P_yxzw$P_yxzw <- P_yxzw

  Cxz <- as_tibble(fake_sub) %>%
    group_by(!!!syms(names(support))) %>%
    summarize(total=n() , .groups = 'drop') %>%
    mutate(total = ifelse(total<=2,NA,total),
           total = total - 2,
           Pxz = total/nrow(datsub),
           y = 1) %>%
    left_join(support_P_yxzw %>% slice(which(row_number() %% 2 == 0)),by=names(fakedat)) %>%
    mutate(P_xw = P_yxzw/Pxz) %>%
    select(-c(Pxz,P_yxzw))


  if (length(exolist)==0){

    support_comb <- support_P_yxzw %>%
      left_join(Cxz,by=names(fakedat)) %>%
      mutate(Pyx_exo = NA)  %>%
      select(c(names(fakedat),P_yxzw,Pyx_exo,P_xw))


    return(list(p_support = support_comb, datsub = datsub))
  }else{

    Cz <- as_tibble(exo_fake_sub) %>%
      group_by(!!!exovars_sym) %>%
      summarize(total=n() , .groups = 'drop') %>%
      left_join(freq_remove,by=exolist,) %>%
      mutate(freq = total-total_remove,
             Pz = freq/nrow(datsub))  # Need to remove counts corresponding to the fake data

    support_comb <- support_P_yxzw %>%
      left_join(Cz,by=exolist) %>%
      mutate(Pyx_exo = P_yxzw/Pz)  %>%
      left_join(Cxz,by=names(fakedat)) %>%
      select(c(names(fakedat),P_yxzw,Pyx_exo,P_xw))

    names(support_comb) <- c(names(fakedat),"P_all","P_exo","P_y1")
    return(list(p_support = support_comb, datsub = datsub))
  }
}

