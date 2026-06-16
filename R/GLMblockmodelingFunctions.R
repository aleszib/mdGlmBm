unlistPar <- function(x) {
  unlist(x, recursive = TRUE, use.names = FALSE)
}


.normalizeOneModeNames<-function(M) {
  if(length(dim(M)) != 2 || nrow(M) != ncol(M)) {
    stop("M must be a square matrix for one-mode GLM blockmodeling.")
  }
  
  n<-ncol(M)
  rn<-rownames(M)
  cn<-colnames(M)
  
  usableNames<-function(x) {
    !is.null(x) && length(x)==n && all(!is.na(x)) && all(nzchar(x)) && anyDuplicated(x)==0
  }
  
  messages<-character(0)
  if(usableNames(cn)) {
    nodeNames<-cn
    if(!identical(rn, cn)) {
      messages<-c(messages, "row names did not match column names")
    }
  } else if(usableNames(rn)) {
    nodeNames<-rn
    messages<-c(messages, "column names were missing, empty, NA or duplicated")
  } else {
    baseNames<-cn
    if(is.null(baseNames) || length(baseNames)!=n) baseNames<-rn
    if(is.null(baseNames) || length(baseNames)!=n) baseNames<-as.character(seq_len(n))
    baseNames<-as.character(baseNames)
    badNames<-is.na(baseNames) | !nzchar(baseNames)
    baseNames[badNames]<-paste0("node", seq_len(n))[badNames]
    nodeNames<-make.unique(baseNames, sep="_")
    messages<-c(messages, "column names were missing, empty, NA or duplicated and row names were not usable")
  }
  
  if(is.null(rn) || is.null(cn) || !identical(rn, nodeNames) || !identical(cn, nodeNames)) {
    rownames(M)<-colnames(M)<-nodeNames
    warning(paste0(
      "Corrected row/column names of M: ",
      paste(messages, collapse="; "),
      ". Both row and column names were set to unique node names."
    ), call.=FALSE)
  }
  
  M
}



#' PPML family for `glm()`
#'
#' Creates a Poisson pseudo-maximum likelihood (PPML) family object for use
#' with [glm()]. The model uses the Poisson log-link mean specification, but
#' avoids the literal Poisson probability mass function in the AIC/deviance
#' components, making it suitable for non-negative non-integer outcomes such as
#' trade flows.
#'
#' PPML assumes the conditional mean
#'
#' \deqn{E(Y_i \mid X_i) = \mu_i = \exp(X_i \beta)}
#'
#' but does not require the outcome to be integer-valued or truly Poisson
#' distributed.
#'
#' @param link Character string specifying the link function. Defaults to
#'   `"log"`, which is the standard PPML specification.
#'
#' @return A family object suitable for use in [glm()].
#'
#' @details
#' This family uses the usual Poisson quasi-deviance contribution
#'
#' \deqn{
#' 2 \left[
#' y_i \log(y_i / \mu_i) - (y_i - \mu_i)
#' \right],
#' }
#'
#' with the convention that the contribution is \eqn{2\mu_i} when
#' \eqn{y_i = 0}.
#'
#' The `aic` component returns a pseudo-likelihood criterion based on
#'
#' \deqn{
#' \sum_i \left[
#' y_i \log(\mu_i) - \mu_i
#' \right],
#' }
#'
#' omitting the problematic \eqn{\log(y_i!)} term, which is not appropriate for
#' non-integer outcomes.
#'
#' @examples
#' fit <- glm(
#'   y ~ x,
#'   family = ppml(),
#'   data = data.frame(y = c(0, 1, 2, 3), x = c(0, 1, 0, 1))
#' )
#'
#' @export
ppml <- function(link = "log") {
  fam <- poisson(link = link)
  
  fam$family <- "ppml"
  attr(fam, "pseudo") <- TRUE
  
  fam$dev.resids <- function(y, mu, wt) {
    2 * wt * ifelse(
      y == 0,
      mu,
      y * log(y / mu) - (y - mu)
    )
  }
  
  fam$aic <- function(y, n, mu, wt, dev) {
    -2 * sum(wt * (ifelse(y == 0, 0, y * log(mu)) - mu))
  }
  
  fam$initialize <- expression({
    if (any(y < 0)) {
      stop("negative values not allowed for 'PPML' family")
    }
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  
  fam
}


#' Fit a GLM-based blockmodel partition
#'
#' @param M A square adjacency matrix.
#' @param clu Initial hard partition.
#' @param formula GLM formula, defaulting to `Freq ~ block`.
#' @param twoStep If `TRUE`, use the two-step offset procedure from the source
#'   package.
#' @param vars Optional list of additional dyadic covariate matrices.
#' @param inDeg,outDeg Include in-degree and out-degree covariates.
#' @param logDeg,degAsOfset Control degree transformations.
#' @param maxIter Maximum local reassignment iterations.
#' @param ignoreDiag Ignore the diagonal when building the dyad table.
#' @param glmFamily Family passed to `glm()`. Use `ppml()` for PPML mode.
#' @param pChange,pChangeRandom Optional change-limit controls.
#' @param verbose Print optimization progress.
#' @param evalMethod Candidate scoring mode.
#' @param maxFastRows Threshold for auto-selecting the fast path.
#' @param cl Optional parallel cluster.
#' @param senderReceiver Add sender/receiver factors instead of block terms.
#'
#' @return A list with `clu`, `fit`, `logLik`, `BIC`, `ICL`, `deviance`,
#'   `objective`, `pseudo`, and supporting diagnostics.
#' @export
optParGlm<-function(M, clu, formula=as.formula(Freq ~ block), twoStep=FALSE, vars=NULL, inDeg=FALSE, outDeg=FALSE, logDeg=TRUE, degAsOfset=TRUE, maxIter=100, ignoreDiag=TRUE, glmFamily=binomial, pChange=0, pChangeRandom=TRUE, verbose=FALSE, evalMethod=c("auto", "fast", "lowMemory", "parallel"), maxFastRows=5e5, cl=NULL, senderReceiver=FALSE) {
  # M: adjacency matrix
  # clu: initial clustering
  # formula: formula for the GLM. The response variable should be "Freq" and the block variable should be "block". The block variable will be created as a combination of the in-cluster and out-cluster of the nodes. For example, if node i is in cluster 1 and node j is in cluster 2, then the block variable for the tie from i to j will be "1#2". The block variable can be omitted and it is added to the model by the function. The formula can include other variables as well. Indegree (inDeg) and outdegree (outDeg) can be included as variables in the model by setting the corresponding arguments to TRUE. The names of the variables will be "inDeg" and "outDeg" respectively.
  # twoStep: Should two-step procedured be used. In the two-step procedure, the model is first fitted without the block variable. The predicted value is used as an offset in the model where the only estimated parameters are associated with the block variable. This means that the block can only explain what is not exaplained by the other variables in the model. The default is FALSE, which means that the model is fitted with all variables at once.
  # vars: additional variables to include in the model as a list of matrices with the same dimensions as M. The names of the list will be used as variable names in the model.
  # inDeg: whether to include indegree as a variable in the model. The name of the variable will be "inDeg".
  # outDeg: whether to include outdegree as a variable in the model. The name of the variable will be "outDeg".
  # senderReceiver: whether to include sender and receiver node effects as factors in the model. The names of the variables will be "sender" and "receiver". This cannot be combined with inDeg or outDeg.
  # logDeg: should the degree variables be log-transformed? The default is TRUE, which means that the degree variables will be log-transformed. If FALSE, the degree variables will be included in the model as they are.
  # degAsOfset: should the degree variables be included as offsets in the model? The default is TRUE, which means that the degree variables will be included as offsets. If FALSE, the degree variables will be included as regular variables in the model.
  # maxIter: number of iterations to perform
  # ignoreDiag: whether to ignore the diagonal of the adjacency matrix (i.e., self ties) when calculating the deviance. If TRUE, the diagonal will be set to NA before calculating the deviance.
  # glmFamily: family to use for the GLM. The default is "binomial", but other families can be used as well (e.g., "poisson", "gaussian", etc.).
  # pChange: the maximal proportion of units that can change cluster membership in each iteration. If the proportion of units that change cluster membership exceeds this threshold. See pChangeRandom for how this selection is made. The default is 0, which means that there is no limit on the proportion of units that can change cluster membership.
  # pChangeRandom = TRUE: Should the units that are allowed to change cluster membership be selected randomly proportionally to the gain in deviance contribution, or should the units with the largest gain be selected? The default is TRUE, which means that the units will be selected randomly. If FALSE, the units with the largest gain will be selected.he 
  # evalMethod: method used to evaluate cluster reassignment scores.
  #   "auto" (default): automatically chooses between "fast" and "lowMemory" based on the estimated size of the temporary data structure.
  #   "fast": uses a vectorized approach by constructing a large temporary data.frame containing all candidate cluster assignments and evaluating them at once. This is typically faster for smaller networks or a small number of clusters but may require substantial memory.
  #   "lowMemory": evaluates cluster reassignment scores unit-by-unit without constructing a large temporary data.frame. This approach is more memory-efficient and suitable for larger networks but may be slower due to reduced vectorization.
  #   "parallel": same as "lowMemory" but computations are distributed across multiple cores using parallel processing. This is useful for larger problems when computational time becomes a bottleneck.
  # maxFastRows: threshold for deciding between "fast" and "lowMemory" when evalMethod = "auto".
  #   It represents the estimated number of rows in the temporary data.frame that would be created by the "fast" method.
  #   The estimate is approximately 2 * nrow(df) * k, where k is the number of clusters.
  #   If the estimated number of rows is less than or equal to maxFastRows, the "fast" method is used; otherwise, "lowMemory" is used.
  #   The default value (e.g., 5e5) should work well in most cases, but it can be adjusted depending on available memory and performance requirements.  
  if(verbose) cat("Start\n")
  if(senderReceiver && (inDeg || outDeg)) {
    stop("senderReceiver=TRUE cannot be used together with inDeg=TRUE or outDeg=TRUE.")
  }

  family_info <- glm_blockmodel_family(glmFamily)
  glmFamily <- family_info$family
  
  formula <- update(formula, .~.+block)
  formula <- update(formula, Freq~.)
  
  devContNewdata <- function(fit, newdata, yvar) {
    dev<-rep(NA, nrow(newdata))
    okBlock <- (newdata$block %in% unique(fit$model$block))
    y <- newdata[okBlock,yvar]
    mu <- predict(fit, newdata[okBlock,], type = "response")
    dev[okBlock]<-family(fit)$dev.resids(y, mu, wt = 1)
    if(any(!okBlock)) {
      warning("Some blocks in newdata were not in the original data. Re-fitting the model without block effects for these ties.")
      modFit<-update(fit, .~.-block)
      y <- newdata[!okBlock,yvar]
      mu <- predict(modFit, newdata[!okBlock,], type = "response")
      dev[!okBlock]<-family(modFit)$dev.resids(y, mu, wt = 1)
    }
    return(dev)
  }
  
  
  evalMethod <- match.arg(evalMethod)
  k<-length(unique(clu))
  M<-.normalizeOneModeNames(M)
  if(ignoreDiag) {
    diag(M)<-NA
  } 
  names(clu)<-colnames(M)
  cluTmp<-clu
  bestClu<-cluTmp
  n<-ncol(M)
  df<-as.data.frame(as.table(M))
  if(!is.null(vars)) {
    for(varName in names(vars)) {
      varMat<-vars[[varName]]
      df[[varName]]<-as.data.frame(as.table(varMat))$Freq
    }
  }
  if(senderReceiver) {
    df$sender<-factor(df$Var1, levels=colnames(M))
    df$receiver<-factor(df$Var2, levels=colnames(M))
    formula<-update(formula, .~.+sender+receiver)
  }
  if(inDeg) {
    inDeg<-colSums(M, na.rm=TRUE)
    names(inDeg)<-colnames(M)
    df$inDeg<-inDeg[df$Var2]
    if(logDeg) df$inDeg<-log(df$inDeg)
    if(degAsOfset) {
      formula<-update(formula, .~.+offset(inDeg))
    } else formula<-update(formula, .~.+inDeg)
  }
  if(outDeg) {
    outDeg<-rowSums(M, na.rm=TRUE)
    names(outDeg)<-colnames(M)
    df$outDeg<-outDeg[df$Var1]
    if(logDeg) df$outDeg<-log(df$outDeg)
    if(degAsOfset) {
      formula<-update(formula, .~.+offset(outDeg))
    } else formula<-update(formula, .~.+outDeg)
  }
  
  df<-df[!is.na(df$Freq), , drop=FALSE]
  
  getScoreByClustersFast <- function(df, fit, cluLevels, pi_k, M) {
    tmpDf<-NULL
    for(i in seq_len(n)) {
      iName<-colnames(M)[i]
      iDf<-df[df$Var1==iName | df$Var2==iName, , drop=FALSE]
      for(ik in cluLevels) {
        iDf$inClu[iDf$Var2==iName]<-ik
        iDf$outClu[iDf$Var1==iName]<-ik
        iDf$block<-paste(iDf$outClu, iDf$inClu, sep="#")
        iDf$unit<-iName
        iDf$clu<-ik
        tmpDf<-rbind(tmpDf, iDf)
      }
    }
    tmpDf<-tmpDf[!is.na(tmpDf$Freq), , drop=FALSE]
    tmpDf$dev<-devContNewdata(fit, tmpDf, "Freq")
    devByClusters<-xtabs(dev ~ unit + clu, data=tmpDf)
    devByClusters<-devByClusters[rownames(M), as.character(cluLevels), drop=FALSE]
    scoreByClusters <- sweep(devByClusters, 2, -2 * log(pi_k), "+")
    return(scoreByClusters)
  }

  oneUnitScores <- function(i, df, fit, cluLevels, pi_k, M) {
    iName <- colnames(M)[i]
    iDf0 <- df[df$Var1 == iName | df$Var2 == iName, , drop = FALSE]
    iDf0 <- iDf0[!is.na(iDf0$Freq), , drop = FALSE]
    isIn <- iDf0$Var2 == iName
    isOut <- iDf0$Var1 == iName
    scores <- numeric(length(cluLevels))

    for (j in seq_along(cluLevels)) {
      ik <- cluLevels[j]
      iDf <- iDf0
      iDf$inClu[isIn] <- ik
      iDf$outClu[isOut] <- ik
      iDf$block <- paste(iDf$outClu, iDf$inClu, sep = "#")
      scores[j] <- sum(devContNewdata(fit, iDf, "Freq")) - 2 * log(pi_k[j])
    }

    return(scores)
  }

  getScoreByClustersLowMemory <- function(df, fit, cluLevels, pi_k, M) {
    scoreList <- lapply(seq_len(n), oneUnitScores, df = df, fit = fit, cluLevels = cluLevels, pi_k = pi_k, M = M)
    scoreByClusters <- do.call(rbind, scoreList)
    rownames(scoreByClusters) <- rownames(M)
    colnames(scoreByClusters) <- as.character(cluLevels)
    return(scoreByClusters)
  }

  getScoreByClustersParallel <- function(df, fit, cluLevels, pi_k, M, cl) {
    parallel::clusterExport(
      cl,
      varlist = c("df", "fit", "cluLevels", "pi_k", "M", "devContNewdata", "oneUnitScores", "n"),
      envir = environment()
    )

    scoreList <- parallel::parLapply(
      cl,
      seq_len(n),
      oneUnitScores,
      df = df, fit = fit, cluLevels = cluLevels, pi_k = pi_k, M = M
    )

    scoreByClusters <- do.call(rbind, scoreList)
    rownames(scoreByClusters) <- rownames(M)
    colnames(scoreByClusters) <- as.character(cluLevels)
    return(scoreByClusters)
  }
  
  resolveEvalMethod <- function(evalMethod, estTmpRows, maxFastRows) {
    if(evalMethod != "auto") return(evalMethod)
    if(is.na(estTmpRows)) estTmpRows <- Inf
    if(estTmpRows <= maxFastRows) "fast" else "lowMemory"
  }

  estTmpRows <- 2 * nrow(df) * k
  evalMethodUsed <- resolveEvalMethod(evalMethod, estTmpRows = estTmpRows, maxFastRows = maxFastRows)
  if(evalMethod == "auto") cat("Evaluation method selected:", evalMethodUsed, "\n")

  if(evalMethodUsed == "parallel" && is.null(cl)) {
    if(verbose) cat("No cluster provided. Creating a local cluster for parallel processing.\n")
    nCores <- max(1L, parallel::detectCores() - 1L)
    cl <- parallel::makeCluster(nCores)
    createdLocalCluster <- TRUE
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  
  df$inClu<-cluTmp[df$Var2]
  df$outClu<-cluTmp[df$Var1]
  df$block<-paste(df$outClu,df$inClu,sep="#")
  
  if(twoStep) {
    orgFormula<-formula
    fitStepOne<-glm(update(formula, .~.-block), data=df, family=glmFamily)
    df$offsetPred<-predict(fitStepOne, type="link", newdata=df)
    formula<-as.formula(Freq ~ offset(offsetPred) + block + 0)
    n_params_step1 <-length(coef(fitStepOne))
  }
  
  if(length(unique(df$block))==1){
    if(twoStep) {
      fit<-fitStepOne
      n_params_step1<-1
    } else {
      tmpFormula<-update(formula, .~.-block)
      fit<-glm(tmpFormula, data=df, family=glmFamily)
    }
  } else {
    fit<-glm(formula, data=df, family=glmFamily)
  }

  nk <- table(clu)
  pi_k <- nk /n

  ll<-rep(NA, maxIter+1)
  ll[1]<-as.numeric(logLik(fit)) + sum(nk * log(nk / n))
  if(verbose) cat("Initial log-likelihood:", ll[1], "\n")
  if(k>1 & maxIter>0) for(iRep in 1:maxIter) {
    cluLevels <- sort(unique(clu))
    nk <- table(factor(clu, levels = cluLevels))
    pi_k <- as.numeric(nk / n)
    names(pi_k) <- as.character(cluLevels)


    scoreByClusters <- switch(
      evalMethodUsed,
      fast = getScoreByClustersFast(df = df, fit = fit, cluLevels = cluLevels, pi_k = pi_k, M = M),
      lowMemory = getScoreByClustersLowMemory(df = df, fit = fit, cluLevels = cluLevels, pi_k = pi_k, M = M),
      parallel = getScoreByClustersParallel(df = df, fit = fit, cluLevels = cluLevels, pi_k = pi_k, M = M, cl = cl),
      stop("Unknown evalMethod: ", evalMethodUsed)
    )

    cluTmp <- apply(scoreByClusters, 1, function(x) {
      clMin <- which(x == min(x))
      if(length(clMin) == 1) return(clMin)
      sample(clMin, 1)
    })
    while(length(unique(cluTmp)) < k) {
      maxClu<-as.integer(names(which.max(table(cluTmp))))
      emptyClu <- setdiff(1:k, unique(cluTmp))[1]
      whichMax <- cluTmp==maxClu
      maxDev<-which.max(apply(scoreByClusters[whichMax,],1,min))
      cluTmp[whichMax][maxDev] <- emptyClu
    }
    cluTmp <- cluLevels[cluTmp]
    #cluTmp<-cluTmp[rownames(M)]
    nChange<-sum(cluTmp!=clu)
    # print timestamp
    if(verbose) cat("Timestamp: ", format(Sys.time()), "\n")
    if(pChange>0 && nChange/n > pChange) {
      # warning(paste("More than", pChange*100, "% of units are changing cluster membership. Limiting the number of changes to", pChange*100, "% of units with the largest gain."))
      oldCluIdx <- match(as.character(clu), colnames(scoreByClusters))
      newCluIdx <- match(as.character(cluTmp), colnames(scoreByClusters))
      gain<- scoreByClusters[cbind(1:n, oldCluIdx)] - scoreByClusters[cbind(1:n, newCluIdx)]
      if(pChangeRandom) {
        change<-sample(n, size = round(pChange*n), prob = gain, replace = FALSE)
      } else change<- which(gain >= sort(gain, decreasing = TRUE)[round(pChange*n)])
      cluTmp[-change]<-clu[-change]
      if(verbose) cat("Iteration ", iRep, ": ", nChange, " units would change cluster membership. Limited to ",length(change)," as pChange = ",pChange,".\n", sep="")
    } else{
      if(verbose) cat("Iteration ", iRep, ": ", nChange, " units changing cluster membership.\n", sep="")
    }
    
    if(all(cluTmp==clu)) {
      break
    } else {
      clu<-cluTmp
      nk <- table(clu)
      pi_k <- nk / n
      df$inClu<-clu[df$Var2]
      df$outClu<-clu[df$Var1]
      df$block<-paste(df$outClu,df$inClu,sep="#")
      fit<-glm(formula, data=df, family=glmFamily)
      ll[iRep+1]<-as.numeric(logLik(fit)) + sum(nk * log(pi_k))
      if(verbose) cat("New log-likelihood:", ll[iRep+1],"\n")
      if(ll[iRep+1]<ll[iRep]) {
        warning("Log-likelihood decreased. Stopping.")
        break
      } else {
        bestClu<-clu
      }
    }
  }
  ll<-na.omit(ll)
  llFinal<-ll[length(ll)]
  n_params <- length(coef(fit))
  if(twoStep) {
    n_params <- n_params + n_params_step1 - 1 # subtract 1 to account for modeling the mean level twice
  }
  n_eff <- length(fit$y)
  ICL <- llFinal - 1/2 * (n_params * log(n_eff) + (k - 1) * log(n))
  if(verbose) cat("Final log-likelihood:", llFinal, "ICL:", ICL, "\n")
  out <- list(clu = bestClu, fit = fit, logLik = llFinal, ICL = ICL, logLikVec = ll, err = fit$deviance)
  meta <- bm_result_metadata(
    fit = fit,
    clu = bestClu,
    ICL = ICL,
    pseudo = family_info$pseudo,
    objective = ICL,
    family = family_info
  )
  out$membership <- meta$membership
  out$BIC <- meta$BIC
  out$deviance <- meta$deviance
  out$objective <- meta$objective
  out$pseudo <- meta$pseudo
  out$family <- meta$family
  out$criterion_note <- meta$criterion_note
  out
}



## This function will likely be removed. It is no longer developed. It does not support twostep procedure. The samefunctionality can be obtained wiht optParGlm with argument maxIter = 0.
glmFit<-function(M, clu, formula=as.formula(Freq ~ block), vars=NULL, inDeg=FALSE, outDeg=FALSE, logDeg=TRUE, degAsOfset=TRUE, ignoreDiag=TRUE, glmFamily=binomial, senderReceiver=FALSE) {
  # M: adjacency matrix
  # clu: initial clustering
  # formula: formula for the GLM. The response variable should be "Freq" and the block variable should be "block". The block variable will be created as a combination of the in-cluster and out-cluster of the nodes. For example, if node i is in cluster 1 and node j is in cluster 2, then the block variable for the tie from i to j will be "1#2". The formula can include other variables as well. Indegree (inDeg) and outdegree (outDeg) can be included as variables in the model by setting the corresponding arguments to TRUE. The names of the variables will be "inDeg" and "outDeg" respectively.
  # vars: additional variables to include in the model as a list of matrices with the same dimensions as M. The names of the list will be used as variable names in the model.
  # inDeg: whether to include indegree as a variable in the model. The name of the variable will be "inDeg".
  # outDeg: whether to include outdegree as a variable in the model. The name of the variable will be "outDeg".
  # senderReceiver: whether to include sender and receiver node effects as factors in the model. The names of the variables will be "sender" and "receiver". This cannot be combined with inDeg or outDeg.
  # logDeg: should the degree variables be log-transformed? The default is TRUE, which means that the degree variables will be log-transformed. If FALSE, the degree variables will be included in the model as they are.
  # degAsOfset: should the degree variables be included as offsets in the model? The default is TRUE, which means that the degree variables will be included as offsets. If FALSE, the degree variables will be included as regular variables in the model.  
  # ignoreDiag: whether to ignore the diagonal of the adjacency matrix (i.e., self ties) when calculating the deviance. If TRUE, the diagonal will be set to NA before calculating the deviance.
  # glmFamily: family to use for the GLM. The default is "binomial", but other families can be used as well (e.g., "poisson", "gaussian", etc.).

  if(senderReceiver && (inDeg || outDeg)) {
    stop("senderReceiver=TRUE cannot be used together with inDeg=TRUE or outDeg=TRUE.")
  }

  family_info <- glm_blockmodel_family(glmFamily)
  glmFamily <- family_info$family
  
  k<-length(unique(clu))
  M<-.normalizeOneModeNames(M)
  if(ignoreDiag) {
    diag(M)<-NA
  } 
  names(clu)<-colnames(M)
  cluTmp<-clu
  n<-ncol(M)
  df<-as.data.frame(as.table(M))
  if(!is.null(vars)) {
    for(varName in names(vars)) {
      varMat<-vars[[varName]]
      df[[varName]]<-as.data.frame(as.table(varMat))$Freq
    }
  }
  if(senderReceiver) {
    df$sender<-factor(df$Var1, levels=colnames(M))
    df$receiver<-factor(df$Var2, levels=colnames(M))
    formula<-update(formula, .~.+sender+receiver)
  }
  if(inDeg) {
    inDeg<-colSums(M, na.rm=TRUE)
    names(inDeg)<-colnames(M)
    df$inDeg<-inDeg[df$Var2]
    if(logDeg) df$inDeg<-log(df$inDeg)
    if(degAsOfset) {
      formula<-update(formula, .~.+offset(inDeg))
    } else formula<-update(formula, .~.+inDeg)
  }
  if(outDeg) {
    outDeg<-rowSums(M, na.rm=TRUE)
    names(outDeg)<-colnames(M)
    df$outDeg<-outDeg[df$Var1]
    if(logDeg) df$outDeg<-log(df$outDeg)
    if(degAsOfset) {
      formula<-update(formula, .~.+offset(outDeg))
    } else formula<-update(formula, .~.+outDeg)
  }

  df$inClu<-cluTmp[df$Var2]
  df$outClu<-cluTmp[df$Var1]
  df$block<-paste(df$outClu,df$inClu,sep="#")
  if(length(unique(df$block))==1){
    formula<-update(formula, .~.-block)
  }
  fit<-glm(formula, data=df, family=glmFamily)
  nk <- table(clu)
  pi_k <- nk / n
  ll<-as.numeric(logLik(fit)) + sum(nk * log(pi_k))
  n_params <- length(coef(fit))
  n_eff <- length(fit$y)
  ICL <- ll - 1/2 * (n_params * log(n_eff) + (k - 1) * log(n))
  
  out <- list(clu = clu, fit = fit, logLik = ll, ICL = ICL, err = fit$deviance)
  meta <- bm_result_metadata(
    fit = fit,
    clu = clu,
    ICL = ICL,
    pseudo = family_info$pseudo,
    objective = ICL,
    family = family_info
  )
  out$membership <- meta$membership
  out$BIC <- meta$BIC
  out$deviance <- meta$deviance
  out$objective <- meta$objective
  out$pseudo <- meta$pseudo
  out$family <- meta$family
  out$criterion_note <- meta$criterion_note
  out
}



#' A function for optimizing multiple random partitions using GLM one-mode blockmodeling. Calls optParGlm for optimizing individual partitions.
#'
#'
#' @param M A square matrix giving the adjaciency relationg between the network's nodes (aka vertexes)
#' @param k The number of clusters used in the generation of partitions.
#' @param rep The number of repetitions/different starting partitions to check.
#' @param save.initial.param Should the inital parameters(\code{approaches}, ...) of using \code{stochBlock} be saved. The default value is \code{TRUE}.
#' @param deleteMs Delete networks/matrices from the results of to save space. Defaults to \code{TRUE}.
#' @param max.iden Maximum number of results that should be saved (in case there are more than \code{max.iden} results with minimal error, only the first \code{max.iden} will be saved).
#' @param return.all If \code{FALSE}, solution for only the best (one or more) partition/s is/are returned.
#' @param return.err Should the error for each optimized partition be returned. Defaults to \code{TRUE}.
#' @param seed Optional. The seed for random generation of partitions.
#' @param parGenFun The function (object) that will generate random partitions. The default function mirrors the imported `genRandomPar()` helper from the `blockmodeling` package. The function has to accept the following parameters: \code{k} (number o of partitions by modes, \code{n} (number of units by modes), \code{seed} (seed value for random generation of partition), \code{addParam} (a list of additional parameters).
#' @param mingr Minimal allowed group size.
#' @param maxgr Maximal allowed group size.
#' @param addParam A list of additional parameters for function specified above. In the usage section they are specified for the default `genRandomPar()` helper.
#' @param maxTriesToFindNewPar The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before - the default value is \code{rep * 1000}.
#' @param skip.par The partitions that are not allowed or were already checked and should therefore be skipped.
#' @param printRep Should some information about each optimization be printed.
#' @param n The number of units by "modes". It is used only for generating random partitions. It has to be set only if there are more than two modes or if there are two modes, but the matrix representing the network is one mode (both modes are in rows and columns).
#' @param nCores Number of cores to be used. Value \code{0} means all available cores. It can also be a cluster object.
#' @param useParLapply Should `parallel::parLapplyLB()` be used. When `FALSE`,
#'   a serial fallback is used. Defaults to `FALSE`.
#' @param cl The cluster to use (if formed beforehand). Defaults to \code{NULL}.
#' @param stopcl Should the cluster be stopped after the function finishes. Defaults to \code{is.null(cl)}.
#' @param \dots Arguments passed to other functions, see `stochBlock`.
#'
#' @return A list of class "opt.more.par" containing:
#'  \item{M}{The one- or multi-mode matrix of the network analyzed}
#'   \item{res}{If \code{return.all = TRUE} - A list of results the same as \code{best} - one \code{best} for each partition optimized.}
#'   \item{best}{A list of results from \code{stochblock}, only without \code{M}.}
#'   \item{err}{If \code{return.err = TRUE} - The vector of errors or inconsistencies = -log-likelihoods.}
#'   \item{ICL}{Integrated classification likelihood for the best partition.}
#'   \item{checked.par}{If selected - A list of checked partitions. If \code{merge.save.skip.par} is \code{TRUE}, this list also includes the partitions in \code{skip.par}.}
#'   \item{call}{The call to this function.}
#'   \item{initial.param}{If selected - The initial parameters are used.}
#'   \item{Random.seed}{.Random.seed at the end of the function.}
#'   \item{cl}{Cluster used for parallel computations if supplied as an input parameter.}
#'   
#' @section Warning:
#' It should be noted that the time needed to optimise the partition depends on the number of units (aka nodes) in the networks as well as the number of clusters
#' due to the underlying algorithm. Hence, partitioning networks with 100 units and large number of blocks (e.g., >5) can take a long time (from 20 minutes to a few hours or even days).
#' 
#' 
#' @author \enc{Aleš, Žiberna}{Ales Ziberna}
#' 
#' @seealso `stochBlock`
#' 
#' @examples
#' \dontrun{
#'# Simple one-mode network
#'library(blockmodeling)
#'k<-2
#'blockSizes<-rep(20,k)
#'IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
#'if(any(dim(IM)!=c(k,k))) stop("invalid dimensions")
#'
#'set.seed(2021)
#'clu<-rep(1:k, times=blockSizes)
#'n<-length(clu)
#'M<-matrix(rbinom(n*n,1,IM[clu,clu]),ncol=n, nrow=n)
#'diag(M)<-0
#'plotMat(M)
#'
#'resORP<-optRandomParGlm(M,k=2, rep=10, return.all = TRUE)
#'resORP$ICL
#'plot(resORP)
#'clu(resORP)
#'
#'
#'# Linked network
#'library(blockmodeling)
#'set.seed(2021)
#'IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
#'clu<-rep(1:2, each=20)
#'n<-length(clu)
#'nClu<-length(unique(clu))
#'M1<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n)
#'M2<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n)
#'M12<-diag(n)
#'nn<-c(n,n)
#'k<-c(2,2)
#'Ml<-matrix(0, nrow=sum(nn),ncol=sum(nn))
#'Ml[1:n,1:n]<-M1
#'Ml[n+1:n,n+1:n]<-M2
#'Ml[n+1:n, 1:n]<-M12
#'plotMat(Ml)
#'
#'resMl<-stochBlockORP(M=Ml, k=k, n=nn, rep=10)
#'resMl$ICL
#'plot(resMl)
#'clu(resMl)
#' }
#'
#' @export
optRandomParGlm<-function(M, #a square matrix
                        k,#number of clusters/groups
                        rep,#number of repetitions/different starting partitions to check
                        save.initial.param=TRUE,  #save the initial parameters of this call
                        deleteMs=TRUE, #delete networks/matrices from results of optParC or optParMultiC to save space
                        max.iden=10, #the maximum number of results that should be saved (in case there are more than max.iden results with minimal error, only the first max.iden will be saved)
                        return.all=FALSE,#if 'FALSE', solution for only the best (one or more) partition/s is/are returned
                        return.err=TRUE,#if 'FALSE', only the results of crit.fun are returned (a list of all (best) solutions including errors), else the result is list
                        seed=NULL,#the seed for random generation of partitions
                        parGenFun = bm_random_partition, #The function that will generate random partitions. It should accept arguments: k (number of partitions by modes, n (number of units by modes), seed (seed value for random generation of partition), addParam (a list of additional parameters)
                        mingr=NULL, #minimal allowed group size (defaults to c(minUnitsRowCluster,minUnitsColCluster) if set, else to 1) - only used for parGenFun function 
                        maxgr=NULL, #maximal allowed group size (default to c(maxUnitsRowCluster,maxUnitsColCluster) if set, else to Inf) - only used for parGenFun function 
                        addParam=list(  #list of additional parameters for generating partitions. Here they are specified for the default function "genRandomPar"
                          genPajekPar = TRUE,     #Should the partitions be generated as in Pajek (the other options is completely random)
                          probGenMech = NULL),    #Here the probabilities for different mechanisms for specifying the partitions are set. If not set this is determined based on the previous parameter.
                        maxTriesToFindNewPar=rep*10,    #The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before 
                        skip.par = NULL, #partitions to be skipped
                        printRep= ifelse(rep<=10,1,round(rep/10)), #should some information about each optimization be printed
                        n=NULL, #the number of units by "modes". It is used only for generating random partitions. It has to be set only if there are more than two modes or if there are two modes, but the matrix representing the network is onemode (both modes are in rows and columns)
                        nCores=1, #number of cores to be used 0 -means all available cores, can also be a cluster object,
                        useParLapply=FALSE, #should parLapply be used instead of foreach
                        cl = NULL, #the cluster to use (if formed beforehand)
                        stopcl = is.null(cl), # should the cluster be stopped
                        ... #paramters to stochBlock
){
  dots<-list(...)
  
  if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters
  
  if(is.null(mingr)){
    if(is.null(dots$minUnitsRowCluster)){
      mingr<-1
    } else {
      mingr<-c(dots$minUnitsRowCluster,dots$minUnitsColCluster)
    }
  }
  
  if(is.null(maxgr)){
    if(is.null(dots$maxUnitsRowCluster)){
      maxgr<-Inf
    } else {
      maxgr<-c(dots$maxUnitsRowCluster,dots$maxUnitsColCluster)
    }
  }
  
  nmode<-length(k)
  
  res<-list(NULL)
  err<-NULL
  dots<-list(...)
  
  if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters
  
  
  if(is.null(n)) if(nmode==1){
    n<-dim(M)[1]
  } else stop("Currently only one-mode networks are implemented!")
  # if(nmode==2){
  #   n<-dim(M)[1:2]
  # } else warning("Number of nodes by modes can not be determined. Parameter 'n' must be supplied!!!")
    
  if(!is.null(seed))set.seed(seed)
  
  
  on.exit({
    res1 <- res[which(err==min(err, na.rm = TRUE))]
    best<-NULL
    best.clu<-NULL
    for(i in 1:length(res1)){
      for(j in 1:length(res1[[i]]$best)){
        if(
          ifelse(is.null(best.clu),
                 TRUE,
                 if(nmode==1){
                   !any(sapply(best.clu,bm_partition_equal,clu2=res1[[i]]$clu)==1)
                 } else {
                   !any(sapply(best.clu,function(x,clu2)bm_partition_equal(unlist(x),clu2),clu2=unlist(res1[[i]]$clu))==1)
                 }
          )
        ){
          best<-c(best,res1[i])
          best.clu<-c(best.clu,list(res1[[i]]$clu))
        }
        
        if(length(best)>=max.iden) {
          warning("Only the first ",max.iden," solutions out of ",length(na.omit(err))," solutions with minimal -loglikelihood will be saved.\n")
          break
        }
        
      }
    }
    
    names(best)<-paste("best",1:length(best),sep="")
    
    if(any(na.omit(err)==-Inf) || ss(na.omit(err))!=0 || length(na.omit(err))==1){
      cat("\n\nOptimization of all partitions completed\n")
      cat(length(best),"solution(s) with minimal deviance =", min(err,na.rm=TRUE), "found.","\n")
    }else {
      cat("\n\nOptimization of all partitions completed\n")
      cat("All",length(na.omit(err)),"solutions have deviance",err[1],"\n")
    }
    
    call<-list(call=match.call())
    best<-list(best=best)
    checked.par<-list(checked.par=skip.par)
    if(return.all) res<-list(res=res) else res<-NULL
    if(return.err) err<-list(err=err) else err<-NULL
    if(!exists("initial.param")){
      initial.param<-NULL
    } else initial.param=list(initial.param)
    
    res<-c(list(M=M),list(ICL=best[[1]][[1]]$ICL),res,best,err,checked.par,call,initial.param=initial.param, list(Random.seed=.Random.seed, cl=cl))
    class(res)<-"opt.more.par"
    return(res)
  })
  
  
  
  if(nCores==1||!requireNamespace('parallel')){
    if(nCores!=1) {
      warning("Only single core is used as package 'parallel' is not available", immediate.=TRUE)
    }
    for(i in 1:rep){
      if(printRep & (i%%printRep==0)) cat("\n\nStarting optimization of the partiton",i,"of",rep,"partitions.\n")
      find.unique.par<-TRUE
      ununiqueParTested=0
      while(find.unique.par){
        temppar<-parGenFun(n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam)
        
        find.unique.par<-
          ifelse(is.null(skip.par),
                 FALSE,
                 if(nmode==1) {
                   any(sapply(skip.par,bm_partition_equal,clu2=temppar)==1)
                 } else any(sapply(skip.par,function(x,clu2)bm_partition_equal(unlist(x),clu2),clu2=unlist(temppar))==1)
          )
        ununiqueParTested=ununiqueParTested+1
        endFun<-ununiqueParTested>=maxTriesToFindNewPar
        if(endFun) {
          break
        } else if(ununiqueParTested%%10==0) cat(ununiqueParTested,"partitions tested for unique partition\n")
      }
      
      if(endFun) break
      
      skip.par<-c(skip.par,list(temppar))
      
      if(printRep==1) cat("Starting partition:",unlistPar(temppar),"\n")
      res[[i]]<-optParGlm(M=M, clu=temppar,  ...)
      if(deleteMs){
        res[[i]]$M<-NULL
      }
      res[[i]]$best<-NULL
      
      err[i]<-res[[i]]$err
      if(printRep==1) cat("Final deviance:",err[i],"\n")
      if(printRep==1) cat("Final partition:   ",unlistPar(res[[i]]$clu),"\n")
    }
  } else {
    oneRep<-function(i,M,n,k,mingr,maxgr,addParam,rep, parGenFun,...){
      temppar<-parGenFun(n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam)
      #skip.par<-c(skip.par,list(temppar))
      
      tres <- try(optParGlm(M=M, clu=temppar,  ...))
      if(inherits(x = tres,what = "try-error")){
        tres<-list("try-error"=tres, err=Inf, startPart=temppar)
      }
      if(deleteMs){
        tres$M<-NULL
      }
      tres$best<-NULL
      return(list(tres))
    }
    
    if(nCores==0){
      nCores<-parallel::detectCores()-1                    
    }
    
    if(useParLapply) {
      if(is.null(cl)) cl<-parallel::makeCluster(nCores)
      parallel::clusterSetRNGStream(cl)
      nC<-nCores
      res<-parallel::parLapplyLB(cl = cl,1:rep, fun = oneRep, M=M,n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam,rep=rep, parGenFun=parGenFun,...)
      if(stopcl) parallel::stopCluster(cl)
      res<-lapply(res,function(x)x[[1]])
    } else {
      res <- lapply(seq_len(rep), function(i) {
        oneRep(
          i = i,
          M = M,
          n = n,
          k = k,
          mingr = mingr,
          maxgr = maxgr,
          addParam = addParam,
          rep = rep,
          parGenFun = parGenFun,
          ...
        )
      })
      res <- lapply(res, function(x) x[[1]])
    }
    err<-sapply(res,function(x)x$err)    
  }
}




#' Compute ICL for a GLM-based blockmodel with hard partition
#'
#' Computes the Integrated Classification Likelihood (ICL) for a fitted
#' \code{glm} object and a corresponding hard partition of units (nodes).
#' The function assumes a conditional model where the partition defines
#' a block variable included in the GLM (e.g., as a factor).
#'
#' The ICL is computed as:
#' \deqn{
#'   \mathrm{ICL} = \log L(\hat{\theta} \mid y) +
#'   \sum_{k=1}^K n_k \log\left(\frac{n_k}{n}\right) -
#'   \frac{\nu}{2} \log(N_{\mathrm{eff}})
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\log L(\hat{\theta} \mid y)} is the log-likelihood of the fitted GLM,
#'   \item \eqn{n_k} are cluster sizes,
#'   \item \eqn{\nu} is the number of free parameters,
#'   \item \eqn{N_{\mathrm{eff}}} is the effective sample size (typically number of dyads).
#' }
#'
#' This formulation corresponds to a BIC approximation to the
#' complete-data likelihood used in classification-based (CEM-like)
#' blockmodeling approaches.
#'
#' @param fit A fitted \code{glm} object.
#' @param clu An integer vector of cluster memberships (length = number of nodes).
#' @param n_eff Optional numeric. Effective sample size. Defaults to
#'   \code{length(fit$y)} (number of observations/dyads).
#' @param include_mixing Logical. If \code{TRUE}, includes \eqn{K-1} parameters
#'   for cluster proportions in the penalty term. Default is \code{TRUE}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{icl}: ICL value
#'   \item \code{logLik}: log-likelihood of the model
#'   \item \code{ll_class}: classification term
#'   \item \code{n_params}: number of parameters used in penalty
#'   \item \code{n_eff}: effective sample size
#'   \item \code{K}: number of clusters
#'   \item \code{nk}: cluster sizes
#' }
#'
#' @details
#' The number of parameters is computed as:
#' \deqn{
#'   \nu = \text{length(coef(fit))} + (K - 1)
#' }
#' if \code{include_mixing = TRUE}, otherwise only the GLM coefficients
#' are counted.
#'
#' This function assumes that the GLM has been fitted conditional on the
#' partition (i.e., the block variable is included as a predictor).
#'
#' @examples
#' \dontrun{
#' fit <- glm(y ~ block + x1, family = poisson(), data = df)
#' clu <- df$cluster
#'
#' icl_res <- icl_from_glm(fit, clu)
#' icl_res$icl
#' }
#'
#' @export
icl_from_glm <- function(
    fit,
    clu,
    n_eff = NULL,
    include_mixing = TRUE
) {
  if (!inherits(fit, "glm")) {
    stop("fit must be a glm object")
  }
  
  nk <- table(clu)
  nk <- nk[nk > 0]
  K <- length(nk)
  n <- sum(nk)
  
  # log-likelihood
  ll <- as.numeric(logLik(fit))
  
  # classification term
  ll_class <- sum(nk * log(nk / n))
  
  # number of parameters
  n_params <- length(coef(fit))
  
  # effective sample size
  if (is.null(n_eff)) {
    if (!is.null(fit$y)) {
      n_eff <- length(fit$y)
    } else {
      stop("Provide n_eff or fit must contain y")
    }
  }
  
  # ICL
  icl <- ll + ll_class - 1/2 * (n_params * log(n_eff) + ifelse(include_mixing, (K - 1) * log(n), 0))
  
  #list(ICL = icl, logLik = ll, ll_class = ll_class, n_params = n_params,  n_eff = n_eff,  K = K, nk = nk)
  return(icl)
}


#' Optimize Random GLM Partitions Over a Range of Block Counts
#'
#' Fits GLM blockmodel partitions for a range of numbers of clusters and selects
#' the solution with the highest ICL value. For each `k` from `minK` to `maxK`,
#' the function runs `optRandomParGlm()`, except when `minK = 1`, where the
#' one-cluster solution is fitted directly with `optParGlm()`.
#'
#' @param M A matrix or network object passed to `optParGlm()` and
#'   `optRandomParGlm()`.
#' @param maxK Integer scalar. Maximum number of clusters to evaluate.
#' @param rep Integer scalar. Number of random repetitions passed to
#'   `optRandomParGlm()` for each value of `k`.
#' @param minK Integer scalar. Minimum number of clusters to evaluate. Defaults
#'   to `1`. If `minK = 1`, the one-cluster solution is fitted using a constant
#'   partition and the random search starts at `k = 2`.
#' @param ... Additional arguments passed to `optParGlm()` and
#'   `optRandomParGlm()`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{bestPart}{The partition of the best-fitting solution.}
#'   \item{bestFit}{The highest ICL value found.}
#'   \item{ICL}{A named vector of ICL values for each evaluated number of clusters.}
#'   \item{bestRes}{The full result object corresponding to the best solution.}
#'   \item{searchHistory}{A list of result objects for all evaluated values of `k`.}
#' }
#' The elapsed computation time is stored as the `"time"` attribute of the
#' returned object.
#'
#' @seealso `optParGlm()`, `optRandomParGlm()`
#'
#' @export
optRandomParRangeGlm<-function(M, maxK,rep,minK=1,...){
  tmpTime<-system.time({
    if(minK==1) {
      minKuse<-2
      initPar<-rep(1, nrow(M))
      res1<-list("1"=optParGlm(M=M, clu=initPar, ...))
    } else{
      minKuse<-minK
      res1<-NULL
    }
    searchHistory<-lapply(minKuse:maxK, \(k){
      cat("k =",k,"\n")
      optRandomParGlm(M=M, k=k, rep=rep,...)
    })
    searchHistory<-c(res1, searchHistory)
    names(searchHistory)<-minK:maxK
    ICL<-sapply(searchHistory, \(x) x$ICL)
    bestInd<-which.max(ICL)
    bestRes<-searchHistory[[bestInd]]
    res<-list(bestPart=bestRes$clu, bestFit=ICL[bestInd], ICL=ICL, bestRes=bestRes, searchHistory=searchHistory)
  })
  attr(res, "time")<-tmpTime
  res
}

optParGlmForMergeSplit<-function(data, initPart, ...){
  res<-optParGlm(M=data, clu=initPart, ...)
  return(list(part=res$clu, fit=res$ICL))
}
