#' Perform Up-and-Down Search Optimization
#'
#' The \code{mergeAndSplit} (alias \code{upAndDownSearch}) function performs an merge and split (increasing and lowering the number of clusters) on a given dataset, starting from an initial partition and fitness value. The function assumes that the higher fitness value is better. Experimental! The \code{stochBlockForUDS} is just a wrapper to the `stochBlock` function to be used within \code{mergeAndSplit} function.
#'
#' @param data A data object (e.g., data frame or matrix) on which the search is performed.
#' @param initPart Initial partition or configuration to start the search from. It can be a vector for simple datasets and a list of vectors for temporal or linked networks
#' @param initFit Initial fitness value associated with the initial partition.
#' @param optimFun A function used for optimization. A function must accept a data and a initial partition and return a list with an element \code{part} holding the final partition and an element \code{fit} holding the fitness value. The function can also return additional elements, which are stored in the search history (if saveAllHist=TRUE, dafult is FALSE) as is the best result. The \code{stochBlockForUDS} function is a wrapper to `stochBlock` that can be used as a default optimization function for stochastic block modeling.
#' @param saveAllHist Logical. Should the full optimization results for all evaluated candidates during the search be saved in the search history? If \code{FALSE}, fit and partitions are saved  Defaults to \code{FALSE}.
#' @param nRep Integer. The number of repetitions for the search algorithm. Defaults to 100.
#' @param minPmove Numeric. The minimum probability of moving to a new partition when the fitness does not improve. Defaults to 0.2. Applicable only to "alternating" and "forwardBackward" strategies when accepting worse solutions is allowed.
#' @param pFitSumTreshMoveBest Numeric. The threshold for the sum of 1 - probabilities of moving to a new partition before reverting to the best partition. Defaults to 5.
#' @param maxTry Integer. The maximum number of attempts to find a valid partition when increasing or decreasing the number of clusters. Defaults to 100.
#' @param minCluSize Integer. The minimum cluster size for a new cluster when splitting a cluster. Defaults to 1. The cluster must be at least twice as large as the minimum cluster size to be considered for splitting.
#' @param minNClust Integer or vector of integers. Minimal allowed number of clusters. For linked networks, a scalar is recycled to all sets and a vector is interpreted per set. Defaults to 1.
#' @param maxNClust Integer, vector of integers, or \code{NULL}. Maximal allowed number of clusters. For linked networks, a scalar is recycled to all sets and a vector is interpreted per set. It must be provided when \code{strategy = "deterministicCycle"}.
#' @param strategy Character. Search strategy. Use \code{"alternating"} for the behavior that alternates split and merge moves by iteration, \code{"forwardBackward"} for a criterion-guided search that evaluates split and merge candidates from the current solution and moves in the better direction, or \code{"deterministicCycle"} for a deterministic search that moves from the current number of clusters to the maximum, then to the minimum, and back again. Defaults to \code{"alternating"}.
#' @param nSplitMergeAttempts Integer. Number of split candidates and number of merge candidates to generate and optimize in each \code{"forwardBackward"} iteration. Defaults to 1.
#' @param maxNoImprove Integer. Number of consecutive \code{"forwardBackward"} iterations without improvement that are still allowed before stopping. Defaults to 0.
#' @param moveToWorseFB Logical. Should \code{"forwardBackward"} use the same Metropolis-style acceptance of worse solutions as \code{"alternating"}? If \code{FALSE}, the search keeps the current solution whenever the best candidate does not improve the fit. Defaults to \code{TRUE}.
#' @param initArgs List of additional arguments to be passed to \code{optimFun} when evaluating the initial partition. This is used only when \code{initFit} is not provided, in which case the initial fitness value is computed by calling \code{optimFun} with the initial partition and the arguments in \code{initArgs}.
#' @param useParallel Logical. Should candidate optimizations be evaluated in
#'   parallel when supported by the selected strategy? This affects
#'   \code{"forwardBackward"} and \code{"deterministicCycle"} only. Defaults to
#'   \code{FALSE}.
#' @param cl Optional cluster object created by \code{\link[parallel]{makeCluster}}.
#'   If not supplied and \code{useParallel = TRUE}, a temporary cluster is
#'   created automatically when at least two workers would be used. For
#'   \code{"forwardBackward"}, the default number of workers is
#'   \code{min(parallel::detectCores() - 1, 2 * nSplitMergeAttempts)} so split
#'   and merge candidates can be evaluated in one parallel batch. For
#'   \code{"deterministicCycle"}, this is used only when
#'   \code{nSplitMergeAttempts > 1}, with default number of workers
#'   \code{min(parallel::detectCores() - 1, nSplitMergeAttempts)}. Parallel
#'   candidate evaluation is skipped when the effective number of workers is 1.
#' @param ... Additional paramters to optimFun.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{data}{The input data.}
#'   \item{finalPart}{The final partition obtained.}
#'   \item{finalFit}{The final fitness value after the search.}
#'   \item{searchHistory}{A list containing the history of partitions and fitness values during the search.}
#'   \item{callUsed}{The call used to invoke the function, capturing the parameters passed.}
#'   \item{initial.param}{A list of initial parameters used in the function call withiut the data.}
#' }
#' @examples
#' \dontrun{
#' # Create a synthetic network matrix
#' set.seed(2022)
#' library(blockmodeling)
#' k<-2 # number of blocks to generate
#' blockSizes<-rep(20,k)
#' IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
#' clu<-rep(1:k, times=blockSizes)
#' n<-length(clu)
#' M<-matrix(rbinom(n*n,1,IM[clu,clu]),ncol=n, nrow=n)
#' initClu<-rep(1, times=n)
#' initFit<-ICLStochBlock(M, initClu) # Initial fitness value
#' # Using up-and-down search to optimise the partition
#' res<-upAndDownSearch(data=M,initPart=initClu, initFit=initFit, optimFun=stochBlockForUDS, nRep=10) 
#' plotMat(res$data, clu=res$bestPart) # Have a look at the optimised parition
#' print(res$bestFit) # Print the final fitness value
#' 
#' # Create a synthetic linked-network matrix
#' set.seed(2022)
#' library(blockmodeling)
#' IM<-matrix(c(0.9,.5,0.1,0.8), nrow=2)
#' clu<-rep(1:2, each=20) # Partition to generate
#' n<-length(clu)
#' nClu<-length(unique(clu)) # Number of clusters to generate
#' M1<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n) # First network
#' M2<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n) # Second network
#' M12<-diag(n) # Linking network
#' nn<-c(n,n)
#' k<-c(2,2)
#' Ml<-matrix(0, nrow=sum(nn),ncol=sum(nn)) 
#' Ml[1:n,1:n]<-M1
#' Ml[n+1:n,n+1:n]<-M2
#' Ml[n+1:n, 1:n]<-M12 
#' plotMat(Ml) # Linked network
#' clu1<-rep(1, n)
#' clu2<-rep(2, n)
#' initClu<-list(clu1, clu2)
#' initFit<-ICLStochBlock(Ml, initClu) # Initial fitness value
#' # Using up-and-down search to optimise the partition
#' res<-upAndDownSearch(data=Ml,initPart=initClu, initFit=initFit, optimFun=stochBlockForUDS, nRep=10)
#' plotMat(res$data, clu=res$bestPart) # Have a look at the optimised parition
#' print(res$bestFit) # Print the final fitness value
#' }
#' @export
upAndDownSearch<-function(data, initPart, initFit, optimFun, saveAllHist= FALSE, nRep=100, minPmove=0.2, pFitSumTreshMoveBest=5, maxTry=100, minCluSize=1, minNClust=1, maxNClust=NULL, strategy=c("alternating", "forwardBackward", "deterministicCycle"), nSplitMergeAttempts=1, maxNoImprove=0, moveToWorseFB=TRUE, initArgs=list(maxIter=0), useParallel=FALSE, cl=NULL,...){
  callUsed <- match.call()
  formal_args <- setdiff(names(formals()), "...")
  evaluated_formals <- lapply(formal_args, function(x) eval.parent(substitute(x)))
  names(evaluated_formals) <- formal_args
  strategy <- match.arg(strategy)
  
  # Capture evaluated ... arguments
  dot_args <- list(...)
  
  # Combine all initial parameters
  initial.param <- c(evaluated_formals, dot_args)
  initial.param$data<-NULL
  startTime <- proc.time()
  
  searchHistory<-list()
  linked<-is.list(initPart)
  if(missing(initFit)){
    initRes <- do.call(optimFun, c(list(data, initPart), dot_args,initArgs))
    initFit <- initRes$fit
    initPart <- initRes$part
    bestRes<-initRes
  } else{
    bestRes<-list(part=initPart, fit=initFit)
  }
  bestPart<-curPart<-initPart
  bestFit<-curFit<-initFit
  nSets <- if(linked) length(initPart) else 1
  if(linked) {
    sets<-sapply(curPart, length)  
  } else sets<-length(curPart)
  pFitSum<-0
  noImproveCount <- 0
  cycleDirection <- "forward"
  createdCluster <- FALSE
  parallelEnabled <- FALSE
  
  normalizeNClustBound <- function(bound, name, defaultValue = NULL){
    if(is.null(bound)){
      if(is.null(defaultValue)) return(NULL)
      return(rep(defaultValue, nSets))
    }
    if(length(bound) == 1) return(rep(bound, nSets))
    if(length(bound) != nSets){
      stop(sprintf("'%s' must have length 1 or %d.", name, nSets))
    }
    bound
  }
  
  minNClust <- normalizeNClustBound(minNClust, "minNClust")
  maxNClust <- normalizeNClustBound(maxNClust, "maxNClust", Inf)
  if(any(minNClust < 1)) stop("'minNClust' must be at least 1.")
  if(any(maxNClust < minNClust)) stop("'maxNClust' must be greater than or equal to 'minNClust'.")
  
  if(strategy == "deterministicCycle"){
    if(linked) stop("The 'deterministicCycle' strategy is not implemented for linked networks.")
    if(any(!is.finite(maxNClust))) stop("The 'deterministicCycle' strategy requires 'maxNClust' to be provided.")
    initK <- length(unique(initPart))
    if(initK > maxNClust[1]){
      warning("The current partition is larger than 'maxNClust'.")
    }
    if(initK >= maxNClust[1]){
      cycleDirection <- "backward"
    }
  }
  
  shouldUseParallel <- useParallel && (
    strategy == "forwardBackward" ||
      (strategy == "deterministicCycle" && nSplitMergeAttempts > 1)
  )
  
  if(shouldUseParallel && is.null(cl)){
    nDetectedCores <- parallel::detectCores()
    nAvailableCores <- if(is.na(nDetectedCores)) 1 else max(1, nDetectedCores - 1)
    nWorkers <- if(strategy == "forwardBackward"){
      min(nAvailableCores, 2 * nSplitMergeAttempts)
    } else{
      min(nAvailableCores, nSplitMergeAttempts)
    }
    if(nWorkers > 1){
      cl <- parallel::makeCluster(nWorkers)
      createdCluster <- TRUE
      parallelEnabled <- TRUE
    }
  } else if(shouldUseParallel && !is.null(cl)){
    parallelEnabled <- length(cl) > 1
  }

  proposeSplit <- function(part){
    if(linked){
      nClust <- sapply(part, function(x) length(unique(x)))
      eligibleSets <- which(nClust < maxNClust)
      if(length(eligibleSets) == 0) return(NULL)
      iSet <- sample(eligibleSets, prob = sets[eligibleSets], size = 1)
      for(i2 in 1:maxTry){
        newPar <- part
        iClu <- sample(newPar[[iSet]], 1)
        maxK <- max(newPar[[iSet]])
        iUnits <- which(newPar[[iSet]] == iClu)
        iCluSize <- length(iUnits)
        if(iCluSize >= minCluSize * 2){
          selUnits <- rbinom(iCluSize, size = 1, prob = 0.5) == 1
          if(sum(selUnits) >= minCluSize && (iCluSize - sum(selUnits)) >= minCluSize){
            newPar[[iSet]][iUnits[selUnits]] <- maxK + 1
            return(newPar)
          }
        }
      }
    } else{
      if(length(unique(part)) >= maxNClust[1]) return(NULL)
      for(i2 in 1:maxTry){
        newPar <- part
        iClu <- sample(newPar, 1)
        maxK <- max(newPar)
        iUnits <- which(newPar == iClu)
        iCluSize <- length(iUnits)
        if(iCluSize >= minCluSize * 2){
          selUnits <- rbinom(iCluSize, size = 1, prob = 0.5) == 1
          if(sum(selUnits) >= minCluSize && (iCluSize - sum(selUnits)) >= minCluSize){
            newPar[iUnits[selUnits]] <- maxK + 1
            return(newPar)
          }
        }
      }
    }
    return(NULL)
  }

  proposeMerge <- function(part){
    if(linked){
      nClust <- sapply(part, function(x) length(unique(x)))
      eligibleSets <- which(nClust > minNClust)
      if(length(eligibleSets) == 0) return(NULL)
      for(i2 in 1:maxTry){
        iSet <- sample(eligibleSets, prob = sets[eligibleSets], size = 1)
        newPar <- part
        clus <- unique(part[[iSet]])
        iClu <- sample(clus, 2)
        newPar[[iSet]][newPar[[iSet]] == iClu[2]] <- iClu[1]
        newPar[[iSet]] <- encodeToZeroIndexed(newPar[[iSet]])
        return(newPar)
      }
    } else{
      if(length(unique(part)) <= minNClust[1]) return(NULL)
      newPar <- part
      clus <- unique(part)
      iClu <- sample(clus, 2)
      newPar[newPar == iClu[2]] <- iClu[1]
      newPar <- encodeToZeroIndexed(newPar)
      return(newPar)
    }
    return(NULL)
  }
  
  makeCandidateTasks <- function(nAttempts, moveType){
    if(moveType == "split"){
      move <- "forward"
      direction <- "split"
    } else if(moveType == "merge"){
      move <- "backward"
      direction <- "merge"
    } else{
      stop(sprintf("Unknown move type '%s'.", moveType))
    }
    lapply(seq_len(nAttempts), function(iAttempt){
      list(
        candidateId = iAttempt,
        moveType = moveType,
        move = move,
        direction = direction
      )
    })
  }

  getWorkerExportNames <- function(args, fun){
    exportNames <- c(
      "data", "dot_args", "optimFun", "proposeSplit", "proposeMerge",
      "encodeToZeroIndexed", "maxTry", "minCluSize", "minNClust",
      "maxNClust", "linked", "sets"
    )

    getObjectGlobals <- function(x){
      if(inherits(x, "formula")){
        objEnv <- environment(x)
        if(is.null(objEnv)) return(character(0))
        candidateNames <- unique(all.names(x))
        return(candidateNames[candidateNames %in% ls(envir = objEnv, all.names = TRUE)])
      }
      if(is.function(x)){
        objEnv <- environment(x)
        if(is.null(objEnv)) return(character(0))
        globals <- tryCatch(
          codetools::findGlobals(x, merge = FALSE),
          error = function(e) list(functions = character(0), variables = character(0))
        )
        candidateNames <- unique(c(globals$variables, globals$functions))
        return(candidateNames[candidateNames %in% ls(envir = objEnv, all.names = TRUE)])
      }
      character(0)
    }

    funEnv <- environment(fun)
    if(!is.null(funEnv)){
      funGlobals <- tryCatch(
        codetools::findGlobals(fun, merge = FALSE),
        error = function(e) list(functions = character(0), variables = character(0))
      )
      funGlobalNames <- unique(c(funGlobals$variables, funGlobals$functions))
      exportNames <- c(
        exportNames,
        funGlobalNames[funGlobalNames %in% ls(envir = funEnv, all.names = TRUE)]
      )
    }

    if(length(args) > 0){
      argGlobals <- unlist(lapply(args, getObjectGlobals), use.names = FALSE)
      exportNames <- c(exportNames, argGlobals)
    }

    unique(exportNames)
  }

  prepareCandidateCluster <- function(cl){
    parallel::clusterEvalQ(cl, {
      NULL
    })
    exportNames <- getWorkerExportNames(dot_args, optimFun)
    if(length(exportNames) > 0){
      parallel::clusterExport(cl, exportNames, envir = environment(proposeSplit))
    }
    invisible(cl)
  }

  if(parallelEnabled){
    prepareCandidateCluster(cl)
  }

  collectCandidates <- function(part, tasks){
    doOneCandidate <- function(task){
      moveFun <- switch(
        task$moveType,
        split = proposeSplit,
        merge = proposeMerge,
        stop(sprintf("Unknown move type '%s'.", task$moveType))
      )
      newPar <- moveFun(part)
      if(is.null(newPar)) return(NULL)
      iRes <- do.call(optimFun, c(list(data, newPar), dot_args))
      iRes$move <- task$move
      iRes$direction <- task$direction
      iRes$candidateId <- task$candidateId
      iRes
    }
    if(parallelEnabled){
      candidates <- parallel::parLapply(cl, tasks, doOneCandidate)
      candidates <- Filter(Negate(is.null), candidates)
    } else{
      candidates <- vector("list", length = length(tasks))
      for(iTask in seq_along(tasks)){
        candidates[[iTask]] <- doOneCandidate(tasks[[iTask]])
      }
      candidates <- Filter(Negate(is.null), candidates)
    }
    candidates
  }
  
  
  on.exit({
    elapsedTime <- proc.time() - startTime
    if(createdCluster) parallel::stopCluster(cl)
    result <- structure(
      list(
        data=data,
        bestPart=bestPart,
        bestFit=bestFit,
        bestRes=bestRes,
        searchHistory=searchHistory,
        callUsed=callUsed,
        initial.param=initial.param
      ),
      class="mergeSplitSearch"
    )
    attr(result, "time") <- elapsedTime
    return(result)
  })
  
  for(i in 1:nRep){
    cat(sprintf("Starting iteration %d/%d\n", i, nRep))
    if(strategy == "alternating"){
      if(i%%2 == 1){ #odd iterations - try increasing number of clusters
        newPar <- proposeSplit(curPart)
      } else { #even iterations - try decreasing number of clusters
        newPar <- proposeMerge(curPart)
      }
      
      if(is.null(newPar)){
        cat("No valid candidate partition could be generated.\n")
        next
      }
      
      iRes<-optimFun(data, newPar, ...)
      cat("Number of clusters:\n")
      cat(numClust(iRes$part))
      cat(sprintf("\nNew fitness %f\n", iRes$fit))
  
      if(iRes$fit>curFit){
        curFit<-iRes$fit
        curPart<-iRes$part
        if(curFit>bestFit){
          cat("Fitness improved!\n")
          bestFit<-curFit
          bestPart<-curPart
          pFitSum<-0
          bestRes <- iRes
        }
      } else{
        pFit<-exp(iRes$fit-curFit)
        p<-max(pFit, minPmove)
        if(runif(1) <= p){
          curFit<-iRes$fit
          curPart<-iRes$part
          cat("Moved to a worse solution.\n")
        }
        cat("pFit:",pFit,"\n")
        pFitSum<-pFitSum + 1 - pFit
        if(pFitSum > pFitSumTreshMoveBest && pFit!=1){
          pFitSum<-0
          curFit<-bestFit
          curPart<-bestPart
          cat("\nFitness did not improve for a while, reverting to best fitness!\n")
          cat("Number of clusters for currently best solution:\n")
          cat(numClust(bestPart))
          cat("\n")
        }
        
      }
      iRes$statePart <- curPart
      iRes$stateFit <- curFit
      iRes$accepted <- identical(curPart, iRes$part) && isTRUE(all.equal(curFit, iRes$fit))
      if(!saveAllHist) iRes <- list(part = iRes$part,fit = iRes$fit)
      searchHistory[[i]]<-iRes
    } else if(strategy == "forwardBackward"){
      candidateTasks <- c(
        makeCandidateTasks(nSplitMergeAttempts, "split"),
        makeCandidateTasks(nSplitMergeAttempts, "merge")
      )
      candidates <- collectCandidates(curPart, candidateTasks)

      if(length(candidates) == 0){
        cat("No valid split or merge candidate could be generated.\n")
        next
      }
      
      candFits <- sapply(candidates, function(x) x$fit)
      iBestCand <- which.max(candFits)
      iRes <- candidates[[iBestCand]]
      cat(sprintf("Best %s move selected (%s, candidate %d).\n", iRes$move, iRes$direction, iRes$candidateId))
      cat("Number of clusters:\n")
      cat(numClust(iRes$part))
      cat(sprintf("\nNew fitness %f\n", iRes$fit))
      accepted <- FALSE
      
      if(iRes$fit > curFit){
        curFit <- iRes$fit
        curPart <- iRes$part
        pFitSum <- 0
        noImproveCount <- 0
        accepted <- TRUE
        if(curFit > bestFit){
          cat("Fitness improved!\n")
          bestFit <- curFit
          bestPart <- curPart
          bestRes <- iRes
        }
      } else{
        if(moveToWorseFB){
          pFit <- exp(iRes$fit - curFit)
          p <- max(pFit, minPmove)
          if(runif(1) <= p){
            curFit <- iRes$fit
            curPart <- iRes$part
            cat("Moved to a worse solution.\n")
            accepted <- TRUE
          } else{
            cat("Kept the current solution.\n")
          }
          cat("pFit:",pFit,"\n")
        } else{
          cat("Kept the current solution.\n")
        }
        noImproveCount <- noImproveCount + 1
        cat(sprintf("Neither forward nor backward move improved the current solution. Consecutive non-improving iterations: %d.\n", noImproveCount))
        if(noImproveCount > maxNoImprove){
          cat("Maximum allowed number of consecutive non-improving iterations reached. Stopping search.\n")
          break
        }
      }
      if(!saveAllHist) iRes <- list(part = iRes$part,fit = iRes$fit)
      searchHistory[[i]] <- c(iRes, list(accepted = accepted))
    } else{
      curK <- length(unique(curPart))
      if(curK >= maxNClust[1]){
        cycleDirection <- "backward"
      } else if(curK <= minNClust[1]){
        cycleDirection <- "forward"
      }
      
      if(cycleDirection == "forward"){
        candidates <- collectCandidates(curPart, makeCandidateTasks(nSplitMergeAttempts, "split"))
      } else{
        candidates <- collectCandidates(curPart, makeCandidateTasks(nSplitMergeAttempts, "merge"))
      }

      if(length(candidates) == 0){
        stop(sprintf("No valid %s candidate could be generated for the 'deterministicCycle' strategy.", cycleDirection))
      }
      
      candFits <- sapply(candidates, function(x) x$fit)
      iBestCand <- which.max(candFits)
      iRes <- candidates[[iBestCand]]
      cat(sprintf("Deterministic %s move selected (%s, candidate %d).\n", iRes$move, iRes$direction, iRes$candidateId))
      cat("Number of clusters:\n")
      cat(numClust(iRes$part))
      cat(sprintf("\nNew fitness %f\n", iRes$fit))
      
      curFit <- iRes$fit
      curPart <- iRes$part
      if(curFit > bestFit){
        cat("Fitness improved!\n")
        bestFit <- curFit
        bestPart <- curPart
        bestRes <- iRes
      }
      if(!saveAllHist) iRes <- list(part = iRes$part,fit = iRes$fit)
      searchHistory[[i]] <- c(iRes, list(accepted = TRUE))
    }
    cat(sprintf("Best fitness %f\n\n", bestFit))
  }
}

#' @rdname upAndDownSearch
#' @export
mergeAndSplit<-upAndDownSearch


#' @rdname upAndDownSearch
#' @export
stochBlockForUDS<-function(data, initPart, ...){
  res<-stochBlock(M=data, clu=initPart, ...)
  return(list(part=res$clu, fit=res$ICL))
}


#' Extract the Best Result for a Given Number of Clusters
#'
#' The \code{getBestNClust} function extracts the best solution encountered by
#' \code{\link{upAndDownSearch}} for a requested number of clusters. For linked
#' networks, the requested number of clusters can be given as a scalar, which is
#' then recycled to all sets, or as a vector with one element per set.
#'
#' @param x Object of class \code{"mergeSplitSearch"} returned by
#'   \code{\link{upAndDownSearch}}.
#' @param nClust Integer or integer vector specifying the requested number of
#'   clusters. For linked networks, a scalar is recycled to all sets.
#'
#' @return A list with elements \code{part}, \code{fit}, \code{iteration}, and
#'   \code{nClust}. Returns \code{NULL} if no matching solution is found.
#' @export
getBestNClust<-function(x, nClust){
  ## Back compatibility fix
  if(is.null(x$searchHistory) && !is.null(x$searchHistroy)){
    x$searchHistory <- x$searchHistroy
    x$searchHistroy <- NULL
  }
  
  if(is.null(x$searchHistory)) stop("Argument 'x' does not appear to be a valid result of 'upAndDownSearch'.")
  
  linked <- is.list(x$bestPart)
  if(linked){
    nSets <- length(x$bestPart)
    if(length(nClust) == 1) nClust <- rep(nClust, nSets)
    if(length(nClust) != nSets) stop(sprintf("'nClust' must have length 1 or %d for linked networks.", nSets))
  } else{
    if(length(nClust) != 1) stop("'nClust' must have length 1 for non-linked networks.")
  }
  
  getEntryState <- function(entry){
    if(is.list(entry) && !is.null(entry$statePart) && !is.null(entry$stateFit)){
      return(list(part=entry$statePart, fit=entry$stateFit))
    }
    if(is.list(entry) && !is.null(entry$part) && !is.null(entry$fit)){
      return(list(part=entry$part, fit=entry$fit))
    }
    if(is.list(entry) && !is.null(entry$selected) && !is.null(entry$selected$part) && !is.null(entry$selected$fit)){
      return(list(part=entry$selected$part, fit=entry$selected$fit))
    }
    return(NULL)
  }
  
  matchesRequestedNClust <- function(part){
    if(linked){
      observedNClust <- sapply(part, function(x) length(unique(x)))
      identical(as.integer(observedNClust), as.integer(nClust))
    } else{
      length(unique(part)) == nClust
    }
  }
  
  ## Back compatibility fix
  if(is.null(x$searchHistory) && !is.null(x$searchHistroy)){
    x$searchHistory <- x$searchHistroy
    x$searchHistroy <- NULL
  }
  history <- x$searchHistory
  history <- history[!vapply(history, is.null, logical(1))]
  
  candidates <- vector("list", length(history))
  iCand <- 0
  for(i in seq_along(history)){
    state <- getEntryState(history[[i]])
    if(is.null(state)) next
    if(matchesRequestedNClust(state$part)){
      iCand <- iCand + 1
      candidates[[iCand]] <- list(
        part=state$part,
        fit=state$fit,
        iteration=i,
        nClust=if(linked) sapply(state$part, function(x) length(unique(x))) else length(unique(state$part))
      )
    }
  }
  candidates <- candidates[seq_len(iCand)]
  
  if(length(candidates) == 0) return(NULL)
  
  bestIdx <- which.max(vapply(candidates, function(x) x$fit, numeric(1)))
  candidates[[bestIdx]]
}


#' Plot Search History from Up-and-Down Search
#'
#' The \code{plot} method for objects of class \code{"mergeSplitSearch"} plots
#' the search history produced by \code{\link{upAndDownSearch}} for non-linked
#' networks. It uses the actual
#' state after each iteration when such information is available in
#' \code{searchHistory}, which is important for recent search strategies where
#' multiple candidates can be tried within one iteration.
#'
#' @param x Object of class \code{"mergeSplitSearch"} returned by
#'   \code{\link{upAndDownSearch}}.
#' @param palette Character vector of colors used to generate the color scale for
#'   the number of clusters. Defaults to \code{c("#440154", "#21908C",
#'   "#FDE725")}.
#' @param what Character. Type of plot to produce. Use \code{"iterations"} for
#'   the search path over iterations, or \code{"best"} to plot only the best fit
#'   found for each number of clusters against \code{K}. Defaults to
#'   \code{"iterations"}.
#' @param useKLabels Logical. If \code{TRUE}, the plot uses labels derived from
#'   the number of clusters. If \code{FALSE}, plain dots are used. Defaults to
#'   \code{TRUE}.
#' @param cexTransform Character or function. Transformation applied to
#'   \code{nchar(as.character(K))} before dividing by its maximum. Supported
#'   character values are \code{"sqrt"}, \code{"identity"}, and
#'   \code{"log"}. Defaults to \code{"sqrt"}.
#' @param dotPch Plotting character to use when \code{useKLabels = FALSE}.
#'   Defaults to \code{16}.
#' @param ... Additional graphical parameters passed to \code{\link{plot}}.
#'
#' @return Invisibly returns a data frame with iteration number, number of
#'   clusters, fit value, plotting color, plotting label, and relative text
#'   size. For \code{what = "best"}, the returned data frame contains one row
#'   for the best fit found at each observed \code{K}.
#' @method plot mergeSplitSearch
#' @export
plot.mergeSplitSearch<-function(x, palette=c("#440154", "#21908C", "#FDE725"), what=c("iterations", "best"), useKLabels=TRUE, cexTransform="sqrt", dotPch=16, ...){
  if(is.null(x$bestPart)) stop("Argument 'x' does not appear to be a valid result of 'upAndDownSearch'.")
  if(is.list(x$bestPart)) stop("The plot method for 'mergeSplitSearch' is implemented only for non-linked networks.")
  if(length(palette) < 2) stop("'palette' must contain at least two colors.")
  what <- match.arg(what)
  if(is.character(cexTransform)){
    cexTransform <- match.arg(cexTransform, c("sqrt", "identity", "log"))
    cexTransform <- switch(cexTransform, sqrt=sqrt, identity=identity, log=function(x)log(x+1))
  } else if(!is.function(cexTransform)){
    stop("'cexTransform' must be either a function or one of 'sqrt', 'identity', and 'log'.")
  }
  
  history <- x$searchHistory
  if(is.null(history))  history <- x$searchHistroy ## back compatibility fix
  history <- history[!vapply(history, is.null, logical(1))]
  if(length(history) == 0) stop("The search history is empty.")
  
  getHistoryState <- function(entry){
    if(is.list(entry) && !is.null(entry$statePart) && !is.null(entry$stateFit)){
      return(c(K=length(unique(entry$statePart)), fit=entry$stateFit))
    }
    if(is.list(entry) && !is.null(entry$part) && !is.null(entry$fit)){
      if(is.list(entry$part)) stop("The plot method for 'mergeSplitSearch' is implemented only for non-linked networks.")
      return(c(K=length(unique(entry$part)), fit=entry$fit))
    }
    if(is.list(entry) && !is.null(entry$selected) && !is.null(entry$selected$part) && !is.null(entry$selected$fit)){
      return(c(K=length(unique(entry$selected$part)), fit=entry$selected$fit))
    }
    try({
      part<-clu(entry)
      if(is.list(entry) && !is.null(part) && !is.null(entry$ICL)){
        if(is.list(part)) stop("The plot method for 'mergeSplitSearch' is implemented only for non-linked networks.")
        return(c(K=length(unique(part)), fit=entry$ICL))
      }
    })
    stop("Unsupported structure in 'searchHistory'.")
  }
  
  resTmpDf <- t(vapply(history, getHistoryState, numeric(2)))
  resTmpDf <- as.data.frame(resTmpDf)
  resTmpDf$iteration <- seq_len(nrow(resTmpDf))
  if(what == "best"){
    resTmpDf <- do.call(rbind, lapply(split(resTmpDf, resTmpDf$K), function(df) df[which.max(df$fit), , drop=FALSE]))
    resTmpDf <- resTmpDf[order(resTmpDf$K), , drop=FALSE]
    rownames(resTmpDf) <- NULL
  }
  resTmpDf$Kcol <- colorRampPalette(palette)(max(resTmpDf$K) - min(resTmpDf$K) + 1)[resTmpDf$K - min(resTmpDf$K) + 1]
  resTmpDf$Kcex <- cexTransform(nchar(as.character(resTmpDf$K)))
  resTmpDf$Kcex <- resTmpDf$Kcex / max(resTmpDf$Kcex)
  resTmpDf$Kchar <- substr(as.character(resTmpDf$K), 1, 1)
  
  if(what == "iterations"){
    plot(resTmpDf$fit ~ resTmpDf$iteration, type="b", pch=if(useKLabels) "" else dotPch, ylab="fit", xlab="iteration", ...)
    if(useKLabels){
      text(resTmpDf$fit ~ resTmpDf$iteration, labels=resTmpDf$Kchar, col=resTmpDf$Kcol, cex=resTmpDf$Kcex)
    }
  } else{
    plot(resTmpDf$fit ~ resTmpDf$K, type="b", pch=if(useKLabels) "" else dotPch, ylab="fit", xlab="number of clusters", ...)
    if(useKLabels){
      text(resTmpDf$fit ~ resTmpDf$K, labels=resTmpDf$Kchar, col=resTmpDf$Kcol, cex=resTmpDf$Kcex)
    }
  }
  
  invisible(resTmpDf)
}
