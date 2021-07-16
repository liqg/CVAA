get_group_end <- function(samples, data){
  # samples must be ordered
  last_col <- samples[, ncol(samples)]# the last col
  sampleGroup <- c(0, cumsum(last_col[-length(last_col)] != last_col[-1]))
  sampleGroupEnd <- c((which(!duplicated(sampleGroup))-1)[-1], nrow(data))
  sampleGroupEnd
}

CVAA_each <- function(id, formula, data, samples, sampleGroupEnd, OR.col=NULL){

  invisible(gc())

  if(missing(sampleGroupEnd))
    sampleGroupEnd = get_group_end(samples=samples, data=data)

  freq <- freq3(data, sampleGroupEnd, id)

  terms <- apply(samples, 2, unique)
  if(ncol(samples)==1){
    terms <- list(a = unique(samples[,1]))
    names(terms) <- colnames(samples)
  }
  dim(freq) <- c(3, rev(sapply(terms, length)), ncol(data)-1)
  terms$compare <- c("smaller", "equal", "greater")
  terms$gene <- 1:(ncol(data)-1)
  dimnames(freq) <- terms[c("compare",rev(colnames(samples)),"gene")]

  a <- MASS::loglm(formula, freq, param=F); lrt <- a$lrt

  res <- c(id, lrt)

  #compute LOD, only for sample_type with two factor levels
  if(!is.null(OR.col)){
    a <- apply(freq, c(OR.col, "compare"), sum)
    lod <- try(log2(as.numeric(a[1,1])*a[2,3]/a[1,3]/a[2,1]))
    if(inherits(lod, "try-error")) lod <- NA
    res <- c(id, lrt, lod)
  }

  res
}

CVAA <- function(formula, data, samples, cores=4, OR.col=NULL, cl.type="PSOCK"){
  #' Croos-Values Association Analysis
  #'
  #' @param formula A formula. "gene" and "compare" are required and used internally.
  #' @param data A numeric matrix with n rows(samples) and m cols(genes/variables).
  #' @param samples A data.frame with sample information.
  #' @param cores How many cores are used.
  #' @param cl.type The cluster type for parallel computation.
  #' @examples
  #' \dontrun{
  #' library(CVAA)
  #' load(system.file("brca.samples.rda", package="CVAA"))
  #' load(system.file("brca.counts.rda", package="CVAA"))
  #' res <- CVAA(
  #' formula = formula("~(compare+sample_type)*cohort*gene"),
  #' data = brca.counts,
  #' samples=brca.samples,
  #' OR.col = "sample_type",
  #' cores = 4)
  #'
  #' }
  #'

  cols_f <- unique(sub(":.+","",labels(terms(formula))))

  stopifnot(all(c("compare", "gene") %in% cols_f))

  cols_f <- setdiff(cols_f, c("compare", "gene"))
  if(!all(cols_f %in% colnames(samples))){
    stop("terms in formula must be in colnames(samples)")
  }

  cols_f <- intersect(colnames(samples), cols_f)

  if(length(cols_f) > 1) {
    samples <- samples[, cols_f]
  }

  # check rownames consistency btw data and samples before sorting
  stopifnot(all(rownames(data) == rownames(samples)))

  data.table::setorderv(samples, colnames(samples))

  if(!all(rownames(data) == rownames(samples)))
    data <- data[rownames(samples), ]

  if(!is.null(OR.col)){
  if(OR.col %in% cols_f){
    if(length(unique(samples[[OR.col]])) != 2){
      stop("the OR.col must has 2 levels when estimating the odd rate")
    }
  }else{
    stop("OR.col must be in the formula")
  }}

  sampleGroupEnd <- get_group_end(samples=samples, data=data)

  # if(cores >1){
  #   cl <- parallel::makePSOCKcluster(cores, cl.type=cl.type)
  #   parallel::clusterExport(cl, c("data","samples","sampleGroupEnd"), envir = environment())
  #
  #   a <- parallel::clusterEvalQ(cl, library(CVAA))
  #
  #   res <- pbapply::pblapply(
  #     1:ncol(data),
  #     CVAA_each,
  #     formula=formula,
  #     data=data,
  #     samples=samples,
  #     sampleGroupEnd=sampleGroupEnd,
  #     OR.col=OR.col,
  #     cl=cl)
  #
  #   parallel::stopCluster(cl)
  # }else{
  #   res <- lapply(
  #     1:ncol(data), CVAA_each,
  #     formula=formula,
  #     data=data,
  #     samples=samples,
  #     sampleGroupEnd=sampleGroupEnd,
  #     OR.col=OR.col)
  # }

  res <- parallel::mclapply(
    1:ncol(data), CVAA_each,
    formula=formula,
    data=data,
    samples=samples,
    sampleGroupEnd=sampleGroupEnd,
    OR.col=OR.col, mc.cores = cores)

  res <- do.call(rbind, res)
  res <- as.data.frame(res)
  rownames(res) <- colnames(data)

  res <- res[order(res$V2, decreasing = T), ]
  res$rank <- 1:nrow(res)

  if(is.null(OR.col)){
    colnames(res) <- c("id","lrt","rank")
  }else{
    colnames(res) <- c("id","lrt", "lod","rank")
  }

  res$id <- NULL

  res
}

