fitted.pheno <- function(object, sub.sample, ...){

    if(object$summary.only){
        stop("fitted cannot be called using a summary.only=TRUE pheno model object.\n")
    }

    fitted <- 1#0 for prediction

    elip.args <- list(...)

    n.omp.threads <- 1
    if("n.omp.threads" %in% names(elip.args)){
        n.omp.threads <-  elip.args[["n.omp.threads"]]
    }
        
    verbose <- 0
    if("verbose" %in% names(elip.args)){
        verbose <-  elip.args[["verbose"]]
    }
    
    n.report <- 100
    if("n.report" %in% names(elip.args)){
        n.report <-  elip.args[["n.report"]]
    }
    
    if(class(object) == "pheno" & object$type[1] == "pheno.mod.0"){
        
        if(missing(sub.sample)){
            sub.sample <- list() #uses default sub.sample start, end, thin
        }        
        
        n.samples <- nrow(object$p.theta.samples)
        start <- ifelse(!"start" %in% names(sub.sample), floor(0.5*n.samples), sub.sample$start)
        end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
        thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)
        if(!is.numeric(start) || start >= n.samples){stop("invalid start")}
        if(!is.numeric(end) || end > n.samples){stop("invalid end")}
        if(!is.numeric(thin) || thin >= n.samples){stop("invalid thin")}
        sub.sample <- list(start=start, end=end, thin=thin)
        s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
        n.samples <- length(s.indx)
        
        ##get samples
        alpha <- t(as.matrix(object$p.theta.samples)[s.indx,grep("alpha.", colnames(object$p.theta.samples))])
        phi <- as.matrix(object$p.theta.samples)[s.indx,"sigma.sq"]
        
        t <- object$t
        n <- length(t)
        family.indx <- object$family.indx        
        t.normal.bounds <- object$t.normal.bounds
        
        storage.mode(t) <- "double"
        storage.mode(n) <- "integer"
        storage.mode(family.indx) <- "integer" ##normal 0, beta 1, truncated normal
        storage.mode(t.normal.bounds) <- "double"
        storage.mode(alpha) <- "double"
        storage.mode(phi) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(fitted) <- "integer"
        storage.mode(n.omp.threads) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        
        ptm <- proc.time()
        
        out <- .Call("phenoMod0Pred", t, n, family.indx, t.normal.bounds,
                     alpha, phi,
                     n.samples, fitted, n.omp.threads, verbose, n.report)
        
        out$run.time <- proc.time() - ptm
        
    }else{
        stop(paste0("fitted not implemented for object class ", class(object),"of model type", object$type[1]))
    }

    if("summary.probs" %in% names(elip.args)){
        probs <- elip.args[["summary.probs"]]
    }else{
        probs <- c(0.025, 0.5, 0.975)
    }

    out$p.fitted.summary <- t(apply(out$p.fitted.samples, 1, function(x){c(mean(x), sd(x), quantile(x, probs=probs))}))
    colnames(out$p.fitted.summary) <- c("mean", "sd", paste("p",100*probs,sep="_"))

    out
}


residuals.pheno <- function(object, sub.sample, ...){

    out <- fitted(object, sub.sample)   
    y <- object$y
        
    residuals.samples <- y-out$p.fitted.samples
    residuals.quants <- t(apply(residuals.samples, 1, function(x) quantile(x, prob=c(0.5, 0.025, 0.975))))
    
    return(list(residuals.samples=residuals.samples, residuals.quantiles=residuals.quants, sub.sample=out$sub.sample))
}

print.pheno <- function(x, ...){
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), "", sep = "\n")
}

plot.pheno <- function(x, ...){
    plot(x$p.theta.samples)
}

summary.pheno <- function(object, sub.sample, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = max(3L, getOption("digits") - 3L), ...){

    if(object$summary.only){
        stop("summary cannot be called using a summary.only=TRUE pheno model object.\n")
    }
    
    n.samples <- nrow(object$p.theta.samples)

    if(missing(sub.sample)){
        sub.sample <- list()
    }
    
    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)
    if(!is.numeric(start) || start >= n.samples){stop("invalid start")}
    if(!is.numeric(end) || end > n.samples){stop("invalid end")}
    if(!is.numeric(thin) || thin >= n.samples){stop("invalid thin")}
    sub.sample <- list(start=start, end=end, thin=thin)
    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))

    cat("Chain sub.sample:\n")
    cat(paste("start = ",start,"\n", sep=""))
    cat(paste("end = ",end,"\n", sep=""))
    cat(paste("thin = ",thin,"\n", sep=""))
    cat(paste("samples size = ",length(s.indx),"\n", sep=""))
    
    ##print(summary(mcmc(object$p.beta.samples[s.indx,]), quantiles=quantiles), digits=digits)
    print(noquote(apply(t(apply(object$p.theta.samples[s.indx,], 2, function(x) quantile(x, prob=quantiles))), 2, function(x) formatC(x, format = "f", digits = digits))))
    
}   

predict.pheno <- function(object, t.0, sub.sample, ...){
    
    if(object$summary.only){
        stop("predict cannot be called using a summary.only=TRUE pheno model object.\n")
    }

    fitted <- 0#1 for fitted
    
    elip.args <- list(...)

    n.omp.threads <- 1
    if("n.omp.threads" %in% names(elip.args)){
        n.omp.threads <-  elip.args[["n.omp.threads"]]
    }
        
    verbose <- 0
    if("verbose" %in% names(elip.args)){
        verbose <-  elip.args[["verbose"]]
    }
    
    n.report <- 100
    if("n.report" %in% names(elip.args)){
        n.report <-  elip.args[["n.report"]]
    }

    ##call
    cl <- match.call()
    
    if(class(object) == "pheno" & object$type[1] == "pheno.mod.0"){
        
        if(missing(sub.sample)){
            sub.sample <- list() #uses default sub.sample start, end, thin
        }        
        
        n.samples <- nrow(object$p.theta.samples)
        start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
        end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
        thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)
        if(!is.numeric(start) || start >= n.samples){stop("invalid start")}
        if(!is.numeric(end) || end > n.samples){stop("invalid end")}
        if(!is.numeric(thin) || thin >= n.samples){stop("invalid thin")}
        sub.sample <- list(start=start, end=end, thin=thin)
        s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
        n.samples <- length(s.indx)
        
        ##get samples
        alpha <- t(as.matrix(object$p.theta.samples)[s.indx,grep("alpha.", colnames(object$p.theta.samples))])
        phi <- as.matrix(object$p.theta.samples)[s.indx,"sigma.sq"]
        
        n <- length(t.0)
        family.indx <- object$family.indx        
        t.normal.bounds <- object$t.normal.bounds
        
        storage.mode(t.0) <- "double"
        storage.mode(n) <- "integer"
        
        storage.mode(family.indx) <- "integer" ##normal 0, beta 1, truncated normal
        storage.mode(t.normal.bounds) <- "double"
        storage.mode(alpha) <- "double"
        storage.mode(phi) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(fitted) <- "integer"
        storage.mode(n.omp.threads) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        
        ptm <- proc.time()
        
        out <- .Call("phenoMod0Pred", t.0, n, family.indx, t.normal.bounds,
                     alpha, phi,
                     n.samples, fitted, n.omp.threads, verbose, n.report)
        
        out$run.time <- proc.time() - ptm
        
    }else{
        stop(paste0("predict not implemented for object class ", class(object),"of model type", object$type[1]))
    }

    if("summary.probs" %in% names(elip.args)){
        probs <- elip.args[["summary.probs"]]
    }else{
        probs <- c(0.025, 0.5, 0.975)
    }

    out$p.predictive.summary <- t(apply(out$p.predictive.samples, 1, function(x){c(mean(x), sd(x), quantile(x, probs=probs))}))
    colnames(out$p.predictive.summary) <- c("mean", "sd", paste("p",100*probs,sep="_"))
    
    out$type <- object$type
    out$call <- cl
    class(out) <- "predict.pheno"
    out
}


print.predict.pheno <- function(x, ...){
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.75)), "", sep = "\n")
  
    cat("Predictions for a model of class ",class(x),", model type ",x$type[2]," and family ", x$type[3],".\n", sep="")
}
