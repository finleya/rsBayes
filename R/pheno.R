pheno <- function(formula, data = parent.frame(), gamma = c(0, 1), family="normal",
                  starting, tuning, priors, n.samples, sub.sample,
                  summary.only = FALSE, fitted = FALSE, verbose = TRUE, n.report=100, ...){

    ####################################################
    ##Check for unused args
    ####################################################
    formal.args <- c(names(formals(sys.function(sys.parent()))), "summary.probs", "t.normal.bounds")

    elip.args <- list(...)
    
    for(i in names(elip.args)){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }

    ##call
    cl <- match.call()

    ####################################################
    ##Formula
    ####################################################
    if(missing(formula)){stop("error: formula must be specified")}

    if(class(formula) == "formula"){

        holder <- parseFormula(formula, data)
        y <- as.vector(holder[[1]])
        t <- as.matrix(holder[[2]])[,2]

    }else{
        stop("error: formula is misspecified")
    }

    n <- length(y)

    ####################################################
    ##gamma (the theoretical range of y)
    ####################################################
    if(!all(is.finite(gamma))){
        stop(paste0("error: values of gamma must be finite."))
    }
    
    if(max(y) > gamma[2] || min(y) < gamma[1]){
     stop(paste0("error: found y values outside the theoretical bounds of y given by the gamma argument."))
    }

    if(gamma[1] >= gamma[2]){
        stop(paste0("error: gamma[2] must be larger than gamma[1]."))
    }

    if(any(gamma[1] < 0, gamma[1] > 1, gamma[2] < 0, gamma[2] > 1) & family == "beta"){
        stop(paste0("error: family = beta is not appropriate for gamma outside the (0,1) bounds. Choose a different family or rescale your y to (0,1)."))
    }

    storage.mode(gamma) <- "double"
    ####################################################
    ##family
    ####################################################
    if(!family %in% c("beta", "normal", "t.normal")){stop(paste0("family=", family, " is not supported."))}

    family.indx <- 0

    if(family == "beta"){
        family.indx <- 0
    }else if(family == "normal"){
        family.indx <- 1
    }else{
        family.indx <- 2
    }
    
    t.normal.bounds <- c(0,1)
    if("t.normal.bounds" %in% names(elip.args)){
        t.normal.bounds <- elip.args[["t.normal.bounds"]]
    }
    
    if(family.indx == 2){
        if(length(t.normal.bounds) != 2 || t.normal.bounds[2] <= t.normal.bounds[1]){
            stop("t.normal.bounds must be a vector of length 2 with element 1 < element 2")
        }
        
        if(any(c(y > t.normal.bounds[2], y < t.normal.bounds[1]))){
            stop("All regression y values must be within t.normal.bounds.")
        }
        storage.mode(t.normal.bounds) <- "double"
    }
            
    storage.mode(family.indx) <- "integer"

    ####################################################
    ##Priors
    ####################################################
    alpha.Unif <- 0
    phi.IG <- 0
    
    if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
    names(priors) <- tolower(names(priors))
    
    if(!"sigma.sq.ig" %in% names(priors)){stop("error: sigmaSq.IG must be specified")}
    phi.IG <- priors[["sigma.sq.ig"]]

    if(!is.vector(phi.IG) || length(phi.IG) != 2){stop("error: sigma.sq.IG must be a vector of length 2")}
    if(any(phi.IG <= 0)){stop("error: sigma.sq.IG must be a positive vector of length 2")}

    ##default alpha priors
    n.alpha <- 7
    alpha.Unif.a <- c(gamma[1],0,0,1,0,0,1)
    alpha.Unif.b <- c(gamma[2],1,1,1,0.001,1,366)
    
    if("alpha" %in% names(priors)){

        alpha <- priors$alpha
        
        for(i in 1:length(alpha)){
            
            alpha.indx <- as.numeric(sub(".*\\.", "",names(alpha[i])))
            if(!alpha.indx %in% 1:n.alpha){
                stop(paste0("error: ",names(alpha)[i])," is not a parameter in the model")
            }
            if(verbose){
                cat(paste0("Modifying default prior for ",names(alpha)[i]),"\n")
            }
            if(!is.vector(alpha[[i]]) || length(alpha[[i]]) != 2){
                stop(paste0("error: ",names(alpha)[i])," must be a vector of length 2")
            }
            
            alpha.Unif.a[alpha.indx] <- alpha[[i]][1]
            alpha.Unif.b[alpha.indx] <- alpha[[i]][2]
        }
    }

    ##make sure alpha.1's prior is within gamma
    if(alpha.Unif.a[1] < gamma[1] || alpha.Unif.b[1] > gamma[2]){
        stop(paste0("error: alpha.1 prior must be within the theoretical bounds of y given by the gamma argument."))
    }
    
    ####################################################
    ##Starting values
    ####################################################
    alphas <- paste0("alpha.",1:n.alpha)
    params <- c("sigma.sq",alphas)

    if(missing(starting)){
        if(verbose){
            cat("Starting values are not provided. Will attempt to generate reasonable ones given the priors.\n")
        }
        starting <- list(phi = 1.0/rgamma(1, shape=phi.IG[1], rate=phi.IG[2]))
        
        starting[["alpha.3"]] <- runif(1, alpha.Unif.a[3], alpha.Unif.b[3])
        starting[["alpha.5"]] <- runif(1, alpha.Unif.a[5], alpha.Unif.b[5])
        starting[["alpha.6"]] <- runif(1, alpha.Unif.a[6], alpha.Unif.b[6])
        
        starting[["alpha.1"]] <- runif(1, alpha.Unif.a[1], alpha.Unif.b[1])
        starting[["alpha.2"]] <- runif(1, alpha.Unif.a[2], gamma[2]-starting[["alpha.1"]])
        starting[["alpha.7"]] <- runif(1, alpha.Unif.a[7], alpha.Unif.b[7])
        starting[["alpha.4"]] <- runif(1, alpha.Unif.a[4], starting[["alpha.7"]])
    }else{
        names(starting) <- tolower(names(starting))
        
        if(!all(params %in% names(starting))){
            stop(paste0("Starting values are not provided for ", paste0(params[!params %in% names(starting)], collapse=", "), "."))
        }
        
        ##check if starting is valid given priors
        for(a in c("alpha.1","alpha.3","alpha.5","alpha.6","alpha.7")){
            a.indx <- which(a == alphas)
            if(starting[[a]] <= alpha.Unif.a[a.indx] || starting[[a]] >= alpha.Unif.b[a.indx]){
                stop(paste0("Starting value for ",a," is not within its prior support of Unif(", alpha.Unif.a[a.indx], ", ", alpha.Unif.b[a.indx],")."))
            }
        }

        if(starting[["alpha.2"]] <= alpha.Unif.a[2] || starting[["alpha.2"]] >= gamma[2]-starting[["alpha.1"]]){
            stop(paste0("Starting value for alpha.2 is not within its prior support of Unif(",
                        alpha.Unif.a[2], ", ", gamma[2]-starting[["alpha.1"]],"), where the upper prior bound is gamma[2]-alpha.1"))
        }
        
        if(starting[["alpha.4"]] <= alpha.Unif.a[4] || starting[["alpha.4"]] >= starting[["alpha.7"]]){
            stop(paste0("Starting value for alpha.4 is not within its prior support of Unif(",
                        alpha.Unif.a[4], ", ", starting[["alpha.7"]],"), where the upper prior bound is alpha.7"))
        }
        
    }
    
    ##update priors based on starting
    alpha.Unif.b[2] <- gamma[2]-starting[["alpha.1"]]
    alpha.Unif.b[4] <- starting[["alpha.7"]]
    
    alpha.Unif.a <- unname(unlist(alpha.Unif.a))
    alpha.Unif.b <- unname(unlist(alpha.Unif.b))
    
    ##get the order right for alpha
    phi.starting <- starting[["sigma.sq"]]
    
    alpha.starting <- starting[order(names(starting[names(starting) %in% alphas]))]
    alpha.starting <- unname(unlist(alpha.starting))

    storage.mode(phi.IG) <- "double"
    storage.mode(alpha.Unif.a) <- "double"
    storage.mode(alpha.Unif.b) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(alpha.starting) <- "double"

    ####################################################
    ##Tuning values
    ####################################################
    if(missing(tuning)){
        stop("Tuning values are not provided.\n")
    }else{
        names(tuning) <- tolower(names(tuning))

        if(!all(params %in% names(tuning))){
            stop(paste0("Tuning values are not provided for ", paste0(params[!params %in% names(tuning)], collapse=", "), "."))
        }
    }

    phi.tuning <- tuning[["sigma.sq"]]
    alpha.tuning <- tuning[order(names(tuning[names(tuning) %in% alphas]))]
    alpha.tuning <- unname(unlist(alpha.tuning))

    storage.mode(phi.tuning) <- "double"
    storage.mode(alpha.tuning) <- "double"

    ####################################################
    ##Other stuff
    ####################################################
    if(missing(sub.sample)){
        sub.sample <- list()
    }

    start <- ifelse(!"start" %in% names(sub.sample), 1, sub.sample$start)
    end <- ifelse(!"end" %in% names(sub.sample), n.samples, sub.sample$end)
    thin <- ifelse(!"thin" %in% names(sub.sample), 1, sub.sample$thin)
    if(!is.numeric(start) || start >= n.samples){stop("error: in sub.sample, invalid start")}
    if(!is.numeric(end) || end > n.samples){stop("error: in sub.sample, invalid end")}
    if(!is.numeric(thin) || thin >= n.samples){stop("error: in sub.sample, invalid thin")}

    s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
    sub.sample <- list(start=start, end=end, thin=thin)

    storage.mode(n.samples) <- "integer"
    storage.mode(fitted) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"

    ####################################################
    ##Pack it up and off it goes
    ####################################################
    ptm <- proc.time()

        out <- .Call("phenoMod0", y, t, n, family.indx, t.normal.bounds,
                     gamma,
                     phi.IG, alpha.Unif.a, alpha.Unif.b,
                     phi.starting, alpha.starting,
                     phi.tuning, alpha.tuning,
                     n.samples, fitted, verbose, n.report)

    out$run.time <- proc.time() - ptm

    rownames(out$p.theta.samples) <- c(alphas, "sigma.sq")
    out$p.theta.samples <- mcmc(t(out$p.theta.samples[,s.indx]))

    if("summary.probs" %in% names(elip.args)){
        probs <- elip.args[["summary.probs"]]
    }else{
        probs <- c(0.025, 0.5, 0.975)
    }
    
    out$p.theta.summary <- apply(out$p.theta.samples, 2, function(x){c(mean(x), sd(x), quantile(x, probs=probs))})
    rownames(out$p.theta.summary) <- c("mean", "sd", paste("p",100*probs,sep="_"))
        
    if(fitted){
        out$p.fitted.samples <- out$p.fitted.samples[,s.indx]
        out$p.fitted.summary <- t(apply(out$p.fitted.samples, 1, function(x){c(mean(x), sd(x), quantile(x, probs=probs))}))
        colnames(out$p.fitted.summary) <- c("mean", "sd", paste("p",100*probs,sep="_"))
    }

    if(summary.only){
        out$p.theta.samples <- NULL
        if(fitted){
            out$p.fitted.samples <- NULL
        }
    }
    
    out$call <- cl
    out$y <- y
    out$t <- t
    out$gamma  <- gamma
    out$priors <- priors
    out$tuning <- tuning
    names(out$MH.acceptance) <- c("overall","last batch")
    out$sub.sample <- sub.sample
    out$summary.only <- summary.only
    out$family.indx <- family.indx
    out$t.normal.bounds <- t.normal.bounds

    out$type <- c("pheno.mod.0", family)
    class(out) <- "pheno"

    out
}
