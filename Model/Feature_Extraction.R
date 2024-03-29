feature_extractor <- function(X, bp = c(0, 22), wl = 2048, threshold = 5, parallel = 1) {
    if(class(X) == "data.frame") {
        if(all(c("sound.files", "selec", "start", "end") %in% colnames(X))) {
            start <- as.numeric(unlist(X$start))
            end <- as.numeric(unlist(X$end))
            sound.files <- as.character(unlist(X$sound.files))
            selec <- as.character(unlist(X$selec))
        } else {
            stop(paste(paste(c("sound.files", "selec", "start", "end")[!(c("sound.files", "selec", "start", "end") %in% colnames(X))], collapse=","), "column(s) not found in the given data frame."))

        }
    } else stop("X is not a data frame.")

    #Stop Cases

    #1. If NA present in 'start' or 'stop'
    if(any(is.na(c(end, start)))) stop("NA found in 'start' and/or 'end'")

    #2. If 'start' or 'end' are not numeric
    if(all(class(end) != "numeric" & class(start) != "numeric")) stop("'end' and 'start' must be numeric")

    #3. If any 'start' is greater than the 'end'
    if(any(end - start<0)) stop(paste("The start is higher than the end in", length(which(end - start<0)), "case(s)"))  

    #4. If any selections are longer than 20 seconds
    if(any(end - start>20)) stop(paste(length(which(end - start>20)), "selection(s) longer than 20 sec"))  
    options( show.error.messages = TRUE)

    #5. If 'bp' is a. not a vector b. length!=2 
    if(!is.vector(bp)) stop("'bp' must be a numeric vector of length 2") 
    else{
        if(!length(bp) == 2) stop("'bp' must be a numeric vector of length 2")
    }

    #6. Return a warning if not all sound files were found
    fs <- list.files(path = getwd(), pattern = ".wav$", ignore.case = TRUE)
    if(length(unique(sound.files[(sound.files %in% fs)])) != length(unique(sound.files))) 
    cat(paste(length(unique(sound.files))-length(unique(sound.files[(sound.files %in% fs)])), 
              ".wav file(s) not found"))
    
    # Step I: Count the number of sound files in working directory
    d <- which(sound.files %in% fs) 
    if(length(d) == 0){
        stop("The .wav files are not in the working directory")
    }  else {
        start <- start[d]
        end <- end[d]
        selec <- selec[d]
        sound.files <- sound.files[d]
    }

    # Step II: 
    
    # A) If parallel is not numeric or if it is not positive
    if(!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
    if(any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
    # B) If parallel greater than 1
    if(parallel > 1) { 
        options(warn = -1)
        if(all(Sys.info()[1] == "Windows",requireNamespace("parallelsugar", quietly = TRUE) == TRUE)) 
        lapp <- function(X, FUN) parallelsugar::mclapply(X, FUN, mc.cores = parallel) else
            if(Sys.info()[1] == "Windows"){ 
            cat("Windows users need to install the 'parallelsugar' package for parallel computing.")
            lapp <- pbapply::pblapply} else lapp <- function(X, FUN) parallel::mclapply(X, FUN, mc.cores = parallel)} else lapp <- pbapply::pblapply
    
    options(warn = 0)

    # C) Parallel = 1 -> Measure the acoustic parameters
    if(parallel == 1) cat("Measuring acoustic parameters:")

    # Step III: Measure the acoustic parameters
    measures <- as.data.frame(lapp(1:length(start), function(i) {
        r <- tuneR::readWave(file.path(getwd(), sound.files[i]), from = start[i], to = end[i], units = "seconds")
        b <- bp
        if(b[2] > ceiling(r@samp.rate/2000) - 1) b[2] <- ceiling(r@samp.rate/2000) - 1

        #Frequency Spectrum Analysis
        songspec <- seewave::spec(r, f = r@samp.rate, plot = FALSE)
        analysis <- seewave::specprop(songspec, f = r@samp.rate, flim - c(0, 280/1000), plot = FALSE)

        # Step IV: Save parameters
        meanfreq <- analysis$mean/1000
        sd <- analysis$sd/1000
        median <- analysis$median/1000
        Q25 <- analysis$Q25/1000
        Q75 <- analysis$Q75/1000
        IQR <- analysis$IQR/1000
        skew <- analysis$skewness
        kurt <- analysis$kurtosis
        sp.ent <- analysis$sh
        sfm <- analysis$sfm
        mode <- analysis$mode/1000
        centroid <- analysis$cent/1000

        # Step V: More Frequency Parameters Calculation
        
        # A) Frequency with amplitude peaks
        peakf <- 0#seewave::fpeaks(songspec, f = r@samp.rate, wl = wl, nmax = 3, plot = FALSE)[1, 1]
        
        # B) Fundamental frequency parameters
        ff <- seewave::fund(r, f = r@samp.rate, ovlp = 50, threshold = threshold, 
                            fmax = 280, ylim=c(0, 280/1000), plot = FALSE, wl = wl)[, 2]
        meanfun<-mean(ff, na.rm = T)
        minfun<-min(ff, na.rm = T)
        maxfun<-max(ff, na.rm = T)
        
        # C) Dominant frecuency parameters
        y <- seewave::dfreq(r, f = r@samp.rate, wl = wl, ylim=c(0, 280/1000), ovlp = 0, plot = F, threshold = threshold, bandpass = b * 1000, fftw = TRUE)[, 2]
        meandom <- mean(y, na.rm = TRUE)
        mindom <- min(y, na.rm = TRUE)
        maxdom <- max(y, na.rm = TRUE)
        dfrange <- (maxdom - mindom)
        duration <- (end[i] - start[i])
        
        # D) Modulation Index Calculation
        changes <- vector()
        for(j in which(!is.na(y))){
        change <- abs(y[j] - y[j + 1])
        changes <- append(changes, change)
        }
        if(mindom==maxdom) modindx<-0 else modindx <- mean(changes, na.rm = T)/dfrange
        
        # E) Save Results
        return(c(duration, meanfreq, sd, median, Q25, Q75, IQR, skew, kurt, sp.ent, sfm, mode, 
                centroid, peakf, meanfun, minfun, maxfun, meandom, mindom, maxdom, dfrange, modindx))
    }))

    # Step VI: Change result names and return the dataframe
    rownames(x) <- c("duration", "meanfreq", "sd", "median", "Q25", "Q75", "IQR", "skew", "kurt", "sp.ent", 
                   "sfm","mode", "centroid", "peakf", "meanfun", "minfun", "maxfun", "meandom", "mindom", "maxdom", "dfrange", "modindx")
    x <- data.frame(sound.files, selec, as.data.frame(t(x)))
    colnames(x)[1:2] <- c("sound.files", "selec")
    rownames(x) <- c(1:nrow(x))
    
    return(x)
}

setwd("C:/Users/Paridhi/Dropbox/PC/Documents/Projects/Gender Recognition using Voice and Speech Analysis/Data/")

data = data.frame("brian.wav", 0, 0, 20)
names(data) <- c('sound.files', 'selec', 'start', 'end')

acoustics <- specan3(data, parallel=1)
write.csv(acoustics, file = "Brian-Acoustics.csv")
