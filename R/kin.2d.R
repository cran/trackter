# halfwave
#' @title Compute half wavelengths from a sine-like waveform
#' @description Computes half wavelengths and their positions and amplitude from a sine-like waveform based on either peak-to-trough or internodal distance.
#'

#' @param x Numeric; x position
#' @param y numeric; y position
#' @param method character; how half waves should be found and classified, where it crosses zero/the internodal length ("zeros") or peak to trough/trough to peak ("p2t"). See Details. 
#' @param zero.begin logical; does wave begin at zero? Default is 'TRUE' and will help find waves beginning at first x,y values if y=0
#' @param fit logical; if 'method="zeros"', should zeros be detected by a fitting operation. See Details.
#' @param smoothing character; the smoothing method when 'fit=TRUE', either 'loess' or 'spline'. See Details.
#' @param dens numeric; factor by which to increase the sample density used in fitting when 'method="zeros"'. See Details.
#'@param smooth numeric; if \code{smoothing} is set to 'loess', 'span' parameter value for \code{\link{loess}}. If \code{smoothing} is set to 'spline' 'spar' parameter value for \code{\link{smooth.spline}}
#'
#' @details If 'method="p2t"', half waves are found using critical points (i.e., local maxima and minima) with \code{\link{features}}. Detected half waves with this method can be either peak to trough or trough to peak.
#' 
#' If 'method="zeros"' and 'fit=TRUE', zero crossings are determined by first increasing the sample density by a factor determined by \code{dens}. A more dense \code{\link{loess}} or \code{\link{smooth.spline}} model is then fit to the data and new y values predicted. Wave positions and lengths are determined based on these predicted values. This option should be useful when the sampling density of the waveform is relatively low and therefor detected wave positions and zero crossings (the internodes) may be rather coarse.
#' @return A list with the following components:
#' 
#' \code{method} the method chosen to find half waves
#' 
#'  \code{names} a data table with columns 'x', 'y', and 'wave' describing the x and y positions of the wave and a numeric name of each half wave detected, resptively.  If 'method="zeros"' and 'fit=TRUE', these values reflect the predicted, more dense data as determined by \code{smoothing},\code{smooth}, and \code{dens}.
#'  
#' \code{dat} a data table describing each have wave detected.
#' \itemize{
#' \item 'zeros': x value where y crosses zero. Returns NA if 'method=p2t'.
#' \item 'wave.begin': x value where each half wave begins.
#' \item 'wave.end': x value where each half wave ends.
#' \item 'begin.index': x index of where each half wave begins.
#' \item 'end.index': x index of where each half wave ends.
#' \item 'wave': numeric name of each wave.
#' \item 'l': the length of each half wave.
#' \item 'amp1':  If method is set to 'p2t' this is the begin amplitude. If "method='zeros'", this is the maximum absolute amplitude between internodes. 
#' \item 'amp2': If method is set to 'p2t', this is the end amplitude. If "method='zeros' value is NA.
#' \item 'pos1': If method is set to 'p2t', the x position of begin amplitude for each half wave and identical to 'begin'. If "method='zeros'", the position of maximum absolute amplitude between the internodes. 
#' \item 'pos2': If method is set to 'p2t', the x position of end amplitude for each half wave and identical to 'end'. If "method='zeros'", value is NA
#' }
#' 
#' If 'method="zeros"' and 'fit=TRUE', these values reflect the predicted, more dense data as determined by \code{smoothing}, \code{smooth}, and \code{dens}.
#' 
#' @import features
#' @import data.table
#' @seealso \code{\link{features}}, \code{\link{loess}}, \code{\link{smooth.spline}}
#' 
#' @examples
#' 
#' require(ggplot2)
#'
#' #Find length of the half waves
#' x <- seq(0,pi,0.01)
#' y <- sin(x^2*pi)
#' qplot(x,y)
#' 
#' #zero method predicting zeros
#' w.z <- halfwave(x,y,method="zeros",fit=TRUE,smoothing="spline")
#'
#' #plot waveform with detected half waves using fitted 'zeros' method
#' p <- ggplot()+geom_point(aes(x=x,y=y))
#' p <- p+geom_line(data=w.z$names,aes(x=x,y=y,col=wave),alpha=0.4,size=3,inherit.aes=FALSE)
#' p+theme_classic()
#' 
#' #plot lambda as it varies with position
#' qplot(data=w.z$dat,x=pos1,y=l)
#' 
#'#peak-to-trough method
#' w.p <- halfwave(x,y,method="p2t")
#'qplot(data=w.p$names,x=x,y=y,col=wave)
#'
#'
#' @export
halfwave <-function(x,y,method = "zeros", zero.begin=TRUE,fit=TRUE,dens=10,smooth=0.1,smoothing="loess") {
  
  wave.begin <- zeros <- l <- begin.index <- end.index <- NULL #to avoid NSE notes in R CMD check
  
  if (!method %in% c("p2t", "zeros"))  stop("method must be set to 'p2t' or 'zeros' (the default)")
  x = c(unlist(x))
  y = c(unlist(y))
  
  added <- FALSE
  #is first value y=0?
  if(y[1]==0){x.0 <- x[1]-diff(x[1:2])
  y.0 <- y[1]-diff(y[1:2])
  x <- c(x.0,x)
  y <- c(y.0,y)
  added <- TRUE
  }
  
  dt <- data.table(x, y)
  
  
  #find peaks and then lengths
  
  
  if (method == "p2t") {
    
    ft <- features(x, y)
    x.pos <- sapply(ft$cpts, function(z) which.min(abs(x - z)))
    y.pk <- y[x.pos]
    x.pk <- x[x.pos]
    
    wave.dat <- data.table(zeros = NA, wave.begin = x.pk)
    
    if(length(x.pk)<=1) stop("method='p2t' and data may have a peak or trough but not both. Try changing method to 'zeros'")
    
    wave.dat[, c("wave.end", "begin.index", "end.index", "wave", "l") := 
               list(dplyr::lead(wave.begin),
                    x.pos,
                    dplyr::lead(x.pos),
                    as.character(1:.N),
                    dplyr::lead(wave.begin) - wave.begin )  ]
  }
  if (method == "zeros") {
    if(fit==TRUE){
      x.pred <- seq(first(x),last(x),length.out = dens*length(x))
      
      if(!smoothing %in% c("spline","loess")) stop("'smoothing' should be set to 'spline' or 'loess'")
      if(smoothing=="loess") y.fit <- predict(loess(y~x,span = smooth),newdata=data.frame(x=x.pred))
      
      if(smoothing=="spline"){ y.sp <- smooth.spline(x,y,spar = smooth)
      y.fit <- predict(y.sp,x=x.pred)$y
      }
      x <- x.pred
      y <- y.fit
      dt <- data.table(x=x.pred,y=y.fit)
    }
    y.signs <- sign(y)
    if(zero.begin) y.signs[y.signs==0] <- -1
    
    z <- which(abs(diff(y.signs))>1)
    
    if(length(z)<=1) stop("method='zeros' and 1 or fewer zero crossings found. Try changing method to 'p2p' or 't2t'")
    names(x) <- NULL
    
    wave.dat <- data.table(zeros = z+1, wave.begin = x[z+1])
    
    
    if (length(z) > 1) {
      wave.dat[, c("wave.end", "begin.index", "end.index", "wave", "l") := list(dplyr::lead(wave.begin)-1,zeros,dplyr::lead(zeros)-1,as.character(1:.N),dplyr::lead(wave.begin) - wave.begin )  ]
    }
  }
  wave.dat <-  wave.dat[!is.na(l)]
  
  if (nrow(wave.dat) != 0) {
    if (all(is.na(wave.dat$l))) {
      wave.dat2 <- dt[, wave := NA]
    } else{
      wave.dat2 <- wave.dat[!is.na(l), list(x = x[begin.index:end.index]), by = list(wave)]
      wave.dat2 <- merge(dt, wave.dat2, all.x = TRUE,by="x")
    }
    
    pks <-wave.dat2[!is.na(wave), list(amp = y[which.max(abs(y))], pos = as.numeric(x[which.max(abs(y))])), by = wave]
    
    if(method=="p2t") pks <-wave.dat2[!is.na(wave), list(amp1 = first(y), amp2=last(y),pos1=first(x),pos2=last(x)), by = wave]
    
    if(method=="zeros") pks <-wave.dat2[!is.na(wave), list(amp1 = max(abs(y)), amp2=NA,pos1=x[which.max(abs(y))],pos2=NA), by = wave]
    
    
    
    if ( nrow(pks)>0){
      wave.dat <- merge(wave.dat, pks, by = "wave")
    }
  } else{
    wave.dat <- data.table(
      zeros = 0,
      wave.begin = 0,
      wave.end = 0,
      begin.index = 0,
      end.index = 0,
      wave = 0,
      l = as.numeric(0),
      pos = 0,
      amp = 0,
      l2 = 0,
      pos2 = 0
    )
    wave.dat2 <- dt[, wave := NA]
  }
  
  #zap dummy data to find wave begining at y=0
  
  if(added){
    wave.dat2 <- wave.dat2[-1,]
    wave.dat2$y <- wave.dat2$y-1e-17
    wave.dat[,zeros:=zeros-1]
    wave.dat[,begin.index:=begin.index-1]
    wave.dat[,end.index:=end.index-1]
  }
  
  return(list(
    method=method,
    names = wave.dat2[!is.na(wave)],
    dat = wave.dat
  ))
}

# wave
#' @title Compute wavelengths from a sine-like waveform
#' @description Computes full wavelengths and their positions and amplitude from a sine-like waveform based on either peak-to-peak, trough-to-trough, or internodal distance.
#'
#' @param x numeric; x position
#' @param y numeric; y position
#' @param method character; how waves should be found and classified, where it crosses zero/the internodal length ("zeros"), peak to peak ("p2p") or trough to trough ("t2t"). See Details. 
#' @param zero.begin logical; does wave begin at zero? Default is 'TRUE' and will help find waves beginning at first x,y values if y=0
#' @param fit logical; if 'method="zeros"', should zeros be detected by a fitting operation. See Details.
#' @param smoothing character; the smoothing method when 'fit=TRUE', either 'loess' or 'spline'. See Details.
#' @param dens numeric; factor by which to increase the sample density used in fitting when 'method="zeros"'. See Details.
#'@param smooth numeric; if \code{smoothing} is set to 'loess', 'span' parameter value for \code{\link{loess}}. If \code{smoothing} is set to 'spline' 'spar' parameter value for \code{\link{smooth.spline}}
#'
#' @details If 'method="p2p"' or 'method="t2t"', full waves are found using critical points (i.e., local maxima, the peaks or minima, the troughs) with \code{\link{features}}.
#' 
#' If 'method="zeros"' and 'fit=TRUE', zero crossings are determined by first increasing the sample density by a factor determined by \code{dens}. A more dense \code{\link{loess}} or \code{\link{smooth.spline}} model is then fit to the data and new y values predicted. Wave positions and lengths are determined based on these predicted values. This option should be useful when the sampling density of the waveform is relatively low and therefor detected wave positions and zero crossings (the internodes) may be rather coarse.
#' 
#' @return A list with the following components:
#' 
#' \code{method} the method chosen to find full waves
#' 
#'  \code{names} a data table with columns 'x', 'y', and 'wave' describing the x and y positions of the wave and a numeric name of each wave detected, respectively. If  'method="zeros"' and 'fit=TRUE', these values reflect the predicted, more dense data as determined by \code{smoothing}, \code{smooth}, and \code{dens}.
#'  
#' \code{dat} a data table describing each wave detected.
#' \itemize{
#' \item 'zeros': x value where y crosses zero. Returns NA if \code{method} is 'p2p' or 't2t', value is NA.
#' \item 'wave.begin': x value where each wave begins.
#' \item 'wave.end': x value where each wave ends.
#' \item 'begin.index': x index of where each wave begins.
#' \item 'end.index': x index of where each  wave ends.
#' \item 'wave': numeric name of each wave.
#' \item 'l': the length of each  wave.
#' \item 'amp1': the peak amplitude of each wave. If method is set to 'p2p' or 't2t' this is the begin amplitude. If "method='zeros'" this is the peak amplitude between internodes. 
#' \item 'amp2': If method is set to 'p2p' or 't2t' this is the end amplitude. If "method='zeros'" this is the minimum amplitude between internodes. 
#' \item 'pos1': If method is set to 'p2p' or 't2t' the x position of begin amplitude for each half wave and identical to 'begin'. If "method='zeros'", the position of peak amplitude between the internodes. 
#' \item 'pos2': If method is set to 'p2p' or 't2t' the x position of end amplitude for each half wave and identical to 'end'. If "method='zeros'", the position of minimum amplitude between the internodes.
#' }
#' 
#' If  'method="zeros"' and 'fit=TRUE', these values reflect the predicted, more dense data as determined by \code{smoothing},\code{smooth}, and \code{dens}.
#' 
#' 
#' @import features
#' @import data.table
#' 
#' 
#' @seealso \code{\link{features}}, \code{\link{loess}}, \code{\link{smooth.spline}}
#' 
#' @examples
#'
#' require(ggplot2)
#' #Find length of the full waves
#' x <- seq(0,pi,0.01)
#' y <- sin(x^2*pi)
#' 
#' #zero method
#' w.z <- wave(x,y,method="zeros",smoothing="spline",smooth=0.1)
#' 
#' #plot wave with detected full waves using fitted 'zeros' method
#'p <- ggplot()+geom_point(aes(x=x,y=y))
#'p <- p+geom_line(data=w.z$names,aes(x=x,y=y,col=wave),alpha=0.4,size=3,inherit.aes=FALSE)
#'p+theme_classic()
#'
#'#plot lambda as it varies with position
#'
#'qplot(data=w.z$dat,x=pos1,y=l)
#'
#'#trough-to-trough method
#' w.p <- wave(x,y,method="t2t")
#' 
#'qplot(data=w.p$names,x=x,y=y,col=wave)
#'
#' @export
wave <-function(x,y,method = "zeros", zero.begin=TRUE,fit=TRUE,dens=10,smooth=0.1,smoothing="loess") {
  
  wave.begin <- zeros <- l <- begin.index <- end.index <- NULL# to avoid NSE erros on R CMD check
  
  if (!method %in% c("p2p", "zeros","t2t"))  stop("method must be set to 'p2p' , 't2t', or 'zeros' (the default)")
  x = c(unlist(x))
  y = c(unlist(y))
  
  added <- FALSE
  #is first value y=0?
  if(y[1]==0){x.0 <- x[1]-diff(x[1:2])
  y.0 <- y[1]-diff(y[1:2])
  x <- c(x.0,x)
  y <- c(y.0,y)
  added <- TRUE
  }
  
  dt <- data.table(x, y)
  
  #find peaks and then lengths
  
  
  if (method == "p2p" | method=="t2t") {
    
    ft <- features(x, y)
    x.pos <- sapply(ft$cpts, function(z) which.min(abs(x - z)))
    
    curve <- ft$curvature
    pt <- sapply(curve,function(x) ifelse(x>=0,"t","p"))
    
    if(method=="ptp") x.pos <- x.pos[pt=="p"]
    if(method=="t2t") x.pos <- x.pos[pt=="t"]
    y.pk <- y[x.pos]
    x.pk <- x[x.pos]
    
    wave.dat <- data.table(zeros = NA, wave.begin = x.pk)
    
    if(length(x.pk)<=1) stop("method='p2p' or 't2t' and data fewer than two peaks or troughs. Try changing method to 'zeros'")
    
    wave.dat[, c("wave.end", "begin.index", "end.index", "wave", "l") := 
               list(dplyr::lead(wave.begin),
                    x.pos,
                    dplyr::lead(x.pos),
                    as.character(1:.N),
                    dplyr::lead(wave.begin) - wave.begin )  ]
  }
  
  if (method == "zeros") {
    if(fit==TRUE){
      x.pred <- seq(first(x),last(x),length.out = dens*length(x))
      
      if(!smoothing %in% c("spline","loess")) stop("'smoothing' should be set to 'spline' or 'loess'")
      if(smoothing=="loess") y.fit <- predict(loess(y~x,span = smooth),newdata=data.frame(x=x.pred))
      if(smoothing=="spline"){ y.sp <- smooth.spline(x,y,spar = smooth)
      y.fit <- predict(y.sp,x=x.pred)$y
      }
      x <- x.pred
      y <- y.fit
      dt <- data.table(x=x.pred,y=y.fit)
    }
    y.signs <- sign(y)
    if(zero.begin) y.signs[y.signs==0] <- -1
    
    z <- which(abs(diff(y.signs))>1)
    z <- z[seq(1,length(z),2)] #every other node
    
    if(length(z)<=1) stop("method='zeros' and 1 or fewer zero crossings found. Try changing method to 'p2p' or 't2t'")
    names(x) <- NULL
    
    wave.dat <- data.table(zeros = z+1, wave.begin = x[z+1])
    
    
    if (length(z) > 1) {
      wave.dat[, c("wave.end", "begin.index", "end.index", "wave", "l") := list(dplyr::lead(wave.begin)-1,zeros,dplyr::lead(zeros)-1,as.character(1:.N),dplyr::lead(wave.begin) - wave.begin )  ]
    }
  }
  
  wave.dat <-  wave.dat[!is.na(l)]
  
  if (nrow(wave.dat) != 0) {
    if (all(is.na(wave.dat$l))) {
      wave.dat2 <- dt[, wave := NA]
    } else{
      wave.dat2 <- wave.dat[!is.na(l), list(x = x[begin.index:end.index]), by = list(wave)]
      wave.dat2 <- merge(dt, wave.dat2, all.x = TRUE,by="x")
      
    }
    
    
    if(method=="p2p" | method=="t2t") pks <-wave.dat2[!is.na(wave), list(amp1 = first(y), amp2=last(y),pos1=first(x),pos2=last(x)), by = wave]
    
    if(method=="zeros") pks <-wave.dat2[!is.na(wave), list(amp1 = max(y), amp2=min(y),pos1=x[which.max(y)],pos2=x[which.min(y)]), by = wave]
    
    
    
    if ( nrow(pks)>0){
      wave.dat <- merge(wave.dat, pks, by = "wave")
    }
  } else{
    wave.dat <- data.table(
      zeros = 0,
      wave.begin = 0,
      wave.end = 0,
      begin.index = 0,
      end.index = 0,
      wave = 0,
      l = as.numeric(0),
      pos = 0,
      amp = 0,
      l2 = 0,
      pos2 = 0
    )
    wave.dat2 <- dt[, wave := NA]
  }
  
  #zap dummy data to find wave begining at y=0
  
  if(added){
    wave.dat2 <- wave.dat2[-1,]
    wave.dat2$y <- wave.dat2$y-1e-17
    wave.dat[,zeros:=zeros-1]
    wave.dat[,begin.index:=begin.index-1]
    wave.dat[,end.index:=end.index-1]
  }
  
  return(list(
    method=method,
    names = wave.dat2[!is.na(wave)],
    dat = wave.dat
  ))
}

#amp.freq
#' @title Computes amplitude and frequency of wave-like data
#' @description Computes amplitude(s) and wavelength(s) of a wave form, amongst other things, based on a sampling frequency
#'
#' @param x Numeric; x position (or sample number)
#' @param y numeric; y position
#' @param sf numeric; sample frequency (i.e., how often was x and y sampled) in Hz
#' @return a list with amplitude "a", frequency "f", amplitude returned from a smoothed sign function "a.f" based on output from \code{features}, signal to noise ratio "snr".
#' @export
#' @import features
#' @seealso \code{\link{features}}
#' @examples
#'\dontrun{
#' #Compute waveform patterns
#' x <- seq(0,pi,0.1)
#' y <- sin(x^1.3*pi)
#' plot(x,y)
#'
#' amp.freq(x=x,y=y)
#' }

amp.freq <- function(x = NULL, y, sf = 100) {
  s <- 1 / sf
  if (is.null(x))
    x <- 1:length(y)
  amp.n <- features(x = x, y = y)
  x.n <- unlist(sapply(amp.n$cpts, function(z)
    which.min(abs(z - x))))
  amp <- abs(diff(y[x.n]) / 2)
  amp.f <- abs(diff(attributes(amp.n)$fits$fn[x.n]) / 2)
  snr <- fget(amp.n)$f["snr"]
  names(snr) <- NULL
  freq <-
    1 / (s * diff(amp.n$cpts[seq(1, length(amp.n$cpts), 2)])) #peak to peak or trough to trough freq
  tail.dat <-
    list(a = amp,
         f = freq,
         a.f = amp.f,
         snr = snr)#a.f is amp according to function
  
  return(tail.dat)
}

