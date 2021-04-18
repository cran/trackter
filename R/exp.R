#' @title  Midline tracking over image sequences with ROI search using LDA

#' @description  Experimental and untested (in the unit-testing sense). Automatically retrieves the midline of a detected ROI in each image of a sequence through thresholding and segmentation. Chose a fish-like ROI class detected through linear discriminate  analysis (LDA) of PCA on elliptical Fourier described shapes. Initial training of ROIs is user defined or with the 'fishshapes' data set loaded with \code{trackter} (see details). For each detected ROI, \code{kin.LDA} finds the y-value midpoint along the x-value array of the ROI and fits a midline according to a chosen smoothing method (loess or spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI. Supported image formats are jpeg, png, and tiff.
#'
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process.
#' @param thr numeric or character ('otsu') threshold to determine binary image. See Details.
#' @param ant.per numeric; left-most percentage of ROI that establishes the horizontal reference for the midline displacement.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position. 
#' @param edges logical, should ROIs on image edges by evaluated. See Details.
#' @param size.min numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame. Default is 0.05.
#' @param enorm logical, should the EFA coefficients from \code{efourier} operations be normalized or not. See \code{details} and \code{\link{efourier}}
#' @param harms numeric, the number of harmonics to use. If missing, \code{Momocs} sets 'nh.b' to 12. Will produce messages.
#'
#'@param rescale logical, should all shapes in PCA be rescaled. Performs best as 'FALSE'.
#'@param train.dat Classified \code{Out} and \code{Coo} outlines that are produced from \code{Momocs}. See Details. 
#'@param retrain numeric, the number of frames on which to retrain the LDA data set. See details.
#'@param after.train character, if set to 'size', LDA will be skipped after \code{retrain} and the ROI with a size closest to the ROI found by the LDA $>=$ will be chosen. This speeds calculations considerably. If 'LDA', the default, LDA will continue using the retraining classifications from frames $<=$ 'train'.
#'
#'@param ties character, how to chose ROI's in any one frame that appear fish-like. See details.
#' @param smoothing character, the midline smoothing method, either 'loess' or 'spline'.
#' @param smooth numeric; if \code{smoothing} is set to 'loess', 'span' parameter value for \code{\link{loess}}. If \code{smoothing} is set to 'spline' 'spar' parameter value for \code{\link{smooth.spline}}
#' @param smooth.points numeric, number of equally spaced points along the ROI midline on which the smoothed midline is computed.
#' @param show.prog logical value indicating if outputted image should be displayed during analysis.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{ant.per} \code{lm()} overlaying original or binary images.
#' @param out.dir character, the directory to which ouputted images should be saved. If NULL, then a subdirectory 'processed_images' in the working directory.
#' @param image.type character; the type of image to be outputted, either 'orig' or 'bin' representing the original or binary images, respectively. Ignored if 'save==FALSE'.
#' @param plot.pml logical, value indicating if outputted images should include the predicted midline (in blue) and the points according to \code{ant.per} used to construct the predicted midline (in green).
#' @param flip logical, indicating if binary should be flipped.
#' 
#' @export
#'
#' @details
#'The algorithm assumes a left-right orientation, i.e., the head of the ROI is positioned left, the tail right. ffmpeg operations or even imageJ can rotate images not in this orientation. The \code{ant.per} value therefor establishes the reference line (theoretical straight midline) based on that portion of the head. The midline is calculated as the midpoints between the y extrema for each x position.  
#'If 'save=TRUE', images are saved as binary or the original with a body midline overlay and, if chosen, with the theoretical midline (based on \code{ant.per}). 
#'
#' Thresholding operations can be performed with an arbitrary (user-defined) numeric value or with Otsu's method ('thr="otsu"'). The latter chooses a threshold value by minimizing the combined intra-class variance. See \code{\link{otsu}}.
#'
#'Before \code{train}, ROIs are chosen according to LDA of a PCA object constructed from \code{efourier} analysis. LDA is trained by a user define 'train.dat' when the frame $<=$ \code{retrain}. LDA will proceed after \code{retrain} if \code{after.train}='LDA', but the LDA will be trained by the contours classified as 'fish' and 'not.fish' found during the chosen training period. 
#'
#'\code{enorm} Normalization of EFA coefficients is often perilous, especially for symmetrical shapes, a conditional met for undulating, bilaterally symmetrical organisms at least some of the time and perhaps for many of the frames included in any analysis. Thus, 'enorm' by default is set to 'FALSE'. 'enorm=TRUE' may produce odd ROI choices and should be used cautiously.
#'
#'\code{train.dat} This should be a \code{Coo} and \code{Out} object produced by \code{efourier} analysis of predefined shapes. A user defined dataset or the \code{fishshapes} dataset in \code{trackter} must be used for training. \code{fishshapes} includes several arbitrary shapes (circles, squares, U-shapes, etc.) as well as several fish shapes: sunfish (genus Lepomis), eel (genus Anguilla), and trout (genus Onchorhynchus) swimming over one tail-beat cycle. A user-defined dataset must have shapes classified with factors identical to the \code{fishshapes} contours, that is by shape, type, and edge. Shape levels should indicate what type of shape is described by the contour (e.g., 'circle', 'L-shape', 'trout', 'eel', etc). The type levels must describe the shape as 'fish' or 'not.fish'. The edge levels must be 'FALSE'. 
#'
#'\code{edges} Set by default to 'FALSE'. It is not advisable to include shapes that are on the edge of any frame and are therefore incomplete.
#'\code{retrain} After this value, the LDA analysis will use the ROIs determined as 'fish' and 'not.fish' in the frames $>=$ \code{retrain} to discriminate fish from non-fish shapes. This speeds up analysis considerably.
#'\code{ties} Determines  how to chose ROIs if more than one fish-like ROI is found in any frame. 'fish' will result in choosing the ROI with shape types in which the best *and* second-best fish-like shape (according to posterior probabilities) match a fish-like shape in the training and/or retraining datasets.'post' will chose the best fish-like shape according the the highest posterior probability from LDA.
#'
#' @return A list with the following components:
#'
#' \code{kin.dat} a data table consisting of frame-by-frame position parameters for the ROI determined by LDA analysis.
#' \itemize{
#' \item the frame number
#'
#' \item 'x' and ''y': the position of the tail (rightmost or posteriormost)
#' \item 'head.x' and 'head.y': the x and y position of the head (leftmost or anteriormost)
#' \item 'amp': the amplitude (\code{amp}) of the tail relative to thr theoretical midline determined by the \code{lm()} predictions from \code{ant.per}
#' \item 'roi': a character indicating the ROI ranked by size ('a' being the largest)
#' \item 'head.pval': p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)}
#'
#' \code{midline} A data table containing, for each frame described by \code{frames}, the following: \itemize{
#' \item 'x' and 'y.m': x and y positions of the midline of the ROI
#' #' \item 'y.min' and 'y.max': min and max y positions ROI's contour used in y.m calculation
#' \item 'mid.pred': the predicted linear midline based on the points/pixels defined by \code{head.per} (green points in the outputted images/video)
#' \item 'y.pred': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
#' \item 'wave.y': midline points 'y.pred' relative to 'mid.pred'
#' \item 'roi': a character indicating ROI size ('a' being the largest)
#' }
#' 
#' \code{cont} A data table containing x and y positions of the contours used to calculate the data in 'kin.dat'. Contains the following: 
#' \itemize{
#' \item 'frame': the frame
#' #' \item 'x' and 'y': the x and y positions of the contours
#' }
#' 
#' \code{all.classes} A data table containing the following for all ROIs detected: 
#' \itemize{
#' \item 'frame': the frame
#' \item 'roi': the name of each ROI found in a frame.
#' \item 'size': the size of each ROI
#' }
#' 
#' \code{dim} the x and y dimensions of the images analyzed
#' 
#' @seealso \code{\link{kin.simple}}, \code{\link{kin.search}}, \code{\link{efourier}} \code{\link{LDA}}, \code{\link{fishshapes}}.
#' 
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess  predict smooth.spline
#' @importFrom utils head setTxtProgressBar tail txtProgressBar
#' @importFrom EBImage ocontour otsu bwlabel
#' @importFrom zoo index
#'
#' 
#' @examples
#' # produce a classic midline waveform plot of swimming 
#' # fish searching a image field with a two fish-like ROIs
#' \dontrun{
#' require(wesanderson)
#' require(ggplot2)
#' require(data.table)
#' require(dplyr)
#' require(EBImage)
#'
#' #download example images and place in 'example' subdirectory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/example.zip?raw=true"
#' download.file(f, "temp.zip")
#' unzip("temp.zip")
#' unlink("temp.zip")
#'
#' #load fishshapes data
#' data(fishshapes)
#' 
#' 
#'kin <- kin.LDA(image.dir = "example",frames=1:20,thr=0.7,
#'               ant.per=.25,enorm=FALSE,show.prog = FALSE,retrain=2,
#'               train.dat = fishshapes,after.train="LDA",edges=FALSE, 
#'               )
#' ml <- kin$midline
#'  #x start at 0
#' ml <-ml[,x2:=x-x[1],by=frame]
#'
#' #compute instantaneous amplitude of tail (last/rightmost point) and wave crest x position  by frame
#' ml2 <-ml[,.(amp.i=abs(last(wave.y))),by=frame]
#'
#' ml <- merge(ml,ml2,by="frame") #merge these
#'
#' pal <- wes_palette("Zissou1", 100, type = "continuous") #"Zissou" color palette
#' p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)+scale_color_gradientn(colours = pal)
#' p <- p+geom_line(aes(group=frame,color=amp.i),
#' stat="smooth",method = "loess", size = 1.5,alpha = 0.5)
#' print(p)
#'
#'
#' ### Make a video of processed frames
#'
#' images.to.video2(image.dir="processed_images",
#' vid.name="trout_test",frame.rate=5,qual=100,raw=FALSE)
#' file.exists("trout_test_red.mp4")
#'
#' }
#' 

kin.LDA <-function(image.dir=NULL,frames=NULL,thr=0.7,ant.per=0.20,tips=0.2,edges=FALSE,train.dat=NULL,rescale=FALSE,harms=15,enorm=TRUE,retrain=5,after.train="LDA",ties="fish",size.min=0.05,show.prog=FALSE,smoothing="loess",smooth=.3,smooth.points=200,save=TRUE,out.dir=NULL,image.type="orig",plot.pml=TRUE,flip=TRUE){

  type <- shape <- post <- type2 <- x <- y.pred <- wave.y <- mid.pred <- roi <- index <-  NULL#to avoid NSE errors in R CMD check
  
  if(is.null(train.dat)) stop("'train.dat' must be specified") ###load training data 
  
  if(save & is.null(out.dir)){
    unlink("processed_images",recursive = TRUE)
    dir.create("processed_images")
    proc.dir <- "processed_images"
  }else{
    if(save){
      proc.dir <- paste0(out.dir,"/processed_images")
      if(dir.exists(proc.dir))  unlink(proc.dir,recursive = TRUE)
      dir.create(proc.dir)
    }
  }
  
  images <- paste0(image.dir,"/",list.files(image.dir)[!grepl("Icon\r",list.files(image.dir))]) #remove pesky Icon\r
  
  if(any(frames>length(images))) stop("variable 'frames' out of range of image sequence")
  if(!is.null(frames)) images <- images[frames]
  
  trial <- gsub("\\.[^.]*$", "", basename(images[1]))
  
  kin.l <- list()
  midline.l<- list()
  classes.l <- list()
  
  lms <- list()
  conts <- list()
  pb = txtProgressBar(min = 0, max = length(images), initial = 0,style=3)
  
  roi.outs <- list() #store the rois for each image
  for(im in images){
    
    frame <- which(im==images)-1
    
    img <- EBImage::readImage(im,all=FALSE) #if don't add package, others use "display"
    
    img.dim <- dim(img)[1:2]
    
    # computes binary mask
    if(thr!="otsu" & !is.numeric(thr)) stop("'thr' must be set to 'otsu' or a numeric value=0-1")
    if(thr=="otsu"){EBImage::colorMode(img)=EBImage::Grayscale
    thr <- EBImage::otsu(img)[1]}
    
    y = img >thr #contrast threshold
    if(flip){#flip binary
      y[y==1] <- 5
      y[y==0] <- 1
      y[y==5] <- 0
    }
    z = EBImage::bwlabel(y)
    rois <- tabulate(z)
    pix <- dim(z[,,1])[1]*dim(z[,,1])[2]
    w <- dim(z[,,1])[1] #width of image
    h <- dim(z[,,1])[2] #height of image
    per <- rois/(w*h) #how big are rois compared to pixel field
    
    c.roi <-  which(per>=size.min) #candidate rois, filtered by size of of pixel field
    
    names(c.roi) <- as.factor(letters[order(rois[c.roi],decreasing = TRUE)])
    
    z.l <- list()
    out.l <- list()
    size.l <- list()
    
    if(frame==retrain){
      kin.train <- do.call(rbind,kin.l)
      out.train <- roi.outs
      size.train <- do.call(rbind,size.l)
      roi.train <-Momocs::combine(roi.outs)
      
    }
    
    if(is.null(after.train)) stop("'after.train' not set to 'size' or 'LDA'.")
    
    if(!after.train %in% c("LDA","size")) stop("after.train not set to 'size' or 'LDA'.")
    
    if(after.train=="size" & frame>retrain){
      sz.m <- mean(kin.train$size)
      size.diff <- abs(rois[c.roi]-sz.m)
      
      for(r in c.roi){
        r.name <- as.character(names(c.roi)[c.roi==r])
        z.r <- z
        z.r[z!=r] <- 0
        z.r[z==r] <- 1
        z.m <- z.r[,,1]
        z.m[1,1] <- 0 #this gets a 1 when
        z.l[[r.name]] <- z.m
        
        z.c <- ocontour(z.m)
        
        wall <- any(z.c[[1]][,1]>dim(z)[1]-2 | z.c[[1]][,1]<2  |z.c[[1]][,2]>dim(z)[2]-2 | z.c[[1]][,2]<2)
        
        
        r.out <- Out(ocontour(z.m))
        if(wall ) edge <- TRUE
        if(!wall) edge <- FALSE
        r.out$fac <- data.frame(shape=paste0("roi-",r.name),type=paste0("roi"),edge=edge)
        out.l[[r.name]] <- r.out
        rois[c.roi[r.name]]
      }
      
      roi.out2 <- Momocs::combine(out.l)
      
      classes <- data.table(roi=gsub("roi-","",roi.out2$fac$shape),type=NA,shape=NA,post=NA,type2=NA,shape2=NA,post2= NA,retrained=frame>retrain,edge=roi.out2$fac$edge,size.diff=size.diff)
      
      classes.l[[paste0(frame)]] <- data.table(frame=frame,classes)
      
      if(edges==FALSE) classes <- classes[edge==FALSE]
      
      z.best <- z.l[[classes[which.min(size.diff)]$roi]]
      r.name <-classes[which.min(size.diff)]$roi
      best.class <- classes[which.min(size.diff)]
    }
    else{
      for(r in c.roi){
        r.name <- as.character(names(c.roi)[c.roi==r])
        z.r <- z
        z.r[z!=r] <- 0
        z.r[z==r] <- 1
        z.m <- z.r[,,1]
        z.m[1,1] <- 0 #this gets a 1 when
        z.l[[r.name]] <- z.m
        
        z.c <- ocontour(z.m)
        
        wall <- any(z.c[[1]][,1]>dim(z)[1]-2 | z.c[[1]][,1]<2  |z.c[[1]][,2]>dim(z)[2]-2 | z.c[[1]][,2]<2)
        
        r.out <- Out(ocontour(z.m))
        if(wall ) edge <- TRUE
        if(!wall) edge <- FALSE
        r.out$fac <- data.frame(shape=paste0("roi-",r.name),type=paste0("roi"),edge=edge)
        out.l[[r.name]] <- r.out
        
        rois[c.roi[r.name]]
      }
      
      ##### efourier and LDA analysis ####
      #train.data from training
      roi.out2 <- Momocs::combine(out.l)
      
      if(frame>retrain){all.outs <- Momocs::combine(roi.train,roi.out2)}else{all.outs <- Momocs::combine(train.dat,roi.out2)}
      #rescale?
      if(rescale) all.outs <- coo_scale(coo_center(all.outs))
      # panel(all.outs)
      if(length(roi.out2)==1) all.outs <- Momocs::combine(all.outs,roi.out2) #must have two roi rows for rePCA to work
      
      shapes.out <- all.outs %>% Momocs::filter(type!="roi")
      
      roi.out <- filter(all.outs,type=="roi") #shape that is roi
      roi.f <- suppressMessages(efourier(roi.out,nb.h = harms,norm=enorm,start=FALSE) )#ef analysis
      
      shapes.f <- suppressMessages(efourier(shapes.out,nb.h = harms,norm=enorm,start = FALSE))
      
      shapes.p <- suppressWarnings(PCA(shapes.f))#PCA of non roi
      #why LDA wouldn't drop levels to prevent warnings is a stumper
      shapes.p$fac$shape <- droplevels(shapes.p$fac$shape)
      
      #if few shapes classified as same, keep a ton of variance
      if(grepl("-",shapes.f$fac$shape[1])){ suppressMessages(shape.l <- LDA(shapes.p,"shape",retain=0.9999999))}else{ 
        suppressMessages(shape.l <- LDA(shapes.p,"shape"))
      }#LDA of non roi
      
      #redo the same PCA and LDA  with roi
      roi.pca <- suppressMessages(rePCA(shapes.p, roi.f))
      
      roi.lda <- suppressMessages(reLDA(roi.pca,shape.l))
      
      roi.shape<- roi.lda$class
      roi.shape2 <- apply(roi.lda$posterior,1, function(x) names(x[order(x,decreasing = TRUE)[2]]))
      
      #retrieve roi type
      roi.type <- sapply(roi.shape,function(x) unique(filter(shapes.p$fac,shape==x)$type))
      
      roi.post <- apply(roi.lda$posterior,1,max)
      roi.post2 <- apply(roi.lda$posterior,1, function(x) x[order(x,decreasing = TRUE)[2]])
      roi.type2 <- sapply(roi.shape2,function(x) unique(filter(shapes.p$fac,shape==x)$type))
      
      classes <- data.table(roi=gsub("roi-","",roi.out2$fac$shape),type=roi.type,shape=roi.shape,post=roi.post,type2=roi.type2,shape2=roi.shape2,post2= roi.post2,retrained=frame>retrain,edge=roi.out2$fac$edge,size.diff=NA)
      
      classes.l[[paste0(frame)]] <- data.table(frame=frame,classes)
      
      if(edges==FALSE) classes <- classes[edge==FALSE]
      
      if(nrow(classes[type=="fish"])<1){stop("No shape type found matching 'fish' or 'edge==FALSE'")}
      
      if(nrow(classes[type=="fish"])>1){
        if(ties=="post"){
          
          z.best <- z.l[[classes[type=="fish"][which.max(post)]$roi]]
          r.name <-classes[type=="fish"][which.max(post)]$roi
          best.class <- classes[type=="fish"][which.max(post)]
          warning("More than one 'fish' found in LDA. Tie broken with highest post. prob.")}
        if(ties=="fish"){
          ###is there just one roi with type and type2 fish?
          if(nrow(classes[type=="fish" & type2=="fish"])==1){
            
            z.best <- z.l[[classes[type=="fish" & type2=="fish"]$roi]]
            r.name <-classes[type=="fish" & type2=="fish"]$roi
            best.class <- classes[type=="fish" & type2=="fish"]
          }
          
          ###are there two rois with type and type2 fish?
          if(nrow(classes[type=="fish" & type2=="fish"])>1){
            
            z.best <- z.l[[classes[type=="fish" & type2=="fish"][which.max(post)]$roi]]
            r.name <-classes[type=="fish" & type2=="fish"][which.max(post)]$roi
            best.class <- classes[type=="fish" & type2=="fish"][which.max(post)]
            warning("'ties=fish' and >1 ROIs have best and second best type as 'fish'; breaking tie with sum of post. prob.")}
          if(nrow(classes[type=="fish" & type2=="fish"])==0){
            classes.l[[paste0(frame)]] <- data.table(frame=frame,classes)
            z.best <- z.l[[classes[type=="fish"][which.max(post)]$roi]]
            r.name <-classes[type=="fish"][which.max(post)]$roi
            best.class <- classes[type=="fish"][which.max(post)]
            warning("'fish' chosen to break tie, but no second best shapes='fish'. Chose 'fish' with highest post. prob")
          }
        }
      }
      else{ 
        r.name <- classes$roi
        z.best <- z.l[[r.name]]
        best.class <- classes
      }
    }  
    if(show.prog) {
      EBImage::display(z.best,method = "raster")
    }
    
    best.cont <- data.table(ocontour(z.best)[[1]])
    colnames(best.cont) <- c("x","y")
    
    conts[[paste0(frame)]] <- data.table(frame=frame,best.cont)
    
    y.df <- best.cont[,list(y.min=min(y),y.max=max(y),y.m=mean(y)),by=list(x)]
    setkey(y.df,"x")
    
    ends <- ceiling(nrow(y.df)*tips)
    tip.y <- mean(tail(y.df$y.m[!is.na(y.df$y.m)],ends))#tip is mean y.m of last 30 pixels
    tip.x <- mean(tail(y.df$x[!is.na(y.df$y.m)],ends))#tip is mean y.m of last 30 pixels
    
    head.y <- mean(head(y.df$y.m[!is.na(y.df$y.m)],ends))#tip is mean y.m of first 30 pixels
    head.x <- mean(head(y.df$x[!is.na(y.df$y.m)],ends))#tip is mean y.m of first 30 pixels
    
    #n midline points
    if(is.null(smooth.points)) smooth.points <- nrow(y.df)
    midline <- y.df[seq(1,nrow(y.df),length.out = smooth.points),] #two hundred points on midline
    
    midline <- midline[complete.cases(midline)]
    midline <- data.table(frame,midline)
    
    
    ####which type of lines to be fitted, spline or loess
    if(!any(c("spline","loess")==smoothing)) stop("'smoothing' must = 'loess' or 'spline'")
    
    if(smoothing=="loess")  ml.pred <- fitted(loess(midline$y.m~midline$x,span=smooth,degree=1))
    if(smoothing=="spline") ml.pred <- smooth.spline(x = midline$x,y=midline$y.m,spar=smooth)$y
    
    midline[,y.pred:=ml.pred]#add smoothed predictions
    
    #head section
    head.dat <- midline[1:(ant.per*smooth.points),]
    head.lm <- lm(y.pred~x,head.dat)
    
    head.p <- summary(head.lm)$r.squared #how well does head lm fit
    
    midline$mid.pred <- predict(head.lm,newdata=midline)#add lm prediction to midline df
    
    midline <- midline[complete.cases(midline),]
    
    midline[,wave.y:=y.pred-mid.pred] #wave y based on midline y and straight head.lm pred points
    
    midline[,roi:=r.name]
    n.roi <- paste0(basename(im),"-",r)
    
    kin.l[[paste(frame)]] <- data.table(frame,x=tip.x,y=tip.y,head.x,head.y,amp=last(midline$wave.y),head.pval=head.p,size=rois[c.roi[r.name]],best.class)
    midline.l[[paste(frame)]] <- midline
    
    if(frame<=retrain){
      roi.out2$fac <- data.frame(shape=roi.shape,type=roi.type,edge=roi.out2$fac$edge)
      #are all shapes the same, if so give them unique levels
      if(all(roi.shape==roi.shape[1]))   roi.out2$fac <-  data.frame(shape=paste0(roi.shape,"-",index(roi.shape)),type=roi.type,edge=roi.out2$fac$edge)
      
      roi.outs[[paste(frame)]] <- roi.out2
    } #save for retraining
    
    
    if(save){
      
      jpeg(paste0(proc.dir,"/",trial,"_",sprintf("%03d",frame),".jpg"),quality = 0.5)
      if(image.type=="bin")EBImage::display(z,method = "raster")
      if(image.type=="orig")EBImage:: display(img,method = "raster")
      
      
      if(plot.pml) lines(predict(lm(mid.pred~x,midline)),x=midline$x,col="blue",lwd=4)
      with(midline,lines(y.pred~x,col="red",lwd=4))
      if(plot.pml) with(midline[1:ceiling(ant.per*smooth.points),],points(x,y.pred,col="green",pch=16,cex=0.75))
      
      dev.off()
    }
    
    
    setTxtProgressBar(pb,which(images==im))
  }
  
  classes.dat <- do.call(rbind,classes.l)
  kin.dat <- do.call(rbind,kin.l)
  midline.dat <- data.table(do.call(rbind,midline.l))
  cont.dat <- do.call(rbind,conts)
  
  return(list(kin.dat=kin.dat,midline=midline.dat,cont=cont.dat,all.classes=classes.dat,dim=img.dim))
  
 
}