#' @title Extracts images from a video file with ffmpeg
#'
#' @description Uses ffmpeg systems calls to extract images from a video.
#'
#' @param vid.path Character; path of video file to be processed.
#' @param out.dir character; directory path in which to store images.
#' @param overwrite logical; should path described by 'out.dir' be overwritten if it exhists. 
#' @param qual numeric; the quality of the jpeg images to be rendered from 1-100\%. Defaults to 50\%.
#' @return Extracts all the images of the video and saves them to an "images" directory with appended number sequence
#' @seealso \code{\link{images.to.video}}
#' @export
#' @examples
#'
#' #make a video with animation package
#' \donttest{
#' require(animation)
#' fun <- function(){
#' y <- sin(1:50)
#' x <- 1:50
#' for(i in 1:50) {
#'   plot(x[i],y[i],col="red",xlim=c(0,50),ylim=range(y))
#'   animation::ani.pause()
#'   }
#' }
#' animation::saveVideo(fun(),video.name=paste0(tempdir(),"/wave.mp4"),interval = 0.2)
#'
#' #create directory in which to store images
#' dir.create(paste0(tempdir(),"/images"))
#' vid.to.images(vid.path=paste0(tempdir(),"/wave.mp4"),
#' out.dir= paste0(tempdir(),"/images"),qual=100)
#'
#' #see the images in the "images" subdirectory
#' list.files( paste0(tempdir(),"/images"))
#' 
#' #clean up
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#' }


vid.to.images <- function(vid.path=NULL,out.dir=NULL,overwrite=FALSE,qual=50)  {
  
  version <-  try(system("ffmpeg -version", intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command 'ffmpeg' is not available in your system. Please install FFmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}
   
  qual <- round(30-30*(qual/100)+1,0)

  out.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",out.dir) #if image.dir has spaces
  
  #delete or create the out directory
  if(is.null(out.dir)|length(out.dir)==0) stop("'out.dir' not specified.")
  if(is.null(vid.path)) stop("'vid.path' not specified.")

  if(file.exists(out.dir) & overwrite==TRUE) unlink(out.dir,recursive = T)
  
  if(!file.exists(out.dir)) stop("Directory specified by 'out.dir' (", paste0(out.dir),") does not exist")
 
  out.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",out.dir) #if image.dir has spaces
  
  if(!file.exists(vid.path)) stop("Path specified by 'vid.path' (", paste0(vid.path),") does not exist")
  
  
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path)
  
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",out.dir)
  
  video.name<- gsub(".avi","",basename(vid.path))
  

  system(paste0("ffmpeg -i ", vid.path, " -q:v ",qual," ", image.dir,"/",video.name,"_%5d.jpg")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
  
}

#' @title Stitches images into a video file
#' @description Stitches images into a video file of type indicated by "vid.ext"
#'
#' @param image.dir character; directory containing images to stitch.
#' @param out.dir character; directory in which to store video.
#' @param vid.name character; file name given to video including extension.  mp4 currently works best.
#' @param qual numeric; the quality of the video rendered from 1-100\%. Defaults to 50\%.
#' @param frame.rate numeric; video frame rate in fps.
#' @param overwrite logical; should path described by vid.name  be overwritten if it exhists. 
#' @param silent logical; should output of \code{system} call for ffmpeg operation be suppressed.
#' @return Outputs a video of name "video.name+vid.ext".
#' @export
#' @details Assumes images are appended with a numeric sequence.
#' @seealso \code{\link{vid.to.images}}
#' @examples
#'
#' #make some images
#' \donttest{
#' dir.create(paste0(tempdir(),"/images")) #make a directory to store images
#'
#' a <- 2
#' b <- 3
#' theta <- seq(0,10*pi,0.01)
#' r <- a + b*theta
#' df <- data.frame(x=r*cos(theta), y=r*sin(theta)) # Cartesian coords
#' every.i <- 30
#' for(i in seq(1,length(theta),30)) {
#'   jpeg(paste0(tempdir(),"/images/image_",sprintf("%03d",which(i==seq(1,length(theta),30))),".jpg"))
#'   with(df[1:i,],plot(x,y,xlim=range(df$x),ylim=range(df$y),col="red"))
#'   dev.off()
#'   }
#'
#'images.to.video(image.dir=paste0(tempdir(),"/images"),
#'vid.name="spiral.mp4",out.dir=tempdir(),
#'frame.rate=5,qual=100,silent=TRUE,overwrite=TRUE)
#'
#'file.exists(paste0(tempdir(),"/spiral.mp4"))
#'
#' #clean up
#' unlink(paste0(tempdir(),"/spiral.mp4"))
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#' }

images.to.video <- function(image.dir=NULL,out.dir=NULL,vid.name=NULL,qual=50,frame.rate=10,overwrite=FALSE,silent=TRUE)  {
  
  version <-  try(system(paste("ffmpeg -version"), intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command 'ffmpeg' is not available in your system. Please install FFmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}
  
  qual <- round(30-30*(qual/100)+1,0)
  
  if(!dir.exists(image.dir)) stop("Directory specified by 'image.dir' (", paste0(image.dir),") does not exist")
 
  if(is.null(out.dir)) stop("'out.dir' not specified.")

  if(!file.exists(out.dir)) stop("Directory specified by 'out.dir' (", paste0(out.dir),") does not exist")
  
  out.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",out.dir) #if image.dir has spaces
  
  vid.path <- paste0(out.dir,"/",vid.name)

  
  if(file.exists(vid.path) & overwrite==FALSE) stop("video with name 'vid.name'  exist in 'out.dir' directory. To save file of this name in this location, 'overwrite' must be 'TRUE'")
  
  if(file.exists(vid.path) & overwrite==TRUE) unlink(vid.path,recursive = T)
  
  images <- paste0(image.dir,"/",list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp"))
  
  
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path) #if vid.path has spaces
  
  images <- paste0(image.dir,"/",list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp",ignore.case = T))
  
  image.name <- gsub("(.+)\\_\\d*\\.\\w+$", "\\1", basename(images[1]))
  
  ext <-    gsub(".*(\\.\\w+$)", "\\1", basename(images[1]))
  
  num <- gsub(".*_(\\d+)\\.\\w+$","\\1",basename(images[1]))
  num.l <- nchar(num)
  num.for <- paste0("_%",num.l,"d",ext)
  
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",image.dir) #if image.dir has spaces

  
  system(paste0("ffmpeg -i ", image.dir,"/", image.name,num.for," -q:v ",qual," -r ", frame.rate," -vcodec mpeg4 ", vid.path),ignore.stderr = silent) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
  
  if(!file.exists(vid.path)) warning("ffmpeg failed to create video. Retry and inspect output of system call with 'silent=FALSE'")
  
  if(file.exists(vid.path)) message("successfully created video \'", paste0(vid.path),"\'")
  
}


#' @title Extracts images from a video file with ffmpeg
#' @description Extract images from video file using ffmpegs flexible video filters and codecs
#'
#' @param vid.path character; path of video file to be processed.
#' @param out.dir character; directory path in which to store images.
#' @param overwrite logical; should path described by 'out.dir' be overwritten if it exhists.
#' @param filt character; video filter that should be applied to ffmpeg operation. See \url{https://ffmpeg.org/ffmpeg-filters.html}
#' @param codec character; video codec to apply in ffmpeg operation
#' @param silent logical; should output of \code{system} call for ffmpeg operation be suppressed.
#' @details Particularly useful for resizing images
#' @return Extracts all the images of the video and saves them to an "images" directory with appended number sequence
#' @seealso \code{\link{images.to.video}}
#' @export
#' @examples
#' 
#' #make a video with animation package
#' \donttest{
#' fun <- function(){
#' y <- sin(1:50)
#' x <- 1:50
#' for(i in 1:50) {
#'   plot(x[i],y[i],col="red",xlim=c(0,50),ylim=range(y))
#'   animation::ani.pause()
#'   }
#' }
#' animation::saveVideo(fun(),video.name=paste0(tempdir(),"/wave.mp4"),interval = 0.2)
#'
#'#reduce the image images to 200 px wide maintaining aspect ratio
#'#notice the spaces at the beginning/end of string
#'filt.red <- " -vf scale=200:-1 "
#'c <- " -c:v libx264 "
#'dir.create(paste0(tempdir(),"/images"))
#' vid.to.images2(vid.path=paste0(tempdir(),"/wave.mp4"),
#' out.dir=paste0(tempdir(),"/images"),filt=filt.red,codec=NULL)
#'
#' #see the images in the "images" directory
#' list.files( paste0(tempdir(),"/images"))
#' 
#' #clean up
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#' }

vid.to.images2 <- function(vid.path=NULL,out.dir=NULL,overwrite=FALSE,filt=NULL,codec=NULL,silent=TRUE)  {
  
  version <-  try(system("ffmpeg -version", intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command 'ffmpeg' is not available in your system. Please install FFmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}
  
  

  #delete or create the out directory
  if(is.null(out.dir)) stop("'out.dir' not specified.")
  if(is.null(vid.path)) stop("'vid.path' not specified.")
  
  out.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",out.dir) #if image.dir has spaces
  
  if(file.exists(out.dir) & overwrite==TRUE) unlink(out.dir,recursive = T)
  
  if(!file.exists(out.dir)) stop("Directory specified by 'out.dir' (", paste0(out.dir),") does not exist")
  
  if(!file.exists(vid.path)) stop("Path specified by 'vid.path' (", paste0(vid.path),") does not exist")
  
  #vid.path <- gsub("Google Drive","\"Google Drive\"",vid.path) ## remove spaces from dir if google drive
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path)
  
  #image.dir <- paste0(gsub("Google Drive","\"Google Drive\"",image.dir),"/") #degooglize path
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",out.dir)
  
  video.name<- gsub(".avi","",basename(vid.path))
  video.name <- gsub(".avi","_red",video.name)
  
  if(is.null(codec)) codec <- " "
  system(paste0("ffmpeg -i ", vid.path, codec, image.dir,"/",video.name,"_%5d.jpg -y"),ignore.stderr = silent) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
  
  if(!is.null(filt)){
    system(paste0("ffmpeg -i ", image.dir,"/",video.name,"_%5d.jpg", filt, image.dir,"/",video.name,"_%5d.jpg  -y"),ignore.stderr = silent) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
    
  }
}

#' @title  Stitches images from video file passing filters to ffmpeg
#'
#' @description Wrapper for ffmpeg video operations. Permits flexible filtering.
#'
#' @param image.dir character; directory containing images to stitch.
#' @param vid.name character; file name to be given video (should not include file extension).
#' @param out.dir character; directory in which to store video.
#' @param qual numeric; the quality of the video rendered from 1-100\%. Defaults to 50\%.
#' @param vid.ext character; video type to output. mp4 currently works best.
#' @param overwrite logical; should path described by vid.name  be overwritten if it exhists. 
#' @param frame.rate numeric; video frame rate in fps.
#' @param raw logical; encodes a raw AVI video with the "rawvideo" codec.
#' @param filt character; video filter that should be applied to ffmpeg operation. See \url{https://ffmpeg.org/ffmpeg-filters.html}.
#' @param silent logical; should output of \code{system} call for ffmpeg operation be suppressed.
#'
#' @return Outputs a video of name "video.name+vid.ext".
#' @export
#' @details Assumes images are appended with a numeric sequence beginning with "_".
#' @seealso \code{\link{vid.to.images2}}
#' @examples
#'
#' #make some spiralled images and video
#'
#'\donttest{
#' dir.create(paste0(tempdir(),"/images")) #make a directory to store images
#'
#' a <- 2
#' b <- 3
#' theta <- seq(0,10*pi,0.01)
#' r <- a + b*theta
#' df <- data.frame(x=r*cos(theta), y=r*sin(theta)) # Cartesian coords
#' every.i <- 30
#' for(i in seq(1,length(theta),30)) {
#'   jpeg(paste0(tempdir(),"/images/image_",sprintf("%03d",which(i==seq(1,length(theta),30))),".jpg"))
#'   with(df[1:i,],plot(x,y,xlim=range(df$x),ylim=range(df$y),col="red"))
#'   dev.off()
#'   }
#'
#'images.to.video2(image.dir=paste0(tempdir(),"/images"),
#'vid.name="spiral",out.dir=tempdir(),
#'frame.rate=5,qual=100,silent=TRUE,overwrite=TRUE)
#'
#'file.exists(paste0(tempdir(),"/spiral.mp4"))
#'
#' #clean up
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#'}

images.to.video2 <- function(image.dir=NULL,out.dir=NULL,vid.name=NULL,overwrite=TRUE,qual=50,vid.ext=".mp4",frame.rate=10,raw=TRUE,filt=NULL,silent=TRUE)  {
    
  if(!raw)   vid.name <- paste0(vid.name,"_red")
  
  version <-  try(system(paste("ffmpeg -version"), intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command 'ffmpeg' is not available on your system. Please install ffmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}
  
  qual <- round(30-30*(qual/100)+1,0)
  
  if(!dir.exists(image.dir)) stop("Directory specified by 'image.dir' (", paste0(image.dir),") does not exist")
  
  if(is.null(out.dir)) stop("'out.dir' not specified.")
  
  if(!file.exists(out.dir)) stop("Directory specified by 'out.dir' (", paste0(out.dir),") does not exist")
  
  if(raw) vid.ext <- ".avi" 
  
  out.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",out.dir) #if image.dir has spaces

  vid.path <- paste0(out.dir,"/",vid.name)
  vid.path2 <- paste0(out.dir,"/",vid.name,vid.ext)

  if(file.exists(vid.path2) & overwrite==FALSE) stop("video with name 'vid.name'  exist in 'out.dir' directory. To save file of this name in this location, 'overwrite' must be 'TRUE'")
  

  if(file.exists(vid.path2) & overwrite==TRUE) unlink(vid.path2,recursive = T)

  
  images <- paste0(image.dir,"/",list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp",ignore.case = T))

  
  image.name <- gsub("(.+)\\_\\d*\\.\\w+$", "\\1", basename(images[1]))
  
  ext <-    gsub(".*(\\.\\w+$)", "\\1", basename(images[1]))
  
  num <- gsub(".*_(\\d+)\\.\\w+$","\\1",basename(images[1]))
  num.l <- nchar(num)
  num.for <- paste0("_%",num.l,"d",ext)
  
  image.dir <- normalizePath(dirname(images[1]))
  
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",image.dir) #if image.dir has spaces
 
  
  if(!raw) system(paste0("ffmpeg -i ", image.dir,"/", image.name,num.for," -q:v ",qual," -r ", frame.rate," -f mp4", filt," -vcodec libx264 -pix_fmt yuv420p ", vid.path,vid.ext, " -y"),ignore.stderr = silent) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
  
 
  if(raw) system(paste0("ffmpeg -i ", image.dir,"/", image.name,num.for," -q:v ",qual," -r ", frame.rate," -f avi -vcodec rawvideo ", vid.path,vid.ext, " -y"),ignore.stderr = silent)
  
  message(paste0("video saved to ", vid.path,vid.ext))
}
