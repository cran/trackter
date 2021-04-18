#' @title Computes distance between two points in Cartesian space.
#'
#' @description Computes distance between two points in Cartesian space using simple trigonometry functions
#'
#' @param x1 Numeric; x position of coordinate 1
#' @param y1 numeric; y position of coordinate 1
#' @param x2 numeric; x position of coordinate 2
#' @param y2 numeric; y position of coordinate 2
#' @return A single value of the distance between p[x1,y1] and p[x2,y2]
#' @export
#' @examples
#' #Find the lengths of the sides of a tringle and print to plot
#' x <- c(0,3,2)
#' y <- c(0,3,0)
#' plot(x,y)
#' lines(x,y)
#' lines(x[c(1,3)],y[c(1,3)])
#' hyp <- dist.2d(x[1],x[2],y[1],y[2])
#' s1 <- dist.2d(x[1],x[3],y[1],y[3])
#' s2 <- dist.2d(x[2],x[3],y[2],y[3])
#' text(mean(x[1:2],mean(y[2:3])),labels=round(hyp,1))
#' text(mean(x[c(1,3)]),y[1]+0.25,labels=round(s1,1))
#' text(mean(x[c(2:3)]),mean(y[2:3]),labels=round(s2,1))

dist.2d <- function(x1, x2, y1, y2) {
  sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2)
}

#' @title Computes angle between two segments sharing a point.
#' @description Computes angle between two segments of a triangle using law of cosines.
#' @param l Numeric; length of segment to the left
#' @param r Numeric; length of segment to the right
#' @param o Numeric; length of opposite segment
#' @return A single value of the angle in radians
#' @export
#' @examples
#'#a right triangle
#'L=3
#'R=3
#'O=sqrt(L^2+R^2)
#'cosine.ang(L,R,O)

cosine.ang <- function(l,r,o) {
  if(o==(l+r)){return(pi)}else{
    a <- acos((o^2-l^2-r^2)/(-2*l*r))
    return(a)}
}

#' @title Converts radians to degrees
#'
#' @param x Numeric; value in radians
#' @return A single value
#' @export
#' @seealso \code{\link{rad}}
deg <- function(x){180/pi*x}


#' convert degrees to radians
#'
#' @param x Numeric; value in degrees
#' @return A single value
#' @export
#' @seealso \code{\link{deg}}
rad <- function(x){pi/180*x}

#' @title Computes the heading between to cartesian points
#' @description Computes the heading (radians counter clockwise from north or vertical)
#'
#' @param x1 Numeric; point A x value
#' @param y1 numeric; point A y value
#' @param x2 numeric; point B x value
#' @param y2 numeric; point B y value
#' @return A single value in radians
#' @export
#' @return a single value in radians
#'
#' @examples
#' #example
#' A <- c(0,0)
#' B <- c(-5,5)
#' thet <- bearing.xy(A[1],B[1],A[2],B[2])
#' deg(thet)
#'
#'
bearing.xy <- function(x1,x2,y1,y2){
  theta <- atan((x2-x1)/(y1-y2))
  return(theta)
}

