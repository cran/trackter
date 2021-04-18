#' Contour data of fish and arbitrary shapes 
#'
#' Contour data constructed with the \code{Momocs} package of several fish shapes and arbitrary polygons. Intended for training of PCA and LDA analysis with \code{LDA}. Shapes are classified with shape, type, and edge factors. 
#' 
#'\itemize{
#'\item 'shape': names the contour
#'\item 'type': classifies the shape as 'fish' or not.fish'
#'\item 'edge': is every coordinate of the contour represented, i.e., not cut off by an edge of the field?  One level of 'FALSE'.
#'} 
#'
#'Shapes of type 'fish' include the follow classification levels:
#'\itemize{
#' \item 'eel' (genus Anguilla) swimming with body undulations
#' \item 'sunfish_BCF' (genus Lepomis) swimming with body-caudal fin propulsion
#' \item 'sunfish_pect' (genus Lepomis) swimming with both body-caudal fin and pectoral fin oscillations
#' \item 'trout' (genus Oncorhynchus) swimming with body-caudal fin propulsion
#' }
#' 
#'   Each of these fish types include contours sampled regularly over one tail-beat cycle. 
#' 
#'Shapes of type 'not.fish' include 'Ellipse', 'Lshape', 'Ushape', Ushape2', 'Ushape3','triangle', 'rectangle', and 'square'. Each 'not.fish' type was resampled 6 times with different \code{efourier} analyses with 'nb.h' values ranging from 5 to 30. This produces shape classes with subtly variable contours.
#'
#'The 'edge' factor is included so as to have the classification factors match those of the \code{efourier} and \code{LDA} analysis in \code{kin.LDA}.
#'
#
#' @docType data
#' @rdname fish-shapes
#' @name fishshapes
#' @usage data(fishshapes)
#' @import data.table
#' @import Momocs
#'
#' @format An object of class \code{"Out"} and \code{"Coo"}; see \code{\link[Momocs]{Out}}.
#' 
#'
#' @keywords datasets
#'
#' @examples
#' library(Momocs)
#' data(fishshapes)
#' panel(fishshapes)
NULL
