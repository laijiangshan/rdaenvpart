#' Hierarchical Partitioning for Canonical Correspondence Analysis and Redundancy Analysis
#'
#' This function calculates the individual explain percentage of each environmental variable for Canonical Correspondence Analysis and Redundancy Analysis,
#' applying the hierarchy algorithm of Chevan and Sutherland (1991) .
#' 
#'
#' @param  Y Community data matrix.
#' @param  X Constraining matrix less than 12 columns, typically of environmental variables.
#' @param  type the Constrained ordination: RDA or CCA, default "RDA"
#' @param  pieplot a pieplot each variable is plotted expressed as percentage of total variation (pieplot="tv") or total explained variation (pieplot="tev").
#' @return a list containing
#' @return \item{R2}{unadjusted R-squared for RDA or CCA with all environmental variables.}
#' @return \item{gfs}{a data frame or vector listing all combinations of individual variables in the first column in ascending order, and the corresponding goodness of fit measure for the model using those variables.}
#' @return \item{IJ.R2}{a data frame of I, the individual and J the joint contribution for each individual variable (based on unadjusted R-squared).}
#' @return \item{I.perc}{a data frame of I as a percentage of total explained variation.}
#' @return \item{adj.R2}{adjusted R-squared for RDA or CCA with all environmental variables.}
#' @return \item{IJ.adjR2}{the individual contribution for each individual variable (based on adjusted R-squared).}
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @references
#' Chevan, A. and Sutherland, M. 1991. Hierarchical Partitioning. The American Statistician 45:90~96
#' @examples
#'require(vegan)
#'data(varespec)
#'data(varechem)
#'rdaenvpart(varespec,varechem[,c("Al","P","K")],pieplot = "tv",type="RDA")
#'rdaenvpart(varespec,varechem[,c("Al","P","K")],pieplot = "tev",type="RDA")
#'rdaenvpart(varespec,varechem[,c("Al","P","K")],pieplot = "tv",type="CCA")
#'rdaenvpart(varespec,varechem[,c("Al","P","K")],pieplot = "tev",type="CCA")


rdaenvpart=function (Y, X,type="RDA", pieplot = "tv")
{
    Env.num <- dim(X)[2]
    if (Env.num > 13)
        stop("Number of variables must be < 13 for current implementation",call. = FALSE)
    else {
        gfs <- allR2(Y, X,type)
        HP <- partition.rda(gfs, Env.num, var.names = names(data.frame(X)))
       
	   if(type=="RDA")
            rsq <- RsquareAdj(rda(Y~., data = X))
       if(type=="CCA")
             rsq <- RsquareAdj(cca(Y~., data = X))
	   if (pieplot=="tv") {
            
           lbls<- c("unexplained",row.names(HP$I.perc)) 
		   pct <- round(c(1-sum(HP$IJ$I),HP$IJ$I)*100,1)
           lbls <- paste(lbls, pct) # add percents to labels
           lbls <- paste(lbls,"%",sep="") # ad % to labels
           pie(pct,labels = lbls, main="individual % on total variation")
        }
		if (pieplot=="tev") {
           lbls<- row.names(HP$I.perc) 
		   pct <- round(HP$I.perc$I,1)
           lbls <- paste(lbls, pct) # add percents to labels
           lbls <- paste(lbls,"%",sep="") # ad % to labels
           pie(pct,labels = lbls,col=rainbow(length(lbls)), main="individual % on total explained variation")
        }
		
         list(R2=rsq$r.squared,gfs = gfs, IJ.R2 = HP$IJ, I.perc = HP$I.perc,adj.R2=rsq$adj.r.squared,IJ.adjR2=HP$I.perc*0.01*rsq$adj.r.squared)
    }
}

