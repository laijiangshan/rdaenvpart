#'Permutation Test for Hierarchical Partitioning for Canonical Correspondence Analysis and Redundancy Analysis
#'
#' This function performs permutation test for the individual explain percentage of each environmental variable for Canonical Correspondence Analysis and Redundancy Analysis,
#' applying the hierarchy algorithm of Chevan and Sutherland (1991) .
#' 
#' @param  object one object from rdaenvpart
#' @param  Y Community data matrix.
#' @param  X Constraining matrix less than 12 columns, typically of environmental variables.
#' @param  type the Constrained ordination: RDA or CCA, default "RDA"
#' @return a list containing
#' @return \item{R2}{a dataframe for unadjusted R-squared for individual environmental variables and p-value.}
#' @return \item{adjR2}{a dataframe for adjusted R-squared for individual environmental variables and p-value.}
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @references
#' Chevan, A. and Sutherland, M. 1991. Hierarchical Partitioning. The American Statistician 45:90~96
#' @examples
#'require(vegan)
#'data(varespec)
#'data(varechem)
#'mod=rdaenvpart(varespec,varechem[,c("Al","P","K")],pieplot = "tv",type="RDA")
#'perm.env(mod,varespec,varechem[,c("Al","P","K")],type="RDA",nper=999)
#'mod2=rdaenvpart(varespec,varechem[,c("Al","P","K")],pieplot = "tv",type="CCA")
#'perm.env(mod2,varespec,varechem[,c("Al","P","K")],type="CCA",nper=999)

perm.env=function(object,Y,X,type="RDA",nper=999)
{require(permute)
obs=object
n=dim(Y)[1]
r2q=obs$IJ.R2["I"]
ar2q=obs$IJ.adjR2
for(i in 1:nper)
{newy=Y[shuffle(n),]
 simu=rdaenvpart(newy,X,type=type,pieplot ="")
 r2q=cbind(r2q,simu$IJ.R2["I"])
 ar2q=cbind(ar2q,simu$IJ.adjR2)
}

Signi=function(x)
{1-ecdf(x)(x[1])+1/(nper+1)}

p.R2=apply(r2q,1,Signi)
p.adjR2=apply(ar2q,1,Signi)
return(list(R2=data.frame(R2=round(obs$IJ.R2["I"],3),Pr=p.R2),adjR2=data.frame(adjR2=round(obs$IJ.adjR2,3),Pr=p.adjR2)))
}




