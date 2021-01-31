# This document has functions that are geared towards the use of odds instead of
# effects. This may be more appropriate for screens that are based on selection,
# such as FACS-based screens, screens with migration or invasion assays, and
# what have you. The use of distributions is quite different from the summed Z
# approach or the growth effect derivation, therefore my hope is that it is more
# robust for these selection screens, in which rate ratios may vary much more
# widely.

#' Calculate odds per guide
#'
#' Returns odds based on rate ratios of guides. The odds are based on the rank
#' of a guide compared to a subset of control guides or the entire library.
#'
#' @param r Numeric vector. Rate ratios of a test versus control arm
#' @param normsubset Integer vector. Specify the indices of features that
#'   function as controls. If omitted, all features are used.
#' @param log Integer or logical. Log-transform odds? If FALSE, odds are not
#'   transformed. If TRUE, odds are transformed by the natural logarithm.
#'   Alternatively, specify the base with which to log-transform. Default = 10
#'
#' @return Returns a numeric vector with (by default log10-transformed) odds, of
#'   the same length as the number of rate ratios provided.
#'
#' @note Odds in this function represent the odds that a feature has a rate
#'   ratio that is this extreme, either relative to the entire data set or to
#'   the normalization subset when provided. These odds are mostly based on
#'   rank. Only if they are more extreme than the most extreme control are they
#'   extrapolated. Even then, extrapolation is quite conservative. There is also
#'   a small interpolation step, which looks at how far the feature is away from
#'   the flanking control features.
#'
#' @author Jos B. Poell
#'
#' @seealso \code{\link{geneodds}}, \code{\link{odds2pq}}
#'
#' @export

oddscores <- function(r, normsubset, log = 10) {
  if (missing(normsubset)) {
    odds <- sapply(r, function(x) {
      sum(x <= r)/sum(x >= r)
    })
    if (log == FALSE) {
      return(odds)
    } else if (log == TRUE) {
      return(log(odds))
    } else {
      return(log(odds, log))
    }
  } else {
    rsub <- sort(r[normsubset])
    rmax <- max(rsub)
    rmin <- min(rsub)
    rmed <- median(rsub)
    n <- length(rsub)
    odds <- sapply(r, function(x) {
      if (x <= rmin) {
        return((rmed-rmin)/((rmed-x)*n))
      } else if (x >= rmax) {
        return(((x-rmed)*n)/(rmax-rmed))
      } else {
        rankup <- sum(rsub <= x)
        rankdown <- sum(rsub >= x)
        mod <- (x-rsub[rankup])/(rsub[rankup+1]-rsub[rankup])
        return((rankup+mod)/(rankdown-mod))
      }
    })
    if (log == FALSE) {
      return(odds)
    } else if (log == TRUE) {
      return(log(odds))
    } else {
      return(log(odds, log))
    }
  }
}

#' Calculate odds per gene
#'
#' geneodds combines the odds of all guides per gene.
#'
#' @param guides Character vector. List of guides
#' @param oddscores Numeric vector. oddscores of guides
#' @param log Logical. Are the oddscores log-transformed or not? Default = TRUE
#'
#' @return Returns a numeric vector of the same length as the number of genes
#'   with the combined odds.
#'
#' @note In a similar fashion as \code{\link{sumZ}} and \code{\link{getdeg}},
#'   this function combines odds of all guides belonging to a gene for all genes
#'   in the data set. Keep in mind that guides have to be sorted by gene name,
#'   and there value is expected to start with the gene name, number or symbol,
#'   followed by an underscore. Odds are combined by multiplying or adding log
#'   odds, after which the odds are corrected for being composed of multiple
#'   odds. The code in the example shows an example with uniformly distributed
#'   p-values.
#'
#' @author Jos B. Poell
#'
#' @seealso \code{\link{oddscores}}, \code{\link{odds2pq}}, \code{\link{sumZ}},
#'   \code{\link{getdeg}}
#'   
#' @examples
#' guides <- paste0("gene", rep(1:2500, each = 4), "_", rep(1:4, 2500))
#' p <- runif(10000)
#' odds <- p/(1-p)
#' hist(log(odds, 10), breaks = 100, xlim = c(-3,3))
#' go <- geneodds(guides, odds, log = FALSE)
#' hist(log(go, 10), breaks = 25, xlim = c(-3,3))
#' pgene <- go/(1+go)
#' hist(pgene, breaks = 25)
#' logodds <- log(odds,10)
#' glo <- geneodds(guides, logodds, log = 10)
#' hist(glo, breaks = 25, xlim = c(-3,3))
#'
#' @export

geneodds <- function(guides, oddscores, log = TRUE) {
  genes <- rle(gsub("_.*", "", guides))$values
  n <- rle(gsub("_.*", "", guides))$lengths
  if (log != FALSE) {
    # Note that we still have to correct for number of guides!
    return(mapply(function(n,gi) {sum(oddscores[(gi-n+1):gi])/sqrt(n)}, 
                  n, cumsum(n)))
  } else {
    return(mapply(function(n,gi) {prod(oddscores[(gi-n+1):gi])^(1/sqrt(n))}, 
                  n, cumsum(n)))
  }
}


# This quirky little function will turn your odds into p-values and q-values.
# But there is a quirk! If your odds are expressed as log-values, your p-values
# and q-values will be expressed in the same log base! There is a tiny cheat
# going on that ensures you can keep working with ridiculous odds. Note that the
# p-values are basically double-sided, meaning that odds of 1 (log odds of 0)
# will return a p-value of 1 instead of 0.5. The q-values are using a
# Benjamini-Hochberg correction, again staying in the same log base if
# applicable.

#' Calculate p-values and q-values for odds
#'
#' Returns a data frame with p-values and q-values for a vector of odds. When
#' using log odds, it keeps the p-values and q-values in the same log base.
#'
#' @param odds Numeric vector. List of odds
#' @param log Integer or logical. Are the odds log-transformed or not? If not,
#'   specify FALSE. If TRUE, the odds are presumed to be transformed by the
#'   natural logarithm. Default = 10
#'
#' @return Returns a data frame with the odds, the p-values and the q-values
#'
#' @note The function basically calculates the probability of something being
#'   more extreme than the corresponding odds. Conversion uses the formula p =
#'   odds/(1+odds) or 1 - odds/(1+odds) if odds are higher than 1. The resulting
#'   value is multiplied by 2 to correct for 2-sided testing. For log odds, the
#'   absolute log value is subtracted from the log of 2 using the corresponding
#'   log base. Though not exact, this is a very close approximation with high or
#'   low odds, which is what this analysis is meant for.
#'
#' @author Jos B. Poell
#'
#' @seealso \code{\link{oddscores}}, \code{\link{geneodds}}
#'
#' @export

odds2pq <- function(odds, log = 10) {
  if (log == TRUE) {
    # I use a bit of a cheat here to convert odds to probabilities, but it
    # converges to the true probability for large or small odds, which are the
    # ones we are going to be interested in anyway.
    p <- pmin(0, log(2)-abs(odds))
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    q <- pmin(0, cummin(log(length(p)/seq(length(p),1)) + p[o]))[ro]
  } else if (log != FALSE) {
    p <- pmin(0, log(2, log)-abs(odds))
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    q <- pmin(0, cummin(log(length(p)/seq(length(p),1),log) + p[o]))[ro]
  } else {
    p <- sapply(odds, function(o) {
      if (o <= 1) {
        (2*o)/(1+o)
      } else {
        2-(2*o)/(1+o)
      }
    })
    q <- p.adjust(p, method = "BH")
  }
  return(data.frame(odds=odds,p=p,q=q))
}
