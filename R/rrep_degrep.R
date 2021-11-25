#' Calculate rate ratios with standard error of count data based on replicates
#'
#' rrep calculates rate ratios and corresponding standard errors of features
#' based on the presence of replicates. The function uses the compound variance
#' of the rate ratio of replicates as well as the variance associated with read
#' depth (i.e. counts). Rate ratio and standard error are expressed as
#' log2-values.
#'
#' @param t1 Matrix or data frame, with rows representing features and columns
#'   representing replicates of measurements at t1 (or treated).
#' @param t0 Matrix or data frame, with rows representing features and columns
#'   representing replicates of measurements at t0 (or untreated).
#' @param paired Logical. Are measurements paired? Default = TRUE
#' @param normfun Character string. Specify with which function to standardize
#'   the data. Default = "sum"
#' @param normsubset Integer vector. Specify the indices of features that are to
#'   be used in standardization
#' @param rstat Character string. Specify whether rate ratios are calculated
#'   over the sum of the counts in replicates or as the median or mean of the
#'   log2 rate ratios. Default = "summed"
#' @param variance Character string. Specify how to calculate variance: using
#'   only \code{count} variance, only \code{replicate} variance, or the both
#'   \code{combined}. Default = "combined"
#'
#' @details This function combines the confidence based on replicate
#'   measurements with the confidence based on counts (e.g. read depth).
#'   Utilizing replicates to assess confidence in point estimates of individual
#'   features is commonplace in many analyses. However, in data sets with many
#'   features just by chance some features will have measurements that lie very
#'   close together. By adding the variance based on the count data, spurious
#'   findings are greatly reduced, especially when counts are low. Variance of
#'   count data on a log2-transformed scale is approximated with the formula
#'   \code{1/(log(2)^2*count)}
#'
#' @return Returns a data frame with the log2-transformed rate ratio and
#'   corresponding standard error of each feature.
#'
#' @note Counts are checked for zeros per feature (row). In case of any zeros, a
#'   pseudocount of 1/replicates is added to all counts in that row. These
#'   pseudocounts are not included in the normalization on the bases of total
#'   counts (of all features are the normalization subset) in an experimental
#'   arm
#'
#' @seealso \code{\link{degrep}}, \code{\link{CRISPRsim}}, \code{\link{ess}},
#'   \code{\link{noness}}
#'
#' @author Jos B. Poell
#'
#' @examples
#'   set.seed(1000)
#'   c0 <- rbind(rpois(3, 60), rpois(3, 5), rpois(3, 10000000))
#'   c1 <- rbind(rpois(3, 20), rpois(3, 20), rpois(3, 10000000))
#'   rrep(c1, c0, paired = FALSE)
#'   c01 <- rbind(rep(60, 3), rep(5, 3), rep(10000000, 3))
#'   c11 <- rbind(rep(20, 3), rep(20, 3), rep(10000000,3))
#'   rrep(c11, c01)
#'   c04 <- 4*c01
#'   c14 <- 4*c11
#'   rrep(c14,c04)
#'
#' @export

rrep <- function(t1, t0, paired = TRUE, normfun = "sum", normsubset, 
                 rstat = "summed", variance = "combined") {
  fun <- get(normfun)
  if (!missing(normsubset)) {
    sr0 <- apply(t0[normsubset,], 2, fun)
    sr1 <- apply(t1[normsubset,], 2, fun)
  } else {
    sr0 <- apply(t0, 2, fun)
    sr1 <- apply(t1, 2, fun)
  }
  reps1 <- ncol(t1)
  reps0 <- ncol(t0)
  if (nrow(t1) != nrow(t0)) {
    stop("Unequal number of features in both arms")
  }
  if (reps1 < 2 || reps0 < 2) {
    stop("All experimental arms need at least 2 replicates for this analysis")
  }
  if (paired==TRUE && reps0 != reps1) {
    stop("In paired analysis both arms need equal number of replicates")
  }
  
  c1 <- apply(cbind(t1,t0), 1, function(x) {
    if (any(x==0)) x[1:reps1]+1/reps1 else x[1:reps1]
  })
  c0 <- apply(cbind(t0,t1), 1, function(x) {
    if (any(x==0)) x[1:reps0]+1/reps0 else x[1:reps0]
  })
  
  if (any(!is.numeric(cbind(c1,c0)))) {
    stop("Make sure all cells of t0 and t1 are numeric")
  }
  
  if (paired == TRUE) {
    if (rstat == "summed") {
      r <- log((apply(c1, 2, sum)/sum(sr1))/(apply(c0, 2, sum)/sum(sr0)), 2)
    } else if (rstat == "median") {
      r <- apply(log((c1/sr1)/(c0/sr0), 2), 2, median)
    } else if (rstat == "mean") {
      r <- apply(log((c1/sr1)/(c0/sr0), 2), 2, mean)
    } else {
      stop ("Invalid rstat: choose from summed, mean or median")
    }
    varcount <- apply(rbind(c1, c0), 2, function(x) sum(1/(log(2)^2*x)))/reps1
    varrep <- apply(log((c1/sr1)/(c0/sr0), 2), 2, var)
    if (variance == "count") {
      se <- sqrt(varcount/reps1)
    } else if (variance == "replicates") {
      se <- sqrt(varrep/reps1)
    } else {
      se <- sqrt((varcount+varrep)/reps1)
    }
  } else {
    if (rstat == "summed") {
      r <- log((apply(c1, 2, sum)/sum(sr1))/(apply(c0, 2, sum)/sum(sr0)), 2)
    } else if (rstat == "median") {
      r <- apply(log((c1/sr1), 2), 2, median) - apply(log((c0/sr0), 2), 2, median)
    } else if (rstat == "mean") {
      r <- apply(log((c1/sr1), 2), 2, mean) - apply(log((c0/sr0), 2), 2, mean)
    } else {
      stop ("Invalid rstat: choose from summed, mean or median")
    }
    varcount <- apply(c1, 2, function(x) mean(1/(log(2)^2*x))) + apply(c0, 2, function(x) mean(1/(log(2)^2*x)))
    varrep <- apply(log(c1/sr1, 2), 2, var) + apply(log(c0/sr0, 2), 2, var)
    if (variance == "count") {
      se <- sqrt(varcount/min(reps0,reps1))
    } else if (variance == "replicates") {
      se <- sqrt(varrep/min(reps0,reps1))
    } else {
      se <- sqrt((varcount+varrep)/min(reps0,reps1))
    }
  }
  
  return(data.frame(r=r,se=se))
}


#' Derive growth-modifying effect of gene knockout in pooled experiments with
#' replicate arms
#'
#' degrep is a variant of getdeg that utilizes confidence measures of rate
#' ratios to find the "best guide". First, the median rate ratio of a group
#' (e.g. a gene) is determined. The best guide has the most extreme rate ratio
#' with the same sign (direction) as the median, after moving two (or another
#' specified number) standard (error) units toward null. See details below. For
#' more context, see \code{\link{getdeg}}
#'
#' @param guides Character vector. Guides are assumed to start with the gene
#'   name, followed by an underscore, followed by a number or sequence unique
#'   within that gene.
#' @param r0 Numeric vector. Log2-transformed rate ratios of features
#'   representing straight lethality.
#' @param se0 Numeric vector. Standard errors corresponding to r0.
#' @param r1 Numeric vector. Log2-transformed rate ratios of features
#'   representing sensitization or synthetic lethality. Optional but required to
#'   calculate e.
#' @param se1 Numeric vector. Standard errors corresponding to r1.
#' @param rt Numeric vector. Log2-transformed rate ratios of features
#'   representing lethality in the test sample. Optional.
#' @param set Numeric vector. Standard errors corresponding to rt.
#' @param a Numeric. Estimated potential population doublings between time
#'   points.
#' @param b Numeric. Estimated potential population doublings between time
#'   points in test sample. Only applicable if r1 is given. If omitted, assumed
#'   equal to a.
#' @param hnull Numeric. Null hypothesis. Growth effects of genes are tested to
#'   be more extreme than this value. Setting hnull can greatly improve the
#'   usefulness of p-values, and can be considered a cutoff for relevance.
#'   Default = 0
#' @param nse Numeric. Number of standard units, used for comparing guides. See
#'   details below. Default = 2
#' @param secondbest Logical. If TRUE, calculate effect sizes based on the
#'   second best guides of each gene as well. Default = TRUE
#' @param correctab Logical. When \code{a != b}, it is be possible (and
#'   necessary?) to mathematically correct for this difference. If you analyze
#'   an experiment with unequal a and b, try both with and without correction.
#'   Default = TRUE
#'
#' @details For more details on basic functionality, see \code{\link{getdeg}}
#'   documentation. The added functionality in this function hinges on the use
#'   of confidence measures of rate ratios. Rate ratios and associated errors
#'   can be derived from other sources, but I recommend using the output of
#'   \code{\link{rrep}}. As in \code{\link{getdeg}}, a single best guide is
#'   determined using all available data. But now, the confidence given by the
#'   standard error is used to help select the best guide. Here is an example to
#'   illustrate. Say guide 1 has a rate ratio of -4 and a standard error of 1.2,
#'   while guide 2 targeting the same gene has a rate ratio of -3 and a standard
#'   error of 0.5 and guide 3 and guide 4 both have little effect. When using
#'   the default \code{nse = 2}, guide 2 scores better than guide 1, and is thus
#'   designated "best guide". However, for the p-value calculation, the lowest
#'   p-value is reported, which is calculated using the rate ratio, its standard
#'   error, and the null hypothesis as determined by \code{hnull} and the number
#'   of population doublings. The p-value reported for a gene does therefore not
#'   necessarily match to the best guide, and can in fact below to an outlier.
#'   The p-values are also not corrected for multiple testing. Instead, p-values
#'   can easily be calculated for all guides using \code{2*pnorm(-abs(r)/se)} or
#'   \code{pnorm((hnull*a-abs(r))/se)}, and corrected for multiple testing using
#'   \code{\link{p.adjust}}.
#'
#' @return Returns a list with the following (depending on input arguments):
#'   \itemize{ \item{genes}{ - list of all gene symbols} \item{n}{ - number of
#'   guides representing the gene} \item{d}{ - gene knockout effects on straight
#'   lethality} \item{d2}{ - gene knockout effects on straight lethality based
#'   on the second-best guide} \item{e}{ - gene knockout effects on
#'   sensitization} \item{e2}{ - gene knockout effects on sensitization based on
#'   the second-best guide} \item{de}{ - gene knockout effects on straight
#'   lethality in the test arm} \item{de2}{ - gene knockout effects on straight
#'   lethality in the test arm based on the second-best guide} \item{g}{ -
#'   estimated guide efficacy} \item{i}{ - within-gene index of the best guide}
#'   \item{j}{ - within-gene index of the second-best guide} \item{pd}{ -
#'   p-value of straight lethality} \item{pe}{ - p-value of sensitization}
#'   \item{pde}{ - p-value of lethality in the test arm} }
#'
#' @note With these analyses, it is important to visually inspect all steps, and
#'   preferentially to analyze a data set with several settings.
#'
#' @seealso \code{\link{getdeg}}, \code{\link{rrep}}
#'
#' @author Jos B. Poell
#'
#' @examples
#' ut1 <- CRISPRsim(1000, 4, a = c(3,3), allseed = 100, t0seed = 10, 
#'                  repseed = 1, perfectseq = TRUE)
#' tr1 <- CRISPRsim(1000, 4, a = c(3,3), e = TRUE, allseed = 100, t0seed = 10, 
#'                  repseed = 2, perfectseq = TRUE)
#' ut2 <- CRISPRsim(1000, 4, a = c(3,3), allseed = 100, t0seed = 20, 
#'                  repseed = 3, perfectseq = TRUE)
#' tr2 <- CRISPRsim(1000, 4, a = c(3,3), e = TRUE, allseed = 100, t0seed = 20, 
#'                  repseed = 4, perfectseq = TRUE)
#' ut3 <- CRISPRsim(1000, 4, a = c(3,3), allseed = 100, t0seed = 30, 
#'                  repseed = 5, perfectseq = TRUE)
#' tr3 <- CRISPRsim(1000, 4, a = c(3,3), e = TRUE, allseed = 100, t0seed = 30, 
#'                  repseed = 6, perfectseq = TRUE)
#' cgi <- tr1$d > -0.05 & tr1$d < 0.025 & tr1$e > -0.05 & tr1$e < 0.025
#' rr0 <- rrep(cbind(ut1$t6, ut2$t6, ut3$t6), cbind(ut1$t0, ut2$t0, ut3$t0), normsubset = cgi)
#' rr1 <- rrep(cbind(tr1$t6, tr2$t6, tr3$t6), cbind(ut1$t6, ut2$t6, ut3$t6), normsubset = cgi)
#' deg <- degrep(ut1$guides, rr0$r, rr0$se, rr1$r, rr1$se, a = 6, b = 6, secondbest = FALSE)
#' reald <- rle(tr1$d)$values
#' reale <- rle(tr1$e)$values
#' plot(reald, deg$d)
#' plot(reale, deg$e)
#'
#' @export

degrep <- function(guides, r0, se0, r1, se1, rt = FALSE, set = FALSE, 
                   a, b, hnull = 0, nse = 2, secondbest = TRUE, correctab = TRUE) {
  
  if (missing(a)) {
    stop("enter the presumed number of population doublings")
  }
  # note that guides are presumed to be named [gene]_[number | sequence]
  genes <- rle(gsub("_.*", "", guides))$values
  n <- rle(gsub("_.*", "", guides))$lengths
  gi <- cumsum(n)
  if (hnull < 0) {
    hnull <- abs(hnull)
    print(paste0("The null hypothesis for negative growth effects is > -", hnull))
    print(paste0("The null hypothesis for positive growth effects is < ", hnull))
  }
  rnulla <- a*hnull
  
  if (missing(r1)) {
    message("no input to calculate effect modifier e")
    dgip <- mapply(function(n, gi) {
      # First, calculate the median log2 rate ratio of all guides targeting a
      # gene. Then, find the rate ratio most se-units from hnull in the same
      # direction as the median. In other words: the median has to be in the
      # correct direction!
      if (median(r0[(gi-n+1):gi]) < 0) {
        i <- which.min(r0[(gi-n+1):gi]+nse*se0[(gi-n+1):gi])
        r <- r0[(gi-n+i)]
      } else {
        i <- which.max(r0[(gi-n+1):gi]-nse*se0[(gi-n+1):gi])
        r <- r0[(gi-n+i)]
      }
      
      d <- r/a
      g <- sapply(seq_len(n), function(i) {
        # The extra if statement is necessary to prevent 0/0
        if (d == 0 && r0[gi-n+i] == 0) {g <- 1} else {g <- (2^(r0[gi-n+i])-1) / (2^(a*d)-1)}
        # Mathematically, g can end up below 0 or above 1 in some instances.
        # We cannot have that ...
        g[g < 0] <- 0
        g[g > 1] <- 1
        return(g)})
      if (r < 0) {
        q <- max((-rnulla-r0[(gi-n+1):gi])/se0[(gi-n+1):gi])
      } else {
        q <- max((r0[(gi-n+1):gi]-rnulla)/se0[(gi-n+1):gi])
      }
      
      if (hnull == 0) {
        pd <- 2*pnorm(-q)
      } else {
        pd <- pnorm(-q)  
      }
      
      return(list(d=d,g=g,i=i,pd=pd))
    }, n, gi)
    if (secondbest == TRUE) {
      # return a list of within-gene indices all the best guides
      i <- unlist(dgip[3,])
      # number of guides per genes is reduced by one
      n2 <- n-1
      gi2 <- cumsum(n2)
      # best guides are excluded to find the second-best guide
      r02 <- r0[-(gi-n+i)]
      
      d2j <- mapply(function(n2, gi2) {
        # if there was only 1 guide to begin with, second best is NaN
        if (n2 == 0) {return(list(d2 = NaN, j = NaN))}
        else {
          if (median(r02[(gi2-n2+1):gi2]) < 0) {
            r <- min(r02[(gi2-n2+1):gi2])
          } else {
            r <- max(r02[(gi2-n2+1):gi2])
          }
          j <- tail(which(r02[(gi2-n2+1):gi2]==r),1)
          # note that guide efficacies and p-values are not calculated again
          return(list(d2=r/a, j = j))}
      }, n2, gi2)
      
      # code below is to find the correct within-gene index of the second-best guide
      j <- unlist(d2j[2,])
      k <- i[!is.na(j)]
      l <- j[!is.na(j)]
      l[l >= k] <- l[l >= k]+1
      j[!is.na(j)] <- l
      
    }
  } else {
    if (missing(b)) {
      warning("b is missing, assuming b equals a")
      b <- a
    }
    rnullb <- b*hnull
    degip <- mapply(function(n, gi) {
      # This time around, most extreme log2 rate ratios are found for all comparisons available
      if (median(r0[(gi-n+1):gi]) < 0) {
        i0 <- which.min(r0[(gi-n+1):gi]+nse*se0[(gi-n+1):gi])
        r00 <- r0[(gi-n+i0)]
      } else {
        i0 <- which.max(r0[(gi-n+1):gi]-nse*se0[(gi-n+1):gi])
        r00 <- r0[(gi-n+i0)]
      }
      if (median(r1[(gi-n+1):gi]) < 0) {
        i1 <- which.min(r1[(gi-n+1):gi]+nse*se1[(gi-n+1):gi])
        r11 <- r1[(gi-n+i1)]
      } else {
        i1 <- which.max(r1[(gi-n+1):gi]-nse*se1[(gi-n+1):gi])
        r11 <- r1[(gi-n+i1)]
      }
      # if (median(r0[(gi-n+1):gi]) < 0) {r00 <- min(r0[(gi-n+1):gi])} else {r00 <- max(r0[(gi-n+1):gi])}
      # if (median(r1[(gi-n+1):gi]) < 0) {r11 <- min(r1[(gi-n+1):gi])} else {r11 <- max(r1[(gi-n+1):gi])}
      # Since the rate ratio of the second arm compared to its own t0 can be
      # given in the argument rt, this part of the function has a lot of extra
      # ifs to find the right parameters... Note that this "rt" is especially
      # relevant in synthetic lethality experimental setups
      if (length(rt) > 1) {
        if (median(rt[(gi-n+1):gi]) < 0) {
          it <- which.min(rt[(gi-n+1):gi]+nse*set[(gi-n+1):gi])
          rtt <- rt[(gi-n+it)]
        } else {
          it <- which.max(rt[(gi-n+1):gi]-nse*set[(gi-n+1):gi])
          rtt <- rt[(gi-n+it)]
        }
        
        #if (median(rt[(gi-n+1):gi]) < 0) {rtt <- min(rt[(gi-n+1):gi])} else {rtt <- max(rt[(gi-n+1):gi])}
        r <- c(abs(r00), abs(r11), abs(rtt))
        w <- which.max(r)
        if (w == 1) {i <- i0}
        else if (w == 2) {i <- i1}
        else {i <- it}
      } else {
        r <- c(abs(r00),abs(r11))
        w <- which.max(r)
        if (w == 1) {i <- i0}
        else {i <- i1}
      }
      
      d <- r0[gi-n+i]/a
      if (length(rt) > 1) {de <- rt[gi-n+i]/b}
      # Code below shows the correction for different number of potential
      # doublings. Although this gets the estimated e-values much closer to
      # the real e-values, I see there is a skew upwards when d decreases
      # (e.g. straight lethality increases).
      if (a != b && correctab == TRUE) {
        e <- (r1[gi-n+i]-r0[gi-n+i]*(b/a-1))/b
      } else {
        e <- r1[gi-n+i]/b
      }
      
      if (r00 < 0) {
        q00 <- max((-rnulla-r0[(gi-n+1):gi])/se0[(gi-n+1):gi])
      } else {
        q00 <- max((r0[(gi-n+1):gi]-rnulla)/se0[(gi-n+1):gi])
      }
      
      if (hnull == 0) {
        pd <- 2*pnorm(-q00)
      } else {
        pd <- pnorm(-q00)
      }
      
      rnullc <- min(rnulla, rnullb)
      if (r11 < 0) {
        q11 <- max((-rnullc-r1[(gi-n+1):gi])/se1[(gi-n+1):gi])
      } else {
        q11 <- max((r1[(gi-n+1):gi]-rnullc)/se1[(gi-n+1):gi])
      }
      
      if (hnull == 0) {
        pe <- 2*pnorm(-q11)
      } else {
        pe <- pnorm(-q11)
      }
      
      if (length(rt) > 1) {
        if (rtt < 0) {
          qtt <- max((-rnullb-rt[(gi-n+1):gi])/set[(gi-n+1):gi])
        } else {
          qtt <- max((rt[(gi-n+1):gi]-rnullb)/set[(gi-n+1):gi])
        }
        
        if (hnull == 0) {
          pde <- 2*pnorm(-qtt)
        } else {
          pde <- pnorm(-qtt)  
        }
        
      }
      
      g <- sapply(seq_len(n), function(i) {
        
        if (w == 1) {
          if (d == 0 && r0[gi-n+i] == 0) {g <- 1} else {g <- (2^(r0[gi-n+i])-1) / (2^(a*d)-1)}
          g[g < 0] <- 0
          g[g > 1] <- 1
          
        } else if (w == 2) {
          if (e == 0 && r1[gi-n+i] == 0) {g <- 1} else {g <- (2^(r1[gi-n+i])-1) / (2^(b*e)-1)}
          g[g < 0] <- 0
          g[g > 1] <- 1
          
        } else {
          if (de == 0 && rt[gi-n+i] == 0) {g <- 1} else {g <- (2^(rt[gi-n+i])-1) / (2^(b*de)-1)}
          g[g < 0] <- 0
          g[g > 1] <- 1
          
        }
        
        return(g)})
      
      if (length(rt) > 1) {
        return(list(d=d,e=e,de=de,g=g,i=i,pd=pd,pe=pe,pde=pde))
      } else {
        return(list(d=d,e=e,g=g,i=i,pd=pd,pe=pe))
      }
    }, n, gi)
    if (secondbest == TRUE) {
      if (length(rt) > 1) {
        i <- unlist(degip[5,])
      } else {
        i <- unlist(degip[4,])  
      }
      n2 <- n-1
      gi2 <- cumsum(n2)
      r02 <- r0[-(gi-n+i)]
      r12 <- r1[-(gi-n+i)]
      if (length(rt) > 1) {rt2 <- rt[-(gi-n+i)]}
      d2e2j <- mapply(function(n2, gi2) {
        
        if (n2 == 0) {
          if (length(rt) > 1) {
            return(list(d2 = NaN, e2 = NaN, de2 = NaN, j = NaN))
          } else {
            return(list(d2 = NaN, e2 = NaN, j = NaN))
          }
        }
        
        else {
          if (median(r02[(gi2-n2+1):gi2]) < 0) {r2 <- min(r02[(gi2-n2+1):gi2])} else {r2 <- max(r02[(gi2-n2+1):gi2])}
          if (median(r12[(gi2-n2+1):gi2]) < 0) {r112 <- min(r12[(gi2-n2+1):gi2])} else {r112 <- max(r12[(gi2-n2+1):gi2])}
          if (length(rt) > 1) {
            if (median(rt2[(gi2-n2+1):gi2]) < 0) {rtt2 <- min(rt2[(gi2-n2+1):gi2])} else {rtt2 <- max(rt2[(gi2-n2+1):gi2])}
            r <- c(abs(r2), abs(r112), abs(rtt2))
            w <- which.max(r)
            if (w == 1) {j <- tail(which(r02[(gi2-n2+1):gi2]==r2),1)}
            else if (w == 2) {j <- tail(which(r12[(gi2-n2+1):gi2]==r112),1)}
            else {j <- tail(which(rt2[(gi2-n2+1):gi2]==rtt2),1)}
          } else {
            r <- c(abs(r2),abs(r112))
            w <- which.max(r)
            if (w == 1) {j <- tail(which(r02[(gi2-n2+1):gi2]==r2),1)}
            else {j <- tail(which(r12[(gi2-n2+1):gi2]==r112),1)}
          }
          
          if (a != b && correctab == TRUE) {
            e2 <- (r12[gi2-n2+j]-r02[gi2-n2+j]*(b/a-1))/b
          } else {e2 <- r12[gi2-n2+j]/b}
          
          if (length(rt) > 1) {
            return(list(d2=r02[gi2-n2+j]/a, e2 = e2, de2 = rt2[gi2-n2+j]/b, j = j))
          } else {
            return(list(d2=r02[gi2-n2+j]/a, e2 = e2, j = j))
          }
          
        }
      }, n2, gi2)
      
      if (length(rt) > 1) {
        j <- unlist(d2e2j[4,])
      } else {
        j <- unlist(d2e2j[3,])
      }
      k <- i[!is.na(j)]
      l <- j[!is.na(j)]
      l[l >= k] <- l[l >= k]+1
      j[!is.na(j)] <- l
      
    }
  }
  
  if (missing(r1)) {
    d <- unlist(dgip[1,])
    g <- unlist(dgip[2,])
    i <- unlist(dgip[3,])
    pd <- unlist(dgip[4,])
    if (secondbest == TRUE) {
      d2 <- unlist(d2j[1,])
      return(list(genes=genes,n=n,d=d,d2=d2,g=g,i=i,j=j,pd=pd))
    } else {
      return(list(genes=genes,n=n,d=d,g=g,i=i,pd=pd))
    }
  } else {
    d <- unlist(degip[1,])
    e <- unlist(degip[2,])
    if (length(rt) > 1) {
      de <- unlist(degip[3,])
      g <- unlist(degip[4,])
      i <- unlist(degip[5,])
      pd <- unlist(degip[6,])
      pe <- unlist(degip[7,])
      pde <- unlist(degip[8,])
    } else {
      g <- unlist(degip[3,])
      i <- unlist(degip[4,])
      pd <- unlist(degip[5,])
      pe <- unlist(degip[6,])
    }
    
    if (secondbest == TRUE) {
      d2 <- unlist(d2e2j[1,])
      e2 <- unlist(d2e2j[2,])
      if (length(rt) > 1) {
        de2 <- unlist(d2e2j[3,])
        return(list(genes=genes,n=n,d=d,d2=d2,e=e,e2=e2,de=de,de2=de2,g=g,i=i,j=j,pd=pd,pe=pe,pde=pde))
      } else {
        return(list(genes=genes,n=n,d=d,d2=d2,e=e,e2=e2,g=g,i=i,j=j,pd=pd,pe=pe))
      }
      
    } else {
      if (length(rt) > 1) {
        return(list(genes=genes,n=n,d=d,e=e,de=de,g=g,i=i,pd=pd,pe=pe,pde=pde))
      } else {
        return(list(genes=genes,n=n,d=d,e=e,g=g,i=i,pd=pd,pe=pe))
      }
      
    }
  }
}


evenlydividecounts <- function(counts, n) {
  m <- sapply(counts, function(x) {
    b <- c()
    for (r in seq_len(n)) {
      b <- append(b, floor(x/n))
      if (x%%n >= r) {b[r] <- b[r]+1}
    }
    return(b)
  })
  return(t(m))
}
