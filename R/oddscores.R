# This document has functions that are geared towards the use of odds instead of
# effects. This may be more appropriate for screens that are based on selection,
# such as FACS-based screens, screens with migration or invasion assays, and
# what have you. The use of distributions is quite different from the summed Z
# approach or the growth effect derivation, therefore my hope is that it is more
# robust for these selection screens, in which rate ratios may vary much more
# widely.

# Having gone this far, I also want a simulator for a selection screen. This is
# almost a copy of CRISPRsim, and actually easier in a lot of respects. The idea
# is as follows. Instead of a growth modifier, each successful knockout now
# changes the odds of being selected (by a factor d). The screen can also be
# done with or without "drug", which represents another modifier (e). Abundance
# and guide efficacy are present as before. This function is more flexible in
# selecting the expected strength and abundance of hits. Depending on the
# baseline probability of being selected, hits are divided. E.g., if baseprob is
# 0.1, 9 out of 10 hits will be positive, and vice versa. The average strength
# of hits is expressed as odds. The final step is just picking cells based on
# their odds!

#' Simulate a selection-based CRISPR-Cas9 pooled screen
#'
#' sortingsim simulates a selection-based CRISPR-Cas9 pooled screen with
#' user-defined parameters. In a way, this is a simplified version of
#' \code{\link{CRISPRsim}}, because there is no simulation of growth and
#' passaging. It is about entering into an assay with a number of cells, and
#' collecting the selected and unselected cells for sequencing. This was mainly
#' created with FACS-based screens in mind, but I suppose one can think of other
#' selections as well. Drug treatment is still an option! Again, "infected"
#' cells either have a successful knockout (with associated effect) or not,
#' based on guide efficacy. Cells are then sorted, which in this simulator means
#' that all cells with a certain knockout have a base probability to be
#' positively selected, modified by a gene-specific score. Note that this
#' function takes the number of sorted cells as input, of which only a fraction
#' is selected. This fraction is roughly equal to the \code{baseprob}. The
#' output of this function is a data frame that contains the guide-relevant
#' parameters and the sequencing coverage per guide for the selected and
#' unselected arms. Simulated screens will aid researchers with their
#' experimental setup. Furthermore, this offers a unique platform for the
#' evaluation of analysis methods for sorting-based pooled gene knockout
#' screens.
#'
#' @param genes Single integer or character vector. Specify how many or which
#'   genes to include in the experiment respectively. Not required when a full
#'   list of guides is given.
#' @param guides Single integer, integer vector or character vector. In case of
#'   single integer, specify by how many guides each gene is represented. In
#'   case of an integer vector, specify per gene by how many guides it is
#'   represented. In case of a character vector, guides are assumed to contain a
#'   gene name, followed by an underscore, followed by an identifier within that
#'   gene (e.g. a number or a nucleotide sequence).
#' @param g Integer vector. Specify guide efficacies per guide. If omitted,
#'   guide efficacies will be sampled from a representative distribution.
#' @param f Integer vector. Specify guide abundance at time of infection per
#'   guide. If omitted, guide abundance will be sampled from a representative
#'   distribution.
#' @param d Integer vector. Specify gene-specific modifier of probability of
#'   selection. This is a modifier of the odds. If omitted, effect of gene
#'   knockout will be sampled from three distributions, depending on base
#'   probability and hit factor. If the length of the vector does not match the
#'   number of genes, values will be randomly sampled from the specified
#'   distribution!
#' @param e Integer vector. Specify treatment-specific selection effect per
#'   gene. If omitted, effects will be sampled from a representative
#'   distribution. If the length of the vector does not match the number of
#'   genes, values will be randomly sampled from the specified distribution!
#' @param baseprob Numeric. Baseline probability of selection. Needs to be
#'   larger than 0 and smaller than 1. Default = 0.1
#' @param hitfraction Numeric. Fraction of genes that affect selection
#'   significantly. Default = 1/200
#' @param hitsup Numeric. Fraction of hits of which the selection probability is
#'   multiplied by \code{hitfactor}. Probability of the other (down) hits are
#'   divided by \code{hitfactor}. Defaults to \code{1-baseprob}
#' @param hitfactor Numeric. Multiplication factor with which hits affect
#'   selection probability on average. Default = 10
#' @param efraction Numeric. Fraction of genes that affects selection
#'   specifically in this treatment arm. Defaults to \code{hitfraction}
#' @param eup Numeric. Same as hitsup, but now relating to treatment effects.
#'   Defaults to 0.5
#' @param efactor Numeric. Multiplication factor for treatment-specific effects.
#'   Defaults to \code{hitfactor}
#' @param sortedcells Integer. Number of cells put through the simulated
#'   selection. Note that this is the sum of the selected and unselected cells!
#' @param seqdepth Integer. Specify the amount of sequencing reads devoted to
#'   each experimental arm. If omitted, depth will default to 500 times the
#'   number of guides
#' @param offtargets Logical or numeric. Specify the fraction of off-targets. If
#'   TRUE, 1 in 1000 guides (0.001) will target a different gene. Default =
#'   FALSE
#' @param allseed Integer. All unspecified seeds default to this plus an
#'   increment of 1 for each different seed. Defaults to NULL, in which case the
#'   unspecified seeds are randomly generated. Default = NULL 
#' @param gseed Integer. Specify seed for guide effiency assignment
#' @param fseed Integer. Specify seed for infectious units assignment, which
#'   dictates a guide's abundance at the start of the experiment
#' @param dseed Integer. Specify seed for straight lethality assignment of genes
#' @param eseed Integer. Specify seed for sensitizer assignment of genes
#' @param oseed Integer. Specify seed for off-target selection
#' @param t0seed Integer. Specify seed for t0, which encompasses assignment of
#'   successful knockout cells versus no knockout cells for each guide
#' @param repseed Integer. Specify the seed after t0
#' @param perfectsampling Logical. If TRUE, all sampling steps are replaced by
#'   simple equations to calculate representation of guides. Useful as null
#'   control to isolate the effect of sampling. Default = FALSE
#' @param perfectseq Logical. If TRUE, sequencing results are a perfect
#'   representation (though still rounded) of guides in the harvested cells.
#'   Applicable to speed up simulations, assuming sequencing is sufficiently
#'   deep. Defaults to \code{perfectsampling}
#' @param returnall Logical. If TRUE, function returns a list with the simulated
#'   data in the guidesdf, summary per gene in the genesdf, and parameters.
#'   Default = FALSE
#' @param outputfile Character string. When used, returned data frame will be
#'   saved as a tab-delimited text to the specified file path
#'
#' @details sortingsim performs a genome-wide (or subsetted) pooled CRISPR
#'   knockout screen ending a binary selection. Perhaps even more so than
#'   growth-based screens, the outcome of such screens can be a massive black
#'   box. As of yet, no specific analysis methods have been published for these
#'   kind of screens, but the simulator below can help assess those. The
#'   parameters are highly customizable, so I sincerely recommend reading the
#'   documentation for all the options. And it is always possible to provide
#'   your own gene-specific selection modifiers or guide efficacies if you are
#'   not happy with the provided distributions. Seeds are relevant if you want
#'   to create replicate screens. You can easily "practice" by simulating some
#'   small experiments (i.e. limit the amount of genes). The basis of the
#'   simulation are as follows. Cells have an a priori probability
#'   \code{baseprob} to be selected. The corresponding odds
#'   \code{baseprop/(1-baseprop)} are multiplied by gene-specific modifier d
#'   (and optionally gene-specific modifier for treatment e). These odds are
#'   converted to the modified probability mod_prob, which is used to determine
#'   how many cells with a specific knockout are selected. Each guide has an
#'   efficacy, which is the chance to create a successful knockout. Only in case
#'   of successful knockout are the modifiers applied. Selected and not-selected
#'   cells are separately sequenced, both to the indicated sequencing depth.
#'
#' @return Returns a data frame with every row representing a single guide.
#'   Contains the pertinent parameters of each guide and the number of
#'   sequencing reads of selected and not selected cells. If the argument
#'   returnall is set to TRUE, the function also returns a data frame with the
#'   true values for the genes, and lists all parameters as well.
#'
#' @note If you specify an inverted hitfactor (e.g. 0.1 instead of 10), your
#'   hits are turned around.
#'
#'   While it also makes sense to be able to specify how many cells are
#'   positively selected (this could be your FACS setup of course), this is not
#'   directly compatible with this simulator. Instead, you can divide the number
#'   of cells you want with \code{baseprob} and use that as input for
#'   \code{sortedcells}. If that does not come close (some wonky parameters
#'   perhaps), or you want to be more precise, you can do a test run with
#'   argument returnall. One of the returned values is selectedcells, which
#'   corresponds to the number of cells used as input for sequencing of the
#'   selected arm. It follows that noteselectedcells equals sortedcells minus
#'   selectedcells
#'
#' @author Jos B. Poell
#'
#' @seealso \code{\link{oddscores}}, \code{\link{CRISPRsim}}
#'
#' @examples
#' sortdf <- sortingsim(18000, 4, e = TRUE, perfectsampling = TRUE)
#' d <- rle(sortdf$d)$values
#' lod <- log(d)
#' e <- rle(sortdf$e)$values
#' loe <- log(e)
#' plot(lod, loe, main = "log odds of selection")
#' enrichment <- log(sortdf$selected+1)-log(sortdf$notselected+1)
#' kocell_logodds <- log(sortdf$mod_prob)-log(1-sortdf$mod_prob)
#' plot(kocell_logodds, enrichment, pch = 16,
#'      cex = 0.75, col = rgb(sortdf$g, 0, 1-sortdf$g))
#'
#' @export

sortingsim <- function(genes, guides, g, f, d, e, baseprob = 0.1, 
                       hitfraction = 1/200, hitsup, hitfactor = 10, 
                       efraction, eup, efactor, sortedcells, 
                       seqdepth, offtargets = FALSE, allseed = NULL, 
                       gseed, fseed, dseed, eseed, oseed, t0seed, 
                       repseed, perfectsampling = FALSE, 
                       perfectseq, returnall = FALSE, 
                       outputfile) {
  
  if(missing(perfectseq)) {perfectseq <- perfectsampling}
  baseodds <- baseprob/(1-baseprob)
  
  if (missing(genes)) {
    if (!missing(guides) && is(guides, "character")) {
      genes <- rle(gsub("_.*", "", guides))$values
      n <- rle(gsub("_.*", "", guides))$lengths
    } else {
      stop("please enter the number of genes you wish to simulate a CRISPR screen for, or provide a vector with gene names")
    }
  }
  if (missing(guides)) {
    stop("please enter the number of guides per gene, a vector with all guides, or a vector with the number of guides for each gene in genes")
  }
  if (length(genes)==1 && (is(genes, "integer") || is(genes, "numeric"))) {
    genes <- paste0("gene", seq_len(genes))
  }
  
  if (length(guides)==1 && (is(guides, "integer") || is(guides, "numeric"))) {
    n <- rep(guides, length(genes))
    guides <- as.vector(sapply(genes, function(x) {paste0(x, "_", seq_len(guides))}))
  } else if (length(guides)==length(genes) && (is(guides, "integer") || is(guides, "numeric"))) {
    n <- guides
    guides <- unlist(mapply(function(guides, genes) {
      paste0(genes, "_", seq_len(guides))
    }, guides, genes))
  } else {
    n <- rle(gsub("_.*", "", guides))$lengths
    if (any(genes != rle(gsub("_.*", "", guides))$values)) {
      warning("guides and genes do not seem to match; using gene list inferred from guides")
      genes <- rle(gsub("_.*", "", guides))$values
    }
  }
  
  if (missing(gseed)) {
    if (is.null(allseed)) {
      gseed <- sample.int(.Machine$integer.max, 1L)
    } else {
      gseed <- allseed
    }
  }
  with_seed(gseed, {
    if (missing(g)) {
      message("guide efficiency is randomly sampled from an empirical distribution")
      bg <- floor(length(guides)/40)
      g <- sample(c(1-runif(length(guides)-bg)^2.5, (runif(bg)/10)^5))
    } else if (length(g)==1) {
      g <- rep(g, length(guides))
    } else if (length(g) != length(guides)) {
      g <- sample(g, length(guides), replace = TRUE)
    } 
  })
  
  if (any(g < 0) || any(g > 1)) {
    warning("g represents guide efficacy and should be a number between 0 and 1; values have been constrained")
    g[g < 0] <- 0
    g[g > 1] <- 1
  }
  
  if (missing(fseed)) {
    if (is.null(allseed)) {
      fseed <- sample.int(.Machine$integer.max, 1L)
    } else {
      fseed <- allseed+1
    }
  }
  with_seed(fseed, {
    if (missing(f)) {
      message("guide abundance (library representation) is randomly assigned")
      # The allseed is shifted to ascertain randomness between g and f
      
      f <- 2^rnorm(length(guides))
    } else if (length(f) == 1) {
      f <- rep(1, length(guides))
      message("assuming equal abundance for all guides")
    } else if (length(f) != length(guides)) {
      f <- sample(f, length(guides), replace = TRUE)
    }
  })
  
  if (missing(dseed)) {
    if (is.null(allseed)) {
      dseed <- sample.int(.Machine$integer.max, 1L)
    } else {
      dseed <- allseed+2
    }
  }
  with_seed(dseed, {
    if (missing(d) || (length(d) == 1 && d == TRUE)) {
      message("effect of gene knockout is randomly sampled from an empirical distribution")
      if (missing(hitsup)) {hitsup <- 1-baseprob}
      oddsup <- floor(length(genes)*hitfraction*hitsup)
      oddsdown <- floor(length(genes)*hitfraction*(1-hitsup))
      lod <- sample(c(rnorm(oddsup, log(hitfactor), abs(log(hitfactor)/3)), 
                      rnorm(oddsdown, -log(hitfactor), abs(log(hitfactor)/3)),
                      rnorm(length(genes)-oddsup-oddsdown, 0, abs(log(hitfactor)/6))))
      d <- exp(lod)
    } else if (length(d)==1) {
      message("assuming equal d for all genes. 
              Note that you can enter custom d-values by providing a vector with length equal to the number of genes, 
              or a vector with unequal size from which will be sampled")
      d <- rep(d, length(genes))
      
    } else if (length(d) != length(genes)) {
      d <- sample(d, length(genes), replace = TRUE)
    }
    unique_d <- d
    d <- rep(d, n)
  })
  
  
  if (any(is.na(d)) || any(is.null(d)) || any(is.character(d))) {
    stop("This is not going to go well, not all d-values are numeric!")
  }
  
  if (!missing(e)) {
    if (missing(eseed)) {
      if (is.null(allseed)) {
        eseed <- sample.int(.Machine$integer.max, 1L)
      } else {
        eseed <- allseed+3
      }
    }
    if (round(dseed) == round(eseed)) {
      warning("d and e sampled using identical seeds: expect weirdly correlated results")
    }
    with_seed(eseed, {
      if (length(e) == 1 && e == TRUE) {
        message("arm-specific effect modifier is randomly sampled from a distribution typical for a dropout screen")
        # The allseed should not be exactly the same between d and e, because it
        # would sample exactly the same order and create huge biases. Therefore I
        # add 1 to the seed!
        
        if (missing(efraction)) {efraction <- hitfraction}
        if (missing(eup)) {eup <- 0.5}
        if (missing(efactor)) {efactor <- hitfactor}
        oddsup <- floor(length(genes)*efraction*eup)
        oddsdown <- floor(length(genes)*efraction*(1-eup))
        loe <- sample(c(rnorm(oddsup, log(efactor), abs(log(efactor)/3)), 
                        rnorm(oddsdown, -log(efactor), abs(log(efactor)/3)),
                        rnorm(length(genes)-oddsup-oddsdown, 0, abs(log(efactor)/6))))
        # Note that e is kept as odds, since it is a modifier of d
        e <- exp(loe)
        
      } else if (length(e) != length(genes)) {
        
        e <- sample(e, length(genes), replace = TRUE)
      }
      unique_e <- e
      e <- rep(e, n)
      if (any(is.na(e)) || any(is.null(e)) || any(is.character(e))) {
        stop("This is not going to go well, not all e-values are numeric!")
      }
    })
  }
  
  offtarget_guides <- rep(0, length(genes))
  
  if (missing(oseed)) {
    if (is.null(allseed)) {
      oseed <- sample.int(.Machine$integer.max, 1L)
    } else {
      oseed <- allseed+4
    }
  }
  if (offtargets != FALSE) {
    with_seed(oseed, {
      if (offtargets == TRUE) {offtargets <- 0.001}
      tochange <- which(runif(length(guides)) < offtargets)
      changeto <- round(runif(length(tochange))*length(genes)+0.5)
      d[tochange] <- unique_d[changeto]
      if (!missing(e)) {
        e[tochange] <- unique_e[changeto]
      }
      off_per_gene <- sapply(tochange, function(x) {
        min(which(x <= cumsum(n)))
      })
      for (o in off_per_gene) {
        offtarget_guides[o] <- offtarget_guides[o]+1
      }
    })
  }
  
  # mod_prob is the probability of selection, potentially modified by e
  if(!missing(e)) {
    mod_prob <- (baseodds*d*e)/(1+baseodds*d*e)
    unique_mod_prob <- (baseodds*unique_d*unique_e)/(1+baseodds*unique_d*unique_e)
  } else {
    mod_prob <- (baseodds*d)/(1+baseodds*d)
    unique_mod_prob <- (baseodds*unique_d)/(1+baseodds*unique_d)
  }
  
  if (!missing(e)) {
    gdf <- data.frame(genes, n, offtarget_guides, d = unique_d, e = unique_e, mod_prob = unique_mod_prob)
  } else {
    gdf <- data.frame(genes, n, offtarget_guides, d = unique_d, mod_prob = unique_mod_prob)
  }
  
  if (missing(sortedcells)) {
    message("assuming mean representation of 200 cells per guide")
    sortedcells <- 200*length(guides)
  } 
  
  if (missing(seqdepth)) {
    message("assuming mean sequencing depth of 500 reads per guide")
    seqdepth <- 500*length(guides)
  }
  
  
  if (missing(t0seed)) {
    if (is.null(allseed)) {
      t0seed <- sample.int(.Machine$integer.max, 1L)
    } else {
      t0seed <- allseed+5
    }
  }
  if (missing(repseed)) {
    if (is.null(allseed)) {
      repseed <- sample.int(.Machine$integer.max, 1L)
    } else {
      repseed <- allseed+6
    }
  }
  
  if (perfectsampling == TRUE) {
    cells <- round(c(f*g*sortedcells, f*(1-g)*sortedcells)/sum(f))
    kocells <- cells[seq_along(guides)]
    nkocells <- cells[seq(length(guides)+1, 2*length(guides))]
    koselected <- round(kocells*mod_prob)
    nkoselected <- round(nkocells*baseprob)
    konotselected <- kocells - koselected
    nkonotselected <- nkocells - nkoselected
    selected <- koselected+nkoselected
    notselected <- konotselected+nkonotselected
    if (perfectseq == TRUE) {
      readsselected <- round(selected*seqdepth/sum(selected))
      readsnotselected <- round(notselected*seqdepth/sum(notselected))
    } else {
      with_seed(repseed, {
        readsselected <- tabulate(c(seq_along(guides), 
                                    sample(seq_along(guides),
                                           seqdepth, replace = TRUE, 
                                           prob = selected)))-1
        readsnotselected <- tabulate(c(seq_along(guides), 
                                       sample(seq_along(guides),
                                              seqdepth, replace = TRUE, 
                                              prob = notselected)))-1
      })
    }
  } else {
    with_seed(t0seed, {
      cells <- tabulate(c(seq_along(c(guides,guides)),
                          sample(seq_along(c(guides,guides)), sortedcells,
                                 replace = TRUE, prob = c(f*g, f*(1-g)))))-1
      kocells <- cells[seq_along(guides)]
      nkocells <- cells[seq(length(guides)+1, 2*length(guides))]
    })
    with_seed(repseed, {
      koselected <- mapply(function(cells, chance) {rbinom(1,cells,chance)},kocells,mod_prob)
      nkoselected <- mapply(function(cells, chance) {rbinom(1,cells,chance)},nkocells,baseprob)
      konotselected <- kocells - koselected
      nkonotselected <- nkocells - nkoselected
      selected <- koselected+nkoselected
      notselected <- konotselected+nkonotselected
      if (perfectseq == TRUE) {
        readsselected <- round(selected*seqdepth/sum(selected))
        readsnotselected <- round(notselected*seqdepth/sum(notselected))
      } else {
        readsselected <- tabulate(c(seq_along(guides), 
                                    sample(seq_along(guides),
                                           seqdepth, replace = TRUE, 
                                           prob = selected)))-1
        readsnotselected <- tabulate(c(seq_along(guides), 
                                       sample(seq_along(guides),
                                              seqdepth, replace = TRUE, 
                                              prob = notselected)))-1
      }
    })
  }
  
  message("wrapping up")
  
  if (missing(e)) {
    df <- data.frame(guides, gene = rep(genes, n), d, mod_prob, f, g, selected = readsselected, 
                     notselected = readsnotselected, stringsAsFactors = FALSE)
  } else {
    df <- data.frame(guides, gene = rep(genes, n), d, e, mod_prob, f, g, selected = readsselected, 
                     notselected = readsnotselected, stringsAsFactors = FALSE)
  }
  
  if (!missing(outputfile)) {
    write.table(df, outputfile, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  if (returnall == TRUE && missing(e)) {
    return(list(
      guidesdf = df,
      genesdf = gdf,
      baseprob = baseprob,
      hitfactor = hitfactor,
      hitsup = hitsup,
      hitfraction = hitfraction,
      offtargets = offtargets,
      sortedcells = sortedcells,
      selectedcells = sum(selected),
      seqdepth = seqdepth,
      perfectsampling = perfectsampling,
      perfectseq = perfectseq,
      allseed = allseed,
      gseed = gseed,
      fseed = fseed,
      dseed = dseed,
      oseed = oseed,
      t0seed = t0seed,
      repseed = repseed
    ))
  } else if (returnall == TRUE && !missing(e)) {
    return(list(
      guidesdf = df,
      genesdf = gdf,
      baseprob = baseprob,
      hitfactor = hitfactor,
      hitsup = hitsup,
      hitfraction = hitfraction,
      efactor = efactor,
      eup = eup,
      efraction = efraction,
      offtargets = offtargets,
      sortedcells = sortedcells, 
      selectedcells = sum(selected),
      seqdepth = seqdepth,
      perfectsampling = perfectsampling,
      perfectseq = perfectseq,
      allseed = allseed,
      gseed = gseed,
      fseed = fseed,
      dseed = dseed,
      eseed = eseed,
      oseed = oseed,
      t0seed = t0seed,
      repseed = repseed
    ))
  } else {
    return(df) 
  }
  
}

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
#'   the same length as the number of rate ratios provided. These are the odds
#'   that the rate ratio is lower than the one observed based on the (control)
#'   population.
#'
#' @note Odds in this function represent the odds that a feature has a rate
#'   ratio that is this low or lower, either relative to the entire data set or to
#'   the normalization subset when provided. These odds are mostly based on
#'   rank. Only if they are more extreme than the most extreme control are they
#'   extrapolated. Even then, extrapolation is quite conservative. There is also
#'   a small interpolation step, which looks at how far the feature is away from
#'   the flanking control features.
#'
#' @author Jos B. Poell
#'
#' @seealso \code{\link{geneodds}}, \code{\link{odds2pq}},
#'   \code{\link{sortingsim}}
#'
#' @export

oddscores <- function(r, normsubset, log = 10) {
  if (missing(normsubset)) {
    odds <- rank(r)/(length(r)-rank(r)+1)
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
#'   and their value is expected to start with the gene name, number or symbol,
#'   followed by an underscore. Odds are combined by multiplying odds or adding
#'   log odds, after which the odds are corrected for being composed of multiple
#'   odds. This correction changes the meaning of the odds: the returned result
#'   does not signify the odds that this specific gene has a lower rate ratio,
#'   but the odds that any gene with this many guide has a lower rate ratio. The
#'   code in the example shows an example with uniformly distributed p-values.
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
#'   low odds, which is what this analysis is meant for. If you require precise
#'   answers or are interested in p-values of weak odds, use non-transformed
#'   odds and specify log = FALSE.
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
    # ones we are going to be interested in anyway. If you are interested in
    # precise answers, use non-transformed odds and specify log = FALSE.
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
    p <- ifelse(odds <= 1, (2*odds)/(1+odds), 2-(2*odds)/(1+odds))
    q <- p.adjust(p, method = "BH")
  }
  return(data.frame(odds=odds,p=p,q=q))
}


# The geteffect functions are a pair of functions that attempt to estimate the
# effect of gene knockout on cell growth by calculating the odds that the gene
# has a lower fold change than observed given a range of potential effect sizes.
# The odds of a lower fold change are calculated for each guide targeting a
# gene, and are then combined. The expected fold changes for certain effect
# sizes are calculated for a range of guide efficacies, which should represent
# the guide efficacies in the screen! At the time of writing I do not have a
# specific function in this package yet. What I have done before, is to
# calculate the guide efficacy of the guides targeting essential genes compared
# to the best guide of that gene, and then exclude the best guides. This
# actually works really well, although the representation of very good guides
# will be underestimated. Your library is better than you think! To do this
# after a CRISPR screen, use getdeg, and check out the estimated g of the
# essential genes. Moving on. Another important consideration is that the odds
# that are calculated for a gene at a specific effect size, are corrected for
# the probability of any gene having that effect size! I call this the empirical
# bayesian correction (is it though?). Here I take a little shortcut. To
# estimate the effect sizes, I take the most extreme fold change of each gene,
# and calculate what the effect size would be given the provided number of cell
# doublings and assuming a perfect guide. Keep in mind, that the fact that I
# assume a perfect guide, actually means that these effect sizes are smaller
# than they really are. Therefore, strong effect sizes are a bit overcorrected.
# I believe the correction is crucial: if you have a mildly lethal gene that
# happens to be targeted by mostly very good guides (yes, that will happen in a
# genome-wide screen), the effect will be greatly overestimated without the
# correction.

#' Approximate effect of gene knockout on growth
#'
#' geteffect functions calculate the combined likelihood that the fold changes
#' observed for all guides of a gene between two time points or conditions is
#' smaller than expected at certain effect size(s).
#'
#' @param guides Character vector. Guides are assumed to start with the gene
#'   name, followed by an underscore, followed by a number or sequence unique
#'   within that gene.
#' @param r Numeric vector. Log2-transformed rate ratios
#' @param rse Numeric vector. Standard error of the rate ratios
#' @param t1 Integer vector. Read counts of the test arm
#' @param t0 Integer vector. Read counts of the control arm
#' @param normfun Character string. Specify with which function to standardize
#'   the data. Default = "sum"
#' @param normsubset Integer vector. Specify the indices of features that
#'   function as controls. If omitted, all features are used.
#' @param a Numeric. Estimated potential population doublings between time
#'   points.
#' @param g Numeric vector. Specify representative guide efficacies. If omitted,
#'   an exponentially decreasing set of guide efficacies with length \code{gl}
#'   will be created.
#' @param gw Numeric vector. Specify the relative prevalence of each guide
#'   efficacy provided with \code{g}. If omitted, a top-heavy distribution will
#'   be created to roughly represent guide efficacy distribution in
#'   \code{\link{CRISPRsim}}
#' @param gl Numeric. Number of guide efficacies to test. Only used when
#'   \code{g} is not specified, in which case it needs to be at least 4. Default
#'   = 11
#' @param subset Integer vector. Specify the indices of features for which to
#'   return output. All features will be used for empirical Bayes correction.
#' @param effectrange Numeric vector. Sequence of effect values to test. If
#'   omitted, a range will be calculated based on the highest and lowest rate
#'   ratios in the data set.
#' @param output Character string. Specify which output to generate. Can be
#'   either "range", "exact", or "both". Default = "range"
#' @param exactci Logical or numeric. Specify the confidence interval of the
#'   exact effect as a fraction between 0 and 1. Only applicable if exact values
#'   are calculated. If FALSE, confidence interval is omitted. Default = FALSE
#' @param semiexact Logical. If TRUE, exact effect values are estimated based on
#'   the likelihoods of effect values flanking a log likelihood of 0 (or the
#'   respective log likelihood of the confidence interval). Note that a finer
#'   resolution of effect values than provided by default are recommended when
#'   applying this approach. Default = FALSE
#' @param ebcfun Character string or logical. Specify the function that is used
#'   to select gene-wise rate ratios for the empirical Bayes correction.
#'   Suggested options are "max", "median", and "mean". Note that "max" selects
#'   the most extreme rate ratio, not the highest. If FALSE, no empirical Bayes
#'   correction is performed, which will generally lead to overestimation of
#'   negative effects. Default = "max"
#' @param minprob Numeric. Cutoff point for the lowest considered likelihood.
#'   Default = 10^-20
#'
#' @details geteffect uses all guides targeting a gene, and incorporates the
#'   distribution of guide efficacies (given as parameters by the user) and
#'   effect values (interpolated from the data). The process is as follows. For
#'   a (predefined) effect value, the expected fold change is calculated at a
#'   range of different guide efficacies. Next, the probability that the
#'   observed fold change is smaller than the expected fold change is calculated
#'   at each guide efficacy. This number is multiplied with the density of each
#'   guide efficacy (i.e. the weight of each guide efficacy, specified by the
#'   user or by default similar to the distribution used in
#'   \code{\link{CRISPRsim}}). This yields the probability that the fold change
#'   observed for a specific guide is smaller than expected given the effect
#'   value. The probability is converted to a log10 odds. The log-odds are
#'   summed for each guide targeting the same gene, and the sum is divided by
#'   the square root of the number of guides (see notes). Finally, the log-odds
#'   are corrected for the prior probability of the effect value, which is
#'   interpolated from the data using the \code{ebcfun}. The best estimate of
#'   effect value for a gene is when the probability of being lower or higher
#'   than that effect value or equal, and therefore the log-odds are 0. Besides
#'   the estimates at the predefined effect values, closer estimations for each
#'   gene (or a subset of the data set) may be obtained in the form of exact or
#'   semiexact values.
#'
#'   Raw count data can be used as input for the \code{geteffect_c} function. If
#'   replicates are present, the function \code{geteffect_r} will take the
#'   precalculated log2 rate ratio (i.e. log2 fold change) and the accompanying
#'   standard error of each rate ratio as input. These values can be obtained by
#'   using the \code{\link{rrep}} function. These may also be derived from other
#'   packages that can analyze high-throughput count data.
#'
#' @return Returns the following (depending on input arguments): \itemize{
#'   \item{range}{ - data frame with log10 odds of effect per gene}
#'   \item{exact}{ - vector with effect estimates, or data frame with effect
#'   estimates including confidence interval} }
#'
#' @note While the calculations on the ranged output are fully vectorized, the
#'   exact output needs to be calculated one gene at a time. And three times
#'   when a confidence interval is included. I have noticed that computation
#'   time does not scale linearly after a certain number of genes, but takes
#'   even longer (this seemed to happen around 500 genes in my case). I guess
#'   that this might mean I ran out of memory, although I do not understand why
#'   this would happen. Therefore, I would not recommend running exact
#'   calculations for a whole genome data set. If you want "more exact" results
#'   than given by the default range, set \code{effectrange} yourself (e.g.
#'   \code{seq(-2, 1, by = 0.01)}). Use \code{semiexact == TRUE} to get closer
#'   approximations and unique effect values without crashing the computer.
#'
#'   Using \code{\link{CRISPRsim}} to validate the performance of this function,
#'   the exact effect values end up very close to the input effect values, but
#'   especially genes that lack an efficacious guide are underestimated. Due to
#'   this, the confidence intervals are a bit off. Specifically, the confidence
#'   limits are too tight on the "more extreme than" side, meaning a gene might
#'   actually have a higher probability of a stronger effect than the confidence
#'   interval would suggest.
#'
#'   A minimum likelihood was introduced with the \code{minprob} option to
#'   prevent Inf and -Inf results. This may prevent errors and help plotting.
#'
#'   The sum of the log-odds is divided by the square root of the number of
#'   guides. Therefore, it is not the odds of this particular gene having a more
#'   negative effect value than the tested value, but the odds of any gene with
#'   this many guides having a more negative effect value. The example code in
#'   \code{\link{geneodds}} visualizes the combination of odds in this fashion.
#'
#' @seealso \code{\link{rrep}}, \code{\link{CRISPRsim}}, \code{\link{odds2pq}}
#'
#' @author Jos B. Poell
#'
#' @examples
#' ut1 <- CRISPRsim(100, 4, a = c(3,3), allseed = 100, t0seed = 10,
#'                  repseed = 1, perfectseq = TRUE)
#' ut2 <- CRISPRsim(100, 4, a = c(3,3), allseed = 100, t0seed = 20,
#'                  repseed = 3, perfectseq = TRUE)
#' ut3 <- CRISPRsim(100, 4, a = c(3,3), allseed = 100, t0seed = 30,
#'                  repseed = 5, perfectseq = TRUE)
#' cgi <- ut1$d > -0.05 & ut1$d < 0.025
#' df <- data.frame(guides = ut1$guides,
#'                  T01 = ut1$t0, T02 = ut2$t0, T03 = ut3$t0,
#'                  UT1 = ut1$t6, UT2 = ut2$t6, UT3 = ut3$t6,
#'                  stringsAsFactors = FALSE)
#' r <- rrep(t1 = df[,5:7], t0 = df[,2:4], normsubset = cgi)
#' results <- geteffect_r(df$guides, r = r$r, rse = r$se, a = 6,
#'                        output = "both", semiexact = TRUE)
#' d <- rle(ut1$d)$values
#' plot(d, results$exact, xlab = "real effect",
#'      ylab = "estimated effect")
#'
#' @name geteffect
#' @aliases geteffect geteffect_c geteffect_r
NULL

#' @rdname geteffect
#' @export
geteffect_c <- function(guides, t1, t0, normfun = "sum", normsubset, 
                        a, g, gw, gl = 11, subset, effectrange, 
                        output = "range", exactci = FALSE, 
                        semiexact = FALSE, ebcfun = "max", 
                        minprob = 10^-20) {
  if (missing(a)) {
    stop("enter the presumed number of population doublings")
  }
  if (missing(g)) {
    if (gl < 4) {
      warning("gl should be at least 4: it is now set to 4")
      gl <- 4
    }
    g <- (11-11^((seq_len(gl)-1)/(gl-1)))/10
  }
  if (missing(gw)) {
    gw <- c(1,3,2,rep(1, length(g)-4),2)
  }
  if (length(g) != length(gw)) {
    warning("Length of g's and weights of g's is unequal")
    if (length(gw) < length(g)) {
      gw <- c(gw, rep(tail(gw, 1), length(g)-length(gw)))
    } else {
      gw <- gw[seq_along(g)]
    }
  }
  # Normalize gene weights so they sum to 1
  gw <- gw/sum(gw)
  
  r <- jar(t1, t0, 1)
  
  if (missing(effectrange)) {
    if (!missing(subset)) {
      effectrange <- round(seq(floor(30*min(r[subset])/a), ceiling(30*max(r[subset])/a))/20, 2)
    } else {
      effectrange <- round(seq(floor(30*min(r)/a), ceiling(30*max(r)/a))/20, 2)
    }
  } else {
    # effectrange needs to be sorted from low to high
    effectrange <- sort(effectrange)
  }
  
  # note that guides are presumed to be named [gene]_[number | sequence]
  genes <- rle(gsub("_.*", "", guides))$values
  nguides <- rle(gsub("_.*", "", guides))$lengths
  genei <- cumsum(nguides)
  
  fcmat <- matrix(nrow = length(effectrange), ncol = length(g))
  for (c in seq_len(ncol(fcmat))) {
    fcmat[,c] <- (g[c]*2^((1+effectrange)*a) + (1-g[c])*2^a)/2^a
  }
  
  fun <- get(normfun)
  if (!missing(normsubset)) {
    sr0 <- fun(t0[normsubset])
    sr1 <- fun(t1[normsubset])
  } else {
    sr0 <- fun(t0)
    sr1 <- fun(t1)
  }
  
  if (ebcfun != FALSE) {
    if (ebcfun == "max") {
      ebcr <- mapply(function(n, gi) {
        r[(gi-n+1):(gi)][which.max(abs(r[(gi-n+1):(gi)]))]
      }, nguides, genei)/a
    } else {
      fun <- get(ebcfun)
      ebcr <- mapply(function(n, gi) {
        fun(r[(gi-n+1):(gi)])
      }, nguides, genei)/a
    }
    effectodds <- oddscores(c(effectrange, ebcr), 
                            normsubset = seq_along(ebcr)+length(effectrange), 
                            log = 10)[seq_along(effectrange)]
  }
  
  if (!missing(subset)) {
    if (!is(subset, "numeric")) {
      subset <- which(gsub("_.*", "", guides) %in% subset)
    }
    guides <- guides[subset]
    genes <- rle(gsub("_.*", "", guides))$values
    nguides <- rle(gsub("_.*", "", guides))$lengths
    genei <- cumsum(nguides)
  }
  
  # I have noted I have to first bind the variables used in foreach...
  n <- gi <- gene <- i <- NULL
  
  if (output == "range" || output == "both") {
    # I imagine it is beneficial to perform the outer loop in parallel using
    # %dopar% if a parallel backend is available
    allodds <- foreach(n=nguides, gi=genei, .combine = cbind) %do% {
      odf <- foreach(i = seq_len(n), .combine = cbind) %do% {
        pmat <- 1-pbinom(t1[gi-n+i], t1[gi-n+i]+t0[gi-n+i], fcmat * sr1 / (fcmat * sr1 + sr0))
        p <- apply(pmat, 1, function(p) {sum(p*gw)})
        o <- log(p/(1-p), 10)
        o[o < log(minprob, 10)] <- log(minprob, 10) 
        o[o > -log(minprob, 10)] <- -log(minprob, 10)
        return(o)
      }
      if (ebcfun != FALSE) {
        if (n==1) {odf+effectodds} else {apply(odf, 1, function(x) {sum(x)/sqrt(length(x))})+effectodds}
      } else {
        if (n==1) {odf} else {apply(odf, 1, function(x) {sum(x)/sqrt(length(x))})}
      }
    }
    allodds <- as.data.frame(t(round(allodds, digits = 3)))
    colnames(allodds) <- effectrange
    rownames(allodds) <- genes
    if (output == "both" || semiexact == TRUE) {
      if (semiexact == TRUE) {
        exactodds <- semiexact(allodds, 0)
      } else {
        exactodds <- foreach(gene=seq_along(genes), .combine = c) %do% {
          lo <- which(allodds[gene,] < 0)
          hi <- which(allodds[gene,] > 0)
          if (length(lo) == 0) {
            highe <- effectrange[min(hi)]
            lowe <- 2*highe
          } else {
            lowe <- effectrange[max(lo)]
            if (length(hi) == 0) {
              highe <- 2*lowe
            } else {
              highe <- effectrange[min(hi)]
            }
          }
          mean(approximate(function(e) {
            fc <- (g*2^((1+e)*a) + (1-g)*2^a)/2^a
            if (ebcfun != FALSE) {
              eodds <- oddscores(c(e, ebcr), 
                                 normsubset = seq_along(ebcr)+1, 
                                 log = 10)[1]
            }
            pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
              index <- genei[gene]-nguides[gene]+i
              1-pbinom(t1[index], t1[index]+t0[index], fc * sr1 / (fc * sr1 + sr0))
            }
            p <- apply(pmat, 2, function(p) {sum(p*gw)})
            o <- log(p/(1-p), 10)
            o[o < log(minprob, 10)] <- log(minprob, 10) 
            o[o > -log(minprob, 10)] <- -log(minprob, 10)
            if (ebcfun != FALSE) {
              sum(o)/sqrt(ncol(pmat))+eodds
            } else {
              sum(o)/sqrt(ncol(pmat))
            }
          }, 0, lowe, highe, 11, 2))
        }
      }
      if (exactci != FALSE) {
        plo <- (1-exactci)/2
        if (semiexact == TRUE) {
          olo <- log(plo/(1-plo),10)
          exactlow <- semiexact(allodds, olo)
          exacthigh <- semiexact(allodds, -olo)
        } else {
          exactlow <- foreach(gene=seq_along(genes), .combine = c) %do% {
            lo <- which(allodds[gene,] < log(plo/(1-plo),10))
            hi <- which(allodds[gene,] > log(plo/(1-plo),10))
            if (length(lo) == 0) {
              highe <- effectrange[min(hi)]
              lowe <- 2*highe
            } else {
              lowe <- effectrange[max(lo)]
              if (length(hi) == 0) {
                highe <- 2*lowe
              } else {
                highe <- effectrange[min(hi)]
              }
            }
            mean(approximate(function(e) {
              fc <- (g*2^((1+e)*a) + (1-g)*2^a)/2^a
              if (ebcfun != FALSE) {
                eodds <- oddscores(c(e, ebcr), 
                                   normsubset = seq_along(ebcr)+1, 
                                   log = 10)[1]
              }
              pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
                index <- genei[gene]-nguides[gene]+i
                1-pbinom(t1[index], t1[index]+t0[index], fc * sr1 / (fc * sr1 + sr0))
              }
              p <- apply(pmat, 2, function(p) {sum(p*gw)})
              o <- log(p/(1-p), 10)
              o[o < log(minprob, 10)] <- log(minprob, 10) 
              o[o > -log(minprob, 10)] <- -log(minprob, 10)
              if (ebcfun != FALSE) {
                sum(o)/sqrt(ncol(pmat))+eodds
              } else {
                sum(o)/sqrt(ncol(pmat))
              }
            }, log(plo/(1-plo),10), lowe, highe, 11, 2))
          }
          exacthigh <- foreach(gene=seq_along(genes), .combine = c) %do% {
            lo <- which(allodds[gene,] < -log(plo/(1-plo),10))
            hi <- which(allodds[gene,] > -log(plo/(1-plo),10))
            if (length(lo) == 0) {
              highe <- effectrange[min(hi)]
              lowe <- 2*highe
            } else {
              lowe <- effectrange[max(lo)]
              if (length(hi) == 0) {
                highe <- 2*lowe
              } else {
                highe <- effectrange[min(hi)]
              }
            }
            mean(approximate(function(e) {
              fc <- (g*2^((1+e)*a) + (1-g)*2^a)/2^a
              if (ebcfun != FALSE) {
                eodds <- oddscores(c(e, ebcr), 
                                   normsubset = seq_along(ebcr)+1, 
                                   log = 10)[1]
              }
              pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
                index <- genei[gene]-nguides[gene]+i
                1-pbinom(t1[index], t1[index]+t0[index], fc * sr1 / (fc * sr1 + sr0))
              }
              p <- apply(pmat, 2, function(p) {sum(p*gw)})
              o <- log(p/(1-p), 10)
              o[o < log(minprob, 10)] <- log(minprob, 10) 
              o[o > -log(minprob, 10)] <- -log(minprob, 10)
              if (ebcfun != FALSE) {
                sum(o)/sqrt(ncol(pmat))+eodds
              } else {
                sum(o)/sqrt(ncol(pmat))
              }
            }, -log(plo/(1-plo),10), lowe, highe, 11, 2))
          }
        }
        exactdf <- data.frame(gene = genes, 
                              exact = exactodds, 
                              lowerci = exactlow,
                              upperci = exacthigh)
        return(list(range=allodds, exact=exactdf))
      } else {
        names(exactodds) <- genes
        return(list(range=allodds, exact=exactodds))
      }
    } else {return(allodds)}
  } else {
    exactodds <- foreach(gene=seq_along(genes), .combine = c) %do% {
      mean(approximate(function(e) {
        fc <- (g*2^((1+e)*a) + (1-g)*2^a)/2^a
        if (ebcfun != FALSE) {
          eodds <- oddscores(c(e, ebcr), 
                             normsubset = seq_along(ebcr)+1, 
                             log = 10)[1]
        }
        pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
          index <- genei[gene]-nguides[gene]+i
          1-pbinom(t1[index], t1[index]+t0[index], fc * sr1 / (fc * sr1 + sr0))
        }
        p <- apply(pmat, 2, function(p) {sum(p*gw)})
        o <- log(p/(1-p), 10)
        o[o < log(minprob, 10)] <- log(minprob, 10) 
        o[o > -log(minprob, 10)] <- -log(minprob, 10)
        if (ebcfun != FALSE) {
          sum(o)/sqrt(ncol(pmat))+eodds
        } else {
          sum(o)/sqrt(ncol(pmat))
        }
      }, 0, min(effectrange), max(effectrange), 11, 3))
    }
    if (exactci != FALSE) {
      plo <- (1-exactci)/2
      exactlow <- foreach(gene=seq_along(genes), .combine = c) %do% {
        mean(approximate(function(e) {
          fc <- (g*2^((1+e)*a) + (1-g)*2^a)/2^a
          if (ebcfun != FALSE) {
            eodds <- oddscores(c(e, ebcr), 
                               normsubset = seq_along(ebcr)+1, 
                               log = 10)[1]
          }
          pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
            index <- genei[gene]-nguides[gene]+i
            1-pbinom(t1[index], t1[index]+t0[index], fc * sr1 / (fc * sr1 + sr0))
          }
          p <- apply(pmat, 2, function(p) {sum(p*gw)})
          o <- log(p/(1-p), 10)
          o[o < log(minprob, 10)] <- log(minprob, 10) 
          o[o > -log(minprob, 10)] <- -log(minprob, 10)
          if (ebcfun != FALSE) {
            sum(o)/sqrt(ncol(pmat))+eodds
          } else {
            sum(o)/sqrt(ncol(pmat))
          }
        }, log(plo/(1-plo),10), min(effectrange), max(effectrange), 11, 3))
      }
      exacthigh <- foreach(gene=seq_along(genes), .combine = c) %do% {
        mean(approximate(function(e) {
          fc <- (g*2^((1+e)*a) + (1-g)*2^a)/2^a
          if (ebcfun != FALSE) {
            eodds <- oddscores(c(e, ebcr), 
                               normsubset = seq_along(ebcr)+1, 
                               log = 10)[1]
          }
          pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
            index <- genei[gene]-nguides[gene]+i
            1-pbinom(t1[index], t1[index]+t0[index], fc * sr1 / (fc * sr1 + sr0))
          }
          p <- apply(pmat, 2, function(p) {sum(p*gw)})
          o <- log(p/(1-p), 10)
          o[o < log(minprob, 10)] <- log(minprob, 10) 
          o[o > -log(minprob, 10)] <- -log(minprob, 10)
          if (ebcfun != FALSE) {
            sum(o)/sqrt(ncol(pmat))+eodds
          } else {
            sum(o)/sqrt(ncol(pmat))
          }
        }, -log(plo/(1-plo),10), min(effectrange), max(effectrange), 11, 3))
      }
      exactdf <- data.frame(gene = genes, 
                            exact = exactodds, 
                            lowerci = exactlow,
                            upperci = exacthigh)
      return(exactdf)
    } else {
      names(exactodds) <- genes
      return(exactodds)
    }
  }
}


#' @rdname geteffect
#' @export
geteffect_r <- function(guides, r, rse, a, g, gw, gl = 11, 
                        subset, effectrange, output = "range", 
                        exactci = FALSE, semiexact = FALSE,  
                        ebcfun = "max", minprob = 10^-20) {
  if (missing(a)) {
    stop("enter the presumed number of population doublings")
  }
  if (missing(g)) {
    if (gl < 4) {
      warning("gl should be at least 4: it is now set to 4")
      gl <- 4
    }
    g <- (11-11^((seq_len(gl)-1)/(gl-1)))/10
  }
  if (missing(gw)) {
    gw <- c(1,3,2,rep(1, length(g)-4),2)
  }
  if (length(g) != length(gw)) {
    warning("Length of g's and weights of g's is unequal")
    if (length(gw) < length(g)) {
      gw <- c(gw, rep(tail(gw, 1), length(g)-length(gw)))
    } else {
      gw <- gw[seq_along(g)]
    }
  }
  # Normalize gene weights so they sum to 1
  gw <- gw/sum(gw)
  
  if (missing(effectrange)) {
    if (!missing(subset)) {
      effectrange <- round(seq(floor(30*min(r[subset])/a), ceiling(30*max(r[subset])/a))/20, 2)
    } else {
      effectrange <- round(seq(floor(30*min(r)/a), ceiling(30*max(r)/a))/20, 2)
    }
  } else {
    # effectrange needs to be sorted from low to high
    effectrange <- sort(effectrange)
  }
  
  # note that guides are presumed to be named [gene]_[number | sequence]
  genes <- rle(gsub("_.*", "", guides))$values
  nguides <- rle(gsub("_.*", "", guides))$lengths
  genei <- cumsum(nguides)
  
  fcmat <- matrix(nrow = length(effectrange), ncol = length(g))
  for (c in seq_len(ncol(fcmat))) {
    fcmat[,c] <- log((g[c]*2^((1+effectrange)*a) + (1-g[c])*2^a)/2^a,2)
  }
  
  if (ebcfun != FALSE) {
    if (ebcfun == "max") {
      ebcr <- mapply(function(n, gi) {
        r[(gi-n+1):(gi)][which.max(abs(r[(gi-n+1):(gi)]))]
      }, nguides, genei)/a
    } else {
      fun <- get(ebcfun)
      ebcr <- mapply(function(n, gi) {
        fun(r[(gi-n+1):(gi)])
      }, nguides, genei)/a
    }
    effectodds <- oddscores(c(effectrange, ebcr), 
                            normsubset = seq_along(ebcr)+length(effectrange), 
                            log = 10)[seq_along(effectrange)]
  }
  
  if (!missing(subset)) {
    if (!is(subset, "numeric")) {
      subset <- which(gsub("_.*", "", guides) %in% subset)
    }
    guides <- guides[subset]
    genes <- rle(gsub("_.*", "", guides))$values
    nguides <- rle(gsub("_.*", "", guides))$lengths
    genei <- cumsum(nguides)
  }
  
  # I have noted I have to first bind the variables used in foreach...
  n <- gi <- gene <- i <- NULL
  
  if (output == "range" || output == "both") {
    # I imagine it is beneficial to perform the outer loop in parallel using
    # %dopar% if a parallel backend is available
    allodds <- foreach(n=nguides, gi=genei, .combine = cbind) %do% {
      odf <- foreach(i = seq_len(n), .combine = cbind) %do% {
        pmat <- pnorm((fcmat-r[gi-n+i])/rse[gi-n+i])
        p <- apply(pmat, 1, function(p) {sum(p*gw)})
        o <- log(p/(1-p), 10)
        o[o < log(minprob, 10)] <- log(minprob, 10) 
        o[o > -log(minprob, 10)] <- -log(minprob, 10)
        return(o)
      }
      if (ebcfun != FALSE) {
        if (n==1) {odf+effectodds} else {apply(odf, 1, function(x) {sum(x)/sqrt(length(x))})+effectodds}
      } else {
        if (n==1) {odf} else {apply(odf, 1, function(x) {sum(x)/sqrt(length(x))})}
      }
      
    }
    allodds <- as.data.frame(t(round(allodds, digits = 3)))
    colnames(allodds) <- effectrange
    rownames(allodds) <- genes
    if (output == "both" || semiexact == TRUE) {
      if (semiexact == TRUE) {
        exactodds <- semiexact(allodds, 0)
      } else {
        exactodds <- foreach(gene=seq_along(genes), .combine = c) %do% {
          lo <- which(allodds[gene,] < 0)
          hi <- which(allodds[gene,] > 0)
          if (length(lo) == 0) {
            highe <- effectrange[min(hi)]
            lowe <- 2*highe
          } else {
            lowe <- effectrange[max(lo)]
            if (length(hi) == 0) {
              highe <- 2*lowe
            } else {
              highe <- effectrange[min(hi)]
            }
          }
          mean(approximate(function(e) {
            fc <- log((g*2^((1+e)*a) + (1-g)*2^a)/2^a,2)
            if (ebcfun != FALSE) {
              eodds <- oddscores(c(e, ebcr), 
                                 normsubset = seq_along(ebcr)+1, 
                                 log = 10)[1]
            }
            pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
              pnorm((fc-r[genei[gene]-nguides[gene]+i])/rse[genei[gene]-nguides[gene]+i])
            }
            p <- apply(pmat, 2, function(p) {sum(p*gw)})
            o <- log(p/(1-p), 10)
            o[o < log(minprob, 10)] <- log(minprob, 10) 
            o[o > -log(minprob, 10)] <- -log(minprob, 10)
            if (ebcfun != FALSE) {
              sum(o)/sqrt(ncol(pmat))+eodds
            } else {
              sum(o)/sqrt(ncol(pmat))
            }
          }, 0, lowe, highe, 11, 2))
        }
      }
      if (exactci != FALSE) {
        plo <- (1-exactci)/2
        if (semiexact == TRUE) {
          olo <- log(plo/(1-plo),10)
          exactlow <- semiexact(allodds, olo)
          exacthigh <- semiexact(allodds, -olo)
        } else {
          exactlow <- foreach(gene=seq_along(genes), .combine = c) %do% {
            lo <- which(allodds[gene,] < log(plo/(1-plo),10))
            hi <- which(allodds[gene,] > log(plo/(1-plo),10))
            if (length(lo) == 0) {
              highe <- effectrange[min(hi)]
              lowe <- 2*highe
            } else {
              lowe <- effectrange[max(lo)]
              if (length(hi) == 0) {
                highe <- 2*lowe
              } else {
                highe <- effectrange[min(hi)]
              }
            }
            mean(approximate(function(e) {
              fc <- log((g*2^((1+e)*a) + (1-g)*2^a)/2^a,2)
              if (ebcfun != FALSE) {
                eodds <- oddscores(c(e, ebcr), 
                                   normsubset = seq_along(ebcr)+1, 
                                   log = 10)[1]
              }
              pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
                pnorm((fc-r[genei[gene]-nguides[gene]+i])/rse[genei[gene]-nguides[gene]+i])
              }
              p <- apply(pmat, 2, function(p) {sum(p*gw)})
              o <- log(p/(1-p), 10)
              o[o < log(minprob, 10)] <- log(minprob, 10) 
              o[o > -log(minprob, 10)] <- -log(minprob, 10)
              if (ebcfun != FALSE) {
                sum(o)/sqrt(ncol(pmat))+eodds
              } else {
                sum(o)/sqrt(ncol(pmat))
              }
            }, log(plo/(1-plo),10), lowe, highe, 11, 2))
          }
          exacthigh <- foreach(gene=seq_along(genes), .combine = c) %do% {
            lo <- which(allodds[gene,] < -log(plo/(1-plo),10))
            hi <- which(allodds[gene,] > -log(plo/(1-plo),10))
            if (length(lo) == 0) {
              highe <- effectrange[min(hi)]
              lowe <- 2*highe
            } else {
              lowe <- effectrange[max(lo)]
              if (length(hi) == 0) {
                highe <- 2*lowe
              } else {
                highe <- effectrange[min(hi)]
              }
            }
            mean(approximate(function(e) {
              fc <- log((g*2^((1+e)*a) + (1-g)*2^a)/2^a,2)
              if (ebcfun != FALSE) {
                eodds <- oddscores(c(e, ebcr), 
                                   normsubset = seq_along(ebcr)+1, 
                                   log = 10)[1]
              }
              pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
                pnorm((fc-r[genei[gene]-nguides[gene]+i])/rse[genei[gene]-nguides[gene]+i])
              }
              p <- apply(pmat, 2, function(p) {sum(p*gw)})
              o <- log(p/(1-p), 10)
              o[o < log(minprob, 10)] <- log(minprob, 10) 
              o[o > -log(minprob, 10)] <- -log(minprob, 10)
              if (ebcfun != FALSE) {
                sum(o)/sqrt(ncol(pmat))+eodds
              } else {
                sum(o)/sqrt(ncol(pmat))
              }
            }, -log(plo/(1-plo),10), lowe, highe, 11, 2))
          }
        }
        exactdf <- data.frame(gene = genes, 
                              exact = exactodds, 
                              lowerci = exactlow,
                              upperci = exacthigh)
        return(list(range=allodds, exact=exactdf))
      } else {
        names(exactodds) <- genes
        return(list(range=allodds, exact=exactodds))
      }
      
    } else {return(allodds)}
  } else {
    exactodds <- foreach(gene=seq_along(genes), .combine = c) %do% {
      mean(approximate(function(e) {
        fc <- log((g*2^((1+e)*a) + (1-g)*2^a)/2^a,2)
        if (ebcfun != FALSE) {
          eodds <- oddscores(c(e, ebcr), 
                             normsubset = seq_along(ebcr)+1, 
                             log = 10)[1]
        }
        pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
          pnorm((fc-r[genei[gene]-nguides[gene]+i])/rse[genei[gene]-nguides[gene]+i])
        }
        p <- apply(pmat, 2, function(p) {sum(p*gw)})
        o <- log(p/(1-p), 10)
        o[o < log(minprob, 10)] <- log(minprob, 10) 
        o[o > -log(minprob, 10)] <- -log(minprob, 10)
        if (ebcfun != FALSE) {
          sum(o)/sqrt(ncol(pmat))+eodds
        } else {
          sum(o)/sqrt(ncol(pmat))
        }
      }, 0, min(effectrange), max(effectrange), 11, 3))
    }
    
    if (exactci != FALSE) {
      plo <- (1-exactci)/2
      exactlow <- foreach(gene=seq_along(genes), .combine = c) %do% {
        mean(approximate(function(e) {
          fc <- log((g*2^((1+e)*a) + (1-g)*2^a)/2^a,2)
          if (ebcfun != FALSE) {
            eodds <- oddscores(c(e, ebcr), 
                               normsubset = seq_along(ebcr)+1, 
                               log = 10)[1]
          }
          pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
            pnorm((fc-r[genei[gene]-nguides[gene]+i])/rse[genei[gene]-nguides[gene]+i])
          }
          p <- apply(pmat, 2, function(p) {sum(p*gw)})
          o <- log(p/(1-p), 10)
          o[o < log(minprob, 10)] <- log(minprob, 10) 
          o[o > -log(minprob, 10)] <- -log(minprob, 10)
          if (ebcfun != FALSE) {
            sum(o)/sqrt(ncol(pmat))+eodds
          } else {
            sum(o)/sqrt(ncol(pmat))
          }
        }, log(plo/(1-plo),10), min(effectrange), max(effectrange), 11, 3))
      }
      exacthigh <- foreach(gene=seq_along(genes), .combine = c) %do% {
        mean(approximate(function(e) {
          fc <- log((g*2^((1+e)*a) + (1-g)*2^a)/2^a,2)
          if (ebcfun != FALSE) {
            eodds <- oddscores(c(e, ebcr), 
                               normsubset = seq_along(ebcr)+1, 
                               log = 10)[1]
          }
          pmat <- foreach(i=seq_len(nguides[gene]), .combine = cbind) %do% {
            pnorm((fc-r[genei[gene]-nguides[gene]+i])/rse[genei[gene]-nguides[gene]+i])
          }
          p <- apply(pmat, 2, function(p) {sum(p*gw)})
          o <- log(p/(1-p), 10)
          o[o < log(minprob, 10)] <- log(minprob, 10) 
          o[o > -log(minprob, 10)] <- -log(minprob, 10)
          if (ebcfun != FALSE) {
            sum(o)/sqrt(ncol(pmat))+eodds
          } else {
            sum(o)/sqrt(ncol(pmat))
          }
        }, -log(plo/(1-plo),10), min(effectrange), max(effectrange), 11, 3))
      }
      exactdf <- data.frame(gene = genes, 
                            exact = exactodds, 
                            lowerci = exactlow,
                            upperci = exacthigh)
      return(exactdf)
    } else {
      names(exactodds) <- genes
      return(exactodds)
    }
  }
}

approximate <- function(fun, solution, from, to, steps, depth) {
  fun <- match.fun(fun)
  for (i in seq_len(depth)) {
    input <- seq(from, to, length.out = steps)
    answers <- sapply(input, fun)
    if (solution %in% answers) {
      return(input[answers == solution])
    } else {
      diff <- answers-solution
      lo <- which(diff < 0)
      if (length(lo)==0) {
        warning("solution out of range")
        return(NaN)
      }
      hi <- which(diff > 0)
      if (length(hi)==0) {
        warning("solution out of range")
        return(NaN)
      }
      if (max(lo) < max(hi)) {
        from <- input[lo[tail(which(diff[lo]==max(diff[lo])),1)]]
        to <- input[hi[head(which(diff[hi]==min(diff[hi])),1)]]
      } else {
        from <- input[lo[head(which(diff[lo]==max(diff[lo])),1)]]
        to <- input[hi[tail(which(diff[hi]==min(diff[hi])),1)]]
      }
    }
  }
  return(c(from, to))
}

semiexact <- function(range, o) {
  as.numeric(apply(range, 1, function(x) {
    l <- tail(which(x < o), 1)
    h <- head(which(x > o), 1)
    lo <- as.numeric(names(x)[l])
    hi <- as.numeric(names(x)[h])
    lo + (hi-lo)*(10^(x[l]+x[h]-2*o) / (1 + 10^(x[l]+x[h]-2*o)))
  }))
}
