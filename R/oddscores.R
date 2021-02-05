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
#' user-defined parameters. In a way, this is simplified version of
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
#' @param allseed Integer. If specified, all unspecified seeds default to this.
#'   Default = NULL
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
  if (!missing(allseed)) {set.seed(allseed)}
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
  
  if (!missing(gseed)) {set.seed(gseed)} else {gseed <- allseed; set.seed(gseed)}
  if (missing(g)) {
    message("guide efficiency is randomly sampled from an empirical distribution")
    bg <- floor(length(guides)/40)
    g <- sample(c(1-runif(length(guides)-bg)^2.5, (runif(bg)/10)^5))
  } else if (length(g)==1) {
    g <- rep(g, length(guides))
  } else if (length(g) != length(guides)) {
    g <- sample(g, length(guides), replace = TRUE)
  } 
  
  if (any(g < 0) || any(g > 1)) {
    warning("g represents guide efficacy and should be a number between 0 and 1; values have been constrained")
    g[g < 0] <- 0
    g[g > 1] <- 1
  }
  
  if (!missing(fseed)) {
    set.seed(fseed)
  } else if (!is.null(gseed)) {
    fseed <- gseed+1
    set.seed(fseed)
  } else {
    # no need in resetting the seed when it's NULL
    fseed <- NULL
  }
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
  
  if (!missing(dseed)) {set.seed(dseed)} else {dseed <- allseed; set.seed(dseed)}
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
  
  if (any(is.na(d)) || any(is.null(d)) || any(is.character(d))) {
    stop("This is not going to go well, not all d-values are numeric!")
  }
  
  if (!missing(e)) {
    if (!missing(eseed)) {
      set.seed(eseed)
    } else if (!is.null(dseed)) {
      eseed <- dseed+1
      set.seed(eseed)
    } else {
      # no need in resetting the seed when it's NULL
      eseed <- NULL
    }
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
  }
  
  offtarget_guides <- rep(0, length(genes))
  
  if (!missing(oseed)) {set.seed(oseed)} else {oseed <- allseed; set.seed(oseed)}
  if (offtargets != FALSE) {
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
  
  if (perfectsampling == TRUE) {
    if (!missing(t0seed)) {set.seed(t0seed)} else {t0seed <- allseed; set.seed(t0seed)}
    cells <- round(c(f*g*sortedcells, f*(1-g)*sortedcells)/sum(f))
    if (!missing(repseed)) {set.seed(repseed)} else {repseed <- allseed; set.seed(repseed)}
    kocells <- cells[seq_along(guides)]
    nkocells <- cells[seq(length(guides)+1, 2*length(guides))]
    koselected <- round(kocells*mod_prob)
    nkoselected <- round(nkocells*baseprob)
  } else {
    if (!missing(t0seed)) {set.seed(t0seed)} else {t0seed <- allseed; set.seed(t0seed)}
    cells <- tabulate(c(seq_along(c(guides,guides)),
                         sample(seq_along(c(guides,guides)), sortedcells,
                                replace = TRUE, prob = c(f*g, f*(1-g)))))-1
    kocells <- cells[seq_along(guides)]
    nkocells <- cells[seq(length(guides)+1, 2*length(guides))]
    if (!missing(repseed)) {set.seed(repseed)} else {repseed <- allseed; set.seed(repseed)}
    koselected <- mapply(function(cells, chance) {rbinom(1,cells,chance)},kocells,mod_prob)
    nkoselected <- mapply(function(cells, chance) {rbinom(1,cells,chance)},nkocells,baseprob)
  }
  
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
  
  print("wrapping up")
  
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
