#' CSSA: Analysis and Simulation Tools for CRISPR-Cas9 Pooled Screens
#'
#' CSSA contains a number of functions to analyze and simulation pooled
#' CRISPR-Cas9-based screens. The simulators, represented by the CRISPRsim and
#' sortingsim functions, offer an easy but also highly customizable tool to
#' create data with guide- and gene-specific variables that can be randomly
#' assigned or specified by the user. The analysis tools currently contain
#' functions to calculate rate ratios, odds based on nonparametric guide
#' distribution, a function to combine genewise Z-values or odds, and a function
#' that calculates effect sizes of gene knockout and guide efficacy based on
#' rate ratios and duration of the experiment.
#'
#' @section CSSA functions:
#' \describe{
#' \item{\code{\link{CRISPRsim}}}{Simulate a CRISPR-Cas9 pooled screen}
#' \item{\code{\link{radjust}}}{Calculate rate ratios restricted by confidence level}
#' \item{\code{\link{nestedradjust}}}{Calculate rate ratios of rate ratios restricted by confidence level}
#' \item{\code{\link{jar}}}{Rate ratios after adding an artificial number to all features}
#' \item{\code{\link{doublejar}}}{Rate ratios of rate ratios after adding an artificial number to all features}
#' \item{\code{\link{r2Z}}}{Convert rate ratios to Z-values}
#' \item{\code{\link{sumZ}}}{Calculate corrected summed Z-values per gene}
#' \item{\code{\link{getdeg}}}{Derive growth-modifying effect of gene knockout in pooled experiments}
#' \item{\code{\link{ess}}}{Get a list of essential genes or corresponding indices in the data set}
#' \item{\code{\link{noness}}}{Get a list of nonessential genes or corresponding indices in the data set}
#' \item{\code{\link{sortingsim}}}{Simulate a selection-based CRISPR-Cas9 pooled screen}
#' \item{\code{\link{oddscores}}}{Calculate odds per guide}
#' \item{\code{\link{geneodds}}}{Calculate odds per gene}
#' \item{\code{\link{odds2pq}}}{Calculate p-values and q-values for odds}
#' \item{\code{\link{geteffect}}}{Approximate effect of gene knockout on growth}
#' \item{\code{\link{rrep}}}{Calculate rate ratios with standard error of count data based on replicates}
#' \item{\code{\link{degrep}}}{Derive growth-modifying effect of gene knockout in pooled experiments with replicate arms}
#' }
#' 
#' @section License:
#' This package is licensed under GPL.
#' 
#' @author Jos B. Poell
#' 
#' @importFrom methods is
#' @importFrom stats median p.adjust pnorm qnorm poisson.test rbinom pbinom rnorm runif sd mad var lm
#' @importFrom utils data head tail write.table
#' @importFrom foreach foreach "%do%"
#' @importFrom withr with_seed
#' @docType package
#' @name CSSA-package
#' @aliases CSSA CSSA-package
NULL

#' Simulate a CRISPR-Cas9 pooled screen
#'
#' CRISPRsim simulates a CRISPR-Cas9 pooled screen with user-defined parameters.
#' These include drug treatment screens! Each "infected" cell expands over time
#' based on effect of gene knockout. Other parameters include the abundance of
#' each guide at the start of the experiment, the efficacy of the guide (chance
#' that it results in a successful gene knockout), and the frequency and depth
#' of sampling. In case of drug treatment, genes are assigned a
#' treatment-specific growth modifier as well. The result is a data frame that
#' contains the guide-relevant parameters and the sequencing coverage per guide
#' for the specified time intervals. Simulated screens will aid researchers with
#' their experimental setup. Furthermore, this offers a unique platform for the
#' evaluation of analysis methods for pooled gene knockout screens.
#'
#' @param genes Single integer or character vector. Specify how many or which
#'   genes to include in the experiment respectively. Not required when a full
#'   list of guides is given.
#' @param guides Single integer, integer vector or character vector. In case of
#'   single integer, specify by how many guides each gene is represented. In
#'   case of an integer vector, specify per gene by how many guides it is
#'   represented. In case of a character vector, guides are assumed to contain a
#'   gene name, followed by an underscore, followed by an identifier within that
#'   gene (e.g. a number or a nucleotide sequence)
#' @param a Numeric. Specify the number of doublings between each "passaging".
#'   For example, in case of an experiment that ends after 12 doublings and was
#'   passaged 3 times, specify a = c(4,4,4)
#' @param g Integer vector. Specify guide efficacies per guide. If omitted,
#'   guide efficacies will be sampled from a representative distribution
#' @param f Integer vector. Specify guide abundance at time of infection per
#'   guide. If omitted, guide abundance will be sampled from a representative
#'   distribution
#' @param d Integer vector. Specify gene-specific growth effect. If omitted,
#'   effect of gene knockout on growth will be sampled from a representative
#'   distribution. If the length of the vector does not match the number of
#'   genes, values will be randomly sampled from the specified distribution!
#' @param e Integer vector. Specify treatment-specific growth effect per gene.
#'   If omitted, effects will be sampled from a representative distribution. If
#'   the length of the vector does not match the number of genes, values will be
#'   randomly sampled from the specified distribution!
#' @param seededcells Integer. Number of cells seeded at the start of each
#'   experimental step. If the length of this argument is smaller than the
#'   number of seedings, all unspecified steps will be assumed equal to the last
#'   specified step! Defaults to 200 times the number of guides
#' @param harvestedcells Integer. Number of cells from which to sample for
#'   subsequent sequencing. This argument is especially useful to restrict the
#'   expected DNA copies present in the PCR reaction. If you wish to do so, make
#'   sure to set harvestall to FALSE. Defaults to be equal to seededcells
#' @param harvestall Logical. If TRUE, all cells are collected and used for
#'   sampling in subsequent sequencing step. Applies to all experimental time
#'   points beyond t0. Default = TRUE
#' @param cellreplace Logical. If FALSE, cells are sampled from the total pool
#'   of cells without replacement. Note that this is the most realistic
#'   simulation of a screen, but then you should also keep realistic passaging
#'   times! It is recommended to keep the total number of cells in the
#'   experiment below 200 million. Setting this to TRUE can dramatically speed
#'   up simulations. Default = FALSE
#' @param treatmentdelay Integer. In case of a treatment experiment, specify
#'   when treatment starts. It is currently only possible to start treatment on
#'   one of the experimental time points. Default = 0
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
#' @param t0seed Integer. Specify seed for t0, which encompasses sampling of the
#'   first seeding and the assignment of successful knockout cells versus no
#'   knockout cells for each guide
#' @param repseed Integer. Specify the seed after t0
#' @param grm Numeric. Growth rate modifier. Specify adjusted growth rate under
#'   treatment conditions. Default = 1
#' @param em Numeric. Effect modifier. Specify how effective treatment is. This
#'   can be used as a proxy for drug concentration. All individual e-values and
#'   grm are modified by this multiplier. Default = 1
#' @param perfectsampling Logical. If TRUE, all sampling steps are replaced by
#'   simple equations to calculate representation of guides. Useful as null
#'   control to isolate the effect of sampling. Default = FALSE
#' @param perfectseq Logical. If TRUE, sequencing results are a perfect
#'   representation (though still rounded) of guides in the harvested cells.
#'   Applicable to speed up simulations, assuming sequencing is sufficiently
#'   deep. Default = FALSE
#' @param returnall Logical. If TRUE, function returns a list with the simulated
#'   data in the guidesdf, summary per gene in the genesdf, and parameters.
#'   Default = FALSE
#' @param outputfile Character string. When used, returned data frame will be
#'   saved as a tab-delimited text to the specified file path
#'
#' @details CRISPRsim performs a genome-wide (or subsetted) pooled CRISPR
#'   knockout screen for you without having to go to the lab and spend
#'   incredible amounts of time and money. This can be a tremendous help if you
#'   want to design an experiment and answer questions such as: how many
#'   replicates do I need, how much coverage, will I pick up genes with x
#'   effect, et cetera. You can give it a spin, but I highly recommend checking
#'   out the documentation for the available parameters! Especially seeds can be
#'   relevant for a proper simulation. You can easily "practice" by simulating
#'   some small experiments. The basis of the simulation is as follows. Between
#'   time points cells with a certain knockout grow according to formula
#'   \code{cellsout = cellsin*2^((grm+d+e)*a))} Each guide has an efficacy,
#'   which is the chance to create a successful knockout. The cellsin is
#'   determined at t0 and depends on guide efficacy and guide abundance. If
#'   there is no successful knockout, d and e are 0. Cells with and without
#'   successful knockout are followed separately throughout the experiment, but
#'   the pairs are pooled in terms of sequencing reads. grm is the growth rate
#'   modifier and is generally 1, but it can be lowered to more properly
#'   simulate resistance screens.
#'
#' @return Returns a data frame with every row representing a single guide.
#'   Contains the pertinent parameters of each guide and the number of
#'   sequencing reads on t0 and all other sampling time points. If the argument
#'   returnall is set to TRUE, the function also returns a data frame with the
#'   true values for the genes, and lists all parameters as well.
#'
#' @note Seeds are set using \code{\link[withr]{with_seed}} from the
#'   \code{\link[withr:withr-package]{withr package}}, thus leaving any
#'   pre-existing seed intact. Avoid using the same seeds for different
#'   arguments. If dseed and eseed are identical, the resulting values for d and
#'   e will have a distinct pattern of correlation. In this case, CRISPRsim will
#'   throw a warning. Given or generated seeds are returned when the option
#'   returnall is set to TRUE.
#'
#' @seealso \code{\link{sortingsim}}, \code{\link{radjust}}, \code{\link{rrep}},
#'   \code{\link{jar}}, \code{\link{nestedradjust}}, \code{\link{doublejar}}
#'
#' @author Jos B. Poell
#'
#' @examples
#' simdf <- CRISPRsim(18000, 4, a = c(3,3), e = TRUE, perfectsampling = TRUE)
#' hist(simdf$g, breaks = 100, main = "distribution of guide efficiencies")
#' d <- rle(simdf$d)$values
#' e <- rle(simdf$e)$values
#' plot(d, e, main = "straight lethality and sensitization")
#'
#' @export

CRISPRsim <- function(genes, guides, a, g, f, d, e, seededcells, harvestedcells,
                      harvestall = TRUE, cellreplace = FALSE, treatmentdelay = 0,
                      seqdepth, offtargets = FALSE, allseed = NULL, gseed,
                      fseed, dseed, eseed, oseed, t0seed, repseed, grm = 1, 
                      em = 1, perfectsampling = FALSE, perfectseq = FALSE, 
                      returnall = FALSE, outputfile) {  
  
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
  if (missing(a)) {
    stop("please enter the number of population doublings before sampling, or a vector with doublings until the next sampling")
  }
  if (cellreplace == FALSE) {
    message("keep in mind sampling without replacement can blow up your computer!")
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
      
      ll <- floor(length(genes)/8)
      ml <- floor(length(genes)/3)
      d <- sample(c(rnorm(ll, -0.8, 0.15), rnorm(ml, -0.25, 0.2),
                    rnorm(length(genes)-ll-ml, -0.05, 0.05)))
    } else if (length(d)==1) {
      message("assuming equal d for all genes. Note that you can enter custom d-values by providing a vector with length equal to the number of genes, or a vector with unequal size from which will be sampled")
      d <- rep(d, length(genes))
      
    } else if (length(d) != length(genes)) {
      d <- sample(d, length(genes), replace = TRUE)
    } 
  })
  
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
        
        ls <- ceiling(length(genes)/200)
        lr <- ceiling(length(genes)/1000)
        e <- sample(c(rnorm(ls, -0.7, 0.2), rnorm(lr, 0.7, 0.2),
                      rnorm(length(genes)-ls-lr, 0, 0.05)))
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
  
  unique_d <- d
  d <- rep(d, n)
  
  if (any(is.na(d)) || any(is.null(d)) || any(is.character(d))) {
    stop("This is not going to go well, not all d-values are numeric!")
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
  
  if (!missing(e)) {
    gdf <- data.frame(genes, n, offtarget_guides, d = unique_d, e = unique_e)
  } else {
    gdf <- data.frame(genes, n, offtarget_guides, d = unique_d)
  }
  
  if (missing(seededcells)) {
    message("assuming mean representation of 200 cells per guide")
    seededcells <- rep(200*length(guides), length(a)+1)
  } else if (length(seededcells) <= length(a)) {
    seededcells <- append(seededcells, rep(tail(seededcells, 1), length(a)+1-length(seededcells)))
  }
  if (missing(harvestedcells)) {
    harvestedcells <- seededcells
  } else if (length(harvestedcells) <= length(a)) {
    harvestedcells <- append(harvestedcells, rep(tail(harvestedcells, 1), length(a)+1-length(harvestedcells)))
  }
  
  if (missing(seqdepth)) {
    message("assuming mean sequencing depth of 500 reads per guide")
    seqdepth <- rep(500*length(guides), length(a)+1)
  } else if (length(seqdepth) <= length(a)) {
    seqdepth <- append(seqdepth, rep(tail(seqdepth, 1), length(a)+1-length(seqdepth)))
  }
  
  if (missing(t0seed)) {
    if (is.null(allseed)) {
      t0seed <- sample.int(.Machine$integer.max, 1L)
    } else {
      t0seed <- allseed+5
    }
  }
  with_seed(t0seed, {
    message("started with t0")
    readmat <- matrix(ncol = length(seqdepth), nrow = length(guides))
    # scells is a vector twice the length of the number of guides. The first half
    # represents the cells with successful knockout of the targeted gene, the
    # second half represents the cells without knockout of the targeted gene
    if (perfectsampling == TRUE) {
      scells <- round(c(f*g*seededcells[1], f*(1-g)*seededcells[1])/sum(f))
      readmat[, 1] <- round(f*seqdepth[1]/sum(f))
    } else {
      scells <- tabulate(c(seq_along(c(guides,guides)),
                           sample(seq_along(c(guides,guides)), seededcells[1],
                                  replace = TRUE, prob = c(f*g, f*(1-g)))))-1
      hcells <- tabulate(c(seq_along(guides),
                           sample(seq_along(guides), harvestedcells[1],
                                  replace = TRUE, prob = f)))-1
      if (perfectseq == TRUE) {
        readmat[, 1] <- round(hcells*seqdepth[1]/sum(hcells))
      } else {
        readmat[, 1] <- tabulate(c(seq_along(guides),
                                   sample(seq_along(guides), seqdepth[1],
                                          replace = TRUE, prob = hcells)))-1
      }
    }
  })
  
  colnames(readmat) <- c("t0", paste0("t", cumsum(a)))
  
  if (missing(repseed)) {
    if (is.null(allseed)) {
      repseed <- sample.int(.Machine$integer.max, 1L)
    } else {
      repseed <- allseed+6
    }
  }
  with_seed(repseed, {
    kocells <- scells[seq_along(guides)]
    nkocells <- scells[seq(length(guides)+1, 2*length(guides))]
    st <- 0
    grme <- 1+grm*em-em
    for (t in seq_along(a)) {
      message("started with t", sum(a[seq(1,t)]))
      
      if (missing(e) || length(e) == 1) {
        cellpool <- mapply(function(stud, dud, d) {
          kocells <- round(stud*2^((1+d)*a[t]))
          nkocells <- round(dud*2^a[t])
          return(list(kocells = kocells, nkocells = nkocells))
        }, kocells, nkocells, d)
      } else {
        if (treatmentdelay < t) {st <- 1}
        cellpool <- mapply(function(stud, dud, d, e) {
          kocells <- round(stud*2^((grme^st+d+e)*a[t]))
          nkocells <- round(dud*2^(grme^st*a[t]))
          return(list(kocells = kocells, nkocells = nkocells))
        }, kocells, nkocells, d, e*st*em)
      }
      
      # Note: the vector hpool is made immediately. Not so the vectors kocells and
      # nkocells. These are reserved for newly sampled batches, to be used in the
      # next iteration of the loop. If the kocells or nkocells are needed, they
      # are unlisted from cellpool.
      hpool <- unlist(cellpool[1,])+unlist(cellpool[2,])
      neededcells <- seededcells[t+1] + harvestedcells[t+1]
      if (neededcells > sum(hpool)) {
        warning(paste0("not enough cells at t", t))
        seededcells[t+1] <- floor(sum(hpool)*seededcells[t+1]/neededcells)
        harvestedcells[t+1] <- floor(sum(hpool)*harvestedcells[t+1]/neededcells)
      }
      if (perfectsampling == TRUE) {
        kocells <- round(unlist(cellpool[1,])*seededcells[t+1]/sum(hpool))
        nkocells <- round(unlist(cellpool[2,])*seededcells[t+1]/sum(hpool))
        readmat[, t+1] <- round(hpool*seqdepth[t+1]/sum(hpool))
      } else {
        
        if (cellreplace == FALSE) {
          
          scells <- tabulate(c(seq_along(c(guides,guides)),
                               sample(rep(seq_along(c(guides,guides)),
                                          c(unlist(cellpool[1,]), unlist(cellpool[2,]))),
                                      seededcells[t+1], replace = FALSE)))-1
          kocells <- scells[seq_along(guides)]
          nkocells <- scells[seq(length(guides)+1, 2*length(guides))]
          if (harvestall == TRUE) {
            hcells <- hpool - kocells - nkocells
          } else {
            temppool <- hpool - kocells - nkocells
            hcells <- tabulate(c(seq_along(guides), sample(rep(seq_along(guides), temppool),
                                                           harvestedcells[t+1], replace = FALSE)))-1
          }
          if (perfectseq == TRUE) {
            readmat[,t+1] <- round(hcells*seqdepth[t+1]/sum(hcells))
          } else {
            readmat[,t+1] <- tabulate(c(seq_along(guides), sample(seq_along(guides),
                                                                  seqdepth[t+1], replace = TRUE, prob = hcells)))-1
          }
          
        } else {
          scells <- tabulate(c(seq_along(c(guides,guides)),
                               sample(seq_along(c(guides,guides)),
                                      seededcells[t+1], replace = TRUE,
                                      prob = c(unlist(cellpool[1,]), unlist(cellpool[2,])))))-1
          kocells <- scells[seq_along(guides)]
          nkocells <- scells[seq(length(guides)+1, 2*length(guides))]
          if (harvestall == TRUE) {
            hcells <- hpool - kocells - nkocells
            hcells[hcells < 0] <- 0
          } else {
            hcells <- tabulate(c(seq_along(guides),
                                 sample(seq_along(guides),
                                        harvestedcells[t+1], replace = TRUE, prob = hpool)))-1
          }
          if (perfectseq == TRUE) {
            readmat[,t+1] <- round(hcells*seqdepth[t+1]/sum(hcells))
          } else {
            readmat[,t+1] <- tabulate(c(seq_along(guides), sample(seq_along(guides),
                                                                  seqdepth[t+1], replace = TRUE, prob = hcells)))-1
          }
          
        }
        
      }
      
    }
  })
  
  message("wrapping up")
  
  df <- as.data.frame(readmat)
  
  if (missing(e)) {
    df <- data.frame(guides, gene = rep(genes, n), d, f, g, df, stringsAsFactors = FALSE)
  } else {
    df <- data.frame(guides, gene = rep(genes, n), d, e, f, g, df, stringsAsFactors = FALSE)
  }
  
  if (!missing(outputfile)) {
    write.table(df, outputfile, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  if (returnall == TRUE && missing(e)) {
    return(list(
      guidesdf = df,
      genesdf = gdf,
      a = a,
      offtargets = offtargets,
      seededcells = seededcells, 
      harvestedcells = harvestedcells,
      harvestall = harvestall,
      seqdepth = seqdepth,
      cellreplace = cellreplace, 
      perfectsampling = perfectsampling,
      perfectseq = perfectseq,
      grm = grm,
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
      a = a,
      offtargets = offtargets,
      treatmentdelay = treatmentdelay,
      seededcells = seededcells, 
      harvestedcells = harvestedcells,
      harvestall = harvestall,
      seqdepth = seqdepth,
      cellreplace = cellreplace, 
      perfectsampling = perfectsampling,
      perfectseq = perfectseq,
      grm = grm,
      em = em,
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

#' Calculate rate ratios restricted by confidence level
#'
#' radjust is designed to calculate rate ratios of sequencing results. It uses
#' the binomial distribution to calculate confidence intervals of rate ratios
#' (and the normal approximation for large counts), and returns the
#' log2-transformed confidence limit that is closest to 0.
#'
#' @param t1 Integer vector or matrix. Raw sequencing reads in test sample
#' @param t0 Integer vector or matrix. Raw sequencing reads in control sample
#' @param conf.level Numeric. Sets the confidence level of the rate ratio. If
#'   FALSE, an unadjusted rate ratio will be returned. Default = 0.8
#' @param normfun Character string. Specify with which function to standardize
#'   the data. Default = "sum"
#' @param normsubset Integer vector. Specify the indices of features that are to
#'   be used in standardization
#' @param log Logical or numeric. Specify whether to log-transform the rate
#'   ratio, and with what base. If TRUE, uses natural logarithm. Default = 2
#' @param belowxreads Logical or numeric. Set threshold above which the
#'   confidence limit is approximated. Default = 300
#'
#' @details The core of this function utilizes the \code{\link{poisson.test}}
#'   function from the stats package, with \code{x = c(t1, t0)} and \code{T =
#'   c(sumreads1, sumreads0)}. Above the belowxreads cutoff, the confidence
#'   interval is approximated as \code{logr +/- Z*sqrt(1/t1 + 1/t0)}. Here logr
#'   is the log-transformed rate ratio and Z is the number of standard units
#'   corresponding with the chosen confidence interval. If conf.level is not
#'   FALSE, the (log-transformed) confidence limit closest to 0 is returned. If
#'   the log-transformed upper and lower limit have opposite signs, i.e. the
#'   null hypothesis lies within the confidence interval, 0 is returned. If the
#'   function is performed on matrices for t0 and t1, then columns of t0 and t1
#'   are expected to be paired. See \code{\link{rrep}} for unpaired data.
#'
#' @return Returns the (log2-transformed) adjusted rate ratio.
#'
#' @note If conf.level is set to 0, it will be interpreted as FALSE. This is
#'   intended, as radjust will then return an unadjusted rate ratio estimate,
#'   instead of giving an error.
#'
#' @seealso \code{\link{rrep}}, \code{\link{nestedradjust}},
#'   \code{\link{getdeg}},\code{\link{oddscores}}, \code{\link{ess}},
#'   \code{\link{noness}}
#'
#' @author Jos B. Poell
#'
#' @examples
#' ut <- CRISPRsim(200, 4, a = c(3,3), allseed = 1, perfectseq = TRUE)
#' tr <- CRISPRsim(200, 4, a = c(3,3), e = TRUE, allseed = 1, perfectseq = TRUE)
#' cgi <- tr$d > -0.05 & tr$d < 0.05 & tr$e > -0.05 & tr$e < 0.05
#' r0 <- radjust(ut$t3, ut$t0, belowxreads = 300, normsubset = cgi)
#' r1 <- radjust(tr$t3, ut$t3, belowxreads = 300, normsubset = cgi)
#' hcr0 <- radjust(ut$t3, ut$t0, conf.level = 0.8, normsubset = cgi)
#' hcr1 <- radjust(tr$t3, ut$t3, conf.level = 0.8, normsubset = cgi)
#' plot(r0, hcr0)
#' abline(0, 1)
#' plot(r1, hcr1)
#' abline(0, 1)
#'
#' @export

radjust <- function(t1, t0, conf.level = 0.8, normfun = "sum", normsubset,
                    log = 2, belowxreads = 300) {
  fun <- get(normfun)
  if (length(ncol(t1))==0 || ncol(t1)==1) {
    if (!missing(normsubset)) {
      sr0 <- fun(t0[normsubset])
      sr1 <- fun(t1[normsubset])
    }
    else {
      sr0 <- fun(t0)
      sr1 <- fun(t1)
    }
  } else if (ncol(t1) != ncol(t0) || nrow(t1) != nrow(t0)) {
    stop("This function requires equal numbers of rows and columns for t0 and t1")
  } else {
    message("Note that this function expects paired columns for t0 and t1")
    if (!missing(normsubset)) {
      if (length(normsubset)==1) {
        sr0 <- rep(t0[normsubset, ], each = nrow(t0))
        sr1 <- rep(t1[normsubset, ], each = nrow(t1))
      } else {
        sr0 <- rep(apply(t0[normsubset, ], 2, fun),  each = nrow(t0))
        sr1 <- rep(apply(t1[normsubset, ], 2, fun),  each = nrow(t0))
      }
    }
    else {
      sr0 <- rep(apply(t0, 2, fun), each = nrow(t0))
      sr1 <- rep(apply(t1, 2, fun), each = nrow(t1))
    }
  } 
  
  Z <- qnorm(1-(1-conf.level)/2)

  r <- mapply(function(t1, t0, sr1, sr0) {
    if (belowxreads != FALSE && min(t1, t0) >= belowxreads) {
      logfc <- log((t1/sr1)/(t0/sr0))
      hi <- logfc + Z*sqrt(1/t1 + 1/t0)
      lo <- logfc - Z*sqrt(1/t1 + 1/t0)
      if (prod(hi,lo) < 0) {
        logr <- 0
      } else {
        logr  <- c(hi,lo)[which.min(abs(c(hi,lo)))]
      }
    } else {
      if (conf.level != FALSE) {
        # log-transformed upper and lower limit of given confidence interval
        logci <- log(poisson.test(x = c(t1, t0),
                                   T = c(sr1, sr0),
                                   conf.level = conf.level)$conf.int)
        # when confidence interval contains 0, log2 rate ratio becomes 0
        if(prod(logci) > 0) {logr  <- logci[which.min(abs(logci))]} else {logr <- 0}
      } else {
        # basically the same as belowxreads
        logr <- log((t1/sr1)/(t0/sr0))
      }
    }
    if(log == FALSE) {
      r <- exp(logr)
      return(r)
    } else if (log == TRUE) {
      return(logr)
    } else {
      return(logr/log(log))
    }
  }, t1, t0, sr1, sr0)
  
  if(length(ncol(t1)) != 0) {
    r <- matrix(r, ncol = ncol(t1))
  }
  
  return(r)
}

#' Calculate rate ratios of rate ratios restricted by confidence level
#'
#' nestedradjust is an extension of radjust. It is specifically applicable for
#' comparison of two samples with an independent time base.
#'
#' @param mt1 Integer vector. Raw sequencing reads in test sample, t1
#' @param wt1 Integer vector. Raw sequencing reads in control sample, t1
#' @param mt0 Integer vector. Raw sequencing reads in test sample, t0
#' @param wt0 Integer vector. Raw sequencing reads in control sample, t0
#' @param conf.level Numeric. Sets the confidence level of the rate ratio. If
#'   FALSE, an unadjusted rate ratio will be returned. Default = 0.8
#' @param normfun Character string. Specify with which function to standardize
#'   the data. Default = "sum"
#' @param normsubset Integer vector. Specify the indices of features that are to
#'   be used in standardization
#' @param log Logical. Specify whether the rate ratio should be log-transformed.
#'   If TRUE, uses natural logarithm. Default = 2
#' @param belowxreads Logical or numeric. Set threshold above which the
#'   unadjusted rate ratio is returned. Default = 300
#'
#' @details The basic functionality is similar to \code{\link{radjust}}. The
#'   major difference lies in the time base factor. Because both samples have an
#'   independent origin, there is added uncertainty of this time base factor.
#'   This is where the nested aspect of this function comes in. Both the upper
#'   and lower confidence limit of the rate ratio at baseline are tested as time
#'   base. From the resulting two rate ratios, the smallest fold change ratio is
#'   returned. It therefore estimates a conservative rate ratio.
#'
#' @note The function will run but give a warning when supplying a matrix per
#'   experimental arm. For most experimental set-ups, replicates between mutant
#'   and wild-type will not be paired, and \code{\link{rrep}} will have to be
#'   used. An example of when you can supply a matrix of data is when you are
#'   analyzing isogenic lines in different parental lines, and each column
#'   represents a different parental line.
#'
#' @return Returns the (log2-transformed) adjusted rate ratio.
#'
#' @seealso \code{\link{radjust}}, \code{\link{doublejar}}, \code{\link{ess}},
#'   \code{\link{noness}}
#'
#' @author Jos B. Poell
#'
#' @examples
#' wt <- CRISPRsim(genes = 10, guides = 4, a = 3, allseed = 1, t0seed = 2, perfectseq = TRUE)
#' mt <- CRISPRsim(genes = 10, guides = 4, a = 3, e = TRUE, f = jitter(wt$f),
#'                 allseed = 1, repseed = 2, perfectseq = TRUE)
#' r <- nestedradjust(mt$t3, wt$t3, mt$t0, wt$t0)
#' plot(mt$e, r)
#'
#' @export

nestedradjust <- function(mt1, wt1, mt0, wt0, conf.level = 0.8, normfun = "sum", 
                          normsubset, log = 2, belowxreads = 300) {
  fun <- get(normfun)
  if (length(ncol(mt1))==0 || ncol(mt1)==1) {
    if (!missing(normsubset)) {
      srmt1 <- fun(mt1[normsubset])
      srwt1 <- fun(wt1[normsubset])
      srmt0 <- fun(mt0[normsubset])
      srwt0 <- fun(wt0[normsubset])
    }
    else {
      srmt1 <- fun(mt1)
      srwt1 <- fun(wt1)
      srmt0 <- fun(mt0)
      srwt0 <- fun(wt0)
    }
  } else if (any(c(dim(mt1) != dim(mt0), dim(wt1) != dim(mt0), dim(wt0) != dim(mt0)))) {
    stop("This function requires equal numbers of rows and columns")
  } else {
    warning("Columns are considered paired for this analysis, even between mt and wt")
    if (!missing(normsubset)) {
      if (length(normsubset)==1) {
        srmt0 <- rep(mt0[normsubset, ], each = nrow(mt0))
        srmt1 <- rep(mt1[normsubset, ], each = nrow(mt1))
        srwt0 <- rep(wt0[normsubset, ], each = nrow(wt0))
        srwt1 <- rep(wt1[normsubset, ], each = nrow(wt1))
      } else {
        srmt0 <- rep(apply(mt0[normsubset, ], 2, fun),  each = nrow(mt0))
        srmt1 <- rep(apply(mt1[normsubset, ], 2, fun),  each = nrow(mt1))
        srwt0 <- rep(apply(wt0[normsubset, ], 2, fun),  each = nrow(wt0))
        srwt1 <- rep(apply(wt1[normsubset, ], 2, fun),  each = nrow(wt1))
      }
    }
    else {
      srmt0 <- rep(apply(mt0, 2, fun), each = nrow(mt0))
      srmt1 <- rep(apply(mt1, 2, fun), each = nrow(mt1))
      srwt0 <- rep(apply(wt0, 2, fun), each = nrow(wt0))
      srwt1 <- rep(apply(wt1, 2, fun), each = nrow(wt1))
    }
  } 
  
  Z <- qnorm(1-(1-conf.level)/2)
  
  r <- mapply(function(mt1, wt1, mt0, wt0, srmt1, srwt1, srmt0, srwt0) {
    # checking whether belowxreads condition is met
    if (belowxreads != FALSE && min(mt1, wt1, mt0, wt0) >= belowxreads) {
      # calculate the sum reads ratio
      
      srr <- (srmt0/srmt1)*(srwt1/srwt0)
      
      logfc <- log(srr*(mt1/wt1)/(mt0/wt0))
      hi <- logfc + Z*sqrt(1/mt1 + 1/mt0 + 1/wt1 + 1/wt0)
      lo <- logfc - Z*sqrt(1/mt1 + 1/mt0 + 1/wt1 + 1/wt0)
      if (prod(hi,lo) < 0) {
        logr <- 0
      } else {
        logr  <- c(hi,lo)[which.min(abs(c(hi,lo)))]
      }
    } else {
      if (conf.level != FALSE) {
        srr1 <- srmt1/srwt1
        # calculate the confidence interval of the rate ratio at t0
        ci0 <- poisson.test(x = c(mt0, wt0),
                            T = c(srmt0, srwt0),
                            conf.level = conf.level)$conf.int
        # rate ratio should not contain 0 or divided by zero
        if (ci0[1] == 0) {ci0[1] <- min(0.01, ci0[2]/100)}
        if (ci0[2] == Inf) {ci0[2] <- max(100, ci0[1]*100)}
        # calculate the confidence interval at t1 with the lower and higher
        # estimates of the rate ratio at t0
        logci1 <- log(poisson.test(x = c(mt1, wt1),
                                   T = c(ci0[1]*srr1, 1),
                                   conf.level = conf.level)$conf.int)
        if(prod(logci1) > 0) {logr1  <- logci1[which.min(abs(logci1))]} else {logr1 <- 0}
        
        logci2 <- log(poisson.test(x = c(mt1, wt1),
                                   T = c(ci0[2]*srr1, 1),
                                   conf.level = conf.level)$conf.int)
        if(prod(logci2) > 0) {logr2  <- logci2[which.min(abs(logci2))]} else {logr2 <- 0}
        
        # return the log-transformed rate ratio that is closest to 0
        logr <- c(logr1, logr2)[which.min(c(abs(logr1), abs(logr2)))]
      } else {
        logr <- log(srr*(mt1/wt1)/(mt0/wt0))
      }
    }

    if(log == FALSE) {
      r <- exp(logr)
      return(r)
    } else if (log == TRUE) {
      return(logr)
    } else {
      return(logr/log(log))
    }
  }, mt1, wt1, mt0, wt0, srmt1, srwt1, srmt0, srwt0)
  
  if(length(ncol(mt1)) != 0) {
    r <- matrix(r, ncol = ncol(mt1))
  }
  
  return(r)
}


#' Rate ratios after adding an artificial number to all features
#'
#' jar adds a specified number to all features of both input vectors, after
#' which it calculates rate ratios for all features. The function was designed
#' to prevent zeros in count data.
#'
#' @param t1 Numeric vector. Feature data in test sample
#' @param t0 Numeric vector. Feature data in control sample
#' @param n Numeric. Specify how much is added to each value before calculating
#'   rate ratios. Default = 5
#' @param log Logical or numeric. Specify whether to log-transform the rate
#'   ratio, and with what base. If TRUE, uses natural logarithm. Default = 2
#' @param normfun Character string. Specify with which function to standardize
#'   the data. Default = "sum"
#' @param normsubset Integer vector. Specify the indices of features that are to
#'   be used in standardization
#'
#' @return Returns a numeric vector of the same length as the input vectors with
#'   the (log2-transformed) rate ratios after adding a specified number to all
#'   features.
#'
#' @seealso \code{\link{radjust}}, \code{\link{doublejar}}, \code{\link{ess}},
#'   \code{\link{noness}}
#'
#' @author Jos B. Poell
#'
#' @export

jar <- function(t1, t0, n = 5, log = 2, normfun = "sum", normsubset) {
  t1 <- t1 + n
  t0 <- t0 + n
  fun <- get(normfun)
  
  if (length(ncol(t1))==0 || ncol(t1)==1) {
    if (!missing(normsubset)) {
      sr0 <- fun(t0[normsubset])
      sr1 <- fun(t1[normsubset])
    }
    else {
      sr0 <- fun(t0)
      sr1 <- fun(t1)
    }
  } else if (ncol(t1) != ncol(t0) || nrow(t1) != nrow(t0)) {
    stop("This function requires equal numbers of rows and columns for t0 and t1")
  } else {
    message("Note that this function expects paired columns for t0 and t1")
    if (!missing(normsubset)) {
      if (length(normsubset)==1) {
        sr0 <- rep(t0[normsubset, ], each = nrow(t0))
        sr1 <- rep(t1[normsubset, ], each = nrow(t1))
      } else {
        sr0 <- rep(apply(t0[normsubset, ], 2, fun),  each = nrow(t0))
        sr1 <- rep(apply(t1[normsubset, ], 2, fun),  each = nrow(t0))
      }
    }
    else {
      sr0 <- rep(apply(t0, 2, fun), each = nrow(t0))
      sr1 <- rep(apply(t1, 2, fun), each = nrow(t1))
    }
  } 
  
  r <- (t1/sr1)/(t0/sr0)
  
  if (log == TRUE) {
    r <- log(r)
  } else if (log != FALSE) {
    r <- log(r, log)
  }
  return(r)
}


#' Rate ratios of rate ratios after adding an artificial number to all features
#'
#' doublejar adds a specified number to all features of all input vectors, after
#' which it calculates rate ratios of the rate ratios for all features. The
#' function was designed to prevent zeros in count data.
#'
#' @param mt1 Numeric vector. Feature data in test sample of line m
#' @param wt1 Numeric vector. Feature data in test sample of line w
#' @param mt0 Numeric vector. Feature data in control sample of line m
#' @param wt0 Numeric vector. Feature data in control sample of line w
#' @param n Numeric. Specify how much is added to each value before calculating
#'   rate ratios of rate ratios. Default = 5
#' @param log Logical or numeric. Specify whether to log-transform the rate
#'   ratio, and with what base. If TRUE, uses natural logarithm. Default = 2
#' @param normfun Character string. Specify with which function to standardize
#'   the data. Default = "sum"
#' @param normsubset Integer vector. Specify the indices of features that are to
#'   be used in standardization
#'
#' @return Returns a numeric vector of the same length as the input vectors with
#'   the (log2-transformed) rate ratios of rate ratios after adding a specified
#'   number to all features.
#'
#' @seealso \code{\link{jar}}, \code{\link{nestedradjust}}, \code{\link{ess}},
#'   \code{\link{noness}}
#'
#' @author Jos B. Poell
#'
#' @export

doublejar <- function(mt1, wt1, mt0, wt0, n = 5, log = 2, normfun = "sum",
                      normsubset) {
  mt1 <- mt1 + n
  wt1 <- wt1 + n
  mt0 <- mt0 + n
  wt0 <- wt0 + n
  fun <- get(normfun)
  
  if (length(ncol(mt1))==0 || ncol(mt1)==1) {
    if (!missing(normsubset)) {
      srmt1 <- fun(mt1[normsubset])
      srwt1 <- fun(wt1[normsubset])
      srmt0 <- fun(mt0[normsubset])
      srwt0 <- fun(wt0[normsubset])
    }
    else {
      srmt1 <- fun(mt1)
      srwt1 <- fun(wt1)
      srmt0 <- fun(mt0)
      srwt0 <- fun(wt0)
    }
  } else if (any(c(dim(mt1) != dim(mt0), dim(wt1) != dim(mt0), dim(wt0) != dim(mt0)))) {
    stop("This function requires equal numbers of rows and columns")
  } else {
    warning("Columns are considered paired for this analysis, even between mt and wt")
    if (!missing(normsubset)) {
      if (length(normsubset)==1) {
        srmt0 <- rep(mt0[normsubset, ], each = nrow(mt0))
        srmt1 <- rep(mt1[normsubset, ], each = nrow(mt1))
        srwt0 <- rep(wt0[normsubset, ], each = nrow(wt0))
        srwt1 <- rep(wt1[normsubset, ], each = nrow(wt1))
      } else {
        srmt0 <- rep(apply(mt0[normsubset, ], 2, fun),  each = nrow(mt0))
        srmt1 <- rep(apply(mt1[normsubset, ], 2, fun),  each = nrow(mt1))
        srwt0 <- rep(apply(wt0[normsubset, ], 2, fun),  each = nrow(wt0))
        srwt1 <- rep(apply(wt1[normsubset, ], 2, fun),  each = nrow(wt1))
      }
    }
    else {
      srmt0 <- rep(apply(mt0, 2, fun), each = nrow(mt0))
      srmt1 <- rep(apply(mt1, 2, fun), each = nrow(mt1))
      srwt0 <- rep(apply(wt0, 2, fun), each = nrow(wt0))
      srwt1 <- rep(apply(wt1, 2, fun), each = nrow(wt1))
    }
  }
  
  r <- (srmt0/srmt1)*(srwt1/srwt0)*(mt1/wt1)/(mt0/wt0)
  
  if (log == TRUE) {
    r <- log(r)
  } else if (log != FALSE) {
    r <- log(r, log)
  }
  return(r)
}

#' Convert rate ratios to Z-values
#'
#' r2Z calculates the Z-values of rate ratios standardized to all rate ratios or
#' a specified normalization subset.
#'
#' @param r Character vector. List of rate ratios
#' @param normsubset Integer vector. Specify the indices of features that are to
#'   be used in standardization
#' @param method Character. Specify what method to use for standardization.
#'   Options are "median" and "mean". Default = "mean"
#'
#' @return Returns a numeric vector of the same length as the number of rate
#'   ratios
#'
#' @details While not strictly required, it is recommended to input
#'   log-transformed rate ratios. Standardization is by default done using the
#'   mean and standard deviation: \code{(r-mean(r))/sd(r)}. By specifying method = "median", it does
#'   so using the median and median absolute deviations: \code{(r-median(r))/mad(r)}. 
#'
#' @note Functions in the CSSA package by default return log-transformed rate
#'   ratios. These can directly be used as input for this function.
#'
#' @author Jos B. Poell
#'
#' @seealso \code{\link{sumZ}}, \code{\link{ess}}, \code{\link{noness}}
#'
#' @export

r2Z <- function(r, normsubset, method = "mean") {
  if (method == "median") {
    if (!missing(normsubset)) {
      Z <- (r-median(r[normsubset]))/mad(r[normsubset])
    } else {
      Z <- (r-median(r))/mad(r) 
    }
  } else {
    if (!missing(normsubset)) {
      Z <- (r-mean(r[normsubset]))/sd(r[normsubset])
    } else {
      Z <- (r-mean(r))/sd(r)
    }
  }
  return(Z)
}


#' Calculate corrected summed Z-values per gene
#'
#' sumZ adds Z-values of all guides per gene, then corrects each value by
#' dividing by the square root of the number of guides representing the gene
#'
#' @param guides Character vector. List of guides
#' @param Z Numeric vector. Z-values of guides
#'
#' @return Returns a numeric vector of the same length as the number of genes
#'   with the corrected summed Z-value.
#'
#' @note Z-values can be based on the entire population, but also on a
#'   population of control genes. See examples below
#'
#' @author Jos B. Poell
#' 
#' @seealso \code{\link{CRISPRsim}}, \code{\link{jar}}, \code{\link{radjust}}, \code{\link{r2Z}}
#'
#' @examples
#' ut <- CRISPRsim(5000, 4, a = c(3,3), allseed = 1, perfectsampling = TRUE)
#' tr <- CRISPRsim(5000, 4, a = c(3,3), e = TRUE, allseed = 1, perfectsampling = TRUE)
#' cgi <- which(tr$d > -0.05 & tr$d < 0.05 & tr$e > -0.05 & tr$e < 0.05)
#' r0 <- jar(ut$t6, ut$t0)
#' r1 <- jar(tr$t6, ut$t6)
#' Z0 <- r2Z(r0, normsubset = cgi)
#' Z1 <- r2Z(r1, normsubset = cgi)
#' sumZ0 <- sumZ(tr$guides, Z0)
#' sumZ1 <- sumZ(tr$guides, Z1)
#' reald <- rle(tr$d)$values
#' reale <- rle(tr$e)$values
#' plot(reald, sumZ0)
#' plot(reale, sumZ1)
#'
#' @export

sumZ <- function(guides, Z) {
  genes <- rle(gsub("_.*", "", guides))$values
  n <- rle(gsub("_.*", "", guides))$lengths
  return(mapply(function(n,gi) {sum(Z[(gi-n+1):gi])/sqrt(n)}, 
                           n, cumsum(n)))
}

# The function below has become quite the behemoth. It started out from the
# apply functions in my older code, that were not generalizable. I wanted to
# rewrite it into a single function that can get the d, e, and g scores for
# genes and guides, with optionally analysis of a second guide. It seemed like
# it should be a walk in the park with all the code already there, but it turned
# out to be very clunky. So I am sorry, it is not well-written, but as far as I
# can see the results are correct. This will then probably become the work horse
# of the analysis pipeline (right after rajdust). The only function missing as
# of this writing, is a function that estimates the values for a and b. Good
# luck!

#' Derive growth-modifying effect of gene knockout in pooled experiments
#'
#' getdeg was specifically designed to derive the effect of gene knockout on
#' cell growth based on results from pooled CRISPR-Cas9 experiments. Using a
#' combination of both rate ratios and (assumed or estimated) maximum population
#' doublings, the straight lethality and optionally sensitization / synthetic
#' lethality are calculated based on the "most efficacious guide targeting the
#' gene", i.e. the feature that shows the most extreme rate ratio change within
#' its group.
#'
#' @param guides Character vector. Guides are assumed to start with the gene
#'   name, followed by an underscore, followed by a number or sequence unique
#'   within that gene.
#' @param r0 Numeric vector. Log2-transformed rate ratios of features
#'   representing straight lethality.
#' @param r1 Numeric vector. Log2-transformed rate ratios of features
#'   representing sensitization or synthetic lethality. Optional but required to
#'   calculate e.
#' @param rt Numeric vector. Log2-transformed rate ratios of features
#'   representing lethality in the test sample. Optional.
#' @param a Numeric. Estimated potential population doublings between time
#'   points.
#' @param b Numeric. Estimated potential population doublings between time
#'   points in test sample. Only applicable if r1 is given. If omitted, assumed
#'   equal to a.
#' @param secondbest Logical. If TRUE, calculate effect sizes based on the
#'   second best guides of each gene as well. Default = TRUE
#' @param skipcutoff Logical or numeric. If specified, do not calculate effect
#'   sizes of genes with maximum absolute rate ratios below this cut-off.
#'   Default = FALSE
#' @param correctab Logical. When \code{a != b}, it is be possible (and
#'   necessary?) to mathematically correct for this difference. If you analyze
#'   an experiment with unequal a and b, try both with and without correction
#'   and read the notes below. Default = TRUE
#'
#' @details getdeg derives gene knockout effect sizes based on rate ratios.
#'   These in turn are derived from sequencing coverage of the features (e.g.
#'   guides in a CRISPR-bases screen). The function expects log2-transformed
#'   rate ratios. It also requires an estimate of potential population
#'   doublings, basically meaning that cells without successful knockout (or
#'   knockout of a gene without any growth-modifying effect) would have divided
#'   this many times. An example: the log2 rate ratio r0 of a guide between t1
#'   and t0 is -3, and the screen encompassed a = 6 doublings. If the guide was
#'   successful in all cells (g = 1), the effect of the corresponding gene
#'   knockout is d = r/a = -0.5. The function assumes the guide with the most
#'   extreme r with the same direction as the median r of all guides targeting
#'   that gene has this efficacy of 1, and then calculates g for the other
#'   guides. Things get more interesting when there is also a treatment effect.
#'   In this case it compares rate ratios of treated versus untreated and t1
#'   versus t0. From these it will decide which is the best guide and calculate
#'   both straight lethal effect d and sensitizing effect e. Although this will
#'   generally improve robustness, it is always a good idea to compare the
#'   results with the analysis of a single rate ratio measure. Optionally, but
#'   by default, the effects based on the second-best guide are also calculated.
#'   This function does not do anything in terms of statistics. It expects
#'   precautions are taken in the calculation of rate ratios!
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
#'   \item{j}{ - within-gene index of the second-best guide} }
#'
#' @note When comparing two experimental arms that have had different numbers of
#'   population doublings, things get quirky. I have put the mathematical
#'   correction in the function, which you can turn off with \code{correctab =
#'   FALSE}. I have noticed that correction gives straight lethal genes
#'   artificially high treatment resistance (positive e). But when I do not
#'   correct, I see a downward skew here. In case of a resistance screen, it may
#'   be more useful to look at the uncorrected variant. If you are interested in
#'   picking up sensitizers, I would recommend correcting.
#'
#' @seealso \code{\link{CRISPRsim}}, \code{\link{jar}}, \code{\link{radjust}}
#'
#' @author Jos B. Poell
#'
#' @examples
#' ut <- CRISPRsim(5000, 4, a = c(3,3), allseed = 1, perfectseq = TRUE)
#' tr <- CRISPRsim(5000, 4, a = c(3,3), e = TRUE, allseed = 1, perfectseq = TRUE)
#' cgi <- tr$d > -0.05 & tr$d < 0.05 & tr$e > -0.05 & tr$e < 0.05
#' r0 <- jar(ut$t6, ut$t0, normsubset = cgi)
#' r1 <- jar(tr$t6, ut$t6, normsubset = cgi)
#' deg <- getdeg(ut$guides, r0, r1, a = 6, b = 6, secondbest = FALSE)
#' reald <- rle(tr$d)$values
#' reale <- rle(tr$e)$values
#' plot(reald, deg$d)
#' plot(reale, deg$e)
#'
#' @export

getdeg <- function(guides, r0, r1, rt = FALSE, a, b, secondbest = TRUE,
                   skipcutoff = FALSE, correctab = TRUE) {
  
  if (missing(a)) {
    stop("enter the presumed number of population doublings")
  }
  # note that guides are presumed to be named [gene]_[number | sequence]
  genes <- rle(gsub("_.*", "", guides))$values
  n <- rle(gsub("_.*", "", guides))$lengths
  gi <- cumsum(n)
  
  if (missing(r1)) {
    message("no input to calculate effect modifier e")
    dgi <- mapply(function(n, gi) {
      # First, calculate the median log2 rate ratio of all guides targeting a
      # gene. Then, find the log2 rate ratio farthest from 0 in the same
      # direction as the median. In other words: the median has to be in the
      # correct direction!
      if (median(r0[(gi-n+1):gi]) < 0) {r <- min(r0[(gi-n+1):gi])} else {r <- max(r0[(gi-n+1):gi])}
      # Note the index of the guide with the selected rate ratio.
      i <- head(which(r0[(gi-n+1):gi]==r),1)
      # If a cutoff was set, the effect has to be more extreme than this cutoff
      # for the guides to be analyzed. If the cutoff is not reached, all guide
      # efficacies of this gene are set to 1.
      if (abs(r) < skipcutoff) {return(list(d=r/a, g = rep(1, n), i = i))} else {
        d <- r/a
        g <- sapply(seq_len(n), function(i) {
          # The extra if statement is necessary to prevent 0/0
          if (d == 0 && r0[gi-n+i] == 0) {g <- 1} else {g <- (2^(r0[gi-n+i])-1) / (2^(a*d)-1)}
          # Mathematically, g can end up below 0 or above 1 in some instances.
          # We cannot have that ...
          g[g < 0] <- 0
          g[g > 1] <- 1
          return(g)})
      }
      return(list(d=d,g=g,i=i))
    }, n, gi)
    if (secondbest == TRUE) {
      # return a list of within-gene indices all the best guides
      i <- unlist(dgi[3,])
      # number of guides per genes is reduced by one
      n2 <- n-1
      gi2 <- cumsum(n2)
      # best guides are excluded to find the second-best guide
      r02 <- r0[-(gi-n+i)]
      
      d2j <- mapply(function(n2, gi2) {
        # if there was only 1 guide to begin with, second best is NaN
        if (n2 == 0) {return(list(d2 = NaN, j = NaN))}
        else {
          if (median(r02[(gi2-n2+1):gi2]) < 0) {r <- min(r02[(gi2-n2+1):gi2])} else {r <- max(r02[(gi2-n2+1):gi2])}
          j <- tail(which(r02[(gi2-n2+1):gi2]==r),1)
          # note that guide efficacies are not calculated again
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
    degi <- mapply(function(n, gi) {
      # This time around, most extreme log2 rate ratios are found for all comparisons available
      if (median(r0[(gi-n+1):gi]) < 0) {r00 <- min(r0[(gi-n+1):gi])} else {r00 <- max(r0[(gi-n+1):gi])}
      if (median(r1[(gi-n+1):gi]) < 0) {r11 <- min(r1[(gi-n+1):gi])} else {r11 <- max(r1[(gi-n+1):gi])}
      # Since the rate ratio of the second arm compared to its own t0 can be
      # given in the argument rt, this part of the function has a lot of extra
      # ifs to find the right parameters... Note that this "rt" is especially
      # relevant in synthetic lethality experimental setups
      if (length(rt) > 1) {
        if (median(rt[(gi-n+1):gi]) < 0) {rtt <- min(rt[(gi-n+1):gi])} else {rtt <- max(rt[(gi-n+1):gi])}
        r <- c(abs(r00), abs(r11), abs(rtt))
        w <- which.max(r)
        if (w == 1) {i <- head(which(r0[(gi-n+1):gi]==r00),1)}
        else if (w == 2) {i <- head(which(r1[(gi-n+1):gi]==r11),1)}
        else {i <- head(which(rt[(gi-n+1):gi]==rtt),1)}
      } else {
        r <- c(abs(r00),abs(r11))
        w <- which.max(r)
        if (w == 1) {i <- head(which(r0[(gi-n+1):gi]==r00),1)}
        else {i <- head(which(r1[(gi-n+1):gi]==r11),1)}
      }
      if (r[w] < skipcutoff) {
        if (a != b && correctab == TRUE) {
          if (length(rt) > 1) {de <- rt[gi-n+i]/b}
          e <- (r1[gi-n+i]-r0[gi-n+i]*(b/a-1))/b
        } else {
          e <- r1[gi-n+i]/b
          if (length(rt) > 1) {de <- rt[gi-n+i]/b}
        }
        if (length(rt) > 1) {
          return(list(d=r0[gi-n+i]/a, e = e, de = de, g = rep(1, n), i = i))
        } else {
          return(list(d=r0[gi-n+i]/a, e = e, g = rep(1, n), i = i))
        }
      } else {
        d <- r0[gi-n+i]/a
        if (length(rt) > 1) {de <- rt[gi-n+i]/b}
        # Code below shows the correction for different number of potential
        # doublings. Although this gets the estimated e-values much closer to
        # the real e-values, I see there is a skew upwards when d decreases
        # (e.g. straight lethality increases).
        if (a != b && correctab == TRUE) {e <- (r1[gi-n+i]-r0[gi-n+i]*(b/a-1))/b}
        else {e <- r1[gi-n+i]/b}
        
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
      }
      if (length(rt) > 1) {
        return(list(d=d,e=e,de=de,g=g,i=i))
      } else {
        return(list(d=d,e=e,g=g,i=i))
      }
    }, n, gi)
    if (secondbest == TRUE) {
      if (length(rt) > 1) {
        i <- unlist(degi[5,])
      } else {
        i <- unlist(degi[4,])  
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
    d <- unlist(dgi[1,])
    g <- unlist(dgi[2,])
    i <- unlist(dgi[3,])
    if (secondbest == TRUE) {
      d2 <- unlist(d2j[1,])
      return(list(genes=genes,n=n,d=d,d2=d2,g=g,i=i,j=j))
    } else {
      return(list(genes=genes,n=n,d=d,g=g,i=i))
    }
  } else {
    d <- unlist(degi[1,])
    e <- unlist(degi[2,])
    if (length(rt) > 1) {
      de <- unlist(degi[3,])
      g <- unlist(degi[4,])
      i <- unlist(degi[5,])  
    } else {
      g <- unlist(degi[3,])
      i <- unlist(degi[4,])
    }
    
    if (secondbest == TRUE) {
      d2 <- unlist(d2e2j[1,])
      e2 <- unlist(d2e2j[2,])
      if (length(rt) > 1) {
        de2 <- unlist(d2e2j[3,])
        return(list(genes=genes,n=n,d=d,d2=d2,e=e,e2=e2,de=de,de2=de2,g=g,i=i,j=j))
      } else {
        return(list(genes=genes,n=n,d=d,d2=d2,e=e,e2=e2,g=g,i=i,j=j))
      }
      
    } else {
      if (length(rt) > 1) {
        return(list(genes=genes,n=n,d=d,e=e,de=de,g=g,i=i))
      } else {
        return(list(genes=genes,n=n,d=d,e=e,g=g,i=i))
      }
      
    }
  }
  
}

#' Get a list of essential genes or corresponding indices in the data set
#' 
#' ess couples reported essential genes to a data set. Without arguments, it 
#' returns the list of essential genes. Otherwise, it returns the indices 
#' corresponding to essential genes.
#' 
#' @param guides Character vector. List of guide names. Names of essential genes
#'   are matched to the string before what is specified in the trnc argument.
#' @param genes Character vector. List of gene names. Names of essential genes 
#'   are matched to the given gene names. Matches must be exact.
#' @param trnc Character. Regular expression pattern to isolate gene name from 
#'   the guide name. Default = "_.*"
#'   
#' @return Returns a character vector with names of essential genes, or a 
#'   numeric vector with indices representing essential genes in the given data 
#'   set.
#'   
#' @note The list of essential genes was retrieved from
#'   http://hart-lab.org/downloads.
#'   
#' @seealso \code{\link{noness}}
#'   
#' @author Jos B. Poell
#'   
#' @export

ess <- function(guides, genes, trnc = "_.*") {
  if (missing(guides) && missing(genes)) {
    return(essentials)
  } else if (!missing(guides)) {
    return(which(gsub(trnc, "", guides) %in% essentials))
  } else {
    return(which(genes %in% essentials))
  }
}

#' Get a list of nonessential genes or corresponding indices in the data set
#' 
#' noness couples reported nonessential genes to a data set. Without arguments, it 
#' returns the list of nonessential genes. Otherwise, it returns the indices 
#' corresponding to nonessential genes.
#' 
#' @param guides Character vector. List of guide names. Names of nonessential genes
#'   are matched to the string before what is specified in the trnc argument.
#' @param genes Character vector. List of gene names. Names of nonessential genes 
#'   are matched to the given gene names. Matches must be exact.
#' @param trnc Character. Regular expression pattern to isolate gene name from 
#'   the guide name. Default = "_.*"
#'   
#' @return Returns a character vector with names of nonessential genes, or a 
#'   numeric vector with indices representing nonessential genes in the given data 
#'   set.
#'   
#' @note The list of nonessential genes was retrieved from
#'   http://hart-lab.org/downloads.
#'   
#' @seealso \code{\link{ess}}
#'   
#' @author Jos B. Poell
#'   
#' @export

noness <- function(guides, genes, trnc = "_.*") {
  if (missing(guides) && missing(genes)) {
    return(nonessentials)
  } else if (!missing(guides)) {
    return(which(gsub(trnc, "", guides) %in% nonessentials))
  } else {
    return(which(genes %in% nonessentials))
  }
}
