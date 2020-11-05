#' Clonal and temporal dissection
#'
#' @param patient patient describes the TCGA identifier
#' @param complete.mutation.table This is the mutation table in matrix format.
#' @param purity ASCAT purity estimate [range, 0.01,1]
#'
#' @return A data.frame.
#' @export
#' @description
#' This function implementents the cancer cell fraction analysis, as well as temporal dissection of mutations
#' @examples
#'
#' @import bootstrap
#' @import boot
#'

earlyORlate <- function(patient, complete.mutation.table, purity) {

  # Get the table ready, with only information for specific patient
  mut.table <- complete.mutation.table[complete.mutation.table$patient == patient, ]

  # and all the other stuff
  cellularity <- as.numeric(purity)
  major.cn <- unlist(mut.table$major_cn)
  abs.cn <- unlist(mut.table$minor_cn) + unlist(mut.table$major_cn)
  depth.t <- unlist(mut.table$ref_counts) + unlist(mut.table$var_counts)
  max.cn <- max(abs.cn)
  VAF <- unlist(mut.table$var_counts) / (unlist(mut.table$var_counts) + unlist(mut.table$ref_counts))

  # Estimate theoretical VAFs for each type of copy number
  
  # types.matrix() is not available in sequenza 3.0.
  # By checking sequenza::baf.types.matrix() and sequenza::mufreq.types.matrix(),
  # It seems like this part is using sequenza::mufreq.types.matrix() to get the data.frame
  # contains "Mt".
  
  types <- sequenza::mufreq.types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 2)
  types.xy <- sequenza::mufreq.types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 1)
  types <- rbind(types, types.xy)
  types <- types[types$Mt >= 1, ]
  types$F <- 0
  for (i in 1:nrow(types)) {
    types$F[i] <- sequenza::theoretical.mufreq(
      cellularity = cellularity,
      CNn = types$CNn[i], CNt = types$CNt[i],
      Mt = types$Mt[i]
    )
  }

  # Let's create some functions that can estimate whether early or late

  get.Mt <- function(F, depth.t, types, CNt, CNn, Mt) {
    types <- types[types$CNn == CNn, ]
    l <- sequenza:::mufreq.dpois(
      mufreq = F, types$F[types$CNt == CNt & types$Mt <= Mt],
      depth.t = depth.t
    )
    l <- l / sum(l)
    L <- data.frame(l = l, Mt = types$Mt[types$CNt == CNt & types$Mt <= Mt])
    return(L)
  }

  get.conf <- function(F, depth.t) {
    conf.int <- data.frame(
      prop.test(round(F * depth.t, 0), depth.t)$conf[1],
      prop.test(round(F * depth.t, 0), depth.t)$conf[2]
    )
    return(conf.int)
  }

  bootstrap.cf <- function(Vaf, cellularity, CNn, CNt, depth.t) {
    # print(i)
    if (Vaf == 1) {
      conf.int <- data.frame(
        prop.test(round(Vaf * depth.t, 0), depth.t)$conf[1],
        prop.test(round(Vaf * depth.t, 0), depth.t)$conf[2]
      )

      lower <- get.mut.mult(Vaf = conf.int[1], cellularity = cellularity, CNt = CNt, CNn = CNn)
      higher <- get.mut.mult(Vaf = conf.int[2], cellularity = cellularity, CNt = CNt, CNn = CNn)
      conf.int <- data.frame(lower, higher)
      return(conf.int)
    }

    x <- c(rep(1, round(Vaf * depth.t, 0)), rep(0, (depth.t - round(Vaf * depth.t, 0))))
    theta <- function(x, i) {
      data <- x[i]
      est <- sum(data) / length(data)
      mut.multi <- (est * 1 / cellularity) * ((cellularity * CNt) + CNn * (1 - cellularity))
      return(mut.multi)
    }

    bt.res <- boot(x, theta, R = 1000)
    bt.ci <- boot.ci(bt.res, type = "norm")
    out <- c(bt.ci$normal[2], bt.ci$normal[3])

    return(out)
  }


  get.mut.mult <- function(CNt, Vaf, cellularity, CNn) {
    return((Vaf * 1 / cellularity) * ((cellularity * CNt) + CNn * (1 - cellularity)))
  }


  get.cancer.cell.fraction <- function(Max.Likelihood, mut.mult) {
    predicted.Mtn <- Max.Likelihood[, "Mt"]
    ccf <- mut.mult / predicted.Mtn
    return(ccf)
  }

  absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number) {
    f.function <- function(c, purity, local.copy.number) {
      return((purity * c) / (2 * (1 - purity) + purity * local.copy.number))
    }
    x <- dbinom(n.alt, depth, prob = sapply(seq(0.01, 1, length.out = 100), f.function, purity, local.copy.number))
    if (min(x) == 0) {
      x[length(x)] <- 1
    }

    names(x) <- seq(0.01, 1, length.out = 100)
    sub.cint <- function(x, prob = 0.95, n.alt, depth) {
      xnorm <- x / sum(x)
      xsort <- sort(xnorm, decreasing = TRUE)
      xcumLik <- cumsum(xsort)
      n <- sum(xcumLik < prob) + 1
      LikThresh <- xsort[n]
      cint <- x[xnorm >= LikThresh]
      all <- as.numeric(names(x))
      cellu <- as.numeric(names(cint))
      l.t <- cellu[1]
      r.t <- cellu[length(cellu)]
      m <- cellu[which.max(cint)]

      prob.subclonal <- sum(xnorm[1:90]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
      prob.clonal <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val

      res.df <- data.frame(left = l.t, est = m, right = r.t, prob.subclonal = prob.subclonal, prob.clonal = prob.clonal)
      return(res.df)
    }


    return(sub.cint(x, n.alt = n.alt, depth = depth))
  }

  absolute.cancer.cell.mtcpn <- function(n.alt, depth, purity, local.copy.number, mtcpn) {
    f.function <- function(c, purity, local.copy.number) {
      return((purity * c) / (2 * (1 - purity) + purity * local.copy.number))
    }
    x <- dbinom(n.alt, depth, prob = sapply(seq(0.01, as.numeric(mtcpn), length.out = (100)), f.function, purity, local.copy.number))
    if (min(x) == 0) {
      x[length(x)] <- 1
    }

    names(x) <- seq(0.01, 1, length.out = 100)
    sub.cint <- function(x, prob = 0.95, n.alt, depth) {
      xnorm <- x / sum(x)
      xsort <- sort(xnorm, decreasing = TRUE)
      xcumLik <- cumsum(xsort)
      n <- sum(xcumLik < prob) + 1
      LikThresh <- xsort[n]
      cint <- x[xnorm >= LikThresh]
      all <- as.numeric(names(x))
      cellu <- as.numeric(names(cint))
      l.t <- cellu[1]
      r.t <- cellu[length(cellu)]
      m <- cellu[which.max(cint)]

      prob.subclonal <- sum(xnorm[1:90]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
      prob.clonal <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val

      res.df <- data.frame(left = l.t, est = m, right = r.t, prob.subclonal = prob.subclonal, prob.clonal = prob.clonal)
      return(res.df)
    }


    return(sub.cint(x, n.alt = n.alt, depth = depth))
  }
  # add an absolute estimate of the cancer cell fraction

  get.all.mut.info <- function(i) {
    # debug:
    # The output returned may contain complex row.names and colnames,
    # which makes the combined data.frame is not aligned.
    FIXED_COLNAMES <- c(
      "obs.VAF",
      "mut.conf.0.05",
      "mut.conf.0.95",
      "mut.multi",
      "mut.multi.0.05",
      "mut.multi.bstr.0.05",
      "mut.multi.0.95",
      "mut.multi.bstr.0.95",
      "Exp.Cpn",
      "Exp.Cpn.Likelihood",
      "ccf",
      "ccf.0.05",
      "ccf.btstr.0.05",
      "ccf.0.95",
      "ccf.btstr.0.95",
      "absolute.ccf",
      "absolute.ccf.0.05",
      "absolute.ccf.0.95",
      "prob.subclonal",
      "prob.clonal",
      "timing"
    )
    # print(i)
    # First estimate the VAF confidence intervals
    obs.VAF <- VAF[i]
    mut.conf.0.05 <- get.conf(F = VAF[i], depth.t = depth.t[i])[1]
    mut.conf.0.95 <- get.conf(F = VAF[i], depth.t = depth.t[i])[2]

    if (abs.cn[i] == 0) {
      output <- data.frame(
        unlist(obs.VAF),
        unlist(mut.conf.0.05),
        unlist(mut.conf.0.95),
        mut.multi = NA,
        mut.multi.0.05 = NA,
        mut.multi.bstr.0.05 = NA,
        mut.multi.0.95 = NA,
        mut.multi.bstr.0.95 = NA,
        Exp.Cpn = NA,
        Exp.Cpn.Likelihood = NA,
        ccf = NA,
        ccf.0.05 = NA,
        ccf.btstr.0.05 = NA,
        ccf.0.95 = NA,
        ccf.btstr.0.95 = NA,
        absolute.ccf = NA,
        absolute.ccf.0.05 = NA,
        absoltue.ccf.0.95 = NA,
        prob.subclonal = NA,
        prob.clonal = NA,
        timing = "Not.Poss",
        row.names = NULL
      )
      colnames(output) <- FIXED_COLNAMES
      return(output)
    }

    # Next estimate the likelihood relating to which copy number the mutation has
    L <- get.Mt(
      F = VAF[i],
      depth.t = depth.t[i], CNt = abs.cn[i],
      types = types, CNn = unlist(mut.table$normal_cn[i]),
      Mt = major.cn[i]
    )

    # Next determine the mut multiplicity
    mut.multi <- get.mut.mult(CNt = abs.cn[i], Vaf = VAF[i], cellularity = cellularity, CNn = unlist(mut.table$normal_cn[i]))
    mut.multi.0.05 <- get.mut.mult(CNt = abs.cn[i], Vaf = mut.conf.0.05, cellularity = cellularity, CNn = unlist(mut.table$normal_cn[i]))
    mut.multi.0.95 <- get.mut.mult(CNt = abs.cn[i], Vaf = mut.conf.0.95, cellularity = cellularity, CNn = unlist(mut.table$normal_cn[i]))
    mut.multi.bstr <- bootstrap.cf(Vaf = VAF[i], cellularity = cellularity, CNn = unlist(mut.table$normal_cn[i]), CNt = abs.cn[i], depth.t = depth.t[i])
    mut.multi.bstr.0.05 <- mut.multi.bstr[1]
    mut.multi.bstr.0.95 <- mut.multi.bstr[2]

    if (is.na(L$l)[1]) {
      output <- data.frame(obs.VAF,
        unlist(mut.conf.0.05),
        unlist(mut.conf.0.95),
        unlist(mut.multi),
        unlist(mut.multi.0.05),
        unlist(mut.multi.bstr.0.05),
        unlist(mut.multi.0.95),
        unlist(mut.multi.bstr.0.95),
        Exp.Cpn = NA,
        Exp.Cpn.Likelihood = NA,
        ccf = NA,
        ccf.0.05 = NA,
        ccf.btstr.0.05 = NA,
        ccf.0.95 = NA,
        ccf.btstr.0.95 = NA,
        absolute.ccf = NA,
        absolute.ccf.0.05 = NA,
        absoltue.ccf.0.95 = NA,
        prob.subclonal = NA,
        prob.clonal = NA,
        timing = "Not.Poss",
        row.names = NULL
      )
      colnames(output) <- FIXED_COLNAMES
      return(output)
    }

    # Now determine which likelihood should be used
    Max.Likelihood <- L[which.max(L$l), ]
    absolute.calc <- absolute.cancer.cell.fraction(n.alt = unlist(mut.table$var_counts)[i], depth = depth.t[i], purity = cellularity, local.copy.number = abs.cn[i])
    absolute.ccf.0.05 <- absolute.calc[1]
    absolute.ccf.0.95 <- absolute.calc[3]
    absolute.ccf <- absolute.calc[2]
    prob.subclonal <- absolute.calc[4]
    prob.clonal <- absolute.calc[5]




    # Next determine the cancer cell fraction
    ccf <- get.cancer.cell.fraction(Max.Likelihood, mut.multi)
    ccf.0.05 <- get.cancer.cell.fraction(Max.Likelihood, mut.multi.0.05)
    ccf.btstr.0.05 <- get.cancer.cell.fraction(Max.Likelihood, mut.multi.bstr.0.05)

    ccf.0.95 <- get.cancer.cell.fraction(Max.Likelihood, mut.multi.0.95)
    ccf.btstr.0.95 <- get.cancer.cell.fraction(Max.Likelihood, mut.multi.bstr.0.95)

    # Next determine the late cancer cell fraction
    # Make sure you also output the theoretical Copy (i.e. what it's closest to using maximum likelihood)
    expected.copy <- Max.Likelihood[2]


    # Finally also make a suggestion about whether the mutation is early late or not possible to tell


    if (Max.Likelihood$Mt > 1) {
      timing <- "early"
    }

    if (Max.Likelihood$Mt <= 1) {
      timing <- "late"
    }

    if (major.cn[i] <= 1) {
      timing <- "Not.Poss"
    }




    # Let's put this all together and output it
    output <- data.frame(
      unlist(obs.VAF),
      unlist(mut.conf.0.05),
      unlist(mut.conf.0.95),
      unlist(mut.multi),
      unlist(mut.multi.0.05),
      unlist(mut.multi.bstr.0.05),
      unlist(mut.multi.0.95),
      unlist(mut.multi.bstr.0.95),
      Exp.Cpn = unlist(Max.Likelihood$Mt),
      Exp.Cpn.Likelihood = unlist(Max.Likelihood$l),
      unlist(ccf),
      unlist(ccf.0.05),
      unlist(ccf.btstr.0.05),
      unlist(ccf.0.95),
      unlist(ccf.btstr.0.95),
      unlist(absolute.ccf),
      unlist(absolute.ccf.0.05),
      unlist(absolute.ccf.0.95),
      unlist(prob.subclonal),
      unlist(prob.clonal),
      unlist(timing),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    colnames(output) <- FIXED_COLNAMES
    return(output)
  }

  # output <- t(sapply(1:nrow(mut.table), get.all.mut.info))
  # output <- data.frame(output, stringsAsFactors = FALSE)

  # get.all.mut.info() returns a data.frame
  # The code above will generate a data.frame of list, which is not best practice.
  output.list <- lapply(1:nrow(mut.table), get.all.mut.info)
  output <- dplyr::bind_rows(output.list)
  output <- as.data.frame(output)

  out <- cbind(mut.table, output)
  return(out)
}
