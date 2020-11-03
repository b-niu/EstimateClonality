#### Functions for cancer cell fraction estimation #####

identify.mut.copy.number.ascat <- function(x,sub.mat.mut, sub.mat.copy)
{
  
  mut                  <- sub.mat.mut[x,,drop=FALSE]
  ww                   <- which(as.numeric(sub.mat.copy$Chr)==as.numeric(mut$Chr)
                                &as.numeric(sub.mat.copy$Start)<=as.numeric(mut$Start_position)
                                &as.numeric(sub.mat.copy$End)>=as.numeric(mut$Start_position))
  copy                 <- sub.mat.copy[ww,,drop=FALSE]
  
  
  mutation_id         <- paste(mut$Patient,mut$Chr,mut$Start_position,mut$Reference,sep=":")
  ref_counts          <- mut$Ref_freq
  var_counts          <- mut$Variant_freq
  normal_cn           <- 2
  patient             <- mut$Patient
  Reference_Base      <- mut$Reference
  Alternate_Base      <- mut$Alternate
  
  if(nrow(copy)!=1)
  {
    minor_cn          <- NA
    major_cn          <- NA
    output            <-  data.frame(mutation_id
                                     ,ref_counts
                                     ,var_counts
                                     ,normal_cn
                                     ,minor_cn
                                     ,major_cn
                                     ,patient
                                     ,Reference_Base
                                     ,Alternate_Base
                                     ,stringsAsFactors=FALSE)
    return(output)
  }
  
  minor_cn            <- min(c(copy$nA,copy$nB))
  major_cn            <- max(c(copy$nA,copy$nB))
  
  
  output              <- data.frame(mutation_id
                                    ,ref_counts
                                    ,var_counts
                                    ,normal_cn
                                    ,minor_cn
                                    ,major_cn
                                    ,patient
                                    ,Reference_Base
                                    ,Alternate_Base
                                    ,stringsAsFactors=FALSE)
  
  return(output)
  
  
  
  
}

#######

earlyORlate <- function(patient,complete.mutation.table,purity){
  
  
  # First of all load the needed packaged
  # 
  suppressPackageStartupMessages(library(sequenza))
  suppressPackageStartupMessages(library(bootstrap))
  suppressPackageStartupMessages(library(boot))
  
  # Get the table ready, with only information for specific patient
  mut.table  <- complete.mutation.table[complete.mutation.table$patient==patient,]
  
  # and all the other stuff
  cellularity <- as.numeric(purity)
  major.cn <- unlist(mut.table$major_cn)
  abs.cn   <- unlist(mut.table$minor_cn) + unlist(mut.table$major_cn)
  depth.t  <- unlist(mut.table$ref_counts) + unlist(mut.table$var_counts)
  max.cn   <- max(abs.cn)
  VAF      <- unlist(mut.table$var_counts)/(unlist(mut.table$var_counts)+unlist(mut.table$ref_counts))
  
  # Estimate theoretical VAFs for each type of copy number
  types    <- types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 2)
  types.xy <- types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 1)
  types    <- rbind(types, types.xy)
  types    <- types[types$Mt >= 1, ]
  types$F  <- 0
  for (i in 1:nrow(types)) {
    types$F[i] <- theoretical.mufreq(cellularity = cellularity,
                                     CNn = types$CNn[i], CNt = types$CNt[i],
                                     Mt = types$Mt[i])
  }
  
  # Let's create some functions that can estimate whether early or late
  
  get.Mt <- function(F, depth.t, types, CNt, CNn, Mt){
    types <- types[types$CNn == CNn, ]
    l <- sequenza:::mufreq.dpois(mufreq = F, types$F[types$CNt== CNt&types$Mt<=Mt],
                                 depth.t = depth.t)
    l <- l/sum(l)
    L <- data.frame(l = l, Mt = types$Mt[types$CNt== CNt&types$Mt<=Mt])
  }
  
  get.conf <- function(F, depth.t){
    conf.int   <- cbind(prop.test(round(F*depth.t,0),depth.t)$conf[1]
                        ,prop.test(round(F*depth.t,0),depth.t)$conf[2])
    return(conf.int)
  }
  
  bootstrap.cf <- function(Vaf, cellularity, CNn, CNt, depth.t)
  {
    #print(i)
    if(Vaf==1)
    {
      conf.int   <- cbind(prop.test(round(Vaf*depth.t,0),depth.t)$conf[1]
                          ,prop.test(round(Vaf*depth.t,0),depth.t)$conf[2])
      
      lower      <- get.mut.mult(Vaf=conf.int[1],cellularity=cellularity,CNt=CNt,CNn=CNn)
      higher     <- get.mut.mult(Vaf=conf.int[2],cellularity=cellularity,CNt=CNt,CNn=CNn)
      conf.int   <- cbind(lower,higher)
      return(conf.int)
      
    }
    
    x          <- c(rep(1,round(Vaf*depth.t,0)),rep(0,(depth.t-round(Vaf*depth.t,0))))
    theta      <- function(x,i)
    {
      data      <- x[i]
      est       <- sum(data)/length(data)
      mut.multi <- (est *1/cellularity)*((cellularity*CNt)+CNn*(1-cellularity))
      return(mut.multi)
    }
    
    bt.res      <- boot(x,theta,R=1000)
    bt.ci       <- boot.ci(bt.res,type='norm')
    out         <- c(bt.ci$normal[2],bt.ci$normal[3])
    
    return(out)
    
  }
  
  
  get.mut.mult <- function(CNt,Vaf,cellularity,CNn)
  {
    
    return((Vaf *1/cellularity)*((cellularity*CNt)+CNn*(1-cellularity)))
    
  }
  
  
  get.cancer.cell.fraction <- function(Max.Likelihood,mut.mult)
  {
    predicted.Mtn   <- Max.Likelihood[,'Mt']
    ccf             <- mut.mult/predicted.Mtn
    return(ccf)
  }
  
  absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number)
  {
    f.function <- function (c,purity,local.copy.number)
    {
      
      return((purity*c) / (2*(1-purity) + purity*local.copy.number))
      
    }
    x              <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number))
    if(min(x)==0)
    {
      x[length(x)] <- 1
    }
    
    names(x)       <- seq(0.01,1,length.out=100)
    sub.cint <- function(x, prob = 0.95,n.alt,depth) {
      xnorm   <- x/sum(x)
      xsort   <- sort(xnorm, decreasing = TRUE)
      xcumLik <- cumsum(xsort)
      n = sum(xcumLik < prob) + 1
      LikThresh <- xsort[n]
      cint  <- x[xnorm >= LikThresh]
      all   <- as.numeric(names(x))
      cellu <- as.numeric(names(cint))
      l.t   <- cellu[1]
      r.t   <- cellu[length(cellu)]
      m     <- cellu[which.max(cint)]
      
      prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
      prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val
      
      data.frame(left = l.t, est = m, right = r.t,prob.subclonal=prob.subclonal,prob.clonal=prob.clonal)
    }
    
    
    return(sub.cint(x,n.alt=n.alt,depth=depth))
    
    
  }
  
  absolute.cancer.cell.mtcpn <- function(n.alt, depth, purity, local.copy.number, mtcpn)
  {
    f.function <- function (c,purity,local.copy.number)
    {
      
      return((purity*c) / (2*(1-purity) + purity*local.copy.number))
      
    }
    x              <- dbinom(n.alt,depth, prob=sapply(seq(0.01,as.numeric(mtcpn),length.out=(100)),f.function,purity,local.copy.number))
    if(min(x)==0)
    {
      x[length(x)] <- 1
    }
    
    names(x)       <- seq(0.01,1,length.out=100)
    sub.cint <- function(x, prob = 0.95,n.alt,depth) {
      xnorm   <- x/sum(x)
      xsort   <- sort(xnorm, decreasing = TRUE)
      xcumLik <- cumsum(xsort)
      n = sum(xcumLik < prob) + 1
      LikThresh <- xsort[n]
      cint  <- x[xnorm >= LikThresh]
      all   <- as.numeric(names(x))
      cellu <- as.numeric(names(cint))
      l.t   <- cellu[1]
      r.t   <- cellu[length(cellu)]
      m     <- cellu[which.max(cint)]
      
      prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
      prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val
      
      data.frame(left = l.t, est = m, right = r.t,prob.subclonal=prob.subclonal,prob.clonal=prob.clonal)
    }
    
    
    return(sub.cint(x,n.alt=n.alt,depth=depth))
    
    
  }
  # add an absolute estimate of the cancer cell fraction
  
  get.all.mut.info <- function(i)
  {
    #print(i)
    # First estimate the VAF confidence intervals
    obs.VAF         <- VAF[i]
    mut.conf.0.05   <- get.conf(F=VAF[i],depth.t=depth.t[i])[1]
    mut.conf.0.95   <- get.conf(F=VAF[i],depth.t=depth.t[i])[2]
    
    if(abs.cn[i]==0)
    {
      output          <- cbind(obs.VAF
                               ,mut.conf.0.05
                               ,mut.conf.0.95
                               ,mut.multi=NA
                               ,mut.multi.0.05=NA
                               ,mut.multi.bstr.0.05=NA
                               ,mut.multi.0.95=NA
                               ,mut.multi.bstr.0.95=NA
                               ,Exp.Cpn =NA
                               ,Exp.Cpn.Likelihood=NA
                               ,ccf=NA
                               ,ccf.0.05=NA
                               ,ccf.btstr.0.05=NA
                               ,ccf.0.95=NA
                               ,ccf.btstr.0.95=NA
                               ,absolute.ccf=NA
                               ,absolute.ccf.0.05=NA
                               ,absoltue.ccf.0.95=NA
                               ,prob.subclonal=NA
                               ,prob.clonal=NA
                               ,timing='Not.Poss')
      return(output)
    }
    
    # Next estimate the likelihood relating to which copy number the mutation has
    L <- get.Mt(F = VAF[i],
                depth.t = depth.t[i], CNt = abs.cn[i],
                types = types, CNn = unlist(mut.table$normal_cn[i])
                ,Mt=major.cn[i])
    
    # Next determine the mut multiplicity
    mut.multi            <- get.mut.mult(CNt=abs.cn[i],Vaf=VAF[i],cellularity=cellularity,CNn=unlist(mut.table$normal_cn[i]))
    mut.multi.0.05       <- get.mut.mult(CNt=abs.cn[i],Vaf=mut.conf.0.05,cellularity=cellularity,CNn=unlist(mut.table$normal_cn[i]))
    mut.multi.0.95       <- get.mut.mult(CNt=abs.cn[i],Vaf=mut.conf.0.95,cellularity=cellularity,CNn=unlist(mut.table$normal_cn[i]))
    mut.multi.bstr       <- bootstrap.cf(Vaf=VAF[i],cellularity=cellularity,CNn=unlist(mut.table$normal_cn[i]),CNt=abs.cn[i],depth.t=depth.t[i])
    mut.multi.bstr.0.05  <- mut.multi.bstr[1]
    mut.multi.bstr.0.95  <- mut.multi.bstr[2]
    
    if(is.na(L$l)[1])
    {
      output          <- cbind(obs.VAF
                               ,mut.conf.0.05
                               ,mut.conf.0.95
                               ,mut.multi
                               ,mut.multi.0.05
                               ,mut.multi.bstr.0.05
                               ,mut.multi.0.95
                               ,mut.multi.bstr.0.95
                               ,Exp.Cpn =NA
                               ,Exp.Cpn.Likelihood=NA
                               ,ccf=NA
                               ,ccf.0.05=NA
                               ,ccf.btstr.0.05=NA
                               ,ccf.0.95=NA
                               ,ccf.btstr.0.95=NA
                               ,absolute.ccf=NA
                               ,absolute.ccf.0.05=NA
                               ,absoltue.ccf.0.95=NA
                               ,prob.subclonal=NA
                               ,prob.clonal=NA
                               ,timing='Not.Poss')
      return(output)
    }
    
    # Now determine which likelihood should be used
    Max.Likelihood   <- L[which.max(L$l),]
    absolute.calc        <- absolute.cancer.cell.fraction(n.alt=unlist(mut.table$var_counts)[i],depth=depth.t[i],purity=cellularity,local.copy.number=abs.cn[i])
    absolute.ccf.0.05    <- absolute.calc[1]
    absolute.ccf.0.95    <- absolute.calc[3]
    absolute.ccf         <- absolute.calc[2]
    prob.subclonal       <- absolute.calc[4]
    prob.clonal          <- absolute.calc[5]
    
    
    
    
    # Next determine the cancer cell fraction
    ccf             <- get.cancer.cell.fraction(Max.Likelihood,mut.multi)
    ccf.0.05        <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.0.05)
    ccf.btstr.0.05  <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.bstr.0.05)
    
    ccf.0.95        <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.0.95)
    ccf.btstr.0.95  <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.bstr.0.95)
    
    # Next determine the late cancer cell fraction
    # Make sure you also output the theoretical Copy (i.e. what it's closest to using maximum likelihood)
    expected.copy     <- Max.Likelihood[2]
    
      
    # Finally also make a suggestion about whether the mutation is early late or not possible to tell
    
    
    if(Max.Likelihood$Mt>1)
    {
      timing    <- 'early'
    }
    
    if(Max.Likelihood$Mt<=1)
    {
      timing    <- 'late'
    }
    
    if(major.cn[i]<=1)
    {
      timing    <- 'Not.Poss'
    }
    
    
    
    
    # Let's put this all together and output it
    output          <- data.frame(obs.VAF
                                  ,mut.conf.0.05
                                  ,mut.conf.0.95
                                  ,mut.multi
                                  ,mut.multi.0.05
                                  ,mut.multi.bstr.0.05
                                  ,mut.multi.0.95
                                  ,mut.multi.bstr.0.95
                                  ,Exp.Cpn =Max.Likelihood$Mt
                                  ,Exp.Cpn.Likelihood=Max.Likelihood$l
                                  ,ccf
                                  ,ccf.0.05
                                  ,ccf.btstr.0.05
                                  ,ccf.0.95
                                  ,ccf.btstr.0.95
                                  ,absolute.ccf
                                  ,absolute.ccf.0.05
                                  ,absolute.ccf.0.95
                                  ,prob.subclonal
                                  ,prob.clonal
                                  ,timing
                                  ,stringsAsFactors=FALSE)
    
    #output           <- data.frame(output,stringsAsFactors=FALSE)
    return(output)
    
  }
  
  output <- t(sapply(1:nrow(mut.table),get.all.mut.info))
  output <- data.frame(output,stringsAsFactors=FALSE)
  
  colnames(output) <-           c('obs.VAF'
                                ,'mut.conf.0.05'
                                ,'mut.conf.0.95'
                                ,'mut.multi'
                                ,'mut.multi.0.05'
                                ,'mut.multi.bstr.0.05'
                                ,'mut.multi.0.95'
                                ,'mut.multi.bstr.0.95'
                                ,'Exp.Cpn' 
                                ,'Exp.Cpn.Likelihood'
                                ,'ccf'
                                ,'ccf.0.05'
                                ,'ccf.btstr.0.05'
                                ,'ccf.0.95'
                                ,'ccf.btstr.0.95'
                                ,'absolute.ccf'
                                ,'absolute.ccf.0.05'
                                ,'absolute.ccf.0.95'
                                ,'prob.subclonal'
                                ,'prob.clonal'
                                ,'timing')
  
  
  out <- cbind(mut.table,output)
  return(out)
  
}

earlyORlate.strict <- function(patient,complete.mutation.table,purity){
  
  # The following function is very similar indeed to the function above
  # However, there is one important difference
  # For this function, copy number states where the minor allele =1 are not classified as late. 
  
  # First of all load the needed packaged
  # 
  suppressPackageStartupMessages(library(sequenza))
  suppressPackageStartupMessages(library(bootstrap))
  suppressPackageStartupMessages(library(boot))
  
  # Get the table ready, with only information for specific patient
  mut.table  <- complete.mutation.table[complete.mutation.table$patient==patient,]
  
  # and all the other stuff
  cellularity <- as.numeric(purity)
  major.cn <- unlist(mut.table$major_cn)
  abs.cn   <- unlist(mut.table$minor_cn) + unlist(mut.table$major_cn)
  minor.cn <- unlist(mut.table$minor_cn)
  depth.t  <- unlist(mut.table$ref_counts) + unlist(mut.table$var_counts)
  max.cn   <- max(abs.cn)
  VAF      <- unlist(mut.table$var_counts)/(unlist(mut.table$var_counts)+unlist(mut.table$ref_counts))
  
  # Use sequenza to estimate theoretical VAFs for each type of copy number
  types    <- types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 2)
  types.xy <- types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 1)
  types    <- rbind(types, types.xy)
  types    <- types[types$Mt >= 1, ]
  types$F  <- 0
  for (i in 1:nrow(types)) {
    types$F[i] <- theoretical.mufreq(cellularity = cellularity,
                                     CNn = types$CNn[i], CNt = types$CNt[i],
                                     Mt = types$Mt[i])
  }
  
  # Let's create some functions that can estimate whether early or late
  
  get.Mt <- function(F, depth.t, types, CNt, CNn, Mt){
    types <- types[types$CNn == CNn, ]
    l <- sequenza:::mufreq.dpois(mufreq = F, types$F[types$CNt== CNt&types$Mt<=Mt],
                                 depth.t = depth.t)
    l <- l/sum(l)
    L <- data.frame(l = l, Mt = types$Mt[types$CNt== CNt&types$Mt<=Mt])
  }
  
  get.conf <- function(F, depth.t){
    conf.int   <- cbind(prop.test(round(F*depth.t,0),depth.t)$conf[1]
                        ,prop.test(round(F*depth.t,0),depth.t)$conf[2])
    return(conf.int)
  }
  
  bootstrap.cf <- function(Vaf, cellularity, CNn, CNt, depth.t)
  {
    #print(i)
    if(Vaf==1)
    {
      conf.int   <- cbind(prop.test(round(Vaf*depth.t,0),depth.t)$conf[1]
                          ,prop.test(round(Vaf*depth.t,0),depth.t)$conf[2])
      
      lower      <- get.mut.mult(Vaf=conf.int[1],cellularity=cellularity,CNt=CNt,CNn=CNn)
      higher     <- get.mut.mult(Vaf=conf.int[2],cellularity=cellularity,CNt=CNt,CNn=CNn)
      conf.int   <- cbind(lower,higher)
      return(conf.int)
      
    }
    
    x          <- c(rep(1,round(Vaf*depth.t,0)),rep(0,(depth.t-round(Vaf*depth.t,0))))
    theta      <- function(x,i)
    {
      data      <- x[i]
      est       <- sum(data)/length(data)
      mut.multi <- (est *1/cellularity)*((cellularity*CNt)+CNn*(1-cellularity))
      return(mut.multi)
    }
    
    bt.res      <- boot(x,theta,R=1000)
    bt.ci       <- boot.ci(bt.res,type='norm')
    out         <- c(bt.ci$normal[2],bt.ci$normal[3])
    
    return(out)
    
  }
  
  
  get.mut.mult <- function(CNt,Vaf,cellularity,CNn)
  {
    
    return((Vaf *1/cellularity)*((cellularity*CNt)+CNn*(1-cellularity)))
    
  }
  
  
  get.cancer.cell.fraction <- function(Max.Likelihood,mut.mult)
  {
    predicted.Mtn   <- Max.Likelihood[,'Mt']
    ccf             <- mut.mult/predicted.Mtn
    return(ccf)
  }
  
  get.all.mut.info <- function(i)
  {
    #print(i)
    # First estimate the VAF confidence intervals
    obs.VAF         <- VAF[i]
    mut.conf.0.05   <- get.conf(F=VAF[i],depth.t=depth.t[i])[1]
    mut.conf.0.95   <- get.conf(F=VAF[i],depth.t=depth.t[i])[2]
    
    if(abs.cn[i]==0)
    {
      output          <- cbind(obs.VAF
                               ,mut.conf.0.05
                               ,mut.conf.0.95
                               ,mut.multi=NA
                               ,mut.multi.0.05=NA
                               ,mut.multi.bstr.0.05=NA
                               ,mut.multi.0.95=NA
                               ,mut.multi.bstr.0.95=NA
                               ,Exp.Cpn =NA
                               ,Exp.Cpn.Likelihood=NA
                               ,ccf=NA
                               ,ccf.0.05=NA
                               ,ccf.btstr.0.05=NA
                               ,ccf.0.95=NA
                               ,ccf.btstr.0.95=NA
                               ,timing='Not.Poss')
      return(output)
    }
    
    # Next estimate the likelihood relating to which copy number the mutation has
    L <- get.Mt(F = VAF[i],
                depth.t = depth.t[i], CNt = abs.cn[i],
                types = types, CNn = unlist(mut.table$normal_cn[i])
                ,Mt=major.cn[i])
    
    # Next determine the mut multiplicity
    mut.multi            <- get.mut.mult(CNt=abs.cn[i],Vaf=VAF[i],cellularity=cellularity,CNn=unlist(mut.table$normal_cn[i]))
    mut.multi.0.05       <- get.mut.mult(CNt=abs.cn[i],Vaf=mut.conf.0.05,cellularity=cellularity,CNn=unlist(mut.table$normal_cn[i]))
    mut.multi.0.95       <- get.mut.mult(CNt=abs.cn[i],Vaf=mut.conf.0.95,cellularity=cellularity,CNn=unlist(mut.table$normal_cn[i]))
    mut.multi.bstr       <- bootstrap.cf(Vaf=VAF[i],cellularity=cellularity,CNn=unlist(mut.table$normal_cn[i]),CNt=abs.cn[i],depth.t=depth.t[i])
    mut.multi.bstr.0.05  <- mut.multi.bstr[1]
    mut.multi.bstr.0.95  <- mut.multi.bstr[2]
    
    
    if(is.na(L$l)[1])
    {
      output          <- cbind(obs.VAF
                               ,mut.conf.0.05
                               ,mut.conf.0.95
                               ,mut.multi
                               ,mut.multi.0.05
                               ,mut.multi.bstr.0.05
                               ,mut.multi.0.95
                               ,mut.multi.bstr.0.95
                               ,Exp.Cpn =NA
                               ,Exp.Cpn.Likelihood=NA
                               ,ccf=NA
                               ,ccf.0.05=NA
                               ,ccf.btstr.0.05=NA
                               ,ccf.0.95=NA
                               ,ccf.btstr.0.95=NA
                               ,timing='Not.Poss')
      return(output)
    }
    
    # Now determine which likelihood should be used
    Max.Likelihood   <- L[which.max(L$l),]
    
    
    
    
    
    # Next determine the cancer cell fraction
    ccf             <- get.cancer.cell.fraction(Max.Likelihood,mut.multi)
    ccf.0.05        <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.0.05)
    ccf.btstr.0.05  <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.bstr.0.05)
    
    ccf.0.95        <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.0.95)
    ccf.btstr.0.95  <- get.cancer.cell.fraction(Max.Likelihood,mut.multi.bstr.0.95)
    
    # Next determine the late cancer cell fraction
    # Make sure you also output the theoretical Copy (i.e. what it's closest to using maximum likelihood)
    expected.copy     <- Max.Likelihood[2]
    
    # Finally also make a suggestion about whether the mutation is early late or not possible to tell
    
    timing      <- 'Not.Poss'
    
    if(Max.Likelihood$Mt>1)
    {
      timing    <- 'early'
    }
    
    if(Max.Likelihood$Mt<=1&minor.cn[i]!=1)
    {
      timing    <- 'late'
    }
    
    if(major.cn[i]<=1)
    {
      timing    <- 'Not.Poss'
    }
    
    
    
    
    # Let's put this all together and output it
    output          <- cbind(obs.VAF
                             ,mut.conf.0.05
                             ,mut.conf.0.95
                             ,mut.multi
                             ,mut.multi.0.05
                             ,mut.multi.bstr.0.05
                             ,mut.multi.0.95
                             ,mut.multi.bstr.0.95
                             ,Exp.Cpn =Max.Likelihood$Mt
                             ,Exp.Cpn.Likelihood=Max.Likelihood$l
                             ,ccf
                             ,ccf.0.05
                             ,ccf.btstr.0.05
                             ,ccf.0.95
                             ,ccf.btstr.0.95
                             ,timing)
    
    return(output)
    
  }
  
  output <- t(sapply(1:nrow(mut.table),get.all.mut.info))
  
  colnames(output) <- c('obs.VAF'
                        ,'mut.conf.0.05'
                        ,'mut.conf.0.95'
                        ,'mut.multi'
                        ,'mut.multi.0.05'
                        ,'mut.multi.bstr.0.05'
                        ,'mut.multi.0.95'
                        ,'mut.multi.bstr.0.95'
                        ,'Exp.Cpn'
                        ,'Exp.Cpn.Likelihood'
                        ,'ccf'
                        ,'ccf.0.05'
                        ,'ccf.btstr.0.05'
                        ,'ccf.0.95'
                        ,'ccf.btstr.0.95'
                        ,'timing')
  
  
  out <- cbind(mut.table,output)
  return(out)
  
}


#### Plotting ######
plot.EarlyOrLate    <- function( seg.mat.patient
                                 , TCGA.earlyLate
                                 , TCGA.purity
                                 , TCGA.barcode
                                 , max.cpn = 5
                                 , min.probes = 10
                                 , sub.clonal = 1
)
{
  
  
  seg.mat.patient$Chromosome <- as.numeric(as.character(seg.mat.patient$Chromosome))
  seg.mat.patient$StartPosition <- as.numeric(as.character(seg.mat.patient$StartPosition))
  seg.mat.patient$EndPosition <- as.numeric(as.character(seg.mat.patient$EndPosition))
  seg.mat.patient$nr.probes <- as.numeric(as.character(seg.mat.patient$nr.probes))
  
  chrom.length.copy   <- fun.chrom.length(seg.mat.patient[,2:4])
  chrom.segs          <- fun.add.chrom(seg.mat.patient[,2:4],chrom.length.copy)
  seg.mat.plot        <- seg.mat.patient
  seg.mat.plot[,2:4]  <- chrom.segs
  major               <- seg.mat.plot$Copy.Number - seg.mat.plot$min.allele
  seg.mat.plot        <- cbind(seg.mat.plot,major)
  seg.mat.plot        <- seg.mat.plot[seg.mat.plot$nr.probes>=min.probes,]
  
  
  # Order correctly
  TCGA.plot          <- TCGA.earlyLate
  Chromosome         <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,2])
  Start_pos          <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,3])
  TCGA.plot          <- cbind(TCGA.plot,Chromosome,Start_pos)
  TCGA.plot          <- data.frame(apply(TCGA.plot,2,unlist),stringsAsFactors=FALSE)
  
  min.x               <- 0
  min.y               <- -0.25
  max.y               <- max.cpn +1
  max.x               <- as.numeric(max(seg.mat.plot$EndPosition))
  
  # Make sure any copy numbers that are too high are still included
  seg.mat.plot$major       <- ifelse(seg.mat.plot$major>max.cpn,max.cpn,seg.mat.plot$major)
  seg.mat.plot$min.allele  <- ifelse(seg.mat.plot$min.allele>max.cpn,max.cpn,seg.mat.plot$min.allele)
  TCGA.plot$mut.multi      <- ifelse(TCGA.plot$mut.multi>max.cpn,max.cpn,TCGA.plot$mut.multi)
  
  #layout(rbind(1,2))
  #par(mar=c(5,5,5,5))
  plot(1
       ,xlim=c(min.x,max.x)
       ,ylim=c(min.y,max.y)
       ,type='n'
       ,xaxt='n'
       ,yaxs='i'
       ,yaxt='n'
       ,xaxs='i'
       ,ylab=""
       ,xlab=""
       ,lwd=2
       ,bty="n")
  
  box(lwd=0.5) 
  
  axis(side=2
       ,at=seq(0,max.cpn,by=1)
       ,labels = c(seq(0,max.cpn-1,by=1),paste('>',max.cpn,sep=""))
       ,las=2
       ,cex.axis=0.7
       ,lwd=0.5)
  
  mtext(text=TCGA.barcode
        ,side=3
        ,line=2
  )
  
  
  
  # Now plot each segment
  fun.plot.segment <- function(x,seg.mat.plot)
  {
    #print(x)
    #start with major allele
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- seg.mat.plot[x,'major']+0.05
    y1 <- y0
  
    
    segments(x0,y0,x1,y1,col='black',lwd=2.4)
    
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- as.numeric(seg.mat.plot[x,'min.allele'])-0.05
    y1 <- y0
    

    segments(x0,y0,x1,y1,col='#009E73',lwd=2.4)
    
  }
  
  sapply(1:nrow(seg.mat.plot),fun.plot.segment,seg.mat.plot)
  # Draw lines to separate the different chromosomes
  for (i in sort(unique(as.numeric(seg.mat.plot$Chromosome))))
  {
    abline(v=max(as.numeric(seg.mat.plot[seg.mat.plot$Chromosome==i,'EndPosition']))
           ,lwd=0.5)
    #,lty='dashed')
  }
  
  # We need to update the 
  TCGA.plot <- TCGA.plot[order(as.numeric(TCGA.plot$Chromosome),as.numeric(TCGA.plot$Start_pos)),]
  
  TCGA.plot$Start_pos <- as.numeric(fun.add.chrom(cbind(TCGA.plot$Chromosome,TCGA.plot$Start_pos,TCGA.plot$Start_pos),chrom.length.copy)[,2])
  
  # Let's add the mutations
  # Start with early mutations
  early.muts      <- TCGA.plot[TCGA.plot$timing%in%c('early'),]
  early.clonal    <- early.muts[early.muts$absoltue.ccf.0.95>=sub.clonal,]
  early.subclonal <- early.muts[early.muts$absoltue.ccf.0.95<sub.clonal,]
  
  late.muts  <- TCGA.plot[TCGA.plot$timing%in%c('late'),]
  late.clonal <- late.muts[late.muts$absoltue.ccf.0.95>=sub.clonal,]
  late.subclonal <- late.muts[late.muts$absoltue.ccf.0.95<sub.clonal,]
  
  nontimed.muts       <- TCGA.plot[!TCGA.plot$mutation_id%in%c(early.muts$mutation_id,late.muts$mutation_id),]
  nontimed.clonal     <- nontimed.muts[nontimed.muts$absoltue.ccf.0.95>=sub.clonal,]
  nontimed.subclonal  <- nontimed.muts[nontimed.muts$absoltue.ccf.0.95<sub.clonal,]
  
  
  # let's determine early or late or not possible
  plot.earlyORlateCol <- function(x,timed.muts,nontimed=FALSE)
  {
    #print(x)
    mut <- timed.muts[x,,drop=FALSE]
    
    #determine cell multiplicity
    mut.multiplicity <- mut$mut.multi
    #mut.multiplicity <- mut.multiplicity*(as.numeric(mut$minor_cn)+as.numeric(mut$major_cn))
    
    if(nontimed)
    {
      if(as.numeric(mut$mut.multi.bstr.0.95)<sub.clonal)
      {
        points(mut$Start_pos
               ,mut.multiplicity
               ,cex=0.7
               ,pch=17
               ,col='#99999965') # cannot be determined
        
      } 
      
      if(as.numeric(mut$mut.multi.bstr.0.95)>=sub.clonal)
      {
        points(mut$Start_pos
               ,mut.multiplicity
               ,cex=0.8
               ,pch=16
               ,col='#99999995') # cannot be determined
        
      } 
      
    }
    
    if(!nontimed)
    {
      # is it an early event?
      if(mut$timing=='early')
      {
        # is it subclonal
        if(as.numeric(mut$mut.multi.bstr.0.95)<sub.clonal)
        {
          points(mut$Start_pos
                 ,mut.multiplicity
                 ,cex=0.7
                 ,pch=17
                 ,col="#0072B245") # early event
        }
        
        # is it clonal
        if(as.numeric(mut$mut.multi.bstr.0.95)>=sub.clonal)
        {
          points(mut$Start_pos
                 ,mut.multiplicity
                 ,cex=0.9
                 ,pch=16
                 ,col="#0072B295") # early event
        }
        
        
      }
      
      if(mut$timing=='late')
      {
        if(as.numeric(mut$mut.multi.bstr.0.95)<sub.clonal)
        {
          points(mut$Start_pos
                 ,mut.multiplicity
                 ,cex=0.7
                 ,pch=17
                 ,col="#D55E0045") # late event
        }
        
        if(as.numeric(mut$mut.multi.bstr.0.95)>=sub.clonal)
        {
          points(mut$Start_pos
                 ,mut.multiplicity
                 ,cex=0.9
                 ,pch=16
                 ,col="#D55E0095") # late event
        }
        
        
        
      }
      
    }
    
  }
  
  
  if(nrow(early.muts)>=1)
  {
    sapply(1:nrow(early.muts),plot.earlyORlateCol,early.muts)  
  }
  if(nrow(late.muts)>=1)
  {
    sapply(1:nrow(late.muts),plot.earlyORlateCol,late.muts)  
  }
  if(nrow(nontimed.muts)>=1)
  {
    sapply(1:nrow(nontimed.muts),plot.earlyORlateCol,nontimed.muts,nontimed=TRUE)
  }
  
  
  
  mtext(side=3
        ,at=fun.chrom.mean(seg.mat.plot[,2:4])
        ,text=sort(unique(seg.mat.plot[,2]))
        ,cex=seq(0.6,0.4,length.out=length(unique(seg.mat.plot[,2])))
        ,line=-1
        #,lwd=0.5
  )
  
  
}




#########################################################################################


fun.chrom.length <- function(seg){
  # This function returns the length of chromosome 1,...,22
  # seg: first three columns of the minimum conistent region matrix
  
  chrom.length <- c()
  for(i in sort(as.numeric(unique(seg[,1])))){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.length <- c(chrom.length, sub.chrom[nrow(sub.chrom), 3])
    
  }
  return(chrom.length)
  
}

#############################################################################################    

fun.add.chrom <- function(seg, chrom.length){
  # This function adds the length of chromosome 1,...,(i - 1) to chromsomes 2,...,22
  # seg: first three columns of the minimum conistent region matrix
  # chrom.length: output of fun.chrom.length(...)
  if(1 %in% seg[, 1])
    seg.tmp <-  subset(seg, seg[, 1] == 1)
  else
    seg.tmp <- matrix(NA, nr = 0, nc = ncol(seg))  
  for(i in 2:max(max(as.numeric(seg[,1])),3)){
    if(!(i %in% seg[, 1]))
      next()
    sub.chrom <- subset(seg, seg[, 1] == i)
    sub.chrom[, 2] <- as.numeric(sub.chrom[, 2]) + sum(as.numeric(chrom.length[1:(i - 1)]))
    sub.chrom[, 3] <- as.numeric(sub.chrom[, 3]) + sum(as.numeric(chrom.length[1:(i - 1)]))
    seg.tmp <- rbind(seg.tmp, sub.chrom)
  }
  return(seg.tmp)
  
} 

####################################################################################################

fun.smooth.score <- function(score, seg){
  # This function applies smoothing (rollmean) to chromosome 1,...,22 separetely
  # score: score to be smoothed (e.g. correlation)
  # seg: Corresponding positions
  # window: window size for rollmean
  # The R package "zoo" is required
  score.tmp <- c()
  for(i in 1:22)
    
    score.tmp <- c(score.tmp, rollmean(score[seg[, 1] == i], k = round(length(score[seg[, 1] == i])/10), na.pad = T))
  
  
  return(score.tmp)
  
  
}  

#######################################################################################################

fun.chrom.mean <- function(seg){
  # This function finds the middle of the positions for each chromosome.
  # seg: Matrix with three columns, containing the positions.
  chrom.mean <- c()
  for(i in sort(as.numeric(unique(seg[,1])))){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.mean <- c(chrom.mean, (min(sub.chrom[, 2]) + max(sub.chrom[, 3]))/2)
    
  }
  return(chrom.mean)
  
}

#######################################################################################################


plot.TCGA.sample.ccfs <- function(TCGA.earlyLate
                                  ,clonal.cut.off=1
                                  ,spec.genes=NULL)
{
  
  
  # get the patient
  cancer.patient    <- data.frame(TCGA.earlyLate)
  
  #par(mfrow=c(2,1))
  layout(rbind(1,2))
  par(mar=c(0.5,5,5,5))
  #density(cancer.patient$ccf,bw='SJ')
  
  plot(density(as.numeric(cancer.patient$ccf),bw='SJ')
       ,xlim=c(0,1.4)
       ,xaxt='n'
       ,yaxt='n'
       #,xaxs='i'
       #,yaxs='i'
       ,xlab=""
       ,ylab='Density (a.u.)'
       ,main=cancer.patient$patient[1]
       ,lwd=3
  )
  par(mar=c(5,5,0.5,5))
  plot(as.numeric(cancer.patient$ccf)
       ,log(as.numeric(cancer.patient$var_counts)+as.numeric(cancer.patient$ref_counts))
       ,ylim=c(log(5),log(1000))
       ,xlim=c(0,1.6)
       ,xaxt='n'
       ,yaxt='n'
       #,xaxs='i'
       #,yaxs='i'
       ,xlab="Adjusted VAF"
       ,ylab="Tumour Coverage"
       ,main=""
       ,lwd=3
       ,cex=1.5
       ,pch=21
       ,col=ifelse(as.numeric(cancer.patient$absolute.ccf.0.95)>=1,'#de2d2699','#0072B299')
       ,bg=ifelse(cancer.patient$absolute.ccf.0.95>=1,'#de2d2699','#0072B299'))
  
  
  #next plot the genes of interest
  if(!is.null(spec.genes))
  {
    for (gene in spec.genes)
    {
      if(!gene%in%cancer.patient$Hugo_Symbol)
      {
        next;
      }
      
      gene.mut <- cancer.patient[cancer.patient$Hugo_Symbol%in%gene,,drop=FALSE]
      gene.mut <- gene.mut[gene.mut$Variant_Classification%in%c('Missense_Mutation'
                                                                ,'Nonsense_Mutation'
                                                                ,'Splice_Site'
                                                                ,'Frame_Shift_Del'
                                                                ,'Translation_Start_Site'
                                                                ,'Nonstop_Mutation'),,drop=FALSE]
      if(nrow(gene.mut)==0)
      {
        next;
      }
      
      points(gene.mut$ccf,log(gene.mut$var_counts+gene.mut$ref_counts)
             ,cex=1.8
             ,lwd=3)
      text(gene.mut$ccf,log(gene.mut$var_counts+gene.mut$ref_counts)
           ,labels=gene
           ,pos=4
           ,cex=1
           ,lwd=3
           ,font=2) 
      
    }
  }
  
  axis(side=1,at=seq(0,1.6,by=0.2))
  cov.seq <- c(5,10,20,50,100,200,500,1000)
  axis(side=2,at=log(cov.seq),labels=cov.seq,las=2)
  
  
}

