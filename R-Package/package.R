require(plyr)
require(dplyr)
require(foreach)
require(Rsubread)
get.retro.enrich.p.value.ip.only <- function(ip.bam.file,rt.annotation.gtf,GTF.featureType="exon", GTF.attrType="gene_id",chr.size.df,number.of.sampled.random.regions=100){
    rs.ip <- featureCounts(files    = ip.bam.file,
                           annot.ext= rt.annotation.gtf, 
                           isGTFAnnotationFile=TRUE,
                           GTF.featureType=GTF.featureType,
                           GTF.attrType=GTF.attrType,  
                           useMetaFeatures=TRUE,
                           allowMultiOverlap=TRUE,
                           primaryOnly = TRUE,
                           isPairedEnd = TRUE,
                           requireBothEndsMapped=TRUE,
                           checkFragLength=FALSE,
                           countChimericFragments=FALSE,
                           nthreads=8,
                           ignoreDup=TRUE
                           )
    rt.ip.read.count        <- c(rs.ip$counts)   
    names(rt.ip.read.count) <- rownames(rs.ip$counts)
    
    annotation            <- rs.ip$annotation 
    shuffle.region.rt.ip.read.count.matrix <- foreach(i = 1:number.of.sampled.random.regions,.combine='cbind') %do% {
            gr                 <- GRanges(seqnames = annotation$Chr,ranges= IRanges(start = annotation$Start,end=annotation$End))
            tmp                <- as.data.frame(randomizeRegions(gr, genome=chr.size.df))
            shuffle.annotation <- data.frame(GeneID= annotation$GeneID, Start = tmp$start,End=tmp$end,Strand=annotation$Strand)
            shuffle.rs <- featureCounts(files    = ip.bam.file,
                                        annot.ext= shuffle.annotation, 
                                        isGTFAnnotationFile=FALSE,
                                        useMetaFeatures=TRUE,
                                        allowMultiOverlap=TRUE,
                                        primaryOnly = TRUE,
                                        isPairedEnd = TRUE,
                                        requireBothEndsMapped=TRUE,
                                        checkFragLength=FALSE,
                                        countChimericFragments=FALSE,
                                        nthreads=8,
                                        ignoreDup=TRUE
                                        )
           shuffle.rs$counts
        }
    lambda.hat <- apply(shuffle.region.rt.ip.read.count.matrix,1,median)        
    p.value    <- ppois(rt.ip.read.count-1,lambda.hat,lower.tail = FALSE,log.p = FALSE)
    df         <- data.frame(rt.name = rownames(rs.ip$counts),ip.read.count = rt.ip.read.count,lambda.hat = lambda.hat,p.value=p.value)    
    df
}




.get.zscore <- function(m,a,n1,n2){ #n1 for IP, n2 for input
    condition.mean <- 0
    tmp            <- 2^(a-log2(sqrt(n1)) - log2(sqrt(n2)))
    condition.var  <- (4*(1-tmp)* (1/log(2)) * (1/log(2))) / ((n1+n2)*tmp)
    z.score        <- (m-condition.mean) / sqrt(condition.var)
    z.score
}


get.retro.enrich.p.value.ip.vs.input <- function(ip.bam.file,input.bam.file,rt.annotation.gtf,GTF.featureType="exon", GTF.attrType="gene_id",p.value.assignment.method='poisson'){
    rs.ip <- featureCounts(files    = ip.bam.file,
                           annot.ext= rt.annotation.gtf, 
                           isGTFAnnotationFile=TRUE,
                           GTF.featureType=GTF.featureType,
                           GTF.attrType=GTF.attrType,  
                           useMetaFeatures=TRUE,
                           allowMultiOverlap=TRUE,
                           primaryOnly = TRUE,
                           isPairedEnd = TRUE,
                           requireBothEndsMapped=TRUE,
                           checkFragLength=FALSE,
                           countChimericFragments=FALSE,
                           nthreads=8,
                           ignoreDup=TRUE
                           )
    rt.ip.read.count        <- c(rs.ip$counts)   
    names(rt.ip.read.count) <- rownames(rs.ip$counts)
    
    rs.input <- featureCounts(files    = input.bam.file,
                              annot.ext= rt.annotation.gtf, 
                              isGTFAnnotationFile=TRUE,
                              GTF.featureType=GTF.featureType,
                              GTF.attrType=GTF.attrType,  
                              useMetaFeatures=TRUE,
                              allowMultiOverlap=TRUE,
                              primaryOnly = TRUE,
                              isPairedEnd = TRUE,
                              requireBothEndsMapped=TRUE,
                              checkFragLength=FALSE,
                              countChimericFragments=FALSE,
                              nthreads=8,
                              ignoreDup=TRUE
                             )   
    rt.input.read.count        <- c(rs.input$counts)   
    names(rt.input.read.count) <- rownames(rs.input$counts)  
    rt.input.read.count        <- rt.input.read.count[names(rt.ip.read.count)]           
    
    ip.lib.size                <- sprintf("samtools view -F 1804  %s | wc -l",ip.bam.file)    %>% system(intern=TRUE,wait=TRUE) %>% as.integer
    input.lib.size             <- sprintf("samtools view -F 1804  %s | wc -l",input.bam.file) %>% system(intern=TRUE,wait=TRUE) %>% as.integer
    
    retrun(list(read.count.df=data.frame(rt.ip.read.count=rt.ip.read.count,rt.input.read.count=rt.input.read.count),ip.lib.size=ip.lib.size,input.lib.size=input.lib.size))
    
    
    if(p.value.assignment.method =='poisson'){
        alpha   <- ip.lib.size /(ip.lib.size  + input.lib.size)
        p.value <- pbinom(rt.ip.read.count-1,rt.ip.read.count + rt.input.read.count,alpha,lower.tail = FALSE, log.p = FALSE)
        df      <- data.frame(rt.name = rownames(rs.ip$counts),ip.read.count = rt.ip.read.count, input.read.count = rt.input.read.count,p.value=p.value,alpha=alpha)    

    }else{
        n1         <-  ip.lib.size
        n2         <-  input.lib.size
        m          <-  (log2(rt.ip.read.count + 1) - log2(rt.input.read.count + 1)) - (log2(n1+1) - log2(n2+1))
        a          <-  (log2(rt.ip.read.count + 1) + log2(rt.input.read.count + 1))/2
        df         <-  data.frame(M=m,A=a)
        loess.fit  <-  loess(formula = M ~ A,data = df,span = 0.3)
        a          <-  loess.fit$x %>% c
        m          <-  loess.fit$residuals
        z.score    <-  .get.zscore(m,a,n1,n2)
        p.value    <-  pnorm(z.score,lower.tail = FALSE)
        names(p.value) <- names(m)
        p.value        <- p.value[rownames(rs.ip$counts)]
        df             <- data.frame(rt.name = rownames(rs.ip$counts),ip.read.count = rt.ip.read.count, input.read.count = rt.input.read.count,p.value=p.value,m=m[rownames(rs.ip$counts)],a=a[rownames(rs.ip$counts)])    

    }
    df
}



rt.annotation.gtf <- 'annotation/mm10.repbase.gtf'
GTF.featureType   <- "exon"
GTF.attrType      <- "repeat"

ip.bam.file       <- 'alignment/SRR5297327.all.rmdup.bam' 

input.bam.file    <- 'alignment/SRR5297331.all.rmdup.bam'
L1                <- get.retro.enrich.p.value.ip.vs.input(ip.bam.file,input.bam.file,rt.annotation.gtf,GTF.featureType,GTF.attrType,p.value.assignment.method='poisson')

input.bam.file    <- 'alignment/SRR5297328.all.rmdup.bam'
L2                <- get.retro.enrich.p.value.ip.vs.input(ip.bam.file,input.bam.file,rt.annotation.gtf,GTF.featureType,GTF.attrType,p.value.assignment.method='poisson')

#df.poisson <- get.retro.enrich.p.value.ip.vs.input(ip.bam.file,input.bam.file,rt.annotation.gtf,GTF.featureType,GTF.attrType,p.value.assignment.method='poisson')
#df.DEGSeq  <- get.retro.enrich.p.value.ip.vs.input(ip.bam.file,input.bam.file,rt.annotation.gtf,GTF.featureType,GTF.attrType,p.value.assignment.method='DEGSeq')
save.image(file='test.RData')

quit()








################################################################################ Trash #######################################################

get.retro.enrich.p.value <- function(ip.bam.file,input.bam.file,rt.annotation.gtf,GTF.featureType="exon", GTF.attrType="gene_id",chr.size.df)
{
    mode <- ifelse(ip.bam.file == input.bam.file,'ip.only','ip.and.input')
    
    rs.ip <- featureCounts(files    = ip.bam.file,
                           annot.ext= rt.annotation.gtf, 
                           isGTFAnnotationFile=TRUE,
                           GTF.featureType=GTF.featureType,
                           GTF.attrType=GTF.attrType,  
                           useMetaFeatures=TRUE,
                           allowMultiOverlap=TRUE,
                           primaryOnly = TRUE,
                           isPairedEnd = TRUE,
                           requireBothEndsMapped=TRUE,
                           checkFragLength=FALSE,
                           countChimericFragments=FALSE,
                           nthreads=8,
                           ignoreDup=TRUE
                           )
    rt.ip.read.count        <- c(rs.ip$counts)   
    names(rt.ip.read.count) <- rownames(rs.ip$counts)
                       
    if(mode == 'ip.only'){
        
    }else{
        rs.input <- featureCounts(files    = input.bam.file,
                                  annot.ext= rt.annotation.gtf, 
                                  isGTFAnnotationFile=TRUE,
                                  GTF.featureType=GTF.featureType,
                                  GTF.attrType=GTF.attrType,  
                                  useMetaFeatures=TRUE,
                                  allowMultiOverlap=TRUE,
                                  primaryOnly = TRUE,
                                  isPairedEnd = TRUE,
                                  requireBothEndsMapped=TRUE,
                                  checkFragLength=FALSE,
                                  countChimericFragments=FALSE,
                                  nthreads=8,
                                  ignoreDup=TRUE
                                )   
        rt.input.read.count        <- c(rs.input$counts)   
        names(rt.input.read.count) <- rownames(rs.input$counts)             
        ip.lib.size                <- sprintf("samtools view -F 1804  %s | wc -l",ip.bam.file)    %>% system %>% as.integer
        input.lib.size             <- sprintf("samtools view -F 1804  %s | wc -l",input.bam.file) %>% system %>% as.integer
        alpha   <- ip.lib.size /(ip.lib.size  + input.lib.size)
        p.value <- pbinom(rt.ip.read.count-1,rt.ip.read.count + rt.input.read.count,alpha,lower.tail = FALSE, log.p = FALSE)
        df      <- data.frame(rt.name = rownames(rs.ip$counts),ip.read.count = rt.ip.read.count, input.read.count = rt.input.read.count,alhpa = alpha,p.value=p.value)    
        df
    }
}

