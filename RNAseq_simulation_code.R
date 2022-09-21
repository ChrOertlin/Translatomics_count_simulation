#### SIMULATE RNASEQ DATA USING MEAN AND DISPERSION PARAMETER ESTIMATES FROM REAL DATA ####
#### 
####  Before using these functions filter your raw data for counts < 44
####  Output will be a dataset with the same amount of genes as in the parent dataset
####  users can specify the amount of replictes in the draw reads function   
####  default number of replicates is set to 4
####
####  functions to call:
####  tcNorm()
####  estimParameters()
####  drawReads_seqDepth()
####
####

tcNorm <- function(data){
    ### get the total number of counts for each lane
    totalPerLane <- apply(data,2,sum)
    ### get the total number of counts across all lanes
    totalOverall <- mean(totalPerLane)
    ### replace each count by (count/ lane total) * overall total
    for(j in 1:dim(data)[[2]]){
        for(i in 1:dim(data)[[1]]){
            data[i,j]<- (data[i,j]/totalPerLane[j])*totalOverall
        }
    }
    return(data)
}




### Pass the library size normalised counts to this function
estimParameters <- function(data,p,iterlim=10000){
    #### log likelihood function
    loglike <- function(dispersion){
        - (sum(lgamma(counts + 1/dispersion)) - length(counts) * lgamma(1/dispersion) - sum(lgamma(counts + 1)) + sum(counts * log(mean * dispersion/(1 + mean * dispersion))) - length(counts)/dispersion * log(1 + mean * dispersion))
    }
    means <- vector("numeric",dim(data)[[1]])
    dsp <- vector("numeric",dim(data)[[1]])
    ####
    #### Estimate mean according to 1/amount of replicates*(sum of counts across replicates)
    #### Estimate per gene disperion by maximasing log- likelihood function, with y set to corrected counts
    ####
    for(i in 1:dim(data)[[1]]){
        mean <- 1/dim(data)[[2]]*sum(data[i,])
        means[i] <- mean
        counts <- data[i,]
        dsp[i] <- nlm(f=loglike,p=p,gradtol=1e-6,iterlim=iterlim)$estimate
    }
    return(list("means"=means, "dispersion"=dsp))
}

### function to draw read counts from a nb distribution with given mean and dispersion
### to simulate DE the data is split into quantiles (6,12.5,25,50,75,87.5,94)
### Genes that should be DE will then be chosen
### control and treatment will then be drawn from a nbin with mean and dispersion
### belonging to different quantiles chosen at random 
###
### indicate proportions of categories
### cat1 : Delta Poly and no delta Cyto
### cat2 : delta poly > delta cyto
### cat3 : no diff (delta poly == delta cyto)

###
### Parameters 
### polyParam - list containing means and dispersion of polysome data
### polyData - Polysome data 
### cytoParam - list containing means and dispersion of cytosolic data
### cytoData - cytosolic data
### percentDE - percentage of genes that shoudl be differential expressed
### percentUpreg - percentage of DE genes that should be upregulated
### percentDownreg - percentage of DE genes that should be downregulated
### nReps - number of replicates that should be simulated
### cat1,cat2,cat3 - percentages of DE genes that should belong to a specified category (explained above)
### noise - percentage of noise that should be added 
### polyParamSeqDetph - list containng means and dispersion of polysome data with a different sequencing depth length has to match nReps*2
### cytoParamSeqDetph - list containng means and dispersion of cytosolic data with a different sequencing depth length has to match nReps*2
### polySeqDepth - vector specifiying what polysome-associated samples should be drawn from distriubtions with different seq depths (0 for not 1 for draw)
### cytoSeqDepth - vector specifiying what cytosolic samples should be drawn from distriubtions with different seq depths (0 for not 1 for draw)  
###
###
###
drawReads_SeqDepth <- function(polyParam=NULL,polyData=NULL,cytoData=NULL,cytoParam=NULL,percentDE=NULL,
                               percentUpreg=NULL,percentDownreg=NULL,
                               nReps=4,quantSeq=NULL,
                               cat1=NULL,cat2=NULL,cat3=NULL,cat4=NULL,noise=NULL,polyVar=NULL,cytoVar=NULL,transVar=FALSE,buffVar=FALSE,
                               polyParamSeqDepth=NULL,
                               cytoParamSeqDepth=NULL,
                               polySeqDepth = c(0,0,0,0,
                                                0,0,0,0),
                               cytoSeqDepth = c(0,0,0,0,
                                                0,0,0,0)
){
    
    polyData <- as.matrix(polyData)
    cytoData <- as.matrix(cytoData)
    selectGenes <- NULL
    ## Normal params
    polyMeans <- round(polyParam$means)
    polyDispersion <- polyParam$dispersion 
    cytoMeans <- round(cytoParam$means)
    cytoDispersion <- cytoParam$dispersion
    ## SeqDepthParams
    polyMeansSeqDepth <- round(polyParamSeqDepth$means)
    polyDispersionSeqDepth <- polyParamSeqDepth$dispersion
    cytoMeansSeqDepth <- round(cytoParamSeqDepth$means)
    cytoDispersionSeqDepth <- cytoParamSeqDepth$dispersion
    
    polyCounts <- matrix(0 ,nrow=dim(polyData)[[1]],ncol=(nReps*2))
    cytoCounts <- matrix(0,nrow=dim(polyData)[[1]],ncol=(nReps*2))
    rownames(polyCounts) <- paste("simGenes", c(1:dim(polyData)[[1]]),sep="_")
    colnames(polyCounts) <- c(paste("control",c(1:nReps),sep="_"),paste("treatment",c(1:nReps),sep="_")) 
    colnames(polyCounts) <- paste(colnames(polyCounts),"poly",sep="_")
    rownames(cytoCounts) <- paste("simGenes", c(1:dim(polyData)[[1]]),sep="_")
    colnames(cytoCounts) <- c(paste("control",c(1:nReps),sep="_"),paste("treatment",c(1:nReps),sep="_")) 
    colnames(cytoCounts) <- paste(colnames(cytoCounts),"cyto",sep="_")
    polyQuantile <- quantile(polyData,probs=c(.5,.95))
    if(!is.null(percentDE)){
        nDE <- round(dim(polyData)[[1]]*percentDE)
        if(is.null(percentUpreg) |is.null(percentDownreg)){
            print("please specify percentages of up and downregulated genes")
            stop()
        }
        upregulatedDE <-round(nDE*percentUpreg)
        downregulatedDE <- round(nDE*percentDownreg)
        regulation <- 1:nrow(polyCounts)
        up <- sample(regulation,upregulatedDE)
        subsetReg <- setdiff(regulation,up)
        down <- sample(subsetReg,downregulatedDE)
        
        regulation[up]<- "UP"
        regulation[down] <- "DOWN"
        for(i in 1:length(regulation)){
            if(regulation[i]!="DOWN" & regulation[i]!="UP"){
                regulation[i] <- "EE" 
            }
        }
        ### second layer of the simulation puts genes into different categories
        ### Translation, Transcription and no difference; based on delta polysome counts
        ### The groups are a proportion of DE genes and proportions are specified in the function
        ###
        ### cat1:translation : Delta Poly and no delta Cyto
        ### cat2: buffering : delta cyto no delta poly
        ### cat3: transcription : (delta poly ><== delta cyto)
        ### 
        #####
        ##### For up and downregulated genes indicate what their categorie 
        ##### of translational or non translation regulation will be 
        ##### by tagging the proportion of DE genes regardless of up or down regulation
        ##### with a categorie flag (cat1,cat2 ... )
        #####
        categories <- vector("numeric",3)
        categories <- c(round(cat1*(upregulatedDE + downregulatedDE)),
                        round(cat2*(upregulatedDE + downregulatedDE)),
                        round(cat3*(upregulatedDE + downregulatedDE)))
        
        if(sum(categories) > nDE){
            diff <- abs(sum(categories)-nDE)
            categories[3] <- categories[3]-diff
        }
        catName <- c("cat1","cat2","cat3")
        for(j in 1:length(categories)){
            if(!categories[j]== 0){
                pos <- which(regulation=="DOWN"| regulation == "UP")
                replace <- base::sample(pos,categories[j])
                for(i in 1:length(replace)){
                    regulation[replace[i]]<- paste(regulation[replace[i]],catName[j],sep="_")
                }  
            }
        }
        ####
        ####
        #### Start drawing reads for the simulated data
        #### polysome data is drawn first
        #### for DE genes based on delta polysome cytosolic reads will be created
        ####
        ####
        
        
        for(i in 1:dim(polyCounts)[[1]]){
            
            effsize <- base::sample(seq(1.5,3,0.2),1)
            
            for(s in 1:length(polySeqDepth)){
                if(polySeqDepth[s] == 0){
                    tmpPolyMeans <- polyMeans
                    tmpPolyDispersion <- polyDispersion
                    polyControl <- tmpPolyMeans[i]
                }
                if(polySeqDepth[s] == 1){
                    tmpPolyMeans <- polyMeansSeqDepth
                    tmpPolyDispersion <- polyDispersionSeqDepth
                    polyControl <- tmpPolyMeans[i]
                }
                if(cytoSeqDepth[s] == 0){
                    tmpCytoMeans <- cytoMeans
                    tmpCytoDispersion <- cytoDispersion
                    cytoControl <- tmpCytoMeans[i]
                }
                if(cytoSeqDepth[s] == 1){
                    tmpCytoMeans <- cytoMeansSeqDepth
                    tmpCytoDispersion <- cytoDispersionSeqDepth
                    cytoControl <- tmpCytoMeans[i]
                }
                ### effsize is used to itntroduce differential expression to the polysome Data
                
                
                #### draw number of EE reads that is counts from the same nbin dist
                #### for both conditions 
                if(regulation[i]=="EE"){
                    polyCounts[i,s] <- rnbinom(n= 1,mu=tmpPolyMeans[i],size=1/tmpPolyDispersion[i])
                    cytoCounts[i,s] <- rnbinom(n=1,mu=tmpCytoMeans[i],size=1/tmpCytoDispersion[i]) 
                    ##### Cytosolic counts should also be EE ...? #####    
                }
                if(grepl("cat1",regulation[i])){
                    cytoCounts[i,s]<- rnbinom(n=1,mu=tmpCytoMeans[i],size=1/tmpCytoDispersion[i])
                }
                if(grepl("cat2",regulation[i])){
                    polyCounts[i,s]<- rnbinom(n=1,mu=tmpPolyMeans[i],size=1/tmpPolyDispersion[i])
                }
                if(grepl("UP",regulation[i])){ 
                    ### base::sample for condition 1
                    if(s <= nReps){
                        polyDispControl <- which.min(abs(tmpPolyMeans-polyControl))
                        if(grepl("UP_cat1",regulation[i])){
                            polyCounts[i,s]<- rnbinom(n=1,mu=polyControl,size=1/tmpPolyDispersion[polyDispControl])
                        }
                        if(grepl("UP_cat2",regulation[i])){
                            cytoDispControl<- which.min(abs(tmpCytoMeans-cytoControl))
                            cytoCounts[i,s]<- rnbinom(n=1,mu=cytoControl,size=1/tmpCytoDispersion[cytoDispControl])
                        }
                        if(grepl("UP_cat3",regulation[i])){
                            
                            polyCounts[i,s]<- rnbinom(n=1,mu=polyControl,size=1/tmpPolyDispersion[polyDispControl])
                            cytoDispControl<- which.min(abs(tmpCytoMeans-cytoControl))
                            cytoCounts[i,s]<- rnbinom(n=1,mu=polyControl,size=1/tmpCytoDispersion[cytoDispControl])
                        }   
                    }
                    ### base::sample for condition 2
                    #μgs(1+μgsϕgs)
                    if(s > nReps){
                        if(grepl("UP_cat1",regulation[i])){
                            polyTreatment <- polyControl*effsize
                            polyDispTreatment <-which.min(abs(tmpPolyMeans-polyTreatment))
                            polyCounts[i,s]<- rnbinom(n=1,mu=polyTreatment,size=1/tmpPolyDispersion[polyDispTreatment])
                        }
                        if(grepl("UP_cat2",regulation[i])){
                            cytoTreatment <- cytoControl*effsize
                            cytoDispTreatment <- which.min(abs(tmpCytoMeans-cytoTreatment))
                            cytoCounts[i,s]<- rnbinom(n=1,mu=cytoTreatment,size=1/tmpCytoDispersion[cytoDispTreatment])
                        }
                        if(grepl("UP_cat3",regulation[i])){
                            polyTreatment <- polyControl*effsize
                            polyDispTreatment <-which.min(abs(tmpPolyMeans-polyTreatment))
                            cytoTreatment <- tmpCytoMeans[which.min(abs(tmpCytoMeans-polyTreatment))]
                            cytoDispTreatment <- which.min(abs(tmpCytoMeans-cytoTreatment))
                            polyCounts[i,s]<- rnbinom(n=1,mu=polyTreatment,size=1/tmpPolyDispersion[polyDispTreatment])
                            cytoCounts[i,s]<- rnbinom(n=1,mu=polyTreatment,size=1/tmpCytoDispersion[cytoDispTreatment])
                        } 
                    }
                }
                if(grepl("DOWN",regulation[i])){
                    if(s <= nReps){
                        polyDispControl <- which.min(abs(tmpPolyMeans-polyControl))
                        if(grepl("DOWN_cat1",regulation[i])){
                            polyCounts[i,s]<- rnbinom(n=1,mu=polyControl,size=1/tmpPolyDispersion[polyDispControl])
                        }
                        if(grepl("DOWN_cat2",regulation[i])){                
                            cytoDispControl<- which.min(abs(tmpCytoMeans-cytoControl))
                            cytoCounts[i,s]<- rnbinom(n=1,mu=cytoControl,size=1/tmpCytoDispersion[cytoDispControl])
                        }
                        if(grepl("DOWN_cat3",regulation[i])){
                            polyCounts[i,s]<- rnbinom(n=1,mu=polyControl,size=1/tmpPolyDispersion[polyDispControl])
                            cytoDispControl<- which.min(abs(tmpCytoMeans-cytoControl))
                            cytoCounts[i,s]<- rnbinom(n=1,mu=polyControl,size=1/tmpCytoDispersion[cytoDispControl])
                        }
                    }
                    ### base::sample for condition 2
                    if(s>nReps){
                        if(grepl("DOWN_cat1",regulation[i])){
                            polyTreatment <- polyControl/effsize
                            polyDispTreatment <- which.min(abs(tmpPolyMeans-polyTreatment))
                            polyCounts[i,s]<- rnbinom(n=1,mu=polyTreatment,size=1/tmpPolyDispersion[polyDispTreatment])
                        }
                        if(grepl("DOWN_cat2",regulation[i])){
                            cytoTreatment <- cytoControl/effsize
                            cytoDispTreatment <- which.min(abs(tmpCytoMeans-cytoTreatment))
                            cytoCounts[i,s]<- rnbinom(n=1,mu=cytoTreatment,size=1/tmpCytoDispersion[cytoDispTreatment])
                        }
                        if(grepl("DOWN_cat3",regulation[i])){
                            polyTreatment <- polyControl/effsize
                            polyDispTreatment <-which.min(abs(tmpPolyMeans-polyTreatment))
                            effsize <- effsize
                            cytoTreatment <- tmpCytoMeans[which.min(abs(tmpCytoMeans-polyTreatment))]
                            cytoDispTreatment <- which.min(abs(tmpCytoMeans-cytoTreatment))
                            polyCounts[i,s]<- rnbinom(n=1,mu=polyTreatment,size=1/tmpPolyDispersion[polyDispTreatment])
                            cytoCounts[i,s]<- rnbinom(n=1,mu=polyTreatment,size=1/tmpCytoDispersion[cytoDispTreatment])
                        }
                    } 
                }
            } 
        }
    }
    
    ### Simulate a dataset with no DE genes that is draw for all base::samples from the same distribution
    ### with the same mean and dispersion. 
    
    
    
    
    if(is.null(percentDE)& is.null(selectGenes)){
        
        for(i in 1:dim(polyData)[1]){
            
            for(s in 1:length(polySeqDepth)){
                if(polySeqDepth[s] == 0){
                    tmpPolyMeans <- polyMeans
                    tmpPolyDispersion <- polyDispersion
                }
                if(polySeqDepth[s] == 1){
                    tmpPolyMeans <- polyMeansSeqDepth
                    tmpPolyDispersion <- polyDispersionSeqDepth
                }
                if(cytoSeqDepth[s] == 0){
                    tmpCytoMeans <- cytoMeans
                    tmpCytoDispersion <- cytoDispersion
                }
                if(cytoSeqDepth[s] == 1){
                    tmpCytoMeans <- cytoMeansSeqDepth
                    tmpCytoDispersion <- cytoDispersionSeqDepth
                }
                polyCounts[i,s] <- rnbinom(n= 1,mu=tmpPolyMeans[i],size=1/tmpPolyDispersion[i])
                cytoCounts[i,s]<- rnbinom(n=1,mu=tmpCytoMeans[i],size=1/tmpCytoDispersion[i])
            }
            
        }
    }  
    
    #### Noise function adds or substracts % as Noise from readCount 
    noisinator <- function(counts,noise){
        for (j in 1:dim(counts)[2]){
            for(i in 1:dim(counts)[1]){
                noiseAmountPos <- counts[i,j] * noise
                noiseAmountNeg <- (counts[i,j] * noise) * -1
                counts[i,j] <- as.integer(counts[i,j] + sample(c(noiseAmountPos,noiseAmountNeg),1))
            } 
        } 
        return(counts)    
    }  
    ######
    if(!is.null(noise)){
        polyCounts <- noisinator(polyCounts,noise)
        cytoCounts <- noisinator(cytoCounts,noise)
    }
    ### IntVar function:
    ### Adds variation to each control and treatment group of a gene 
    ### by adding a percentage of the standard deviation of the gene condition
    ### To ensure a spread along the axis, counts first get sorted and the lowest count 
    ### will be assigned the smallest change in variation
    intVar <- function(counts,minVar,maxVar,steps=nReps){
        varRange <- seq(from=minVar,to=maxVar,by= (maxVar-minVar)/(steps-1))
        varRange <- varRange
        for(i in 1:dim(counts)[1]){
            controlCounts <- sort(counts[i,1:steps],decreasing = FALSE)
            treatmentCounts <- sort(counts[i,(steps+1):(steps*2)],decreasing = FALSE)
            controlSD <- sd(controlCounts)
            treatmentSD <- sd(treatmentCounts)
            for(j in 1:steps){
                controlCounts[j] <- controlCounts[j] + (controlSD * varRange[j])
                treatmentCounts[j] <- treatmentCounts[j] + (treatmentSD * varRange[j])
            }
            counts[i,1:steps] <- controlCounts
            counts[i,(steps+1):(steps*2)] <- treatmentCounts
        }
        
        return(list(as.integer(counts)))
    }
    
    if(!is.null(percentDE)|!is.null(selectGenes)){
        deGenes <- rownames(polyCounts)[which(grepl("UP",regulation) | grepl("DOWN",regulation))]
        cat1Genes <- rownames(polyCounts)[which(grepl("cat1",regulation))]
        cat2Genes <- rownames(polyCounts)[which(grepl("cat2",regulation))]
        cat3Genes <- rownames(polyCounts)[which(grepl("cat3",regulation))]
        
        ### add noise in specified count set only to the DE genes to induce variation
        ### at only cytosolic or polysome associated mRNA levels
        if(transVar==TRUE){
            deGenesVar <- cat1Genes
        }
        
        if(buffVar==TRUE){
            deGenesVar <- cat2Genes
        }
        
        if(!is.null(polyVar)){
            varCounts <- intVar(polyCounts[deGenesVar,],polyVar[1],polyVar[2])
            polyCounts[deGenesVar,] <- varCounts[[1]]
        }
        if(!is.null(cytoVar)){
            varCounts <- intVar(cytoCounts[deGenesVar,],cytoVar[1],cytoVar[2])
            cytoCounts[deGenesVar,] <- varCounts[[1]]  
        }
        
        #### combine Poly and Cyto set into one big dataset
        combinedCounts <- cbind(polyCounts,cytoCounts)
        regulationMatrix <- cbind(rownames(polyCounts),regulation)
        rownames(regulationMatrix) <- regulationMatrix[,1]
        
        return(list("polyCounts"= polyCounts,"cytoCounts"=cytoCounts ,"combined"=combinedCounts, "deGenes"= deGenes,"translation"=cat1Genes,"buffering"=cat2Genes,"mRNAAbundance"=cat3Genes,"regulation"=regulationMatrix))
    }
    if(is.null(percentDE) & is.null(selectGenes)){
        deGenes <- vector("character",1)
        #### combine Poly and Cyto set into one big dataset
        combinedCounts <- cbind(polyCounts,cytoCounts)
        return(list("polyCounts"= polyCounts,"cytoCounts"=cytoCounts, "combined"=combinedCounts, "deGenes"= deGenes))
    }
} 

mySampleFun <- function(x, names, sampSize=1e6){
    cat ("Starting sampling\n")
    data <- x
    ##The idea here is to sample a fixed number of reads to enable a more comparable analysis
    ##The input is a matrix so we will use an apply function
    tmpData <- rep(names, data)
    tmpSamp <- sample(tmpData, sampSize)
    tmpSum <- table(tmpSamp)
    tmpSum2 <- as.vector(tmpSum)
    names(tmpSum2) <- names(tmpSum)
    tmpOut <- c(rep(0, length(names)))
    names(tmpOut) <- names
    tmpOut[names(tmpSum2)] <- tmpSum2[names(tmpSum2)]
    return(tmpOut)
}

resample <- function(simCountsObj,samSize){
    mySampleFun <- function(x, names, sampSize=1e6){
        cat ("Starting sampling\n")
        data <- x
        ##The idea here is to sample a fixed number of reads to enable a more comparable analysis
        ##The input is a matrix so we will use an apply function
        tmpData <- rep(names, data)
        tmpSamp <- sample(tmpData, sampSize)
        tmpSum <- table(tmpSamp)
        tmpSum2 <- as.vector(tmpSum)
        names(tmpSum2) <- names(tmpSum)
        tmpOut <- c(rep(0, length(names)))
        names(tmpOut) <- names
        tmpOut[names(tmpSum2)] <- tmpSum2[names(tmpSum2)]
        return(tmpOut)
    }
    
    simCountsObj$combined <- apply(simCountsObj$combined,2,mySampleFun,
                                   names=rownames(simCountsObj$combined),
                                   sampSize=samSize)
    simCountsObj$polyCounts <- simCountsObj$combined[,grepl("poly",colnames(simCountsObj$combined))]
    simCountsObj$cytoCounts <- simCountsObj$combined[,grepl("cyto",colnames(simCountsObj$combined))]                               
    
    return(simCountsObj)
    
}




### Example code


#tcNormDat <- tcNorm(<Your RNAseq data>)

#dataP <- tcNormDat[,<YOUR Polysome-assocaited mRNA samples>]
#dataT <- tcNormDat[,<your total mRNA samples>]

## get parameter estimates for simulations --
#pDataParam <- estimParameters(dataP,10,25)
#cDataParam <- estimParameters(dataT,10,25)

##
## get parameter estimates for simulations with reduced sample size
#dataP2.5M <- apply(dataP,2, mySampleFun,names=rownames(dataP),sampSize = 2.5e6)
#dataP5M <-apply(dataP,2, mySampleFun,names=rownames(dataP),sampSize = 5e6)

## generate a simulated data set with 4 samples of a reduced sequencing depth (i.e. 2 poly and 2 total)

#simCounts <- drawReads_SeqDepth(polyParam=pDataParam5M,polyData=pData,cytoData=cData,
#                                                 cytoParam=cDataParam5M,percentDE=0.15,
#                                                 percentUpreg=.5,percentDownreg=.5,
#                                                 nReps=4,cat1=1/3,cat2=1/3,cat3=1/3,noise=NULL,
#                                                 polyParamSeqDepth = pDataParam2.5M,
#                                                cytoParamSeqDepth = cDataParam2.5M,
#                                               polySeqDepth = c(1,0,0,0,
#                                                               1,0,0,0)
#                                             cytoSeqDepth = c(1,0,0,0,
#                                                             1,0,0,0))

