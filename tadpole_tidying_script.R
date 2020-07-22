
library(FRESA.CAD)

clean_tadpole <- function(path){
  #load ds
    DatosANDI <- read.delim(path,na.strings=c("NA","","#N/A","#DIV/0!"),sep = ",")
    #ordenar por RID
    DatosANDI <- DatosANDI[with(DatosANDI,order(RID,EXAMDATE)),]
    #numero de pacientes (no observaciones)
    length(unique(DatosANDI$RID))
    
    
    #create the subsets of type of patients
    baseline <- subset(DatosANDI,VISCODE=="bl" & DX=="MCI")
    #control
    cn <- subset(DatosANDI,VISCODE=="bl" & DX_bl=="CN")
    conEvento <- subset(DatosANDI,DXCHANGE==5)
    
    #repetidos
    conEvento<-conEvento[!duplicated(conEvento$RID),]
    baselineDeEvento <- baseline[which(baseline$RID %in% conEvento$RID),]
    baseline$status <- 0
    # changing status
    baseline[which(baseline$RID %in% conEvento$RID),]$status<-1
    
    #load the csv with the important MRI variables
    varnames <- read.delim("./varnames.csv",sep=",")
    colnames(varnames)<-c("names","description")
    baseline<-as.data.frame(baseline);
    baselineSoloMRI <- baseline[,c("RID","status",as.character(varnames$names))]
    
    #normalizo datos de MRI
    mriIzquierdas <- varnames[grep("Left",varnames$description),]
    mriIzquierdas$description <- gsub("Left","",mriIzquierdas$description)
    mriDerechas <- varnames[grep("Right",varnames$description),]
    mriDerechas$description <- gsub("Right","",mriDerechas$description)
    
    mergeMRI <- merge(mriIzquierdas, mriDerechas, by="description")
    varnames[which((!varnames$names %in% mergeMRI$names.x) & (!varnames$names %in% mergeMRI$names.y)),]
    
    #transformacion de variables de izquierdas y derechas
    sujetos <- baselineSoloMRI
    for(i in 1:nrow(mergeMRI)) {
      column = mergeMRI[i,]
      diffLeftRight <- abs(sujetos[,as.character(column$names.x)]-sujetos[,as.character(column$names.y)])
      sumLeftRight <- sujetos[,as.character(column$names.x)]+sujetos[,as.character(column$names.y)]
      sujetos[,as.character(column$names.x)] <- diffLeftRight;
      sujetos[,as.character(column$names.y)] <- sumLeftRight;
    }
    
    sujetosMCItoAD <- sujetos;
    #diferencias y sumas de cn
    sujetos <- as.data.frame(cn)
    for(i in 1:nrow(mergeMRI)) {
      column = mergeMRI[i,]
      diffLeftRight <- abs(sujetos[,as.character(column$names.x)]-sujetos[,as.character(column$names.y)])
      sumLeftRight <- sujetos[,as.character(column$names.x)]+sujetos[,as.character(column$names.y)]
      sujetos[,as.character(column$names.x)] <- diffLeftRight;
      sujetos[,as.character(column$names.y)] <- sumLeftRight;
    }
    sujetosNC <- sujetos;
    
    
    #scaling preprosess
    featuresToScale <- as.character(varnames$names)
    featuresToScale <- featuresToScale[-c(1,2)]
    MCIToADToScale <- sujetosMCItoAD[,featuresToScale]
    NCToScale <- sujetosNC[,featuresToScale]
    
    #saling using fresa.cad
    scale <- FRESAScale(data=MCIToADToScale,refFrame = NCToScale,method = "RankInv")
    scaledData <- scale$scaledData
    
    #imputing  preprocess
    normData <- cbind(sujetosMCItoAD[,c("RID","status","PTGENDER","APOE4")],scaledData)
    withOutNorm <- sujetosMCItoAD;
    
    baselineSoloMRISinColumnasNA <- normData[,colSums(is.na(normData))<(0.5*nrow(normData))]
    baselineSoloMRISinColumnasNAWN <- withOutNorm[,colSums(is.na(withOutNorm))<(0.5*nrow(withOutNorm))]
    baselineSoloMRISinColumnasNAWNWAR <- baselineSoloMRI[,colSums(is.na(baselineSoloMRI))<(0.5*nrow(baselineSoloMRI))]
    
    #elimino sujetos que tienen todo NA
    sujetosSinTodoNA<-baselineSoloMRISinColumnasNA[rowSums(is.na(baselineSoloMRISinColumnasNA))<309,]
    sujetosSinTodoNAWN<-baselineSoloMRISinColumnasNAWN[rowSums(is.na(baselineSoloMRISinColumnasNAWN))<309,]
    sujetosSinTodoNAWNAR<-baselineSoloMRISinColumnasNAWNWAR[rowSums(is.na(baselineSoloMRISinColumnasNAWNWAR))<309,]
    
    #cambio genero para impute
    sujetosSinTodoNA$PTGENDER<- gsub("Male",1,sujetosSinTodoNA$PTGENDER)
    sujetosSinTodoNA$PTGENDER<- gsub("Female",2,sujetosSinTodoNA$PTGENDER)
    sujetosSinTodoNA$PTGENDER<-as.numeric(sujetosSinTodoNA$PTGENDER)
    
    sujetosSinTodoNAWN$PTGENDER<- gsub("Male",1,sujetosSinTodoNAWN$PTGENDER)
    sujetosSinTodoNAWN$PTGENDER<- gsub("Female",2,sujetosSinTodoNAWN$PTGENDER)
    sujetosSinTodoNAWN$PTGENDER<-as.numeric(sujetosSinTodoNAWN$PTGENDER)
    
    sujetosSinTodoNAWNAR$PTGENDER<- gsub("Male",1,sujetosSinTodoNAWNAR$PTGENDER)
    sujetosSinTodoNAWNAR$PTGENDER<- gsub("Female",2,sujetosSinTodoNAWNAR$PTGENDER)
    sujetosSinTodoNAWNAR$PTGENDER<-as.numeric(sujetosSinTodoNAWNAR$PTGENDER)  
    
    
    #imputing NA with knn
    # library(FRESA.CAD)
    baselineImputed <-  nearestNeighborImpute(sujetosSinTodoNA)
    baselineImputed <- as.data.frame(baselineImputed)
    RIDS <- baselineImputed[,1]
    rownames(baselineImputed)<-RIDS
    baselineImputed <- baselineImputed[,-c(1)]
    fnames <- colnames(baselineImputed)[-c(1,2)]
    columnsNotCorralated <- correlated_Remove(baselineImputed, fnames=fnames, thr=0.7)
    baselineImputedNoCorrelated<-baselineImputed[,c("status",as.character(columnsNotCorralated))]
    #save(baselineImputedNoCorrelated,file="baselineImputedNoCorrelated.RDATA")
    return(baselineImputedNoCorrelated)
    # WN
#   baselineImputedWN <-  nearestNeighborImpute(sujetosSinTodoNAWN)
#   baselineImputedWN <- as.data.frame(baselineImputedWN)
#   RIDS <- baselineImputedWN[,1]
#   rownames(baselineImputedWN)<-RIDS
#   baselineImputedWN <- baselineImputedWN[,-c(1)]
#   #save(baselineImputedWN,file="baselineImputedWN.RDATA")
#   #load("baselineImputedWN.RDATA")
#   
#   #WNAR
#   baselineImputedWNAR <-  nearestNeighborImpute(sujetosSinTodoNAWNAR)
#   baselineImputedWNAR <- as.data.frame(baselineImputedWNAR)
#   RIDS <- baselineImputedWNAR[,1]
#  rownames(baselineImputedWNAR)<-RIDS
#  baselineImputedWNAR <- baselineImputedWNAR[,-c(1)]
#  
#  fnames <- colnames(baselineImputedWNAR)[-c(1,2)]
#  columnsNotCorralated <- correlated_Remove(baselineImputedWNAR, fnames=fnames, thr=0.7)
#  baselineImputedWNAR<-baselineImputedWNAR[,c("status",as.character(columnsNotCorralated))]
#  
#  #save(baselineImputedWNAR,file="baselineImputedWNAR.RDATA")
#  #load("baselineImputedWNAR.RDATA")
#  
#  
  }
  