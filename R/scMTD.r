#' A statistical multidimensional imputation method for single-cell RNA-seq data.
#' @param data, Gene expression matrix (gene by cell).
#' @param do.nor, Normalized step.Default is TRUE.
#' @param do.log, logarithm transformation step.Default is TRUE.
#' @param cores, A integer specifying the number of cores used for parallel computation.
#' @export
#' @import foreach
#' @import parallel
#' @import doParallel
#' @import TSCAN
#' @return An imputation matrix
#' @examples
#' library("scMTD")
#' data(data)
#' imputed_data<-scMTD(data,do.nor=TRUE,do.log=TRUE,cores=2)



scMTD<-function(data,do.nor=TRUE,do.log=TRUE,cores){
  #### Preprocess data
  message("Data preprocessing...")
  #normalization
  nor <- function(x) {
    if (x == TRUE) {
      exmatrix <- t(t(data)/colSums(data)) * 1e+06
    }
    else {
      exmatrix = data
    }
  }
  exmatrix <- nor(do.nor)
  #logarithm transformation
  logt <- function(x) {
    if (x == TRUE) {
      exmatrix <- log2(exmatrix + 1)
    }
    else {
      exmatrix = exmatrix
    }
  }
  exmatrix <- logt(do.log)




  # order cells
  message("Order cells...")
  Ordercell <- function(exmatrix) {
    warnings("off")
    options(warn = -1)
    procdata <- preprocess(exmatrix,takelog = FALSE)
    lpsmclust <- exprmclust(procdata)
    cellorder<-TSCANorder(lpsmclust,MSTorder = c(1:length(table(lpsmclust$clusterid))),orderonly = T,flip = F)
    return(cellorder)
  }

  cellorder<-Ordercell(exmatrix)
  exmatrix<-exmatrix[,cellorder]



  #pseudo cell matrix
  message("Pseudo cell matrix...")
  mnn=5*(dim(exmatrix)[2]%/%1000+1)
  pseudocellm<- function(exmatrix,mnn){
    if((dim(exmatrix)[2]/mnn)==(dim(exmatrix)[2]%/%mnn)){
      n<-dim(exmatrix)[2]%/%mnn}else{
        n<-dim(exmatrix)[2]%/%mnn+1
      }
    supmatix<-matrix(0,ncol =mnn*n-dim(exmatrix)[2] ,nrow = dim(exmatrix)[1])
    supmatix<-cbind(exmatrix,supmatix)
    pseudomatrix<-matrix(0,ncol =n,nrow = dim(exmatrix)[1])
    for (i in 1:n) {
      for (j in 1:mnn) {
        pseudomatrix[,i]=pseudomatrix[,i]+supmatix[,(i-1)*mnn+j]
      }
    }
    rownames(pseudomatrix)<-rownames(supmatix)
    colnames(pseudomatrix)<- paste("pseudocell_", 1:dim(pseudomatrix)[2], sep = "")
    colnames(pseudomatrix)<- c(1:dim(pseudomatrix)[2])
    return(pseudomatrix)
  }

  pseudomatrix<-pseudocellm(exmatrix,mnn)


  # Dropout  probabilities
  message("Dropout probabilities...")
  p<-function(exmatrix,mnn){
    if((dim(exmatrix)[2]/mnn)==(dim(exmatrix)[2]%/%mnn)){
      n<-dim(exmatrix)[2]%/%mnn}else{
        n<-dim(exmatrix)[2]%/%mnn+1
      }
    aedata<-matrix(0,ncol = n,nrow = dim(exmatrix)[1])
    for (i in 1:(n-1)) {
      aedata[,i]=rowMeans(exmatrix[,(1+(i-1)*mnn):(i*mnn)])
    }
    aedata[,n]=rowMeans(exmatrix[,(1+(i-1)*mnn):dim(exmatrix)[2]])
    zcm<-matrix(0,ncol = n,nrow = dim(exmatrix)[1])
    for (j in 1:dim(exmatrix)[1]){
      for (i in 1:(n-1)) {
        zcm[j,i]=length(which(exmatrix[j,(1+(i-1)*mnn):(i*mnn)]==0))/mnn
      }
      zcm[j,n]=length(which(exmatrix[j,(1+(i-1)*mnn):dim(exmatrix)[2]]==0))/(dim(exmatrix)[2]-(i-1)*mnn)
    }
    data0<-cbind(zcm,aedata)
    data0<-as.data.frame(data0)
    e<-matrix(0,ncol =n,nrow = 2)
    for (c in 1: n) {
      anes=glm(zcm[,c]~aedata[,c],family=quasibinomial(link = "logit"),data=data0,control=list(maxit=100))
      e[,c]<-anes$coefficients
    }
    lable<-c(rep(0,dim(exmatrix)[2]))
    for (i in 1:(n-1)) {
      lable[(1+(i-1)*mnn):(i*mnn)]<-c(rep(i,mnn))
    }
    lable[(1+(n-1)*mnn):dim(exmatrix)[2]]<-c(rep(n,(dim(exmatrix)[2]-(i-1)*mnn)))
    pm<-exmatrix
    for(c in 1:n){
      for (j in 1:dim(exmatrix)[2]){
        for (i in 1:dim(exmatrix)[1]){
          if (lable[j] == c)
            pm[i,j]=exp(e[1,c]+e[2,c]*exmatrix[i,j])/(1+exp(e[1,c]+e[2,c]*exmatrix[i,j]))
        }
      }
    }

    return(pm)
  }

  dropoutp<-p(exmatrix,mnn)

  #cell state network
  message("cell-state specific network...")
  cogene<-function(pseudomatrix){
    log2_FC<-matrix(0,nrow = nrow(pseudomatrix),ncol =(ncol(pseudomatrix)-1))
    rownames(log2_FC)<-rownames(exmatrix)
    for(i in 1:nrow(log2_FC)){
      for (j in 1:ncol(log2_FC)) {
        log2_FC[i,j]<-log2((pseudomatrix[i,j]+0.001)/(pseudomatrix[i,(j+1)]+0.001))
      }
    }
    warnings("off")
    options(warn = -1)
    cormatrix<-cor(t(prcomp(log2_FC)$x[,1:20]))
    Mcogene<-list()
    for (i in 1:dim(pseudomatrix)[1]){
      Mcogene[[i]]<-c(which(cormatrix[i,]>=0.8))
    }
    res<-Mcogene
    return(res)
  }

  Mcogene<-cogene(pseudomatrix)

  corgene<-function(i,m){
    edgem<-list()
    edgem0<-list()
    geneid<-c(as.numeric(Mcogene[[i]]))
    dx<-abs(pseudomatrix[i,]-pseudomatrix[i,m])
    ss<-sort(dx)
    for (j in 1:length(geneid)) {
      edgem[[geneid[j]]]<-list()
      edgem0[[geneid[j]]]<-list()

      dy<-abs(pseudomatrix[geneid[j],]-pseudomatrix[geneid[j],m])
      ss1<-sort(dy)

      for (mm in 1:5) {
        a=0.1*(mm)
        cellidx<-ss[1:round(a*dim(pseudomatrix)[2])]
        cellidy<-ss1[1:round(a*dim(pseudomatrix)[2])]
        p<-length(intersect(as.character(names(cellidx)),as.character(names(cellidy))))/ dim(pseudomatrix)[2]
        if((p-a^2)<0.01){
          edgem0[[geneid[j]]][[mm]]<-0
        }else{
          edgem0[[geneid[j]]][[mm]]<-1
        }
      }
      if(sum(as.numeric(edgem0[[geneid[j]]]))==5){

        edgem[[geneid[j]]]<-1 }else{edgem[[geneid[j]]]<-"NULL"}
    }

    return(which(edgem!="NULL"))
  }

  log2_FC<-matrix(0,nrow = nrow(pseudomatrix),ncol =(ncol(pseudomatrix)-1))
  rownames(log2_FC)<-rownames(exmatrix)
  for(i in 1:nrow(log2_FC)){
    for (j in 1:ncol(log2_FC)) {
      log2_FC[i,j]<-log2((pseudomatrix[i,j]+0.001)/(pseudomatrix[i,(j+1)]+0.001))
    }
  }

  cormatrix<-cor(t(prcomp(log2_FC)$x[,1:20]))

  ####
  message("gene-level imputation...")
  imputedmatrix<-exmatrix
  dd<-20
  prexgene<-function(i,m){
    if((dim(exmatrix)[2]/mnn)==(dim(exmatrix)[2]%/%mnn)){
      n=dim(exmatrix)[2]%/%mnn
      ll=corgene(i,m)
      if(length(ll)>dd){
        lll=ll[order(-as.numeric(cormatrix[i,ll]))[1:dd]]
        a=exmatrix[lll,(1+(m-1)*mnn):(m*mnn)]
        b=as.numeric(cormatrix[i,lll])/sum(as.numeric(cormatrix[i,lll]))
        imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=b%*%a

      }else if(length(ll)>1&length(ll)<=dd){
        a=exmatrix[ll,(1+(m-1)*mnn):(m*mnn)]
        b=as.numeric(cormatrix[i,ll]) /sum(as.numeric(cormatrix[i,ll]))
        imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=b%*%a

      }
      else{
        imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=exmatrix[i,(1+(m-1)*mnn):(m*mnn)]
      }
      return(imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)])
    }


    else{
      n<-dim(exmatrix)[2]%/%mnn+1
      if(m<n){
        ll<-corgene(i,m)
        if(length(ll)>dd){
          lll=ll[order(-as.numeric(cormatrix[i,ll]))[1:dd]]
          a=exmatrix[lll,(1+(m-1)*mnn):(m*mnn)]
          b=as.numeric(cormatrix[i,lll])/sum(as.numeric(cormatrix[i,lll]))
          imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=b%*%a


        }else if(length(ll)>1&length(ll)<=dd){
          a=exmatrix[ll,(1+(m-1)*mnn):(m*mnn)]
          b=as.numeric(cormatrix[i,ll]) /sum(as.numeric(cormatrix[i,ll]))
          imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=b%*%a
        }
        else{
          imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=exmatrix[i,(1+(m-1)*mnn):(m*mnn)]
        }
        return(imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)])

      }else{

        ll<-corgene(i,m)
        if(length(ll)>dd){

          lll=ll[order(-as.numeric(cormatrix[i,ll]))[1:dd]]
          a=exmatrix[lll,(1+(m-1)*mnn):dim(exmatrix)[2]]
          b=as.numeric(cormatrix[i,lll])/sum(as.numeric(cormatrix[i,lll]))
          imputedmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]]==b%*%a

        }else if(length(ll)>1&length(ll)<=dd){
          a=exmatrix[ll,(1+(m-1)*mnn):dim(exmatrix)[2]]
          b=as.numeric(cormatrix[i,ll]) /sum(as.numeric(cormatrix[i,ll]))
          imputedmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]]=b%*%a
        }
        else{
          imputedmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]]=exmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]]
        }
        return(imputedmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]])
      }
    }
  }



  message("Gaussian kernel coefficient matrix...")
  ### Gaussian kernel coefficient matrix
  zerone<-exmatrix
  zerone[zerone>=0]=1
  gausskmatrix<-function(smatrix){
    g<-prcomp(t(smatrix))$x[,1:2]
    dw<-as.matrix(dist(g))
    tknn<-c(rep(1,dim(smatrix)[2]))
    for (i in 1:dim(smatrix)[2]){
      tknn[i]=mean(dw[i,])
    }
    gaussk=matrix(0,ncol = dim(smatrix)[2],nrow = dim(smatrix)[2])
    for (i in 1:dim(smatrix)[2]) {
      for (j in 1:dim(smatrix)[2]) {
        gaussk[i,j]=exp(-dw[i,j]/tknn[i])
      }
    }
    for (m in 1:dim(smatrix)[2]) {
      gaussk[,m]<-gaussk[,m]/sum(gaussk[,m])
    }
    return(as.matrix(gaussk))
  }

  gausscoeff<-function(ve,ma) {
    if(sum(ve)>0){
      for (i in 1:dim(ma)[2]) {
        ma[,i]=ve*ma[,i]
        ma[,i]=ma[,i]/sum(ma[,i])
      }
    }else{
      ma=ma
    }
    return(ma)
  }

  message("Cell-level imputation...")
  ### Cell-level imputation
  precellmatrix<-exmatrix
  if((dim(exmatrix)[2]/mnn)==(dim(exmatrix)[2]%/%mnn)){
    n<-dim(exmatrix)[2]%/%mnn
    for (m in 1:n) {
      gkmatrix=gausskmatrix(imputedmatrix[,(1+(m-1)*mnn):(m*mnn)])

      for (i in 1:dim(exmatrix)[2]) {
        precellmatrix[i,(1+(m-1)*mnn):(m*mnn)]=exmatrix[i,(1+(m-1)*mnn):(m*mnn)]%*%gausscoeff(zerone[i,(1+(m-1)*mnn):(m*mnn)],gkmatrix)
      }
    }
  }else{

    n<-dim(exmatrix)[2]%/%mnn+1
    for (m in 1:(n-1)) {
      gkmatrix=gausskmatrix(imputedmatrix[,(1+(m-1)*mnn):(m*mnn)])
      for (i in 1:dim(exmatrix)[2]) {
        precellmatrix[i,(1+(m-1)*mnn):(m*mnn)]=exmatrix[i,(1+(m-1)*mnn):(m*mnn)]%*%gausscoeff(zerone[i,(1+(m-1)*mnn):(m*mnn)],gkmatrix)
      }
    }
    m=n
    gkmatrix=gausskmatrix(imputedmatrix[,(1+(m-1)*mnn):dim(exmatrix)[2]])
    for (i in 1:dim(exmatrix)[2]) {
      precellmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]]=exmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]]%*%gausscoeff(zerone[i,(1+(m-1)*mnn):dim(exmatrix)[2]],gkmatrix)
    }

  }




  prexcell<-function(i,m){
    if((dim(exmatrix)[2]/mnn)==(dim(exmatrix)[2]%/%mnn)){
      n<-dim(exmatrix)[2]%/%mnn
      imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=precellmatrix[i,(1+(m-1)*mnn):(m*mnn)]
      return(imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)])
    }

    else{
      n<-dim(exmatrix)[2]%/%mnn+1
      if(m<n){
        imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=precellmatrix[i,(1+(m-1)*mnn):(m*mnn)]
        return(imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)])

      }else{
        imputedmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]]= precellmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]]
        return(imputedmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]])
      }
    }
  }


  message("Multidimensional imputation...")
  ### Multidimensional imputation
  dropp<-function(i,m){
    if((dim(exmatrix)[2]/mnn)==(dim(exmatrix)[2]%/%mnn)){
      n<-dim(exmatrix)[2]%/%mnn
      return(dropoutp[i,(1+(m-1)*mnn):(m*mnn)])
    }

    else{
      n<-dim(exmatrix)[2]%/%mnn+1
      if(m<n){
        return(dropoutp[i,(1+(m-1)*mnn):(m*mnn)])

      }else{
        return(dropoutp[i,(1+(m-1)*mnn):dim(exmatrix)[2]])
      }
    }
  }

  obzex<-function(i,m){
    if((dim(exmatrix)[2]/mnn)==(dim(exmatrix)[2]%/%mnn)){
      n<-dim(exmatrix)[2]%/%mnn
      return(exmatrix[i,(1+(m-1)*mnn):(m*mnn)])
    }

    else{
      n<-dim(exmatrix)[2]%/%mnn+1
      if(m<n){
        return(exmatrix[i,(1+(m-1)*mnn):(m*mnn)])

      }else{
        return(exmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]])
      }
    }
  }


  prex<-function(i,m){
    prexc<-prexcell(i,m)
    prexg<-prexgene(i,m)
    sdc<-sd(prexc)
    sdg<-sd(prexg)
    if(sdc*sdg!=0){
      prexp<-(prexg*(sdc/(sdc+sdg))+prexc*(sdg/(sdc+sdg)))*dropp(i,m)+(1-dropp(i,m))*obzex(i,m)

    }else{
      prexp<-((prexg+prexc)/2)*dropp(i,m)+(1-dropp(i,m))*obzex(i,m)
    }
    return(prexp)
  }

  cl<-makeCluster(cores)
  registerDoParallel(cl, cores=cores)
  if((dim(exmatrix)[2]/mnn)==(dim(exmatrix)[2]%/%mnn)){
    n<-dim(exmatrix)[2]%/%mnn
    res <- foreach(m=1:n, .combine='cbind') %dopar% {
      for (i in 1:dim(exmatrix)[1]) {
        imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=prex(i,m)
      }
      return(imputedmatrix[,(1+(m-1)*mnn):(m*mnn)])
    }
  }else{
    n<-dim(exmatrix)[2]%/%mnn+1
    res1 <- foreach(m=1:(n-1), .combine='cbind') %dopar% {
      for (i in 1:dim(exmatrix)[1]) {
        imputedmatrix[i,(1+(m-1)*mnn):(m*mnn)]=prex(i,m)
      }
      return(imputedmatrix[,(1+(m-1)*mnn):(m*mnn)])
    }
    m=n
    for (i in 1:dim(exmatrix)[1]) {
      imputedmatrix[i,(1+(m-1)*mnn):dim(exmatrix)[2]]=prex(i,m)
    }
    res<-cbind(res1,imputedmatrix[,(1+(m-1)*mnn):dim(exmatrix)[2]])
  }

  res<-res[,colnames(data)]
  stopImplicitCluster()
  stopCluster(cl)
  message("Done!")
  return(res)
}
