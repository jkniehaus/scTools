#' Binomial test
#'
#' This function tests the probability of a difference in the proportion of non-zero values b/w groups; used by other functions. Adopted from Shekhar 2016
#' @param object Seurat object
#' @param cells.1 Group 1; a working group within seurat object levels/Idents
#' @param cells.2 Group 2; if NULL = all other cells
#' @param effect.size Filter DEGs by effect size
#' @keywords binom
#' @export
#' @importFrom Seurat GetAssayData
#' @examples
#' binomcount.test()
binomcount.test=function(object, cells.1,cells.2, effect.size) {
  x=GetAssayData(object=object,layer='counts')
  
  #Test for enrichments in cluster #1
  m=rowSums(x[,cells.2]>0) #Number of cells expressing marker in cluster 2
  m1 = m; m1[m==0]=1; # Regularization. Add a pseudo count of 1 to unexpressed markers to avoid false positives
  n=rowSums(x[,cells.1]>0) #number of cells expressing marker in cluster 1
  #Find the probability of finding n or more cells +ve for marker in cluster 1 given the fraction in cluster 2
  pv1 = stats::pbinom(n, length(cells.1), m1/length(cells.2), lower.tail = FALSE) + dbinom(n, length(cells.1), m1/length(cells.2))
  
  log_fold_express = log((n+1)*length(cells.2)/((m+1)*length(cells.1)),base=2) #log10 proportion of expressing cells
  
  d1 <- data.frame(log.effect=log_fold_express,pval=pv1)
  if(!is.na(effect.size)) {
    d1 <- subset(d1, log.effect >= effect.size)
  }
  d1 <- d1[order(d1$pval,decreasing=FALSE),]
  
  #Enrichments in cells.2
  n1 = n; n1[n==0]=1; # Regularization.  Add a pseudo count of 1 to unexpressed markers to avoid false positives
  #Find the probability of finding m or more cells +ve for marker in cluster 2 given the fraction in cluster 1
  pv2 = stats::pbinom(m, length(cells.2), n1/length(cells.1), lower.tail=FALSE) + dbinom(m, length(cells.2), n1/length(cells.1))
  d2 <- data.frame(log.effect=log_fold_express,pval=pv2)
  if(!is.na(effect.size)) {
    d2 <- subset(d2, log.effect <= -effect.size)
  }            
  d2 <- d2[order(d2$pval,decreasing=FALSE),]
  
  d = rbind(d1, d2);
  d = d[order(d$pval, decreasing=FALSE),]
  return(d)
} 




###identify marker genes based on presence/absence binomial test (can be done vs single, muliple, or all other clusters. Adopted from Shekhar 2016
###Note, use of seurat_clusters as cluster IDs.
###' markers.binom
#'
#' This function uses the binomial test function to identify differentially expressed genes
#' @param object Seurat object
#' @param clust.1 Group 1; a working group within seurat object levels/Idents
#' @param clust.2 Group 2; if NULL = all other cells
#' @param effect.size Filter DEGs by effect size
#' @param assay Seurat object assay
#' @param posFrac Include only genes expressed by at least posFrac proportion of cells in either group; default = 25%
#' @keywords markers
#' @export
#' @examples
#' markers.binom()
markers.binom = function(object, clust.1,clust.2=NULL,effect.size,assay,posFrac=0.25,posOnly=FALSE) {
  genes.use=rownames(object)
  clust.use=seurat::Idents(object)
  cells.1=names(clust.use[which(clust.use%in%clust.1)])
  
  if (is.null(clust.2)) {
    clust.2="rest"
    cells.2=names(clust.use)
    cells.2=cells.2[!(cells.2%in%cells.1)]
  } else {
    cells.2=names(clust.use[which(clust.use%in%clust.2)])
  }
  seurat::DefaultAssay(object)=assay
  Count.mat = object[[assay]]$counts
  TPM.mat = object[[assay]]$data
  result=binomcount.test(object, cells.1, cells.2, effect.size)
  if (nrow(result)>1){
    #ensure countmat = dataframe; throws error for only 1 DEG
    #get prop posFrac per marker gene
    #posFrac.1 = apply(Count.mat[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2)) 
    posFrac.1 <- round(rowSums(Count.mat[rownames(result), cells.1] > 0) / length(cells.1),digits = 2)
    #posFrac.2 = apply(Count.mat[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))
    posFrac.2 <- round(rowSums(Count.mat[rownames(result), cells.2] > 0) / length(cells.2),digits= 2)
  } else if (nrow(result)==1) {
    posFrac.1 = length(which((Count.mat[rownames(result),cells.1]>0)))/length(cells.1)
    posFrac.2 = length(which((Count.mat[rownames(result),cells.2]>0)))/length(cells.2)
  }
  
  if (clust.2=="rest" && posOnly==FALSE){
    genes.include = posFrac.1 >= posFrac
  } else{
    genes.include = (posFrac.1 >= posFrac) | (posFrac.2 >= posFrac)
  }
  
  result = result[genes.include,]
  result = result[order(abs(result$log.effect), decreasing=TRUE),]
  if (nrow(result)>0){
    if (nrow(result)>1){
      #Mean number of transcripts per cell
      nTrans.1 = apply(TPM.mat[rownames(result), cells.1], 1, function(x) round(mean(x),3))
      nTrans.2 = apply(TPM.mat[rownames(result), cells.2], 1, function(x) round(mean(x),3))
    } else if (nrow(result)==1) {
      nTrans.1=mean(TPM.mat[rownames(result),cells.1])
      nTrans.2=mean(TPM.mat[rownames(result),cells.2])
    }
    result[,paste0("nTrans_", clust.1)] = nTrans.1
    result[, paste0("nTrans_", clust.2)] = nTrans.2
    #THIS IF STATEMENT NEEDS TO BE RUN, OR IT WILL FAIL WHEN NO DEGS EXIST BETWEEN CLUSTERS
    result$ClusterID=clust.1 
  }           
  return(result)
} 



####determine cluster centroids... This should be easy. they're using PCs.
#####' Cluster Centrois
#'
#' Identify cluster centroids in PC space; used by other functions
#' @param object Seurat object
#' @param reduction.use dimensionality reduction to use
#' @param pcs.use Number of PCs to use
#' @param cells.use Group of cells to calculate centroids; default = all
#' @keywords centroids
#' @export
#' @examples
#' ClusterCentroids()
ClusterCentroids = function(object,reduction.use="pca", pcs.use = 1:10, cells.use=NULL) {
  cells.use = colnames(object)
  group.use=seurat::Idents(object)[cells.use]
  if (reduction.use == "pca"){
    data.use = object@reductions$pca@cell.embeddings[cells.use,pcs.use]
  }
  
  centroids = c()
  for (i in levels(group.use)){
    cells.in.cluster = names(Idents(object))[Idents(object) %in% i]
    cells.in.cluster = cells.in.cluster[cells.in.cluster %in% cells.use]
    centroids = rbind(centroids, colMeans(data.use[cells.in.cluster,]))
  }
  centroids = as.data.frame(centroids)
  colnames(centroids) = colnames(data.use)
  rownames(centroids) = as.numeric(levels(Idents(object)))
  
  return(centroids)
}





## Centroid distances
#'
#' This function calculates next closest cluster for each cell; used for merging
#' @param object Seurat object
#' @param reduction.use Dimenionality reduction used to calculate distances
#' @param dist.type How to calculate distances; 'nn' or 'centroid' 
#' @param pcs.use PC dims used to calculate distances
#' @param cells.use subset of cells to calc distances; default NULL = all cells
#' @keywords binom
#' @export
#' @examples
#' ComputeClusterDistances()
ComputeClusterDistances=function(object,reduction.use="pca",dist.type="nn", pcs.use = 1:10, cells.use=NULL) {
  cells.use =  colnames(object)
  group.use=seurat::Idents(object)[cells.use]
  if (reduction.use == "pca"){
    data.use = object@reductions$pca@cell.embeddings[cells.use,pcs.use]
    centroids = ClusterCentroids(object, reduction.use="pca", pcs.use=pcs.use, cells.use=cells.use)
  }
  
  if (dist.type=="centroid"){
    clust.dists = as.matrix(dist(centroids, upper=TRUE))
    diag(clust.dists) = 1e6
  }
  
  num.clust = length(levels(group.use))
  
  
  if (dist.type == "nn"){
    clust.dists = matrix(0, nrow=num.clust, ncol=num.clust)
    diag(clust.dists) = 1e6
    rownames(clust.dists) = levels(group.use)
    colnames(clust.dists) = rownames(clust.dists)
    for (i in 1:nrow(clust.dists)){
      for(j in ((i+1):ncol(clust.dists))){
        if (j>nrow(clust.dists)) break
        cells.in.cluster_i = names(seurat::Idents(object))[seurat::Idents(object) %in% i]
        cells.in.cluster_i = cells.in.cluster_i[cells.in.cluster_i %in% cells.use]
        cells.in.cluster_j = names(seurat::Idents(object))[seurat::Idents(object) %in% j]
        cells.in.cluster_j = cells.in.cluster_j[cells.in.cluster_j %in% cells.use]
        
        nnA = RANN::nn2(data.use[cells.in.cluster_i,], query = centroids[j,], k=1)
        nnB = RANN::nn2(data.use[cells.in.cluster_j,], query = centroids[i,],k=1)
        clust.dists[i,j] = min(c(nnA$nn.dists, nnB$nn.dists))
        
        clust.dists[j,i] = clust.dists[i,j]
      }
    }
  }
  
  colnames(clust.dists) = c(1:ncol(clust.dists))
  rownames(clust.dists) = colnames(clust.dists)
  return(clust.dists)
}




#' Merge clusters
#'
#' This function merges clusters that don't have enough DEGs to delineate them
#' @param object Seurat object
#' @param min.de.genes Min # DEGs calculated by binom to keep clusters distinct
#' @param effect.size Min effect size to call DEG in binom test
#' @param pval.cutoff Max pvalue to call DEG in binom test
#' @param pcs.use If merging clusters, PCs used to reassign cells to nearest cluster
#' @param assay Seurat assay used
#' @keywords binom
#' @export
#' @examples
#' merge.clusters.DE()
merge.clusters.DE=function(object,min.de.genes, effect.size,pval.cutoff, pcs.use,assay='RNA') {
  
  genes.use = rownames(object)
  clust.test = as.numeric(levels(seurat::Idents(object)))
  # if (is.null(tag)){
  #   filename = "CLUSTER_PAIRWISE_MARKERS.txt"
  # } else {
  #   filename = paste0("CLUSTER_PAIRWISE_MARKERS_",tag,".txt")
  # }
  # zz = file(filename,open="wt")
  
  
  num.clust=length(clust.test) 
  print(paste0("Starting with ", num.clust, " clusters"))
  
  pass.thresh=1e6*data.frame(diag(length(levels(seurat::Idents(object))))); 
  
  for (i in setdiff(as.numeric(levels(seurat::Idents(object))), clust.test)){
    pass.thresh[i,]=1e6; pass.thresh[,i]=1e6;
  } 
  
  dist.clust = pass.thresh
  
  #Find the number of differentially expressed genes between every pair of clusters
  for(k in 1:num.clust) {
    i=clust.test[k]
    print(paste0("Testing Cluster ", i))
    for(m in ((k+1):num.clust)) {
      j=clust.test[m]
      print(j)
      if (m>num.clust) break
      if (pass.thresh[i,j]==0) {
        marker=markers.binom(object=object, clust.1=i,clust.2=j,effect.size=effect.size,assay=assay)
        paste0('found markers for',i,' vs ',j)
        
        
        num.de.genes = nrow(marker)
        pass.thresh[i,j]=num.de.genes; pass.thresh[j,i]=pass.thresh[i,j];
        
      }
      
    }
    
    print(pass.thresh[i,])
    
  }
  
  colnames(pass.thresh) = levels(Idents(object))
  rownames(pass.thresh) = levels(Idents(object))
  
  write.table(pass.thresh, file=paste0("DE_genes_matrix_2.txt"), sep="\t", quote=FALSE)
  
  #iteratively merge clusters
  min.val = min(pass.thresh)
  min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
  min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
  min.val.ind$val = min(pass.thresh)
  rownames(min.val.ind) = 1:nrow(min.val.ind)
  
  merge.ind=-1
  while(min.val <= min.de.genes & length(unique(Idents(object)))>2) {
    merge.ind=merge.ind+1
    
    #In case of ties, merge clusters that are closest in PC space
    print('ComputingClusterDistances')
    clust.dists = ComputeClusterDistances(object, reduction.use="pca", dist.type="centroid", pcs.use=pcs.use) 
    ind.min=which.min(clust.dists[cbind(min.val.ind$row, min.val.ind$col)])
    test.1 = min.val.ind[ind.min,]$row; test.2 = min.val.ind[ind.min,]$col
    
    if (pass.thresh[test.1,test.2]<= min.de.genes) {
      seurat::Idents(object)[which(seurat::Idents(object)==test.2)]=test.1
      pass.thresh = pass.thresh[-test.2,]; pass.thresh = pass.thresh[,-test.2]
      old.group.levels = as.numeric(levels(Idents(object)))
      old.group.levels = setdiff(old.group.levels, test.2)
      clust.test = setdiff(clust.test, test.2)
      
      seurat::Idents(object) = droplevels(seurat::Idents(object))
      levels(seurat::Idents(object)) = c(1:length(levels(seurat::Idents(object))))
      object$seurat_clusters = seurat::Idents(object)
      
      new.group.levels = as.numeric(levels(seurat::Idents(object)))
      names(new.group.levels) = as.character(old.group.levels)
      clust.test = new.group.levels[as.character(clust.test)]
      
      
      
      #Recompute pairwise markers for merged cluster
      print(paste0("Recomputing pairwise markers for new clust ", test.1))
      for (i in setdiff(clust.test, test.1)){
        print(i)
        marker = markers.binom(object,test.1,i,effect.size = effect.size,assay = assay)
        P2 = stats::p.adjust(marker$pval, method="fdr")
        marker$pval = P2
        marker.pass=subset(marker,pval<pval.cutoff)
        pass.thresh[test.1,i]=2*min(nrow(subset(marker.pass, log.effect>0)),nrow(subset(marker.pass, log.effect<0))); pass.thresh[i,test.1]=pass.thresh[test.1,i];
        #pass.thresh[test.1,i]=nrow(marker.pass); 
        pass.thresh[i,test.1]=pass.thresh[test.1,i];
        
      }
      
    }
    colnames(pass.thresh) = 1:length(levels(Idents(object)))
    rownames(pass.thresh) = colnames(pass.thresh)
    
    min.val = min(pass.thresh)
    min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
    min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
    min.val.ind$val = min(pass.thresh)
    rownames(min.val.ind) = 1:nrow(min.val.ind)
    
  }
  return(object)
}
#' Find doublets
#'
#' This function uses the DoubletFinder package to identify and remove doublets; parallel computing supported
#' @param object Seurat object
#' @param cores Number of working cores for parallelization
#' @param pcs Number of PCs used in doubletfinder
#' @param assay Seurat assay to normalize/scale data on. value 'SCT' performs 'SCTransform' and increases processing time significantly
#' @keywords doublets
#' @export
#' @examples
#' findDoublets()
findDoublets <- function(object,cores = 1, pcs = 30,assay='RNA') {
  if (assay != 'SCT') {
  seurat::DefaultAssay(object) = assay
  print('Normalizing data')
  object <- seurat::NormalizeData(object)
  object <- seurat::FindVariableFeatures(object)
  object <- seurat::ScaleData(object,vars.to.regress='percent.mt')
  } else if (assay == 'SCT') {
  object <- seurat::SCTransform(object, vst.flavor = 'v2', verbose = FALSE, vars.to.regress = 'percent.mt')
  }
  object <- seurat::RunPCA(object, verbose = FALSE, npcs = pcs)
  #object <- RunUMAP(object, dims = 1:pcs)
  print('making artificial doublets')
  if (assay == 'SCT') {
  sweep.res.list <- DoubletFinder::paramSweep(object, PCs = 1:pcs, sct = TRUE, num.cores = cores)
  } else {
  sweep.res.list <- DoubletFinder::paramSweep(object, PCs = 1:pcs, num.cores = cores)
  }
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  pk.use <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  print('runnding doubletFinder')
  if (assay == 'SCT') {
  object <- DoubletFinder::doubletFinder(object, PCs = 1:pcs, pN = 0.25, pK = 0.06, nExp = round(0.08 * ncol(object)), sct = TRUE)
  } else {
  object <- DoubletFinder::doubletFinder(object, PCs = 1:pcs, pN = 0.25, pK = 0.06, nExp = round(0.08 * ncol(object)))
  }
  dfc=grep('classifications',colnames(object@meta.data))
  pnn=grep('pANN',colnames(object@meta.data))
  colnames(object@meta.data)[dfc]='doublet_finder'
  colnames(object@meta.data)[pnn]='pANN'
  s=subset(object,subset=doublet_finder=='Singlet')
  s[['pANN']]=NULL
  s[['doublet_finder']]=NULL
  return(s)
}

#' Modify Violin
#'
#' This function is used in StackedVlnPlot for formatting
#' @param object Seurat object
#' @param feature gene
#' @param pt.size point size
#' @param plot.margin Margins
#' @keywords  violin
#' @export
#' @examples
#' modify_vlnplot()()
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- seurat::VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin )
  return(p)
}

#' ExtractMax
#'
#' This function is used in StackedVlnPlot to get max Y-axis values
#' @param p plot
#' @keywords extract
#' @export
#' @examples
#' extract_max()
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


#' Stacked Violin
#'
#' This function creates a stacked violin plot from the seurat VlnPlot function
#' @param obj Seurat object
#' @param features gene names
#' @param plot.margin Plot Margins
#' @keywords stacked violin
#' @export
#' @examples
#' StackedVlnPlot()
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, hjust=1,vjust=1), axis.ticks.x = element_line()) 
  
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



###identify marker genes based on presence/absence binomial test (can be done vs single, muliple, or all other clusters.
###Note, use of seurat_clusters as cluster IDs.
#' Metagroup Binomial test
#'
#' This function applies the markers binom function within specific metadata groups
#' @param object Seurat object
#' @param metagroup colname in object metadata
#' @param clust.1 Group 1
#' @param clust.2 Group 2; if NULL = all other cells
#' @param effect.size Filter DEGs by effect size
#' @param assay seurat object assay
#' @keywords metagroup
#' @export
#' @importFrom Seurat Idents DefaultAssay GetAssayData 
#' @examples
#' markers.binom.metagroup()
markers.binom.metagroup = function(object,metagroup, clust.1,clust.2=NULL,effect.size,assay='RNA',posFrac=0.25) {
  DefaultAssay(object)=assay
  genes.use=rownames(object)
  clust.use=unlist(object[[metagroup]])
  names(clust.use)=colnames(object)
  cells.1=names(clust.use[which(clust.use%in%clust.1)])
  
  if (is.null(clust.2)) {
    clust.2="rest"
    cells.2=names(clust.use)
    cells.2=cells.2[!(cells.2%in%cells.1)]
  } else {
    cells.2=names(clust.use[which(clust.use%in%clust.2)])
  }
  
  Count.mat = GetAssayData(object,assay=assay,layer='counts')
  TPM.mat=GetAssayData(object,assay=assay,layer='data')
  result=binomcount.test(object, cells.1, cells.2, effect.size)
  
  
  posFrac.1 = apply(Count.mat[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
  posFrac.2 = apply(Count.mat[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))
  
  if (clust.2=="rest"){
    genes.include = posFrac.1 >= posFrac
  } else{
    genes.include = (posFrac.1 >= posFrac) | (posFrac.2 >= posFrac)
  }
  
  result = result[genes.include,]
  result = result[order(abs(result$log.effect), decreasing=TRUE),]
  
  #Mean number of transcripts per cell
  if (nrow(result)>0){
    nTrans.1 = apply(TPM.mat[rownames(result), cells.1], 1, function(x) round(mean(x),3))
    nTrans.2 = apply(TPM.mat[rownames(result), cells.2], 1, function(x) round(mean(x),3))
    result[,paste0("nTrans_", clust.1)] = nTrans.1
    result[, paste0("nTrans_", clust.2)] = nTrans.2
    #THIS IF STATEMENT NEEDS TO BE RUN, OR IT WILL FAIL WHEN NO DEGS EXIST BETWEEN CLUSTERS
    result$ClusterID=clust.1
  }
  return(result)
}

#' Euler plots
#'
#' This function generates euler plots from a seurat object using the eulerr package
#' @param object Seurat object
#' @param genes Gene names
#' @param metacol Column in metadata to subset cells by
#' @param metavar group in metadata column to subset cells by
#' @param title Title in plot
#' @param assay Seurat object assay
#' @keywords euler
#' @export
#' @importFrom Seurat Idents DefaultAssay GetAssayData 
#' @importFrom eulerr euler
#' @examples
#' euler_plot()
euler_plot <- function(object, genes,metacol=NULL, metavar=NULL, title=NULL,assay='RNA') {
  if (!is.null(metacol)) {
    cells.use=which(object[[metacol]]==metavar)
    object=subset(object,cells=cells.use)
  }
  opdf=data.frame(row.names=colnames(object[[assay]]),a=rep(TRUE,ncol(object[[assay]])))
  for (i in 1:length(genes)){
    opdf[,i+1]=object[[assay]]$counts[genes[i],]>0
  }
  colnames(opdf)=c('Remaining',paste0(genes,'+'))
  if(!is.null(metacol)) {
    p=plot(euler(opdf,shape='ellipse'),counts=TRUE,quantities=list(type = c('counts','percent')),main=paste0(metacol,' = ',metavar))
  } else if (!is.null(title)) {
    p=plot(euler(opdf,shape='ellipse'),counts=TRUE,quantities=list(type = c('counts','percent')),main=paste0(title))
  } else {
    p=plot(euler(opdf,shape='ellipse'),counts=TRUE,quantities=list(type = c('counts','percent')))
  }
  return(p)
}

#' Euler plot v2
#'
#' This function generates euler plots with eulerr package; can be used across multiple assays
#' @param object Seurat object
#' @param genes Gene names with assay info; e.g. c('rna_Snap25','exons_Snap25')
#' @param metacol Column in metadata to subset cells by
#' @param metavar group in metadata column to subset cells by
#' @param title Title in plot
#' @keywords euler
#' @importFrom eulerr euler
#' @importFrom Seurat FetchData
#' @export
#' @examples
#' euler_plotv2()
euler_plotv2 <- function(object, genes,metacol=NULL, metavar=NULL, title=NULL) {
  #Note all assays need identical cells (none missing or renamed)
  if (!is.null(metacol)) {
    cells.use=which(object[[metacol]]==metavar)
    object=subset(object,cells=cells.use)
  }
  opdf=data.frame(row.names=colnames(object),a=rep(TRUE,ncol(object)))
  for (i in 1:length(genes)){
    opdf[,i+1] <- FetchData(object = object, vars = genes[i],layer = 'counts')[,1]>0
  }
  colnames(opdf)=c('Remaining',paste0(genes,'+'))
  if(!is.null(metacol)) {
    p=plot(euler(opdf,shape='ellipse'),counts=TRUE,quantities=list(type = c('counts','percent')),main=paste0(metacol,' = ',metavar))
  } else if (!is.null(title)) {
    p=plot(euler(opdf,shape='ellipse'),counts=TRUE,quantities=list(type = c('counts','percent')),main=paste0(title))
  } else {
    p=plot(euler(opdf,shape='ellipse'),counts=TRUE,quantities=list(type = c('counts','percent')))
  }
  return(p)
}

#' upset plot
#'
#' This function generates upset plots with UpSetR package
#' @param object Seurat object
#' @param genes Gene names
#' @param metacol Column in metadata to subset cells by
#' @param metavar group in metadata column to subset cells by
#' @keywords upset
#' @export
#' @importFrom UpSetR upset
#' @importFrom Seurat FetchData
#' @examples
#' upset_plot()
upset_plot <- function(object, genes, metacol=NULL, metavar=NULL) { 
  if (!is.null(metacol)) {
    cells.use=which(s[[metacol]]==metavar) 
    object=subset(object,cells=cells.use)
  }
  exp<-FetchData(object,vars=genes,layer='counts',assay='RNA')#dataframe 
  exp[exp>0]=1 
  if (length(which(colSums(exp)>0))>1) { 
    print(upset(exp,order.by='freq',empty.intersections='on')) 
    grid::grid.text(paste0(nrow(exp),' cells total'),x = 0.65, y=0.95, gp=gpar(fontsize=20))
  } else {
    print('Genes not expressed')
  }
}

#' Percent Expressing
#'
#' This function calculates the percentage of cells expressing a given gene across groups and pads non-expressing groups with zeros.
#' @param object Seurat object
#' @param features Gene names
#' @param assay Seurat object assay
#' @param group_by metadata variable to group cells by. defaults to Idents
#' @param split_by Further subset data by metadata variable
#' @keywords PercExp
#' @export
#' @importFrom Seurat DefaultAssay Idents FetchData 
#' @examples
#' PercExp()
PercExp <- function(object, features, assay=NULL, group_by='NULL', split_by=NULL) {
  if (is.null(assay) == FALSE) {
    DefaultAssay(object)=assay
  }
  if (is.null(group_by) == FALSE) {
    Idents(object) = group_by
  }
  
  dat <- FetchData(object,vars = features, layer = 'counts')
  if (is.null(split_by) == FALSE) {
  mat <- matrix(nrow = length(levels(object)), ncol = length(table(object[[split_by]])))
  for (i in 1:nrow(mat)){
    j=levels(object)[i]
    d2 <- dat[Idents(object) == levels(object)[i],] #subset data by cluster
    names(d2) <- rownames(dat)[Idents(object) == levels(object)[i]]
    for (ii in 1:ncol(mat)){
      nn <- colnames(object)[which(object[[split_by]] == names(table(object[[split_by]]))[ii])] #all cell IDs by orig ident
      if (length(nn) == 0) {
        mat[i,ii]=0
      } else {
        d3 <- d2[names(d2) %in% nn]
        if (length(d3) == 0) {
          mat[i,ii]=0
        } else {
        mat[i,ii] <- length(which(d3>0)) / length(d3) * 100
      }
    }
  }
  }
  rownames(mat) <-levels(object)
  colnames(mat) <- names(table(object[[split_by]]))
  return(mat)
  } else {
    mat=matrix(nrow = 1, ncol = length(levels(object)))
    for (i in 1:nrow(mat)){
      j=levels(object)[i]
      d2 <- dat[Idents(object) == levels(object)[i],] #subset data by cluster
      mat[1,i] <- length(which(d2 > 0)) / length(d2) * 100
    }
  }
}

#' Sex differences
#'
#' test whether there are differences in the proportion of GeneX+ cells between sexes by sample; generates plots of significant differences (p val <.05 or abs(effect size) > 1) by cluster. metadata column 'Sex' with values 'Male' or 'Female' must exist.
#' @param object Seurat object
#' @param gene Gene name
#' @param prefix File name prefix
#' @param assay seurat object assay
#' @param plot Whether to generate box plots split by sex
#' @keywords Sex differences
#' @export
#' @importFrom Seurat DefaultAssay 
#' @examples
#' sexDif()
sexDif <- function(object, gene,prefix, assay, plot = FALSE) {
  DefaultAssay(object)=assay
  male <- unique(object$orig.ident[which(object$Sex == 'Male')])
  maleindex <- which(names(table(object$orig.ident)) %in% male)
  female <- unique(object$orig.ident[which(object$Sex == 'Female')])
  femaleindex <- which(names(table(object$orig.ident)) %in% female)
  #create a dataframe of all samples + clusters
  l1=levels(object)
  l2=names(table(object$orig.ident))
  l3 <- unlist(lapply(l1, function(x) sapply(l2, function(y) paste(x, y, sep = "_"))))
  names(l3)=NULL
  
  percs <- rep(0,length(l3))
  names(percs)=l3
  perc <- PercExp(object, features = gene, group_by='seurat_clusters',split_by = 'orig.ident') #corrected for depth..
  percmat <- as.matrix(perc)
  
  percsub=percmat
  tt <- c()
  lt <- c()
  lowexp=c()
  es <- c()
  for (i in 1:nrow(percmat)) {
    print(i)
    dat <- percmat[i, ]
    dat <- data.frame(PercentExpressing=unlist(percmat[i,]),Sex=rep(0,ncol(percmat)),row.names=names(table(object$orig.ident)))
    dat[maleindex,2]='Male'
    dat[femaleindex,2]='Female'
    mdat <- dat[maleindex,1]
    fdat <- dat[femaleindex,1]
    #must be expressed in average of 5% of cells in either male or female samples
    if(mean(mdat)<10 && mean(fdat)<10 ){
      lowexp=c(lowexp,i)
      pval <- 1
      effect <- 0
    } else if (all(mdat==100) && all(fdat==100)){
      pval <- 1
      effect <- 0
    }else {
      tstat <- matrixTests::row_t_equalvar(x = mdat, y = fdat)
      if(all(mdat==0)) {
        effect=-mean(fdat)
      } else if(all(fdat==0)){
        effect=mean(mdat)
      } else {
        effect <- effectsize::cohens_d(mdat,fdat)
        effect=effect$Cohens_d
      }
      
      pval <- tstat$pvalue
    }
    pval <- ifelse(is.na(pval) == TRUE,1,pval)
    if (pval < 0.05 || abs(effect) > 1) {
      p=ggpubr::ggdotplot(dat,x="Sex",y="PercentExpressing",
                  add='boxplot',color="Sex",fill='Sex',
                  palette =c('red','blue')) + 
        ylim(0,max(dat[,1]+5)) + 
        NoLegend() +
        ggtitle(label = paste0(levels(object)[i],': ',gene))
      if (plot == TRUE) {
      pdf(paste0(prefix,'SexProp_Sig_',assay,levels(object)[i],'_',gene,'.pdf'),width=2,height=3)
      print(p)
      dev.off()
      }
    } else {
      p=ggpubr::ggdotplot(dat,x="Sex",y="PercentExpressing",
                  add='boxplot',color="Sex",fill='Sex',
                  palette =c('red','blue')) + 
        ylim(0,max(dat[,1]+5)) + 
        NoLegend() +
        ggtitle(label = paste0(levels(object)[i],': ',gene))
      if (plot == TRUE) {
      pdf(paste0(prefix,'SexProp_',assay,levels(object)[i],'_',gene,'.pdf'),width=2,height=3)
      print(p)
      dev.off() 
      }
    }
    
    tt[i] <- pval
    es[i] <- effect
    lt[i] <- ifelse(mean(fdat) > mean(mdat), -log10(pval), log10(pval))
  }
  percmat=as.data.frame(percmat)
  if (length(maleindex)>1){
    percmat$avgMale=rowMeans(percmat[,maleindex])
  } else {
    percmat$avgMale=percmat[,maleindex]
  }
  if (length(femaleindex)>1){
    percmat$avgFemale=rowMeans(percmat[,femaleindex])
  } else {
    percmat$avgFemale = percmat[,femaleindex]
  }
  percmat$SexEffect <- tt
  percmat$log10Pval <- lt
  percmat$effectSize <- es
  colnames(percmat)[1:length(unique(object$orig.ident))] <- names(table(object$orig.ident))
  rownames(percmat) <- paste0(prefix,'_',levels(object))
  #if (!is.null(lowexp)){
  #percsub=percmat[-lowexp,]
  #}
  df <- apply(percmat,2,as.character)
  rownames(df) <- rownames(percmat)
  
  write.table(df, paste0(prefix,assay,gene, 'SexProps.txt'), sep = '\t', quote = F)
  return(df)
}



#' Prop Expressing Test
#'
#' Test to determine whether assays (alignment) have different proportions of GeneX+ cells per sample by cluster. Also tests for sex effect 
#' @param object Seurat object
#' @param gene Gene name
#' @param prefix Output file name prefix
#' @param group_by metadata variable to group cells by. defaults to Idents
#' @param split_by Further subset data by metadata variable
#' @keywords PercExp
#' @export
#' @importFrom Seurat DefaultAssay
#' @examples
#' difAlignment()
difAlignment <- function(object, gene,prefix, plot = FALSE,assay1, assay2) {
  #create a dataframe of all samples + clusters
  male <- unique(object$orig.ident[which(object$Sex == 'Male')])
  maleindex <- which(names(table(object$orig.ident)) %in% male)
  female <- unique(object$orig.ident[which(object$Sex == 'Female')])
  femaleindex <- which(names(table(object$orig.ident)) %in% female)
  xx <- rep('Female',length(table(object$orig.ident)))
  xx[maleindex] <- 'Male'
  xx <- c(xx, xx)
  l1=levels(object)
  l2=names(table(object$orig.ident))
  l3 <- unlist(lapply(l1, function(x) sapply(l2, function(y) paste(x, y, sep = "_"))))
  names(l3)=NULL
  
  percs <- rep(0,length(l3))
  names(percs)=l3
  DefaultAssay(object)=assay1
  #perc1 <- Percent_Expressing(object, features = gene, split_by = 'orig.ident',group_by = 'seurat_clusters') #corrected for depth? need to sctransform on joined layers. 
  perc1 <- PercExp(object, features='Oprm1',assay=assay1,group_by='seurat_clusters',split_by='orig.ident')
  #note: not much difference in corrected vs raw counts; just use raw)
  percmat <-as.matrix(perc1)
  DefaultAssay(object)=assay2
  #perc2 <- Percent_Expressing(object, features = gene, split_by = 'orig.ident')
  perc2 <- PercExp(object,features = 'Oprm1',assay=assay2,group_by='seurat_clusters',split_by='orig.ident')
  percmat2 <- as.matrix(perc2)
  
  tt <- c()
  lt <- c()
  lowexp=c()
  es <- c()
  sex <- c()
  for (i in 1:nrow(percmat)) {
    print(i)
    d1 <- unlist(percmat[i,])
    d2 <- unlist(percmat2[i,])
    d1m <- d1[maleindex]
    d2m <- d2[maleindex]
    d1f <- d1[femaleindex]
    d2f <- d2[femaleindex]
    #must be expressed in average of 5% of cells in either male or female samples
    if(mean(d1)<10 && mean(d2)<10 ){
      lowexp=c(lowexp,i)
      pval <- 1
      effect <- 0
      sexef <- 0
    } else if (all(d1==100) && all(d2==100)){
      pval <- 1
      effect <- 0
      sexef <- 0
    }else {
      tstat <- matrixTests::row_t_equalvar(x = d2, y = d1)
      sex1 <- matrixTests::row_t_equalvar(x = d1m, y = d1f)
      psex1 <- c('Female','Male',formatC(sex1$pvalue,format = 'e',digits=2))
      names(psex1)=c('group1','group2','p')
      psex1 = as.data.frame(t(psex1))
      sex2 <- matrixTests::row_t_equalvar(x = d2m, y = d2f)
      psex2 <- c('Female','Male',formatC(sex2$pvalue,format='e',digits=2))
      names(psex2) <- c('group1','group2','p')
      psex2 <- as.data.frame(t(psex2))
      tf <-c(psex1$p < .05, psex2$p < .05)
      sexef <- 'No'
      if (length(unique(tf)) > 1){
        sexef <- 'Yes'
      }
      if(all(d1==0)) {
        effect=-mean(d2)
      } else if(all(d2==0)){
        effect=mean(d1)
      } else {
        effect <- effectsize::cohen.d(d1,d2)#negative = down in adjusted
        effect=effect$estimate
        
      }
      
      pval <- tstat$pvalue
    }
    pval <- ifelse(is.na(pval) == TRUE,1,pval)
    dat=data.frame(PercentExpressing=c(d2,d1),Alignment = c(rep(assay2,length(d2)),rep(assay1,length(d1))),Sex = xx)
    if (pval < 0.05 || abs(effect) > 1) {
      p=ggpubr::ggdotplot(dat,x="Alignment",y="PercentExpressing",
                  add='boxplot',color="Alignment",fill='Alignment',
                  palette =c('red','blue')) + 
        ylim(0,max(dat[,1]+5)) + 
        NoLegend() +
        ggtitle(label = paste0(levels(object)[i],': ',gene))
      if (plot == TRUE) {
      pdf(paste0(prefix,'AlignerProp_',levels(object)[i],'_',gene,'.pdf'),width=2,height=3)
      print(p)
      dev.off()
      }
      if (sexef == 'Yes') {
        dat2=dat[which(dat$Alignment==assay2),]
        p2=ggpubr::ggdotplot(dat2,x = 'Sex', y = 'PercentExpressing',
                     add='boxplot',color='Sex',fill='Sex',
                     palette = c('green','magenta')) +
          ylim(0,max(dat2[,1]) + (max(dat2[,1])/10)) +
          NoLegend()
        p2.1 = p2 + stat_pvalue_manual(psex2, label = 'p',y.position = max(d2)+(max(d2)/12))
        dat1 <- dat[which(dat$Alignment==assay1),]
        p1=ggpubr::ggdotplot(dat1,x = 'Sex', y = 'PercentExpressing',
                     add='boxplot',color='Sex',fill='Sex',
                     palette = c('green','magenta')) +
          ylim(0,max(dat1[,1] + (max(dat1[,1])/10))) +
          NoLegend()
        p1.1 = p1 + stat_pvalue_manual(psex1, label = 'p',y.position = max(d1)+ (max(d1)/12))
        if (plot == TRUE) {
        pdf(paste0(prefix,'AlignerProp_bySex_',levels(object)[i],'_',gene,'.pdf'),width=4.2,height=3.9)
        pt <- ggpubr::ggarrange(p2.1, p1.1, ncol=2,nrow=1,common.legend = TRUE, legend = "bottom")
        
        pt + plot_annotation(
          title = "My main title",
          theme = theme(
            plot.title = element_text(hjust = 0.5)
          )
        )
        print(annotate_figure(pt, top = text_grob(paste0(levels(object)[i],': ',gene))))
        dev.off()
        }
      }
    }
    
    
    tt[i] <- pval
    es[i] <- effect
    lt[i] <- ifelse(mean(d2) > mean(d1), -log10(pval), log10(pval))
    
  }
  
  percmat <- apply(percmat, c(1, 2), as.numeric)
  percmat2 <- apply(percmat2, c(1, 2), as.numeric)
  
  pmat <- percmat2 - percmat #differences
  pmat <- as.data.frame(pmat)
  if (length(maleindex)>1){
    pmat$avgMale=rowMeans(pmat[,maleindex])
  } else {
    pmat$avgMale=pmat[,maleindex]
  }
  if (length(femaleindex)>1){
    pmat$avgFemale=rowMeans(pmat[,femaleindex])
  } else {
    pmat$avgFemale = pmat[,femaleindex]
  }
  pmat$alignerPval <- tt
  pmat$log10Pval <- lt
  pmat$effectSize <- es
  colnames(pmat)[1:length(unique(object$orig.ident))] <- names(table(object$orig.ident))
  rownames(pmat) <- paste0(prefix,'_',levels(object))
  #if (!is.null(lowexp)){
  #percsub=percmat[-lowexp,]
  #}
  df <- apply(pmat,2,as.character)
  rownames(df) <- rownames(pmat)
  
  write.table(df, paste0(prefix, gene, 'DifferentAlignmentProps.txt'), sep = '\t', quote = F)
  dfl <- c()
  dfl[[1]] <- percmat
  dfl[[2]] <- percmat2
  dfl[[3]] <- df
  return(dfl)
  
}


###################
#incision vs morphine IEG profiles in control (SS) acute pain (IS) and acute analgesic (IM) samples
##################
#' Incision Morphine analysis
#' 
#' Function for all acute pain vs morphine data analysis and plotting? too much likely
#' @param object Seurat object
#' @param module module score metadata column name
#' @param prefix File name prefix
#' @param gene_list list of genes in module score
#' @param Sex Whether to test for sex as a variable
#' @keywords Pain Morphine Analysis
#' @export
#' @examples
#' IncMorphDif()
IncMorphDif <- function(object = seurat, module = 'signature_1DEGIEG',prefix,gene_list, Sex=FALSE) {
  pvaldf=as.data.frame(matrix(nrow = length(levels(object)),ncol =length(module)*3))
  cnames=c()
  for (i in 1:length(module)){
    j=paste0(module[i],'_SSIS')
    k=paste0(module[i],'_ISIM')
    l=paste0(module[i],'_SSIM')
    cnames[[i]]=c(j,k,l)
  }
  rownames(pvaldf) = levels(object)
  colnames(pvaldf)=unlist(cnames)
  esdf=as.data.frame(matrix(nrow = length(levels(object)),ncol =length(module)*3))
  rownames(esdf)=levels(object)
  colnames(esdf)=unlist(cnames)
  lpdf=as.data.frame(matrix(nrow = length(levels(object)), ncol = length(module) *3))
  rownames(lpdf) = levels(object)
  colnames(lpdf) = unlist(cnames)
  pm <- pvaldf
  pf <- pvaldf
  em <- esdf
  ef <- esdf
  lpm <-  lpdf
  lpf <- lpdf
  k = object@meta.data[,c(module,'Condition','Sex','seurat_clusters')]
  for (i in 1:length(levels(object))){
    j=levels(object)[i]
    print(paste0('Working on ',j))
    ksub <- k[which(k$seurat_clusters==j),]
    #k=subset(object,idents = j) #STOP SUBSETTING EACH TIME! SPLIT BY CLUSTER AND CALL EACH ELEMENT!!!
    #check if cluster is represented by at least 2 cells from all conditions
      ps=c()
      es=c()
      lp=c()
      for (ii in 1:length((module))){
        ss=ksub[,ii][which(ksub$Condition=='SS')]
        is=ksub[,ii][which(ksub$Condition=='IS')]
        im=ksub[,ii][which(ksub$Condition=='IM')]
        if (length(ss) >= 3 && length(is) >= 3) {
          ssisp=wilcox.test(is,ss,exact = F)
          ssise=rank_biserial(is,ss)
        } else {
          ssisp <- c()
          ssisp$p.value = 1
          ssise <- c()
          ssise$r_rank_biserial = 0
        }
        if (length(is) >= 3 && length(im) >= 3) {
          isimp=wilcox.test(im,is,exact = F)
          isime=rank_biserial(im,is)
        } else {
          isimp <- c()
          isimp$p.value = 1
          isime <- c()
          isime$r_rank_biserial = 0
        }
        if (length(ss) >= 3 && length(im) >= 3) {
          ssimp=wilcox.test(im,ss,exact = F)
          ssime=rank_biserial(im,ss)
        } else {
          ssimp <- c()
          ssimp$p.value = 1
          ssime <- c()
          ssime$r_rank_biserial = 0
        }
        ps[[ii]]=c(ssisp$p.value,isimp$p.value,ssimp$p.value)
        es[[ii]]=c(ssise$r_rank_biserial,isime$r_rank_biserial,ssime$r_rank_biserial)
        lp[[ii]]=ifelse(es[[ii]]< 0,log10(ps[[ii]]),-log10(ps[[ii]]))
      }
      pvaldf[i,]=unlist(ps)
      esdf[i,]=unlist(es)
      lpdf[i,]=unlist(lp)
  }
  
  boop=which(is.na(pvaldf),arr.ind = TRUE) #check for NaN
  pvaldf[boop]=1
  esdf[boop]=0
  lpdf[boop]=0
  pval_melt=melt(as.matrix(pvaldf))
  es_melt=melt(as.matrix(esdf))
  
  pssis=pvaldf[,grep('SSIS',colnames(pvaldf))]
  pisim=pvaldf[,grep('ISIM',colnames(pvaldf))]
  pssim=pvaldf[,grep('SSIM',colnames(pvaldf))]
  plist=list(pssis,pisim,pssim)
  
  essis=esdf[,grep('SSIS',colnames(esdf))]
  eisim=esdf[,grep('ISIM',colnames(esdf))]
  essim=esdf[,grep('SSIM',colnames(esdf))]
  elist=list(essis,eisim,essim)
  
  lssis=lpdf[,grep('SSIS',colnames(lpdf))]
  lisim=lpdf[,grep('ISIM',colnames(lpdf))]
  lssim=lpdf[,grep('SSIM',colnames(lpdf))]
  llist=list(lssis,lisim,lssim)
  
  
  #make a df of effect size and pvalue
  op=Percent_Expressing(object,features='Oprm1',assay='RNA',layer='counts')
  df1=data.frame(Cluster=levels(object),Pvalue=pssis,Effect.Size=essis,Oprm1Perc=unlist(op))
  df2=data.frame(Cluster=levels(object),Pvalue=pisim,Effect.Size=eisim,Oprm1Perc=unlist(op))
  df3=data.frame(Cluster=levels(object),Pvalue=pssim,Effect.Size=essim,Oprm1Perc=unlist(op))
  #indicate significant clusters w/ '*'
  df1$Cluster[which(df1$Pvalue<.05)]=paste0(df1$Cluster[which(df1$Pvalue<.05)],'*')
  df2$Cluster[which(df2$Pvalue<.05)]=paste0(df2$Cluster[which(df2$Pvalue<.05)],'*')
  df3$Cluster[which(df3$Pvalue<.05)]=paste0(df3$Cluster[which(df3$Pvalue<.05)],'*')
  #bold clusters w/ >25% oprm1 expression?
  #op=Percent_Expressing(s,features='Oprm1')
  #df1$Oprm1Perc=unlist(op[grep('_SS',names(op))])
  df1$Bold='plain'
  df1$Bold[which(df1$Oprm1Perc>50)]='bold'
  df2$Bold='plain'
  df2$Bold[which(df2$Oprm1Perc>50)]='bold'
  df3$Bold='plain'
  df3$Bold[which(df3$Oprm1Perc>50)]='bold'
  #make sure order is consistent
  levels(df1$Cluster) <- df1$Cluster
  df1$Cluster <- factor(df1$Cluster, levels = df1$Cluster)
  levels(df2$Cluster) <- df2$Cluster
  df2$Cluster <- factor(df2$Cluster, levels = df2$Cluster)
  levels(df3$Cluster) <- df3$Cluster
  df3$Cluster <- factor(df3$Cluster, levels = df3$Cluster)
  
  dfall <- data.frame(Cluster = rownames(df1),
                      SSIS_p = pssis, SSIS_l = lssis, SSIS_e = essis, 
                      ISIM_p = pisim, ISIM_l = lisim, ISIM_e = eisim,
                      SSIM_p = pssim, SSIM_l = lssim, SSIM_e = essim
  )
  dir.create('data')
  write.table(dfall,paste0('data/',prefix,'_',module,'_data.txt'),quote=F,sep='\t',col.names=NA)
  #draw bar
  #subset df by incisionmorphine only
  pain <- ggplot(data=df1, aes(x=Cluster,y=Effect.Size,fill=Pvalue)) +
    geom_bar(stat='identity',colour='black') + 
    scale_fill_gradientn(limits = c(0,.05), colours = c('red','white'),na.value = 'white')+ coord_flip() + theme_classic() + labs(title='Sham vs Incision') + theme(axis.text.y=element_text(face=df1$Bold)) + theme(axis.text.x=element_text(angle=45,hjust=1))
  
  analgesia <- ggplot(data=df2, aes(x=Cluster,y=Effect.Size,fill=Pvalue)) +
    geom_bar(stat='identity',colour='black') + 
    scale_fill_gradientn(limits = c(0,.05), colours = c('red','white'),na.value = 'white')+ coord_flip() + theme_classic() + labs(title='Incision vs Morphine') + theme(axis.text.y=element_text(face=df2$Bold)) + theme(axis.text.x=element_text(angle=45,hjust=1))
  
  opioid <- ggplot(data=df3, aes(x=Cluster,y=Effect.Size,fill=Pvalue)) +
    geom_bar(stat='identity',colour='black') + 
    scale_fill_gradientn(limits = c(0,.05), colours = c('red','white'),na.value = 'white')+ coord_flip() + theme_classic() + labs(title='Sham vs Morphine') + theme(axis.text.y=element_text(face=df3$Bold)) + theme(axis.text.x=element_text(angle=45,hjust=1))
  dname <- paste0('plots/',prefix)
  pdf(paste0(dname,'_',module,'_bars.pdf'),width=9.5, height=4)
  print(plot_grid(pain,analgesia,opioid,nrow=1))
  dev.off()
  
  pdf(paste0(dname,'_',module,'_VlnBox.pdf'),width=6,height=length(levels(s))*2/3)
  print(VlnPlot(object,features=module,pt.size=0,split.by = 'Condition',cols=c('blue','red','grey')) + coord_flip()  + ggtitle('Module score') +  geom_boxplot(position=position_dodge(1)))
  dev.off()
  
  gene_list = unlist(gene_list)
  gene_frequencies <- table(gene_list)
  sorted_genes <- sort(gene_frequencies, decreasing = TRUE)
  genes <- names(sorted_genes)[1:12]
  genes <- genes[!is.na(genes)]
  # Display the sorted list of genes
  object <- NormalizeData(object)
  pdf(paste0(dname,'_',module,'_dot.pdf'),width=6.5,height=length(levels(object))*.75)
  print(DotPlot(object,features=genes,split.by='Condition',cols=c('blue','red','grey'),assay='RNA') + 
          theme(axis.text.x=element_text(angle=45,hjust=1)))
  dev.off()
  
  if (Sex == TRUE){
    #subset by sex
  
    psm=c()
    llpm=c()
    esm=c()
    psf=c()
    llpf=c()
    esf=c()
    k = object@meta.data[,c(module,'Condition','Sex','seurat_clusters')]
    for (i in 1:length(levels(object))){
      j=levels(object)[i]
      print(paste0('Working on ',j))
      ksub=k[which(k$seurat_clusters == j),]
      mm <- ksub[which(ksub$Sex == 'Male'),]
      ff <- ksub[which(ksub$Sex == 'Female'),]
      for (ii in 1:length((module))){
      mscore <- mm[,ii]
      fscore <- ff[,ii]
      #male scores
      ss=mscore[which(mm$Condition=='SS')]
      is=mscore[which(mm$Condition=='IS')]
      im=mscore[which(mm$Condition=='IM')]
      # need at least 3 cells per group per cluster
      if (length(ss) >= 3 && length(is) >= 3) {
        ssisp=wilcox.test(is,ss,exact = F)
        ssise=rank_biserial(is,ss)
      } else {
        ssisp <- c()
        ssisp$p.value = 1
        ssise <- c()
        ssise$r_rank_biserial = 0
      }
      if (length(is) >= 3 && length(im) >= 3) {
        isimp=wilcox.test(im,is,exact = F)
        isime=rank_biserial(im,is)
      } else {
        isimp <- c()
        isimp$p.value = 1
        isime <- c()
        isime$r_rank_biserial = 0
      }
      if (length(ss) >= 3 && length(im) >= 3) {
        ssimp=wilcox.test(im,ss,exact = F)
        ssime=rank_biserial(im,ss)
      } else {
        ssimp <- c()
        ssimp$p.value = 1
        ssime <- c()
        ssime$r_rank_biserial = 0
      }
      psm[[i]]=c(ssisp$p.value,isimp$p.value,ssimp$p.value)
      esm[[i]]=c(ssise$r_rank_biserial,isime$r_rank_biserial,ssime$r_rank_biserial)
      llpm[[i]]=ifelse(esm[[i]] < 0,log10(psm[[i]]),-log10(psm[[i]]))
      #female scores
      ss=fscore[which(ff$Condition=='SS')]
      is=fscore[which(ff$Condition=='IS')]
      im=fscore[which(ff$Condition=='IM')]
      # need at least 3 cells per group per cluster
      if (length(ss) >= 3 && length(is) >= 3) {
        ssisp=wilcox.test(is,ss,exact = F)
        ssise=rank_biserial(is,ss)
      } else {
        ssisp <- c()
        ssisp$p.value = 1
        ssise <- c()
        ssise$r_rank_biserial = 0
      }
      if (length(is) >= 3 && length(im) >= 3) {
        isimp=wilcox.test(im,is,exact = F)
        isime=rank_biserial(im,is)
      } else {
        isimp <- c()
        isimp$p.value = 1
        isime <- c()
        isime$r_rank_biserial = 0
      }
      if (length(ss) >= 3 && length(im) >= 3) {
        ssimp=wilcox.test(im,ss,exact = F)
        ssime=rank_biserial(im,ss)
      } else {
        ssimp <- c()
        ssimp$p.value = 1
        ssime <- c()
        ssime$r_rank_biserial = 0
      }
      psf[[i]]=c(ssisp$p.value,isimp$p.value,ssimp$p.value)
      esf[[i]]=c(ssise$r_rank_biserial,isime$r_rank_biserial,ssime$r_rank_biserial)
      llpf[[i]]=ifelse(esf[[i]]< 0,log10(psf[[i]]),-log10(psf[[i]]))
      
      pm[i,]=unlist(psm[[i]])
      em[i,]=unlist(esm[[i]])
      pf[i,]=unlist(psf[[i]])
      ef[i,]=unlist(esf[[i]])
      lpm[i,]=unlist(llpm[[i]])
      lpf[i,]=unlist(llpf[[i]])
      }
    }
    #process male data
    boop=which(is.na(pm),arr.ind = TRUE) #check for NaN
    pm[boop]=1
    em[boop]=0
    lpm[boop]=0
    pval_melt=melt(as.matrix(pm))
    es_melt=melt(as.matrix(em))
    
    pssis=pm[,grep('SSIS',colnames(pm))]
    pisim=pm[,grep('ISIM',colnames(pm))]
    pssim=pm[,grep('SSIM',colnames(pm))]
    plist=list(pssis,pisim,pssim)
    
    essis=em[,grep('SSIS',colnames(em))]
    eisim=em[,grep('ISIM',colnames(em))]
    essim=em[,grep('SSIM',colnames(em))]
    elist=list(essis,eisim,essim)
    
    lssis=lpm[,grep('SSIS',colnames(lpm))]
    lisim=lpm[,grep('ISIM',colnames(lpm))]
    lssim=lpm[,grep('SSIM',colnames(lpm))]
    llist=list(lssis,lisim,lssim)
    
    
    #make a df of effect size and pvalue
    op=Percent_Expressing(s,features='Oprm1',assay='RNA',layer='counts')
    df1=data.frame(Cluster=levels(object),Pvalue=pssis,Effect.Size=essis,Oprm1Perc=unlist(op))
    df2=data.frame(Cluster=levels(object),Pvalue=pisim,Effect.Size=eisim,Oprm1Perc=unlist(op))
    df3=data.frame(Cluster=levels(object),Pvalue=pssim,Effect.Size=essim,Oprm1Perc=unlist(op))
    #indicate significant clusters w/ '*'
    df1$Cluster[which(df1$Pvalue<.05)]=paste0(df1$Cluster[which(df1$Pvalue<.05)],'*')
    df2$Cluster[which(df2$Pvalue<.05)]=paste0(df2$Cluster[which(df2$Pvalue<.05)],'*')
    df3$Cluster[which(df3$Pvalue<.05)]=paste0(df3$Cluster[which(df3$Pvalue<.05)],'*')
    #bold clusters w/ >25% oprm1 expression?
    #op=Percent_Expressing(s,features='Oprm1')
    #df1$Oprm1Perc=unlist(op[grep('_SS',names(op))])
    df1$Bold='plain'
    df1$Bold[which(df1$Oprm1Perc>50)]='bold'
    df2$Bold='plain'
    df2$Bold[which(df2$Oprm1Perc>50)]='bold'
    df3$Bold='plain'
    df3$Bold[which(df3$Oprm1Perc>50)]='bold'
    #make sure order is consistent
    levels(df1$Cluster) <- df1$Cluster
    df1$Cluster <- factor(df1$Cluster, levels = df1$Cluster)
    levels(df2$Cluster) <- df2$Cluster
    df2$Cluster <- factor(df2$Cluster, levels = df2$Cluster)
    levels(df3$Cluster) <- df3$Cluster
    df3$Cluster <- factor(df3$Cluster, levels = df3$Cluster)
    
    dfall <- data.frame(Cluster = rownames(df1),
                        SSIS_p = pssis, SSIS_l = lssis, SSIS_e = essis, 
                        ISIM_p = pisim, ISIM_l = lisim, ISIM_e = eisim,
                        SSIM_p = pssim, SSIM_l = lssim, SSIM_e = essim
    )
    dir.create('data')
    write.table(dfall,paste0('data/',prefix,'_',module,'_male_data.txt'),quote=F,sep='\t',col.names=NA)
    #draw bar
    #subset df by incisionmorphine only
    pain <- ggplot(data=df1, aes(x=Cluster,y=Effect.Size,fill=Pvalue)) +
      geom_bar(stat='identity',colour='black') + 
      scale_fill_gradientn(limits = c(0,.05), colours = c('red','white'),na.value = 'white')+ coord_flip() + theme_classic() + labs(title='Sham vs Incision') + theme(axis.text.y=element_text(face=df1$Bold)) + theme(axis.text.x=element_text(angle=45,hjust=1))
    
    analgesia <- ggplot(data=df2, aes(x=Cluster,y=Effect.Size,fill=Pvalue)) +
      geom_bar(stat='identity',colour='black') + 
      scale_fill_gradientn(limits = c(0,.05), colours = c('red','white'),na.value = 'white')+ coord_flip() + theme_classic() + labs(title='Incision vs Morphine') + theme(axis.text.y=element_text(face=df2$Bold)) + theme(axis.text.x=element_text(angle=45,hjust=1))
    
    opioid <- ggplot(data=df3, aes(x=Cluster,y=Effect.Size,fill=Pvalue)) +
      geom_bar(stat='identity',colour='black') + 
      scale_fill_gradientn(limits = c(0,.05), colours = c('red','white'),na.value = 'white')+ coord_flip() + theme_classic() + labs(title='Sham vs Morphine') + theme(axis.text.y=element_text(face=df3$Bold)) + theme(axis.text.x=element_text(angle=45,hjust=1))
    dname <- paste0('plots/',prefix)
    pdf(paste0(dname,'_',module,'_male_bars.pdf'),width=9.5, height=4)
    print(plot_grid(pain,analgesia,opioid,nrow=1))
    dev.off()
    
    mm = subset(object,subset=Sex=='Male')
    mm = NormalizeData(mm)
    
    pdf(paste0(dname,'_',module,'_male_VlnBox.pdf'),width=6,height=length(levels(s))*2/3)
    print(VlnPlot(mm,features=module,pt.size=0,split.by = 'Condition',cols=c('blue','red','grey')) + coord_flip()  + ggtitle('Module score') +  geom_boxplot(position=position_dodge(1)))
    dev.off()
    
    gene_list = unlist(gene_list)
    gene_frequencies <- table(gene_list)
    sorted_genes <- sort(gene_frequencies, decreasing = TRUE)
    genes <- names(sorted_genes)[1:12]
    genes <- genes[!is.na(genes)]
    # Display the sorted list of genes
    pdf(paste0(dname,'_',module,'_male_dot.pdf'),width=6.5,height=length(levels(object))*.75)
    print(DotPlot(mm,features=genes,split.by='Condition',cols=c('blue','red','grey'),assay='RNA') + 
            theme(axis.text.x=element_text(angle=45,hjust=1)))
    dev.off()
    
    #process female data
    boop=which(is.na(pf),arr.ind = TRUE) #check for NaN
    pm[boop]=1
    em[boop]=0
    lpm[boop]=0
    pval_melt=melt(as.matrix(pf))
    es_melt=melt(as.matrix(ef))
    
    pssis=pf[,grep('SSIS',colnames(pf))]
    pisim=pf[,grep('ISIM',colnames(pf))]
    pssim=pf[,grep('SSIM',colnames(pf))]
    plist=list(pssis,pisim,pssim)
    
    essis=ef[,grep('SSIS',colnames(ef))]
    eisim=ef[,grep('ISIM',colnames(ef))]
    essim=ef[,grep('SSIM',colnames(ef))]
    elist=list(essis,eisim,essim)
    
    lssis=lpf[,grep('SSIS',colnames(lpf))]
    lisim=lpf[,grep('ISIM',colnames(lpf))]
    lssim=lpf[,grep('SSIM',colnames(lpf))]
    llist=list(lssis,lisim,lssim)
    
    
    #make a df of effect size and pvalue
    op=Percent_Expressing(s,features='Oprm1',assay='RNA',layer='counts')
    df1=data.frame(Cluster=levels(object),Pvalue=pssis,Effect.Size=essis,Oprm1Perc=unlist(op))
    df2=data.frame(Cluster=levels(object),Pvalue=pisim,Effect.Size=eisim,Oprm1Perc=unlist(op))
    df3=data.frame(Cluster=levels(object),Pvalue=pssim,Effect.Size=essim,Oprm1Perc=unlist(op))
    #indicate significant clusters w/ '*'
    df1$Cluster[which(df1$Pvalue<.05)]=paste0(df1$Cluster[which(df1$Pvalue<.05)],'*')
    df2$Cluster[which(df2$Pvalue<.05)]=paste0(df2$Cluster[which(df2$Pvalue<.05)],'*')
    df3$Cluster[which(df3$Pvalue<.05)]=paste0(df3$Cluster[which(df3$Pvalue<.05)],'*')
    #bold clusters w/ >25% oprm1 expression?
    #op=Percent_Expressing(s,features='Oprm1')
    #df1$Oprm1Perc=unlist(op[grep('_SS',names(op))])
    df1$Bold='plain'
    df1$Bold[which(df1$Oprm1Perc>50)]='bold'
    df2$Bold='plain'
    df2$Bold[which(df2$Oprm1Perc>50)]='bold'
    df3$Bold='plain'
    df3$Bold[which(df3$Oprm1Perc>50)]='bold'
    #make sure order is consistent
    levels(df1$Cluster) <- df1$Cluster
    df1$Cluster <- factor(df1$Cluster, levels = df1$Cluster)
    levels(df2$Cluster) <- df2$Cluster
    df2$Cluster <- factor(df2$Cluster, levels = df2$Cluster)
    levels(df3$Cluster) <- df3$Cluster
    df3$Cluster <- factor(df3$Cluster, levels = df3$Cluster)
    
    dfall <- data.frame(Cluster = rownames(df1),
                        SSIS_p = pssis, SSIS_l = lssis, SSIS_e = essis, 
                        ISIM_p = pisim, ISIM_l = lisim, ISIM_e = eisim,
                        SSIM_p = pssim, SSIM_l = lssim, SSIM_e = essim
    )
    dir.create('data')
    write.table(dfall,paste0('data/',prefix,'_',module,'_female_data.txt'),quote=F,sep='\t',col.names=NA)
    #draw bar
    #subset df by incisionmorphine only
    pain <- ggplot(data=df1, aes(x=Cluster,y=Effect.Size,fill=Pvalue)) +
      geom_bar(stat='identity',colour='black') + 
      scale_fill_gradientn(limits = c(0,.05), colours = c('red','white'),na.value = 'white')+ coord_flip() + theme_classic() + labs(title='Sham vs Incision') + theme(axis.text.y=element_text(face=df1$Bold)) + theme(axis.text.x=element_text(angle=45,hjust=1))
    
    analgesia <- ggplot(data=df2, aes(x=Cluster,y=Effect.Size,fill=Pvalue)) +
      geom_bar(stat='identity',colour='black') + 
      scale_fill_gradientn(limits = c(0,.05), colours = c('red','white'),na.value = 'white')+ coord_flip() + theme_classic() + labs(title='Incision vs Morphine') + theme(axis.text.y=element_text(face=df2$Bold)) + theme(axis.text.x=element_text(angle=45,hjust=1))
    
    opioid <- ggplot(data=df3, aes(x=Cluster,y=Effect.Size,fill=Pvalue)) +
      geom_bar(stat='identity',colour='black') + 
      scale_fill_gradientn(limits = c(0,.05), colours = c('red','white'),na.value = 'white')+ coord_flip() + theme_classic() + labs(title='Sham vs Morphine') + theme(axis.text.y=element_text(face=df3$Bold)) + theme(axis.text.x=element_text(angle=45,hjust=1))
    dname <- paste0('plots/',prefix)
    pdf(paste0(dname,'_',module,'_female_bars.pdf'),width=9.5, height=4)
    print(plot_grid(pain,analgesia,opioid,nrow=1))
    dev.off()
    ff=subset(object,subset=Sex=='Female')
    ff=NormalizeData(ff)
    pdf(paste0(dname,'_',module,'_female_VlnBox.pdf'),width=6,height=length(levels(s))*2/3)
    print(VlnPlot(ff,features=module,pt.size=0,split.by = 'Condition',cols=c('blue','red','grey')) + coord_flip()  + ggtitle('Module score') +  geom_boxplot(position=position_dodge(1)))
    dev.off()
    
    gene_list = unlist(gene_list)
    gene_frequencies <- table(gene_list)
    sorted_genes <- sort(gene_frequencies, decreasing = TRUE)
    genes <- names(sorted_genes)[1:12]
    genes <- genes[!is.na(genes)]
    # Display the sorted list of genes
    pdf(paste0(dname,'_',module,'_female_dot.pdf'),width=6.5,height=length(levels(object))*.75)
    print(DotPlot(ff,features=genes,split.by='Condition',cols=c('blue','red','grey'),assay='RNA') + 
            theme(axis.text.x=element_text(angle=45,hjust=1)))
    dev.off()
  }
  
  gc()
}

