##############################
###Carey MYC activity score###
##############################

MYC.Clusters <- function(i,x,y,c1){
  #i: Matrix with transformed gene expression values (z-score), genes in lines and samples in columns.
  #x: Named Factor. Class for high and low MYC IHC expression for training dataset (e.g. Carey dataset).
  #y: Named Factor. Class for high and low MYC IHC expression for test dataset (e.g. HOVON-84).
  #c1: Character vector. Name of the class with high MYC IHC expression. (e.g. "HIGH")

  #Define a model with the training dataset
  elastic_net_model <- train(x= t(i),
                             y = x,
                             method = "glmnet",tuneGrid = newgrid, trControl = train_control, standardize=T)
  
  #extract stats
  mypred <- extractProb(list(elastic_net_model)) 
  rownames(mypred) <- colnames(i)
  mypred <- mypred[order(mypred[,1]),]
  summary.stat.train <- confusionMatrix(data=mypred[,4], reference=mypred[,3] , positive=c1 )$byClass
  
  #extract gene importance
  mypredImp <- varImp( elastic_net_model, scale = T)
  gene.imp <- rownames(subset(mypredImp$importance, mypredImp$importance > 0 ))
  
  #predict on test dataset the myc activity score
  pred_elasticnet <- predict(elastic_net_model, newdata= t(i)[,gene.imp]  )

  #extract stats
  mypred2 <- extractProb(list(elastic_net_model), unkX=t(normalized.data.3.alt.scale)[,gene.imp] ) 
  rownames(mypred2) <- colnames(i)
  mypred2 <- mypred2[order(mypred2[,1]),]
  summary.stat.test <-  confusionMatrix(data=pred_elasticnet2, reference=y , positive=c1 )$byClass
  
  #Final MYC activity score
  test.class <- y[rownames(mypred2)]
  
  output <- list(summary.stat.train,summary.stat.test,gene.imp,test.class)
  
  return(output)
  
}

#################################
###Monti consensus clustering###
#################################

Monti.Clusters <- function(i,x,y){
  #i: Matrix with transformed gene expression values (z-score), genes in lines and samples in columns.
  #x: Numerical value for dimension x of your expression matrix. Used 10. For more info see: kohonen::somgrid
  #y: Numerical value for dimension y of your expression matrix. Used 14. For more info see: kohonen::somgrid
  
  #Apply the three algorithms (CC, HC and PC)
  EMclustering <- Mclust(t(i), G = 2:9 ) 
  HCclustering <- ConsensusClusterPlus(as.matrix(i),clusterAlg="hc",distance="euclidean",maxK=9,pItem=0.8,reps=200,plot="png",title="./CC_HC_clustering")
  mysom <- function(this_dist,k){
    sommap = kohonen::som(t(as.matrix(this_dist)), grid = somgrid( x,y, topo= "hexagonal", toroidal=T))
    assignment = cutree(hclust(dist(sommap$codes[[1]],method = "euclidean"),method="average"), k)
    return(assignment)  
  }
  SOMclustering <- ConsensusClusterPlus(as.matrix(i),clusterAlg="mysom",distance="euclidean",maxK=9,pItem=0.8,reps=200,plot="png",title="./CC_SOM_clustering")

  #subset samples from same cluster for the three algorithms
  ccomb_class <- data.frame(HC=HCclustering[[2]][["consensusClass"]], SOM=SOMclustering[[2]][["consensusClass"]] , PC=EMclustering$classification)
  ccomb <- ccomb_class[which(ccomb_class[,1]==ccomb_class[,3] & ccomb_class[,2]==ccomb_class[,3] ),]
  
  #define classes for unclassified samples
  model = train(t(i[,rownames(ccomb)]), factor(ccomb[,1]) ,'nb',trControl=trainControl(method='cv',number=10))
  non.ccomb.pred <- predict(model$finalModel, t(i[, -which(colnames(i) %in% rownames(ccomb)) ]))$class
  
  #final classes
  monti.classes <- i[complete.cases(i), c(rownames(ccomb),names(non.ccomb.pred)) ]
  return(monti.classes)
  
}
