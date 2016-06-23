#############################################################
PlacePt <- function(x, y, verts, norms, mesh, start){

	temp <- rgl.user2window(x = verts[, 1], y = verts[, 2],
				z = verts[, 3], projection = start$projection)
  	# conversion window coordinates in 2D screen coordinates
  	X <- temp[, 1] * start$viewport[3]
  	Y <- (1 - temp[, 2]) * start$viewport[4]

  	# 3D window coordinates of cliked point
  	X1 <- x / start$viewport[3]
  	Y1 <- 1 - y / start$viewport[4]

  # d?termination des vertices qui se situent dans un carr? de s?lection
  # mX, mY : matrices contenant les x et y (coordonn?es 2D ?cran) des 3 sommets de chaque face retenue
  mX <- matrix(X[mesh$it], nrow = 3, byrow = FALSE)
  mY <- matrix(Y[mesh$it], nrow = 3, byrow = FALSE)
  # calcul des longueurs de chaque arr?te
  d <- mX
  d[1,] <- (mX[1,] - mX[2,])^2 + (mY[1,] - mY[2,])^2
  d[2,] <- (mX[1,] - mX[3,])^2 + (mY[1,] - mY[3,])^2
  d[3,] <- (mX[2,] - mX[3,])^2 + (mY[2,] - mY[3,])^2
  dm <- sqrt(max(d))
  # indices des vertices dans carr? de s?lection
  sqIdx <- (X >= (x-dm)) & (X <= (x+dm)) & (Y >= (y-dm)) & (Y <= (y+dm))

  # extraction du sous mesh dans carr? de s?lection (sans utiliser ma fonction submesh, mais avec rmVertex de Morpho)
  #selecMesh <- rmVertex(mesh, which(sqIdx), keep=TRUE)
  selecMesh <- subset.mesh(mesh, sqIdx)

  # d?finition du point cliqu? comme un mesh3d
  Q <- rgl.window2user(X1, Y1, 0)
  normQ <- rgl.window2user(X1, Y1, 1) - Q
  normQ <- normQ/sqrt(sum(normQ^2))
  lQ <- list(vb = matrix(c(Q, 1), 4, 1), normals = matrix(c(normQ, 1), 4, 1))
  class(lQ) <- "mesh3d"

   # Search for the intersection
   int <- vcgRaySearch(lQ, mesh)
   cas <- 1
   # if there is no, then we need to search the nearest point
   # Use the angle to filter out possible points
   # However this return
   if (int$quality == 0){
        temp2 <- rgl.user2window(x = verts[, 1] + norms[, 1],
     				y = verts[,2] + norms[, 2],
     				z = verts[, 3] + norms[, 3],
     				projection = start$projection)
        normals <- temp2 - temp

        u <- par3d()$observer
     alpha <- acos((t(u) %*% t(normals)) / (sqrt(rowSums(normals^2)) * sqrt(sum(u^2))))

     Idx <- alpha > pi/2
     # sometimes the angle failed then took all vertices
     if (sum(is.na(Idx))>0 || length(Idx)==0) {
         Idx <- rep(TRUE, length(X))
     }
     Xs <- X[Idx]
     Ys <- Y[Idx]

     Dist <- sqrt((Xs - x)^2 + (Ys - y)^2)
      idx <- which.min(Dist)
     int <- subset.mesh(mesh, which(Idx)==idx, select = "vb")
     str(int)
     if (length(int$vb)==0) {
         int$vb <- c(0,0,0,1) # put the center as closest point if the submesh is empty...
     }
     cas <- 2
   }

  if (length(selecMesh$vb) == 0 | !is.matrix(selecMesh$it)) {
    # clic dans une zone de vide ou une seule face triangulaire (fait planter vcgIsolated) => zone de zoom pas d?finissable
    visibles <- idx <-NULL
  } else {
    isol <- vcgIsolated(selecMesh, split=TRUE)
    vd <- matrix(0, length(isol), 2)
    for (i in 1:length(isol)){
      temp <- sqrt(colMeans((isol[[i]]$vb - matrix(int$vb, 4, dim(isol[[i]]$vb)[2]))^2))
      vd[i, 1] <- min(temp)
      vd[i, 2] <- which.min(temp)
    }
    idx <- which(vd[, 1] == min(vd[, 1]))
    visibles <- cbind(isol[[idx]]$vb[1,], isol[[idx]]$vb[2,], isol[[idx]]$vb[3,])
    idx <- vd[idx, 2]
  }
  return(list(visibles=visibles, idx=idx))
}


#############################################################
rgl.select2<-function (button = c("left", "middle", "right"), verts, norms, mesh,
                       modify=FALSE, A=NULL, IdxPts=NULL, grDev) {

  # Fonction pour placer un point ? l'aide la souris
  #
  # D?finition de 3 sous-fonctions Begin(), Update() et End() qui d?terminent ce qu'il faut faire quand l'utilisateur commence ? appuyer
  # sur le bouton gauche de la souris (Begin()), se d?place avec la souris tout en maintenant le bouton appuy? (Update()), et le rel?che (End())

  start <- list()
  Sp<-Tx<-idx<-Idx<-NULL
  firstTime<-TRUE

  # Fonction de d?but de clic
  Begin <- function(x, y) {

    # r?cup?ration d'infos sur la projection actuelle du mesh (variable globale pour utilisation dans Update())
    start$viewport <<- par3d("viewport")
    start$projection <<- rgl.projection()



    # d?termination du point du mesh le plus proche du clic de l'utilisateur (x,y)
      temp<-PlacePt(x,y,verts,norms,mesh,start)

      visibles<<-temp$visibles
      idx<<-temp$idx



    if (modify){
      # quand l'utilisateur en est ? l'?tape o? il peut modifier les points d?j? plac?s :
      # il faut d?terminer quel point existant il faut modifier

      # calcul rapide sur la distance : que le point soit visible ou non => l'utilisateur doit cliquer ? proximit? imm?diate du point ? d?placer pour ?viter qu'il y ait des probl?mes...
      distPt<-sqrt(rowSums((A-matrix(visibles[idx,],dim(A)[1],dim(A)[2],byrow=TRUE))^2))
      Idx<<-which(distPt==min(distPt))
      IdxPts<<-Idx
      if (firstTime){
        # effacement du point/label ? modifier : ? faire seulement une fois...
        rgl.pop("shapes",grDev$vSp[Idx])
        rgl.pop("shapes",grDev$vTx[Idx])
        firstTime<<-FALSE
      }else{
        # on efface le point/label plac? lors du pr?c?dent clic non valid? par ?chappe
        rgl.pop("shapes",Sp)
        rgl.pop("shapes",Tx)
      }
    }else{
      # on efface le point/label plac? lors du pr?c?dent clic non valid? par ?chappe
      if (!is.null(Sp)){
        rgl.pop("shapes",Sp)
        rgl.pop("shapes",Tx)
      }
      Idx<<-NULL
    }

    # plot du point/label
    Sp<<-spheres3d(x=visibles[idx,1],y=visibles[idx,2],z=visibles[idx,3],alpha=0.5,radius=grDev$spradius)
    Tx<<-text3d(x=visibles[idx,1],y=visibles[idx,2],z=visibles[idx,3],texts=as.character(IdxPts),cex=grDev$tcex)
  }

  Update <- function(x, y) {

   # d?termination du point du mesh le plus proche du clic de l'utilisateur (x,y)
    temp<-PlacePt(x,y,verts,norms,mesh,start)
    visibles<<-temp$visibles
    idx<<-temp$idx

    # on efface le point/label plac? lors du pr?c?dent clic non valid? par ?chappe
    if (!is.null(Sp)){
      rgl.pop("shapes",Sp)
      rgl.pop("shapes",Tx)
    }

    # plot du point/label
    Sp<<-spheres3d(x=visibles[idx,1],y=visibles[idx,2],z=visibles[idx,3],alpha=0.5,radius=grDev$spradius)
    Tx<<-text3d(x=visibles[idx,1],y=visibles[idx,2],z=visibles[idx,3],texts=as.character(IdxPts),cex=grDev$tcex)
  }

  End <- function(x,y){
    # bouton de la souris rel?ch?...
  }

  # lignes suivantes : reprises de rgl.select()
  button <- match.arg(button)
  newhandler <- par3d("mouseMode")
  newhandler[button] <- "selecting"
  oldhandler <- par3d(mouseMode = newhandler)
  on.exit(par3d(mouseMode = oldhandler))

  # modification de l'action d'interface quand l'utilisateur clique sur le bouton droit de la souris :
  # avant zoom, maintenant pla?age d'un point par magn?tisme (utilisation des sous-fonctions Begin, Update et End)
  rMul<-rgl.setMouseCallbacks(2, Begin, Update, End)

  # il faut maintenant suspendre l'?xecution du code et maintenir cette d?finition du clic droit tant que l'utilisateur n'a pas appuy? sur la touche ?chappe
  dev<-rgl.cur()
  while (dev==rgl.cur()){
    if (!is.null(idx)){ # v?rification que clic pas dans le vide ou que trop peu de faces triangulaires dans le sous-mesh
      # r?cup?ration soit des coordonn?es de position du pointeur soit le cas ?ch?ant de l'action de la touche ?chappe
      result <- rgl:::rgl.selectstate()
      # si ?chappe : on sort d'ici
      if (result$state >= rgl:::msDONE)
      {break}
    }
  }

  # si la fen?tre a ?t? ferm?e (lorsqu'on on souhaite arr?t? de modifier des points) : on sort de cette fonction
  if (dev!=rgl.cur())
  {
    if (modify){
      isDone<-TRUE
      return(list(isDone=isDone, isClosed=TRUE))
    }
  }

  # Sinon
  # on red?finit l'action du clic droit par l'action par d?faut comme ?tant le zoom
  rgl.setMouseCallbacks(2)
  # ligne de rgl.select :
  rgl:::rgl.setselectstate("none")

  # Exports
  isDone<-TRUE
  if (result$state == rgl:::msDONE) isDone<-FALSE

    return(list(isDone=isDone, coords=visibles[idx,], sp=Sp, Idx=Idx, isClosed=FALSE, tx=Tx))

}


SelectPoints3d<-function (mesh, modify=FALSE, A=NULL, IdxPts=NULL, grDev) {
    # Do the transpose only once
    verts <- t(mesh$vb[1:3, ])
    norms <- t(mesh$normals[1:3, ])
    StopPts <- 0
    # Stop when landmark is validated by ESC or when the window is closed
    while (StopPts==0) {
        temp <- rgl.select2(button="right",verts = verts, norms = norms, mesh = mesh,
                            modify = modify, A = A, IdxPts = IdxPts, grDev = grDev)
        if (temp$isClosed | temp$isDone){
          if (temp$isClosed){
            # Need to close twice becaus of rgl...
            rgl.close()
          }
          break
        }
    }
    return(temp)
}



#############################################################
SetPtZoom <- function(dd, specFull, Trans, Pt, IdxPts=NULL, orthoplanes,
                      percDist=0.15, modify=FALSE, A=NULL, grDev) {

  # Fonction pour effectuer "un zoom" autour du point s?lectionn? sur le mesh d?cim? :
  # ouvre une 2?me fen?tre graphique avec la partie du mesh non d?cim? correspondante et pla?age du point

  # les points ? conserver se situent ? moins de 10% de la distance max observ?e entre le point s?lectionn? et l'ensmble des autres vertices du mesh
  keep <- dd < (percDist * max(dd)) # distances inf?rieures ? 10% distance max

  # extraction du sous-mesh (Use rmVertex instead ? actually subset.mesh performs quicker)
  specFull2 <- subset.mesh(specFull, keep)

  # si le sous-mesh contient plusieurs meshs isol?s : on ne conserve que le sous-mesh ? proximit? imm?diate du point cliqu?
  temp <- vcgIsolated(specFull2, split = TRUE)
  Min <- +Inf
  for (ii in 1:length(temp)) {
    vb <- temp[[ii]]$vb[1:3, ]
    cs <- sqrt(colSums((vb - (Pt + Trans))^2)) #!problem
    MinT <- cs[which(cs==min(cs))]
    if (MinT < Min){
      specFull2 <- temp[[ii]]
      Min <- MinT
    }
  }

  # center the vertices of the submesh
  Trans2 <- apply(specFull2$vb[1:3, ], 1, mean)
  specFull2$vb[1:3,] <- specFull2$vb[1:3, ] - Trans2
  # plot
  param3d <- par3d()
  d2 <- open3d()
  par3d(windowRect = grDev$windowRect[2, ])
  ids2 <- plot3d(specFull2$vb[1, ], specFull2$vb[2, ], specFull2$vb[3, ],
                 size = grDev$ptSize, aspect = FALSE,
                 axes = F, box = F, xlab="", ylab="", zlab="", main = paste("Land - ", IdxPts))
  if (is.null(specFull2$material))
      specFull2$material <- "gray"
  shade3d(specFull2)

  # draws eventually the intersections visibles on the submesh
  if (!is.null(orthoplanes$vInter)){
    for (i in 1:3){
      # Only intersection points closed to the intersection planes
      ddi <- sqrt(colSums((t(orthoplanes$vInter[[i]])-Trans-Pt)^2))
      keep <- ddi < (percDist*max(ddi))
      if (sum(keep)> 0){
        inter <- sweep(orthoplanes$vInter[[i]][keep, ], 2, Trans2)
        lines3d(inter, col="red", lwd=2)
      }
    }
  }

  # Adjust the orientation of the zoomed mesh to correspond to the one of the decim mesh
  rgl.viewpoint(userMatrix = param3d$userMatrix)

  # Add the point on the zoomed mesh
  if (is.null(grDev$spradius)) {
     tmp <- apply(specFull2$vb[1:3,],1,range)
     tmp <- tmp[2,]-tmp[1,]
     grDev$spradius <- (1/50)*min(tmp)
  }
  res2 <- SelectPoints3d(specFull2, modify, A, IdxPts, grDev)

  rgl.close()
  return(list(coords = res2$coords, sp = res2$sp, Trans2 = Trans2, tx = res2$tx))
}




