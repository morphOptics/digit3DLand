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
  selecMesh <- rmVertex(mesh, which(sqIdx), keep=TRUE)

  # d?finition du point cliqu? comme un mesh3d
  Q <- rgl.window2user(X1, Y1, 0)
  normQ <- rgl.window2user(X1, Y1, 1) - Q
  normQ <- normQ/sqrt(sum(normQ^2))
  lQ <- list(vb = matrix(c(Q, 1), 4, 1), normals = matrix(c(normQ, 1), 4, 1))
  class(lQ) <- "mesh3d"

   # recherche de l'intersection
   int <- vcgRaySearch(lQ, mesh)

#   spheres3d(int$vb[1],int$vb[2],int$vb[3],col="cyan",radius=0.2,alpha=0.1)
   # si rayon n'intersecte pas le mesh : vcgRaySearch ne retourne pas un r?sultat appropri?
   # => n?cessit? de rechercher le point le plus proche
   cas <- 1
   if (int$quality == 0){
     # normales du mesh (coordonn?es fen?tre 3D)
     temp2 <- rgl.user2window(x = verts[, 1]+norms[, 1],
     				y = verts[,2]+norms[, 2],
     				z = verts[, 3] + norms[, 3],
     				projection = start$projection)
     normals <- temp2 - temp

     u <- par3d()$observer
     alpha <- acos((t(u) %*% t(normals)) / (sqrt(rowSums(normals^2)) * sqrt(sum(u^2))))
     Idx <- alpha > pi/2

     Xs <- X[Idx]
     Ys <- Y[Idx]

     Dist <- sqrt((Xs - x)^2 + (Ys - y)^2)
     idx <- which.min(Dist)
     int <- rmVertex(mesh, which(Idx)[idx], keep=TRUE)
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


#############################################################
selectPlanes<-function (button = c("left", "middle", "right"), mesh, planes,
                        Pl=NULL, pl=NULL, dev = rgl.cur(), uplanes=NULL) {

# R?cup?ration de la fonction mouseTrackball() de rgl (via les binary) utilis?e pour effectuer la rotation manuel du mesh
# compl?t?e ici pour qu'ne fonction de cette rotation (via le bouton droit de la souris), les plans orthogonaux restent "fixes" par rapport ? cette rotation.
# En r?alit?, c'est bien le mesh qui reste fixe en coordonn?es absolues (car avec le mesh ,les axes tournent), et ce sont des nouveaux plans qui sont calcul?s ? chaque fois.

    # initialisations r?cup?r?es de mouseTrackball()
    width <- height <- rotBase <- NULL
    userMatrix <- list()
    cur <- rgl.cur()

    vInter<-NULL
    vPlanes<-NULL

    # fonctions utilitaires r?cup?r?es de mouseTrackball()
    screenToVector <- function(x, y) {
      # doit permettre ? partir d'un point 2D ?cran de calculer un vecteur
      # entre 2 positions successives de la souris, fonction appel?e 2 fois, et les deux vecteurs permettent de d?terminer un angle et un axe de rotation 3D
      radius <- max(width, height)/2
      centre <- c(width, height)/2
      pt <- (c(x, y) - centre)/radius
      len <- vlen(pt)

      if (len > 1.e-6) pt <- pt/len

      maxlen <- sqrt(2)
      angle <- (maxlen - len)/maxlen*pi/2
      z <- sin(angle)
      len <- sqrt(1 - z^2)
      pt <- pt * len
      return (c(pt, z))
    }

    # fonction de d?but de clic (droit)
    trackballBegin <- function(x, y) {

        # lignes inchang?es de mouseTrackball()
        vp <- par3d("viewport")
        width <<- vp[3]
        height <<- vp[4]
        cur <<- rgl.cur()
        for (i in dev) {
            if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
            else userMatrix[[i]] <<- par3d("userMatrix")
        }
        rgl.set(cur, TRUE)
        rotBase <<- screenToVector(x, height - y)

        # initialisations rajout?es pour les plans orthogonaux
        if (is.null(uplanes)){
          # premi?re rotation effectu? : on initialise uplanes (updated planes) ? planes initial (array contenant pour chaque plan des points contenus dans le plan le long de 2 vecteurs orthogonaux)
          uplanes<<-planes
          # vInter et vPlanes contiendront les donn?es n?cessaires pour dessiner les plans et intersections
          vInter<<-list()
          vPlanes<<-matrix(NA,3,4)
        }else{
          # ni?me rotation effectu? : planes est initialis? ? uplanes (derni?re rotation effectu?e)
          planes<<-uplanes
        }
    }

    # fonction de maintien du clic (droit)
    trackballUpdate <- function(x,y) {

        # lignes inchang?es de mouseTrackball()
        rotCurrent <- screenToVector(x, height - y)
        angle <- angle(rotBase, rotCurrent) # d?termination de l'angle...
        axis <- xprod(rotBase, rotCurrent)  # ...et de l'axe de rotation

        mouseMatrix <- rotationMatrix(angle, axis[1], axis[2], axis[3]) # calcul de la matrice de rotation correspondante
        for (i in dev) {
            if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
            else par3d(userMatrix = mouseMatrix %*% userMatrix[[i]]) # en fait ce qui tourne : c'est le point de vue utilisateur (userMatrix) par la matrice de rotation pr?c?dente
            # => les coordonn?es r?elles du mesh restent inchang?es car les axes tournent en m?me temps
        }
        rgl.set(cur, TRUE)

        # lignes rajout?es

        # mouseMatrix affectant la userMatrix, je suppose donc que les coordonn?es d'axe ne sont pas en coordonn?es absolues du mesh, mais ont d?j? subit une partie du traitement des coordonn?es 3D (cf aide de par3d())
        # Plus pr?cis?ment, on doit en ?tre ? l'?tape o? les coordonn?es sont multipli?es par userMatrix
        # Si on reprend les ?tapes :
        # - coordonn?es absolues mesh : vector (x,y,z)
        # - v = (x,y,z,1)
        # - v'=v*scale : scale d?pendant de scaling des axes, pas d'importance ici...
        # - v''= userMatrix*v*scale
        # donc :
        # axis = userMatrix%*%axis_abs_coord
        # inv(userMatrix)%*%axis = axis_abs_coord
        axis<-c(axis,1)
        axis<-solve(par3d("userMatrix"))%*%axis
        axis<-axis[1:3]
        # d?termination de la matrice de rotation ? appliquer aux coordonn?es absolues (axis modifi?, angle restant le m?me)
        uMat <- rotationMatrix(angle, axis[1], axis[2], axis[3])

        # boucle pour d?terminer et dessiner les nouvelles intersections
        for (i in 1:3){

          # plan consid?r? : plan PC2-PC3 si i=1, PC1-PC3 si i=2, PC1-PC2 si i=3
          ci<-setdiff(1:3,i)

          # extraction des valeurs du plan consid?r?
          A<-p1<-planes[1,,i]
          p2<-planes[2,,i]
          p3<-planes[3,,i]

          # calcul des vecteurs du plan
          p2<-p2-p1
          p3<-p3-p1

          # rotation par uMat
          p2<-c(t(p2)%*%uMat[1:3,1:3])+p1
          p3<-c(t(p3)%*%uMat[1:3,1:3])+p1

          # stcokage pour utilisation ult?rieure d'un nouvel appel de la fonction
          uplanes[,,i]<<-rbind(p1,p2,p3)

          # intersection plan/mesh
          inter<- meshPlaneIntersect(mesh,p1,p2,p3)

          # les points fournis par meshPlane Intersect ne sont pas tri?s (pas g?nant quand on utilise points3d(), mais ?a l'est avec lines3d()...)
          # => on trie les points en fonction de l'angle qu'ils forment dans le plan de l'intersection
          # centrage des points
          int<-inter-matrix(A,nrow=dim(inter)[1],ncol=3,byrow=TRUE)
          # rotation pour que les points d'intersection soient dans un plan
          int2<-(int%*%svd(t(int)%*%int)$u)[,1:2]
          # d?termination du rang de l'angle form? par chaque vecteur point d'intersection (via les complexes)
          idx<-order(Arg(int2[,1]+1i*int2[,2]))
          # tri des points
          inter<-inter[idx,]

          # effacement des plans et intersections pr?c?dents
          rgl.pop("shapes",Pl[[i]])
          rgl.pop("shapes",pl[[i]])

          # plot des nouvelles intersections et stockage
          Pl[[i]]<<-lines3d(inter,col="red",lwd=2)
          vInter[[i]]<<-inter

        }

        # boucle pour dessiner les plans
        for (i in 1:3){
          A <- uplanes[1,,1]
          if (i==1){
            # vect norm ? plan 2-3
            n <- uplanes[2,,3] - uplanes[1,,3]
          } else {
            n <- uplanes[i,,1] - A
          }

          # calcul des param?tres du plan
          d<- -t(n)%*%A
          # plot des nouveaux plans et stockage
          pl[[i]]<<-planes3d(n[1],n[2],n[3],d,alpha=0.7,col="cyan")
          vPlanes[i,]<<-c(n[1:3],d)
       }
    }

  # lignes suivantes : reprises de rgl.select()
  button <- match.arg(button)
  newhandler <- par3d("mouseMode")
  newhandler[button] <- "selecting"
  oldhandler <- par3d(mouseMode = newhandler)
  on.exit(par3d(mouseMode = oldhandler))

  # modification de l'action d'interface quand l'utilisateur clique sur le bouton droit de la souris :
  # avant zoom, maintenant pla?age d'un point par magn?tisme (utilisation des sous-fonctions Begin, Update et End)
  rMul<-rgl.setMouseCallbacks(2, begin = trackballBegin, update = trackballUpdate, end = NULL)

  # il faut maintenant suspendre l'?xecution du code et maintenir cette d?finition du clic droit tant que l'utilisateur n'a pas appuy? sur la touche ?chappe
  dev<-rgl.cur()

  while (dev==rgl.cur()){
      result <- rgl:::rgl.selectstate()
      # if ESC -> get-out
      if (result$state >= rgl:::msDONE) break
  }

  # si la fen?tre a ?t? ferm?e (?a ne devrait pas arriver) : on sort de cette fonction
  if (dev!=rgl.cur())
  {
    if (modify){
      isDone<-TRUE
      return(list(isDone=isDone,isClosed=TRUE))
    }
  }

  # Sinon
  # on red?finit l'action du clic droit par l'action par d?faut comme ?tant le zoom
  rgl.setMouseCallbacks(2)
  # ligne de rgl.select :
  rgl:::rgl.setselectstate("none")

  # Exports
  isDone <- TRUE
  if (result$state == rgl:::msDONE) isDone <- FALSE
  return(list(isDone=isDone,isClosed=FALSE,vInter=vInter,vPlanes=vPlanes))
}

#############################################################
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
RotateMeshPlane3d <- function(mesh, planes, Pl, pl) {
    Stop <- 0
    # stop when the rotation is validated by ESC
    while (Stop==0) {
        temp <- selectPlanes(button="right", mesh=mesh, planes=planes, Pl=Pl, pl=pl)
        if (temp$isDone){
          # because of rgl.... need to close twice the window
          if (temp$isClosed){
            rgl.close()
          }
          break
        }
    }
    return(temp)
}

#############################################################
submesh<-function(spec,keep){

  # Fonction qui construit un sous mesh ? partir de l'objet mesh3d "mesh" selon les points ? conserver contenus dans le vecteur logique keep
  # Retourne un objet mesh3d

  spec2<-list()

  # extraction des vertices
  spec2$vb<-spec$vb[,keep]

  # extraction des edges
  idxV<-is.element(spec$it,which(keep))
  idxV<-matrix(idxV,dim(spec$it)[1],dim(spec$it)[2],byrow=FALSE)
  M<-spec$it[,which(colSums(idxV)==3)]
  spec2$it<-matrix(match(M,which(keep)),dim(M)[1],dim(M)[2])

  # extraction des normales
  spec2$normals<-spec$normals[,keep]

  # extraction du "material"
  spec2$material<-spec$material
  if (length(spec$material)>0){
    spec2$material$color<-spec$material$color[,which(colSums(idxV)==3)]
  }

  # d?finition de la class de spec2
  class(spec2) <- "mesh3d"

  return(spec2)
}


#############################################################
SetPtZoom <- function(dd, specFull, Trans, Pt, IdxPts=NULL, orthoplanes,
                      percDist=0.15, modify=FALSE, A=NULL, grDev) {

  # Fonction pour effectuer "un zoom" autour du point s?lectionn? sur le mesh d?cim? :
  # ouvre une 2?me fen?tre graphique avec la partie du mesh non d?cim? correspondante et pla?age du point

  # les points ? conserver se situent ? moins de 10% de la distance max observ?e entre le point s?lectionn? et l'ensmble des autres vertices du mesh
  keep <- dd < (percDist * max(dd)) # distances inf?rieures ? 10% distance max

  # extraction du sous-mesh (Use rmVertex instead ?)
  specFull2 <- submesh(specFull, keep)

  # si le sous-mesh contient plusieurs meshs isol?s : on ne conserve que le sous-mesh ? proximit? imm?diate du point cliqu?
  temp <- vcgIsolated(specFull2, split = TRUE)
  Min <- +Inf
  for (ii in 1:length(temp)) {
    vb <- temp[[ii]]$vb[1:3, ]
    cs <- sqrt(colSums((vb - (Pt + Trans))^2))
    MinT <- cs[which(cs==min(cs))]
    if (MinT < Min){
      specFull2 <- temp[[ii]]
      Min <- MinT
    }
  }

  # plot
  Trans2 <- rowMeans(specFull2$vb[1:3, ])
  specFull2$vb[1:3,] <- specFull2$vb[1:3, ] - Trans2

  param3d <- par3d()
  d2 <- open3d()
  par3d(windowRect = grDev$windowRect[2, ])
  ids2 <- plot3d(specFull2$vb[1, ], specFull2$vb[2, ], specFull2$vb[3, ],
                 size = grDev$ptSize, aspect = FALSE,
                 axes = F, box = F, xlab="", ylab="", zlab="", main = paste("Land - ", IdxPts))
  if (is.null(specFull2$material)) {
      specFull2$material$color <- matrix("gray", 3, dim(specFull2$it)[2])
  }
  shade3d(specFull2) # pb de trous parfois... ??????????????????????????

  # dessin des ?ventuelles intersections qui seraient visibles sur le sous-mesh
  if (!is.null(orthoplanes$vInter)){
    for (i in 1:3){
      # on ne conserve que les points d'intersections ? proximit? imm?diate du point cliqu?...
      ddi <- sqrt(colSums((t(orthoplanes$vInter[[i]])-Trans-Pt)^2))
      keep <- ddi < (percDist*max(ddi))

      if (sum(keep)> 0){
        # ... s'il y en a, on les dessine
        inter <- sweep(orthoplanes$vInter[[i]][keep, ], 2, Trans2)
            # - matrix(Trans2, nrow = dim(inter)[1], ncol=3, byrow=TRUE)
        lines3d(inter, col="red", lwd=2)
      }
    }
  }


  # r?glage de l'orientation du mesh de zoom pour qu'elle corresponde ? l'orientation du mesh d?cim?
  rgl.viewpoint(userMatrix = param3d$userMatrix)

  # pla?age du point sur le mesh de zoom
  if (is.null(grDev$spradius)) {
     tmp<-apply(specFull2$vb[1:3,],1,range)
     tmp<-tmp[2,]-tmp[1,]
     grDev$spradius<-(1/50)*min(tmp)
  }
  res2 <- SelectPoints3d(specFull2, modify, A, IdxPts, grDev)

  rgl.close()
  return(list(coords=res2$coords,sp=res2$sp,Trans2=Trans2,tx=res2$tx))
}


#############################################################
#' @title digitFixed
#' @description XXXX
#' @details
#' @param specFull full resolution mesh3d object
#' @param specDecim decimated mesh3d object
#' @param decim 0-1 number setting the amount of reduction relative to existing face number to decimate the full resolution mesh. Used for multiresolution and only if specDecim is NULL. See vcg:::vcgQEdecim
#' @param fixed number of landmarks to digitalize
#' @param templateFile template of 3D coordinates. Order of landmarks must be the
#' same than the one specified in index
#' @param idxPtsTemplate Indices of landmarks used to fit the template
#' @param index Digitalization sequence of landmarks
#' @param center logical whether or not centering is required
#' @param orthoplane logical whether or not orthoplane are draw
#' @param percDist
#' @param grDev list with ptSize the point size, windowRect xxxx, spradius and tcex
#' @export
#' @author Remi Laffont and Nicolas Navarro
#' @return
#' @examples
#'
DigitFixed <- function (specFull, specDecim = NULL, decim = 0.5, fixed, templateFile = NULL, idxPtsTemplate,
                        index = 1:fixed, center = TRUE, orthoplane = TRUE,
                        percDist = 0.15, grDev = list(windowRect = rbind(c(0,50,830,904), c(840,50,1672,904)),
                                                      ptSize = 1, spradius = NULL, tcex=2)) {

    if (class(specFull) != "mesh3d") stop("specFull must have class \"mesh3d\".")
    if (is.null(specFull$material$color)) {
        specFull$material$color <- matrix("gray", 3, dim(specFull$it)[2])
    }
    if (is.null(specDecim)) {
        print("Mesh decimation for multiresolution view")
        specDecim <- vcgQEdecim(specFull, percent = decim)
    }
    if (is.null(specDecim$material$color)) {
            specDecim$material$color <- matrix("gray", 3, dim(specDecim$it)[2])
    }
    if (missing(fixed)) {
        stop("missing number of landmarks to digitalize")
    } else if (!is.numeric(fixed) || length(fixed) > 1) stop("fixed must be a single number")

    # Use or not of a template
    if (is.null(templateFile)){
        idxPtsTemplate <- 1:fixed
    } else {
      if (missing(idxPtsTemplate)) {
          idxPtsTemplate <- 1:3
          warning("idxPtsTemplate was missing. First 3 landmarks will be used to align the template")
      } else
          if (length(idxPtsTemplate) < 3) stop("idxPtsTemplate must contain at least 3 landmarks")
      p1 <- length(idxPtsTemplate)
      template <- list()
      template$M <- templateFile
    }
    # Define default values for graphics interactivity
    if (is.null(grDev$windowRect))
        grDev$windowRect[1, ] <- c(0,50,830,904)
    if (length(grDev$windowRect) == 4)
        grDev$windowRect <- rbind(grDev$windowRect, c(840,50,1672,904))
    grDev$vSp <- grDev$vTx <- Sp <- Tx <- rep(NA, fixed)

    # Centering
  if (center == TRUE) {
      tmp <- scale(t(specDecim$vb[-4, ]), scale = FALSE)
      specDecim$vb <- rbind(t(tmp), 1)
      Trans <- attr(specDecim,"scaled:center") <- attr(tmp,"scaled:center")
  } else Trans <- attr(specDecim,"scaled:center") <- rep(0, 3)


  # plot decimated mesh
  d1 <- Clear3d()
  par3d(windowRect = grDev$windowRect[1, ])
  ids1 <- plot3d(specDecim$vb[1, ], specDecim$vb[2, ], specDecim$vb[3, ],
                 size = grDev$ptSize, aspect = FALSE,
                 axes = F, box = F, xlab="", ylab="", zlab="")
  shade3d(specDecim, add = TRUE)
  grDev$dev <- rgl.cur()

  if (is.null(grDev$spradius)) {
     tmp <- apply(specDecim$vb[1:3,], 1, range)
     tmp <- tmp[2,] - tmp[1,]
     grDev$spradius <- (1/50) * min(tmp)
  }

  # plot of orthogonal planes: they are initialized as major axes of the mesh
  orthoplanes <- list(vInter=NULL, vPlanes = NULL)
  if (orthoplane) orthoplanes <- DrawOrthoplanes(specDecim)

  # Landmark selection - A is the individual configuration matrix [k x 3]
  A <- Adeci <- matrix(NA, fixed, 3)
  rownames(A) <- 1:fixed
  colnames(A) <- c("x","y","z")
  attr(A, which="spec.name") <- deparse(substitute(specFull))

  Idx <- setdiff(1:fixed, index[idxPtsTemplate])
  for (i in 1:fixed){
    if (i <= length(idxPtsTemplate)){
      # Place 1st points require to adjust the template if it exist otherwise take all pts
      idx_pts <- index[idxPtsTemplate[i]]
      res <- SelectPoints3d(mesh=specDecim, A=A, IdxPts=idx_pts, grDev=grDev)
      Sp[idx_pts] <- res$sp
      Tx[idx_pts] <- res$tx
      Pt <- res$coords

      # distances full resolution mesh to template landmark
      dd <- sqrt(colSums((specFull$vb[1:3,] - Trans - Pt)^2))

      # zoom on full resolution mesh around the selected landmark
      res2 <- SetPtZoom(dd, specFull=specFull, Trans=Trans, Pt=Pt, IdxPts = idx_pts,
                        orthoplanes = orthoplanes, percDist = percDist, grDev=grDev)
      # A coordinates
      A[idx_pts, ] <- res2$coords + res2$Trans2 - Trans
      # Projection of landmarks on decimated mesh for graphics
      Adeci[idx_pts,] <- project(t(A[idx_pts, ]), specDecim,sign=FALSE)$vb[1:3]      # sign ????????????????????
      # plot
      grDev <- plotLandmarks(d1, Adeci[idx_pts,], Sp, Tx, idx_pts, grDev)

      if(!is.null(templateFile) & i==length(idxPtsTemplate)){
        # tous les points de r?f?rencement du template sont plac?s => calculs pour ajuster le template au mesh
        # d?finition de 3 matrices de configurations :
        # - configA : points plac?s sur le mesh
        # - configB : points contenus dans templateFile
        # - configC : points contenus dans templateFile comparables ? configA
        configA<-as.matrix(A[!is.na(A[,1]), ])
        configB<-as.matrix(template$M)
        p2<-dim(configB)[1]
        configC<-configB[idxPtsTemplate,]

        # transation & scaling de configA & configC
        transA<-colMeans(configA)
        AA<-configA-matrix(transA,p1,3,byrow=TRUE)
        scaleA<-1/sqrt(sum(AA^2))
        AA<-AA*scaleA
        transB<-colMeans(configC)
        BB<-configC-matrix(transB,p1,3,byrow=TRUE)
        scaleB<-1/sqrt(sum(BB^2))
        BB<-BB*scaleB

        # rotation de configC sur configA
        sv<-svd(t(AA)%*%BB)
        U<-sv$v
        V<-sv$u
        sig<-sign(det(t(AA)%*%BB))
        V[,3]<-sig * V[,3]
        rot<- U%*%t(V)

        # applications des transformations calcul?es ? partir de configC sur configB
        BB<-configB
        BB<-BB-matrix(transB,p2,3,byrow=TRUE)
        BB<-BB*scaleB
        BB<-BB%*%rot

        # mise ? l'?chelle et translation pour que configB retrouve une taille et une position comparable ? configA
        BB<-BB/scaleA
        BB<-BB+matrix(transA,p2,3,byrow=TRUE)

        # remise des valeurs des coord des points d?j? plac?s au pr?alable contenues dans A vers B
        B <- BB
        B[!is.na(A[,1])] <- A[!is.na(A[,1])]
        ptsB <- project(B, specDecim,sign=FALSE) # sign ??????????????????????????

        # plot des points/labels de B pas plac?s au pr?alables
        vv<-index[Idx]
        for (ii in 1:length(vv)){
          grDev$vSp[vv[ii]]<-spheres3d(t(ptsB$vb[1:3,vv[ii]]),color = "blue", alpha=0.5,radius=grDev$spradius)
          grDev$vTx[vv[ii]]<-text3d(t(ptsB$vb[1:3,vv[ii]]),texts=as.character(vv[ii]),col="red",cex=grDev$tcex)
        }
      }

    }else{
      # Selection of remaining landmarks (if any)
      idx_pts <- index[Idx[i-length(idxPtsTemplate)]]
      # distances full resolution mesh to automatic template landmark
      Pt <- B[idx_pts,] + Trans
      dd <- sqrt(colSums((specFull$vb[1:3,] - Pt)^2))

      # zoom on full resolution mesh around the selected landmark
      res2 <- SetPtZoom(dd, specFull=specFull, Trans=Trans, Pt=Pt, IdxPts = idx_pts,
                        orthoplanes = orthoplanes, percDist = percDist, grDev=grDev)
      # A coordinates
      A[idx_pts, ] <- res2$coords + res2$Trans2 - Trans
      # Projection of landmarks on decimated mesh for graphics
      Adeci[idx_pts,] <- project(t(A[idx_pts, ]), specDecim, sign=FALSE)$vb[1:3]
      # plot
      grDev <- plotLandmarks(d1, Adeci[idx_pts,], grDev$vSp, grDev$vTx, idx_pts, grDev)
    }
  }

  # Now, wait if the user want changed any landmark. Stop when the graphics is closed
  Stop <- 0
  grDev$dev <- d1
  while((Stop==0)){
    # on place le point sur le mesh d?cim?
    res <- SelectPoints3d(mesh = specDecim, modify = TRUE, A = Adeci, grDev = grDev)

    if (res$isClosed) break

    idx_pts <- res$Idx
    Sp[idx_pts] <- res$sp
    Tx[idx_pts] <- res$tx
    Pt <- res$coords

    # calcul des distances mesh complet/ point plac?
    dd <- sqrt(colSums((specFull$vb[1:3,] - Trans - Pt)^2))

    # zoom sur le mesh complet autour du point plac? + pla?age du point plus pr?cis
    res2 <- SetPtZoom(dd, specFull=specFull, Trans = Trans, Pt = Pt, IdxPts=idx_pts,
                      orthoplanes = orthoplanes, percDist = percDist, grDev =  grDev)
    # A coordinates
    A[idx_pts,] <- res2$coords + res2$Trans2 - Trans
    # Projection of landmarks on decimated mesh for graphics
    Adeci[idx_pts,] <- project(t(A[idx_pts,]), specDecim, sign=FALSE)$vb[1:3]      # sign ??????????????????????
    # plot
    grDev <- plotLandmarks(d1, Adeci[idx_pts,], Sp, Tx, idx_pts, grDev)

  }
  return(A)
}
plotLandmarks <- function(d1, pts, Sp, Tx, idx_pts, grDev){
    rgl.set(d1)
    rgl.pop("shapes", Sp[idx_pts])
    rgl.pop("shapes", Tx[idx_pts])
    grDev$vSp[idx_pts]<-spheres3d(pts,color = "green",alpha=0.5,radius=grDev$spradius)
    grDev$vTx[idx_pts]<-text3d(pts,texts=as.character(idx_pts),col="red",cex=grDev$tcex)

    return(grDev)
}
DrawOrthoplanes <- function(specDecim) {
        # centering of vertex coordinates
        pts <- t(specDecim$vb[1:3, ])
        A <- apply(pts, 2, mean)
        pts <- sweep(pts, 2, A)
        # Planes
        sv <- svd(crossprod(pts))

        # Compute and draw the intersection mesh/planes
        ptsPlanes <- array(NA, c(3,3,3))
        Pl <- pl <- list()
        for (i in 1:3){
            # plane PC2-PC3 when i=1, PC1-PC3 when i=2, PC1-PC2 when i=3
            ci <- setdiff(1:3, i)
            # intersection planes/mesh
            inter <- meshPlaneIntersect(specDecim, A, A+sv$u[,ci[1]], A+sv$u[,ci[2]])
            # Centering
            int <- sweep(inter, 2, A)
            # Rotation in order that the intersection pts lied in the plan
            int2 <- (int %*% svd(crossprod(int))$u)[, 1:2]
            # Order according to the angle of each vector (needed before to used lines3d())
            idx <- order(Arg(int2[,1] + 1i * int2[,2]))
            inter <- inter[idx, ]
            # plot of the intersection
            Pl[[i]] <- lines3d(inter, col = "red", lwd=2)
            # stockage
            #vInter[[i]] <- inter
            ptsPlanes[, , i] <- rbind(A, A+sv$u[,ci[1]], A+sv$u[,ci[2]])
            # for each plane, contains 3 points in the plane (needed just after)
        }

        # Draw orthogonal planes via planes3d()
        # coefficients a,b,c,d of the plane equation can be determined by a point of the plane
        # and a normal vector to this plane
        # https://fr.wikipedia.org/wiki/Plan_%28math%C3%A9matiques%29#D.C3.A9finition_par_un_vecteur_normal_et_un_point
        A <- ptsPlanes[1, , 1]
        for (i in 1:3) {
            if (i == 1) {
                # vect norm to plane 2-3
                n <- ptsPlanes[2, , 3] - ptsPlanes[1, , 3] #Nico ??? why not A and 3
            } else {
                # vect norm to plan 1-2 and plan 1-3
                n <- ptsPlanes[i, , 1] - A
            }
            d <- -t(n) %*% A
            pl[[i]] <- planes3d(n[1], n[2], n[3], d, alpha = 0.7, col="cyan")
        }
        # User interaction: manual rotation of the mesh (the orthogonal planes being fixed)
        # until the wanted orientation
        return(RotateMeshPlane3d(specDecim, planes = ptsPlanes, Pl=Pl, pl=pl))
}

