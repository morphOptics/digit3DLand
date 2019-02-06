###########################
PlacePt <- function(x, y, verts, norms, mesh, start){

    # conversion window coordinates in 2D screen coordinates
    # mesh -----------------------------------
    temp <- rgl.user2window(x = verts[, 1], y = verts[, 2],
                            z = verts[, 3], projection = start$projection)
    X <- temp[, 1] * start$viewport[3]
    Y <- (1 - temp[, 2]) * start$viewport[4]

    # click point ----------------------------
    subS<-rgl.ids(type="subscene",subscene=0)
    if (dim(subS)[1]>1 & tail(subS[,1],1)==currentSubscene3d()){ # a single window, set on the second subscene
        width<-par3d()$windowRect[3]-par3d()$windowRect[1] # screen width => width/2: subscene width
        x<-x+width/2
    }

    # 3D window coordinates of click point
    X1 <- x / start$viewport[3]
    Y1 <- 1 - y / start$viewport[4]

    # Detrmine which vertices are in the selected square
    # mX, mY : matrices with xy-coordinates (2D screen) of the 3 vertices of each retained faces
    mX <- matrix(X[mesh$it], nrow = 3, byrow = FALSE)
    mY <- matrix(Y[mesh$it], nrow = 3, byrow = FALSE)
    # Length of each edge
    d <- mX
    d[1,] <- (mX[1,] - mX[2,])^2 + (mY[1,] - mY[2,])^2
    d[2,] <- (mX[1,] - mX[3,])^2 + (mY[1,] - mY[3,])^2
    d[3,] <- (mX[2,] - mX[3,])^2 + (mY[2,] - mY[3,])^2
    dm <- sqrt(max(d))
    # indices within the selection
    sqIdx <- (X >= (x-dm)) & (X <= (x+dm)) & (Y >= (y-dm)) & (Y <= (y+dm))
    # Extract the submesh
    selecMesh <- subset(mesh, subset = sqIdx)

    # defined the clicked point as a mesh3d
    Q <- rgl.window2user(X1, Y1, 0)
    normQ <- rgl.window2user(X1, Y1, 1) - Q
    normQ <- normQ / sqrt(sum(normQ^2))
    lQ <- list(vb = matrix(c(Q, 1), 4, 1), normals = matrix(c(normQ, 1), 4, 1))
    class(lQ) <- "mesh3d"

    # Search for the intersection ----------------------------
    int <- vcgRaySearch(lQ, mesh)
    # if there is no, then we need to search the nearest point
    # Use the angle to filter out possible points
    # However this return
    if (int$quality == 0){
         idx <- which.min((X - x)^2 + (Y - y)^2)
         int <- subset(mesh, subset = (1:dim(mesh$vb)[2])==idx, select = "vb")
    }
    # Defined visibles ----------------------------
    if (length(selecMesh$vb) == 0 || !is.matrix(selecMesh$it) || length(int$vb)==0 || attr(try(meshintegrity(selecMesh,facecheck=TRUE),silent=TRUE),"class")=="try-error" ) {
        # clic in empty zone or only 1 face or empty submesh -> undefined visibles
        visibles <- idx <-NULL
    } else {
        isol <- vcgIsolated(selecMesh, split=TRUE,silent=TRUE)
        vd <- lapply(isol, distMin, int)
        vd <- matrix(unlist(vd), length(isol), 2, byrow=TRUE)
        idx <- which.min(vd[, 1])
        visibles <- cbind(isol[[idx]]$vb[1,], isol[[idx]]$vb[2,], isol[[idx]]$vb[3,])
        idx <- vd[idx, 2]
    }
    return(list(visibles=visibles, idx=idx))
}

###########################
distMin <- function(x, y){
    tmp <- sqrt(apply(sweep(x$vb, 1, y$vb)^2, 2, sum)) # To Remi why it was colMeans?
    return(c(min(tmp), which.min(tmp)))
}

#############################################################
rgl.select2<-function (button = c("left", "middle", "right"), verts, norms, mesh,
                       modify=FALSE, A=NULL, IdxPts=NULL, grDev, whichMesh, prePlaced = NULL) {

    # Defines 3 sub-fonctions Begin(), Update() et End() that determine what
    # should be done when the user starts to click
    # Left click -> (Begin())
    # Move the mouse maintaining the button clicked -> (Update())
    # stop -> (End())

    start <- list()
    Sp <- Tx <- idx <- Idx <- NULL
    firstTime <- TRUE

    # Beginning of click -------------------------------
    Begin <- function(x, y) {
        # Infos on the projection of the mesh (global variable to use in Update())
        start$viewport <<- par3d("viewport")
        start$projection <<- rgl.projection()
        # Determine the mesh vertex the closest of the user click
        tmp <- PlacePt(x, y, verts, norms, mesh, start)
        visibles <<- tmp$visibles
        idx <<- tmp$idx
        if (modify) {
            if (firstTime) {
                # When the user can modified the points already placed, we need to
                # determine which existing point is modified
                Idx <<- which.min(sqrt(apply(sweep(A, 2, visibles[idx, ])^2, 1, sum)))
                IdxPts <<- Idx

                # remove the point and label to modify (needed only once)
                rgl.pop("shapes", grDev$vSp[Idx])
                rgl.pop("shapes", grDev$vTx[Idx])
                firstTime <<- FALSE
            } else {
                rgl.pop("shapes", Sp)
                rgl.pop("shapes", Tx)
            }
        } else {
            if (!is.null(Sp)){
                rgl.pop("shapes", Sp)
                rgl.pop("shapes", Tx)
            }
            Idx <<- NULL
        }
        if (grDev$zoomOptions$zoomSeeLm & whichMesh==2){
            Ids<-rgl.ids()
            rgl.pop("shapes",Ids[Ids[,2]=="spheres",1])
        }
        # plot point/label
        Sp <<- spheres3d(visibles[idx, ], alpha = grDev$spheresOptions$spheresAlpha[1,whichMesh],
                         radius=grDev$spradius[1,whichMesh],
                         col=grDev$spheresOptions$spheresColor[1,whichMesh])
        Tx <<- text3d(visibles[idx, ], texts = as.character(IdxPts),
                      cex=grDev$labelOptions$labelCex[1,whichMesh],
                      col= grDev$labelOptions$labelColor[1,whichMesh],
                      adj = grDev$labadj[1,whichMesh,])
    }
    # Updating the position --------------------------
    Update <- function(x, y) {
        temp <- PlacePt(x, y, verts, norms, mesh, start)
        visibles <<- temp$visibles
        idx <<- temp$idx
        if (!is.null(Sp)){
            rgl.pop("shapes", Sp)
            rgl.pop("shapes", Tx)
        }
        Sp <<- spheres3d(visibles[idx, ], alpha = grDev$spheresOptions$spheresAlpha[1,whichMesh],
                         radius=grDev$spradius[1,whichMesh],
                         col=grDev$spheresOptions$spheresColor[1,whichMesh])
        Tx <<- text3d(visibles[idx, ], texts = as.character(IdxPts),
                      cex=grDev$labelOptions$labelCex[1,whichMesh],
                      col= grDev$labelOptions$labelColor[1,whichMesh],
                      adj = grDev$labadj[1,whichMesh,])
    }
    # Finalizing ----------------------------
    End <- function(x,y){ }

    # Codes based on rgl.select()
    button <- match.arg(button)
    newhandler <- par3d("mouseMode")
    newhandler[button] <- "selecting"
    oldhandler <- par3d(mouseMode = newhandler)
    on.exit(par3d(mouseMode = oldhandler))

    if (!is.null(prePlaced)){
        # In the case where the landmark placed on the decimated mesh is pre-positionned on the full mesh and is
        # directly validated by user (via esc key) without any manual change, the outputs returned by the mouse action
        # need to be returned anyway.
        # The solution choosen here is to call the PlacePt() function with the coordinates of the prePlaced landmark,
        # to obtain those outputs for this unmodified landmark. Because PlacePt() requires 2D screen coordinates, we
        # first need to convert the absolute 3D coordinates in prePlaced into 2D screen coordinates:

        start$viewport <- par3d("viewport")
        start$projection <- rgl.projection()

        # conversion true 3D coordinates in 3D window coordinates
        temp <- rgl.user2window(prePlaced[1], prePlaced[2], prePlaced[3])
        # conversion window coordinates in 2D screen coordinates
        X <- temp[, 1] * start$viewport[3]
        Y <- (1 - temp[, 2]) * start$viewport[4]
        subS<-rgl.ids(type="subscene",subscene=0)
        if (dim(subS)[1]>1 & tail(subS[,1],1)==currentSubscene3d()){ # a single window, set on the second subscene
            width<-par3d()$windowRect[3]-par3d()$windowRect[1] # screen width => width/2: subscene width
            X<-X-width/2
        }

        # call of PlacePt
        tmp <- PlacePt(X, Y, verts, norms, mesh, start)
        visibles <- tmp$visibles
        idx <- tmp$idx

        # plot
        Sp <- spheres3d(visibles[idx, ], alpha = grDev$spheresOptions$spheresAlpha[1,whichMesh],
                         radius=grDev$spradius[1,whichMesh],
                         col=grDev$spheresOptions$spheresColor[1,whichMesh])
        Tx <- text3d(visibles[idx, ], texts = as.character(IdxPts),
                      cex=grDev$labelOptions$labelCex[1,whichMesh],
                      col= grDev$labelOptions$labelColor[1,whichMesh],
                      adj = grDev$labadj[1,whichMesh,])
    }

    # Modification of the user action when the user use right click
    # before it was a zoom, now track the surface by magnetism
    # (use of sub-functions Begin, Update and End)
    rMul <- rgl.setMouseCallbacks(2, Begin, Update, End)

    # Execution of the code is waiting until the user press ESC
    if (grDev$winOptions$winNb==2){
        dev <- rgl.cur()
        while (dev == rgl.cur()) {
            if (!is.null(idx)) {
                result <- rgl:::rgl.selectstate()
                # if ESC -> get out
                if (result$state >= rgl:::msDONE)
                    break
            }
        }
    }else{
        Dev <- rgl.cur()
        dev <- currentSubscene3d()
        while (dev == currentSubscene3d()) {
            if (!is.null(idx)) {
                result <- rgl:::rgl.selectstate()
                # if ESC -> get out
                if (result$state >= rgl:::msDONE)
                    break
            }
        }
        dev<-Dev
    }

    # if the window has been closed -> get out
    if (dev != rgl.cur()) {
        if (modify){
            isDone <- TRUE
            return(list(isDone = isDone, isClosed = TRUE))
        }
    }

    # Otherwise, redefined right click by default zoom
    rgl.setMouseCallbacks(2)
    rgl:::rgl.setselectstate("none")

    # Exports
    isDone <- TRUE
    if (result$state == rgl:::msDONE) isDone <- FALSE

    return(list(coords=visibles[idx, ], Idx=Idx, isDone=isDone, isClosed=FALSE, sp=Sp, tx=Tx))
}

###########################
SelectPoints3d<-function (mesh, modify=FALSE, A=NULL, IdxPts=NULL, grDev, whichMesh, prePlaced = NULL) {
    # Do the transpose only once
    verts <- t(mesh$vb[1:3, ])
    norms <- t(mesh$normals[1:3, ])
    StopPts <- 0
    # Stop when landmark is validated by ESC or when the window is closed
    while (StopPts==0) {
        temp <- rgl.select2(button="right", verts = verts, norms = norms, mesh = mesh,
                            modify = modify, A = A, IdxPts = IdxPts, grDev = grDev, whichMesh=whichMesh,
                            prePlaced=prePlaced)
        if (temp$isClosed || temp$isDone){
            if (temp$isClosed){
                # Need to close twice because of rgl...
                rgl.close()
            }
            break
        }
    }
    return(temp)
}

###########################
SetPtZoom <- function(specFull, Pt, IdxPts=NULL, orthoplanes, idxPlanes,
                      modify=FALSE, A=NULL, grDev) {

    if (missing(grDev)) stop("grDev missing without default value. See 'DigitFixed' ")
    # Conserved only vertices at some distances (eg 15%) maximum of the cliked point
    dd <- sqrt(apply(sweep(specFull$vb[1:3,], 1, Pt)^2, 2, sum))
    maxRad<-(grDev$zoomOptions$zoomPercDist * max(dd))
    keep <- dd < maxRad
    specFull2 <- subset(specFull, subset = keep)

    # if the submesh contains several isolated meshes, we keep the closest to the clicked point
    tmp <- vcgIsolated(specFull2, split = TRUE, silent=TRUE)
    vd <- lapply(tmp, distMin, list(vb = matrix(c(Pt, 1), 4, 1)))
    vd <- matrix(unlist(vd), length(tmp), 2, byrow=TRUE)
    specFull2 <- tmp[[which.min(vd[, 1])]]

    # transfers vertex colors (if any)
    specFull2 <- colTransf(specFull, specFull2)

    # center the vertices of the submesh
    Trans2 <- apply(specFull2$vb[1:3, ], 1, mean)
    specFull2$vb[1:3,] <- specFull2$vb[1:3, ] - Trans2

    # plot
    param3d <- par3d()
    vb<-specFull2$vb[1:3,]+Trans2
    if (grDev$zoomOptions$zoomPtsDraw){
        # projection of zoom extent on the decimated mesh
        points3d(vb[1, ], vb[2, ], vb[3, ], col=grDev$zoomOptions$zoomPtsCol)
    }
    if (grDev$winOptions$winNb==1){
        next3d()
    }else{
        d2 <- open3d()
        par3d(windowRect = grDev$winOptions$winSize[2,])
    }

    # plot zoomed mesh
    if (grDev$meshOptions$meshVertCol){
        Col <- NULL
        ColPts <- specFull2$material$color[match(1:ncol(specFull2$vb),specFull2$it)]
    }else{
        Col <- ColPts <- grDev$meshOptions$meshColor[2]
    }
    if (grDev$meshOptions$meshShade[2]){
        shade3d(specFull2, col=Col, alpha=grDev$meshOptions$meshAlpha[2])
    }
    if (grDev$meshOptions$meshWire[2]){
        wire3d(specFull2, col=Col, alpha=grDev$meshOptions$meshAlpha[2])
    }
    if (grDev$meshOptions$meshPoints[2]){
        points3d(t(specFull2$vb[1:3,]), col=ColPts, alpha=grDev$meshOptions$meshAlpha[2])
    }

    # draws eventually the intersections visibles on the submesh
    if (!is.null(idxPlanes)){
        cpt<-0
        for (i in idxPlanes){
            cpt<-cpt+1
            # Only intersection points closed to the intersection planes
            ddi <- sqrt(apply(sweep(t(orthoplanes$vInter[[i]]), 1, Pt)^2, 2, sum))
            keep <- ddi < maxRad
            if (sum(keep, na.rm=TRUE)> 0){
                inter <- sweep(orthoplanes$vInter[[i]][keep, ], 2, Trans2)
                if (grDev$intersectOptions$intersectLines[cpt]){
                    lines3d(inter, col = grDev$intersectOptions$intersectColor[cpt])
                }
                if (grDev$intersectOptions$intersectPoints[cpt]){
                    points3d(inter, col = grDev$intersectOptions$intersectColor[cpt])
                }
            }
        }
    }

    # Adjust the orientation of the zoomed mesh to correspond to the one of the decim mesh
    rgl.viewpoint(userMatrix = param3d$userMatrix)

    # set sphere radius
    tmp <- diff(apply(specFull2$vb[1:3,], 1, range))
    grDev$spradius[,2] <- grDev$spheresOptions$spheresRad[,2] * mean(tmp)

    # Add the point on the zoomed mesh (if asked)
    if (!is.null(IdxPts) & grDev$zoomOptions$zoomSeeLm) {
        res2 <- SelectPoints3d(specFull2, modify, A, IdxPts, grDev, whichMesh = 2, prePlaced = Pt - Trans2)
    }else{
        res2 <- SelectPoints3d(specFull2, modify, A, IdxPts, grDev, whichMesh = 2)
    }

    res2$coords <- matrix(res2$coords + Trans2, 1, 3)
    if (grDev$winOptions$winNb>1){
        # Adjust the orientation of the decimated mesh to correspond to the one of zoomed mesh
        par3d(dev=grDev$dev, userMatrix = par3d(dev=d2)$userMatrix)
        grDev$winOptions$winSize[2, ]<-par3d()$windowRect
        rgl.close()
    }else{
        tmp<-rgl.ids()
        delFromSubscene3d(ids=tmp[,1])
        subS<-rgl.ids(type="subscene",subscene=0)
        useSubscene3d(subS[2,1])
    }

    # pop the projection of zoom extent on the decimated mesh (taking care to don't delete possible points for mesh
    # plotting)
    if (grDev$zoomOptions$zoomPtsDraw){
        Ids<-rgl.ids()
        idxIdsPts<-which(Ids[,2]=="points")
        li<-length(idxIdsPts)
        if (li>1){
            idxIdsPts<-idxIdsPts[li]
            rgl.pop("shapes",Ids[idxIdsPts,1])
        }else{
            if (li==1 & !grDev$meshOptions$meshPoints[1]){
                rgl.pop("shapes",Ids[idxIdsPts,1])
            }
        }
    }

    return(list(coords = res2$coords, sp = res2$sp, tx = res2$tx, grDev=grDev))
}

# DrawSpheres <- function(Pt){
#     alpha <- matrix(seq(0, 2*pi, by=pi/8), ncol=1)
#     phi <- seq(-pi/2, pi/2, by=pi/50)
#     r <- grDev$spradius * 5
#
#     circl <- function(alpha, phi, r, Pt){
#         x <- r * cos(alpha) * cos(phi) + Pt[1]
#         y <- r * sin(alpha) * cos(phi)
#         z <- r * sin(phi)
#         lines3d(x, z + Pt[2], y + Pt[3], col="lightskyblue2")
#         lines3d(x, y + Pt[2], z + Pt[3], col="lightskyblue2")
#     }
#     apply(alpha, 1, circl, phi, r, Pt)
# }


