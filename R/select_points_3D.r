# Cross-product between two 3D vectors
CrossProd <- function(v1, v2){

    v <- v1
    v[1] <- v1[2]*v2[3] - v1[3]*v2[2]
    v[2] <- v1[3]*v2[1] - v1[1]*v2[3]
    v[3] <- v1[1]*v2[2] - v1[2]*v2[1]

    return(v)
}

# Simple R translation of the codes given in Möller & Trumbore 1997 to compute the intersection
# (if any) between a 3D triangle and a 3D line. The computation of 3D coordinates of this intersection,
# rather tan the barycentric coordinates, was added at the end.
IntersectionLineTriangle <- function(lineOri, lineDir, vert0, vert1, vert2, eps = 1e-6, cull = TRUE){

    # find vectors for two edges sharing vert0
    edge1 <- vert1 - vert0
    edge2 <- vert2 - vert0

    # begin calculating determinant - also used to calculate U parameter
    pvec <-  CrossProd(lineDir, edge2)

    # if determinant is near zero, line lies in plane of triangle
    det <-  as.numeric(edge1%*%pvec)

    if (cull){

        if (det < eps){
            return(0)
        }

        # calculate distance from vert0 to ray origin
        tvec <- lineOri - vert0

        # calculate U parameter and test bounds
        u <-  as.numeric(tvec%*%pvec)
        if ( u < 0 |  u > det){
            return(0)
        }

        # prepare to test V parameter
        qvec <- CrossProd(tvec, edge1)

        # calculate V parameter and test bounds
        v <-  as.numeric(lineDir%*%qvec)
        if ( v < 0 |  u +  v > det){
            return(0)
        }

        # calculate t, scale parameters, ray intersects triangle
        inv_det <- 1 / det
        t <- as.numeric(edge2%*%qvec)  * inv_det
        u <- u * inv_det
        v <- v * inv_det

    } else {
        # the non-culling branch

        if (det > -eps & det < eps){
            return(0)
        }
        inv_det <- 1 / det

        #calculate distance from vert0 to ray origin
        tvec <- lineOri - vert0

        # calculate U parameter and test bounds
        u <-  as.numeric(tvec%*%pvec) * inv_det
        if ( u < 0 |  u > 1){
            return(0)
        }

        # prepare to test V parameter
        qvec <- CrossProd(tvec, edge1)

        # calculate V parameter and test bounds
        v <-  as.numeric(lineDir%*%qvec) * inv_det
        if ( v < 0 |  u +  v > 1){
            return(0)
        }

        # calculate t, ray intersects triangle
        t <- as.numeric(edge2%*%qvec) * inv_det
    }

    #########
    # addition from the original algorithm:
    # computation of 3D coordinates for the intersection following the formulas in Möller & Trumbore
    #########

    coords <- c(0, u, v)

    tmp <- matrix(0,3,3)
    tmp[1,] <- lineDir*-1
    tmp[2,] <- edge1
    tmp[3,] <- edge2

    v1 <- tmp[,1]
    v2 <- tmp[,2]
    v3 <- tmp[,3]

    Tmp <- v1
    Tmp[1] <- as.numeric(v1%*%coords)
    Tmp[2] <- as.numeric(v2%*%coords)
    Tmp[3] <- as.numeric(v3%*%coords)

    coords <- Tmp

    coords <- coords + vert0

    return(coords)

}

whichIntersect <- function(Msign){
    # given a set of triangles centred on a given horizontal of vertical line, it identifies
    # which ones are overlapping this line, basing on the signs of the triangle vertices

    ss <- colSums(Msign)
    idx <- abs(ss) < 3

    return(idx)
}

rgl.user2window2 <- function(points, idata, projection = rgl.projection()){
    # slight simplification of the rgl.user2window function from rgl package

    ret <- .C(rgl:::rgl_user2window, success = FALSE, idata, as.double(points),
              window = double(length(points)), model = as.double(projection$model),
              proj = as.double(projection$proj), view = as.integer(projection$view))
    if (!ret$success)
        stop("'rgl_user2window2' failed")
    return(matrix(ret$window, ncol(points), 3, byrow = TRUE))

}

###########################
PlacePt <- function(x, y, verts, norms, it, start){

    # Given 2D coordinates of a point resulting from a user click on the screen, a mesh (though
    # its vertex coordinates, its normals and its triangular faces, a particular projection
    # matrix linked to the plot of the mesh, and the window dimensions, this function computes
    # the 3D coordinates of the intersections between a ray starting from the 2D clicked point, and
    # orthogonal to the screen, and the mesh, and keeps the closest to the virtual position of the
    # observer. The closest mesh vertex to this intersection is then returned as the intersection.
    #
    # The 2 main steps of this function are:
    # - identifying which triangles are possibly intersected: in a 2D screen projection of mesh
    # triangles, those triangles should lay both on a vertical and a horizontal line intersecting
    # themselves on the clicked point (necessary but not sufficient condition).
    # - computing the intersections between those triangles and the 3D ray orthogonal to the screen
    # and starting at the clicked point using Möller & Trumbire algorithm . Only the closest intersection
    # to the observer will be kept. A magnetism to the closest triangle vertex is then done.
    #
    # If no intersection is found, returns the closest mesh vertex to the clicked point.

    # conversion window coordinates in 2D screen coordinates
    # mesh -----------------------------------
    #temp <- rgl.user2window2(t(verts), as.integer(nrow(verts)), projection = start$projection)
    temp <- rgl.user2window(verts[,1], verts[,2], verts[,3], projection = start$projection)
    X <- temp[, 1] * start$viewport[3]
    Y <- (1 - temp[, 2]) * start$viewport[4]

    # 3D window coordinates of click point
    X1 <- x / start$viewport[3]
    Y1 <- 1 - y / start$viewport[4]

    # click point ----------------------------
    subS <- rgl.ids(type = "subscene", subscene = 0)
    if (dim(subS)[1] > 1 & tail(subS[, 1], 1) == currentSubscene3d()) {
        # a single window, set on the second subscene
        width <- par3d()$windowRect[3] - par3d()$windowRect[1]
        # screen width => width/2: subscene width
        x <- x + width/2
    }

    # centering 2D vertex point cloud on clicked point
    X_ <- X - x
    Y_ <- Y - y

    # getting signs
    sX <- sign(X_)
    sY <- sign(Y_)

    # finding which triangles intersect a vertical line passing on the clicked point
    id1 <- whichIntersect(matrix(sX[it], nrow = 3))

    noInter <- FALSE
    if (any(id1)){

        # finding, among those triangles, which ones intersect an horizontal line passing on the clicked point
        id2 <- whichIntersect(matrix(sY[it[, id1]], nrow = 3))

        if (any(id2)){
            # at least one triangle is intersected

            # index of triangles intersecting both lines (among those will be the one(s), if any, really
            # intersected by the ray between the observer and the clicked point)
            idTri <- which(id1)[which(id2)]

            # 3D coordinates of the given triangles
            idxVb <- it[, idTri, drop = FALSE]
            Coords <- matrix(t(verts[idxVb, ]), nrow=3)

            # conversion from 2D to 3D coordinates of the clicked point, and the associated ray
            Q <- as.numeric(rgl.window2user(X1, Y1, 0))
            normQ <- as.numeric(rgl.window2user(X1, Y1, 1)) - Q

            # computing 3D intersections between the ray and the given triangles
            tt <- sapply(1:length(idTri),
                         function(i){IntersectionLineTriangle(Q, normQ,Coords[,3*(i-1)+1], Coords[,3*(i-1)+2],
                                                              Coords[,3*(i-1)+3], cull = FALSE)})

            # keeping only the intersection data
            idx_1 <- which(sapply(tt,length)>1)
            tt <- tt[idx_1]
            Mt <- matrix(unlist(tt), nrow=3)

            # find the closest intersection to the observer
            closestInter <- project2(Q, Mt, idx = TRUE)
            theInter <- Mt[1, closestInter]
            the_Tri <- idTri[idx_1[closestInter]]

            # vertex indexes of the intersected triangle
            the_vb <- idxVb[,idx_1[closestInter]]

            # get shortest distance among the intersected point and the triangle vertex
            coo <- rbind(X_[the_vb], Y_[the_vb])
            vd <- sqrt(colSums(coo^2))

            visibles <- verts[the_vb,]
            idx <- which.min(vd)
            the_idx <- the_vb[idx]

        }else{
            # no triangle is intersected
            noInter <- TRUE
        }

    }else{
        # no triangle is intersected
        noInter <- TRUE
    }

    if (noInter){
        # find the closest vertex
        the_idx <- which.min((X - x)^2 + (Y - y)^2)
        visibles <- verts[the_idx,, drop=FALSE]
        idx <- 1
    }

    return(list(visibles = visibles, idx = idx, the_idx = the_idx))
}

###########################
distMin <- function(x, y) {
    tmp <- sqrt(colSums((x$vb - c(y$vb))^2))
    return(c(min(tmp), which.min(tmp)))
}

#############################################################
rgl.select2 <- function (button = c("left", "middle", "right"), verts, norms, it, #, mesh
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
        tmp <- PlacePt(x, y, verts, norms, it, start)
        visibles <<- tmp$visibles
        idx <<- tmp$idx
        the_idx <<- tmp$the_idx
        if (modify) {
            if (firstTime) {
                # When the user can modified the points already placed, we need to
                # determine which existing point is modified
                #Idx <<- which.min(sqrt(colSums((t(A)-c(visibles[idx, ]))^2)))
                Idx <<- project2(visibles[idx, ], t(A), idx = TRUE)
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
        if (grDev$zoomOptions$zoomSeeLm & whichMesh == 2){
            Ids <- rgl.ids()
            rgl.pop("shapes", Ids[Ids[,2] == "spheres",1])
        }
        # plot point/label
        Sp <<- spheres3d(visibles[idx, ], alpha = grDev$spheresOptions$spheresAlpha[1, whichMesh],
                         radius = grDev$spradius[1, whichMesh],
                         col = grDev$spheresOptions$spheresColor[1, whichMesh])
        Tx <<- text3d(visibles[idx, ], texts = as.character(IdxPts),
                      cex = grDev$labelOptions$labelCex[1, whichMesh],
                      col = grDev$labelOptions$labelColor[1, whichMesh],
                      adj = grDev$labadj[1, whichMesh,])
    }
    # Updating the position --------------------------
    Update <- function(x, y) {
        temp <- PlacePt(x, y, verts, norms, it, start)
        visibles <<- temp$visibles
        idx <<- temp$idx
        the_idx <<- temp$the_idx
        if (!is.null(Sp)){
            rgl.pop("shapes", Sp)
            rgl.pop("shapes", Tx)
        }
        Sp <<- spheres3d(visibles[idx, ], alpha = grDev$spheresOptions$spheresAlpha[1, whichMesh],
                         radius = grDev$spradius[1, whichMesh],
                         col = grDev$spheresOptions$spheresColor[1, whichMesh])
        Tx <<- text3d(visibles[idx, ], texts = as.character(IdxPts),
                      cex = grDev$labelOptions$labelCex[1, whichMesh],
                      col = grDev$labelOptions$labelColor[1, whichMesh],
                      adj = grDev$labadj[1, whichMesh,])
    }
    # Finalizing ----------------------------
    End <- function(x, y){ }

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
        subS <- rgl.ids(type = "subscene", subscene = 0)
        if (dim(subS)[1] > 1 & tail(subS[, 1], 1) == currentSubscene3d()) {
            # a single window, set on the second subscene
            width <- par3d()$windowRect[3]-par3d()$windowRect[1] # screen width => width/2: subscene width
            X <- X - width/2
        }

        # call of PlacePt
        tmp <- PlacePt(X, Y, verts, norms, it, start)
        visibles <- tmp$visibles
        idx <- tmp$idx
        the_idx <- tmp$the_idx

        # plot
        Sp <- spheres3d(visibles[idx, ], alpha = grDev$spheresOptions$spheresAlpha[1, whichMesh],
                         radius = grDev$spradius[1, whichMesh],
                         col = grDev$spheresOptions$spheresColor[1, whichMesh])
        Tx <- text3d(visibles[idx, ], texts = as.character(IdxPts),
                      cex = grDev$labelOptions$labelCex[1, whichMesh],
                      col = grDev$labelOptions$labelColor[1, whichMesh],
                      adj = grDev$labadj[1, whichMesh,])
    }

    # Modification of the user action when the user use right click
    # before it was a zoom, now track the surface by magnetism
    # (use of sub-functions Begin, Update and End)
    rMul <- rgl.setMouseCallbacks(2, Begin, Update, End)

    # Execution of the code is waiting until the user press ESC
    if (grDev$winOptions$winNb == 2){
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
        dev <- Dev
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

    return(list(coords = visibles[idx, ], Idx = Idx, isDone = isDone,
                isClosed = FALSE, sp = Sp, tx = Tx, the_idx = the_idx))
}

###########################
SelectPoints3d <- function (verts, it, norms, modify = FALSE, A = NULL, IdxPts = NULL, grDev,
                          whichMesh, prePlaced = NULL) {

    StopPts <- 0
    # Stop when landmark is validated by ESC or when the window is closed
    while (StopPts == 0) {
        temp <- rgl.select2(button="right", verts = verts, norms = norms,
                            it = it, modify = modify, A = A, IdxPts = IdxPts, grDev = grDev,
                            whichMesh = whichMesh, prePlaced = prePlaced)
        if (temp$isClosed || temp$isDone) {
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
SetPtZoom <- function(specFull, verts, it, norms, Pt, IdxPts = NULL, orthoplanes, idxPlanes,
                      modify = FALSE, A = NULL, grDev) {

    if (missing(grDev)) {
        stop("grDev missing without default value. See 'DigitFixed' ")
    }

    # Conserved only vertices at some distances (eg 15%) maximum of the cliked point
    dd <- sqrt(colSums((t(verts)-Pt)^2))
    maxRad <- grDev$zoomOptions$zoomPercDist * max(dd)
    IdxVert <- 1:nrow(verts)
    keep <- dd < maxRad
    specFull2 <- subset(specFull, subset = keep)
    IdxVert <- IdxVert[keep]

    # if the submesh contains several isolated meshes, we keep the closest to the clicked point
    tmp <- vcgIsolatedIndex(specFull2)
    # list of vertex coordinates split by isolated components
    isol <- tmp[[1]]
    # list of vertex indexes split by isolated components
    idx_vb <- tmp[[2]]
    # minimal distances (& associated index) intersection point/mesh by isoltated components
    vd <- lapply(isol, distMin, list(vb = matrix(c(Pt, 1), 4, 1)))
    vd <- matrix(unlist(vd), nrow = length(isol), ncol = 2, byrow = TRUE)
    # index of isolated component the closest to the intersection point
    idxL <- which.min(vd[, 1])
    # 'visible' submesh
    specFull2 <- isol[[idxL]]
    # 'absolute' indexes of the vertexes composing the 'visible' submesh
    IdxVert <- IdxVert[idx_vb[[idxL]]]

    # looking for possibly already digitized landmarks in the zoomed mesh
    prevLm <- NULL
    if (!is.null(A)){
        # indexes of the vertexes composing the 'visible' submesh
        v1 <- as.vector(IdxVert)
        # indexes of vertexes correponding to the already digitized landmarks
        v2 <- as.vector(attr(A, "vertex.idx"))
        # intersection of both index vectors
        inters_idx <- match(v1, v2, 0L)
        inters <- unique(v2[match(v1, v2, 0L)])
        # common vertexes (if any), ie landmarks visible on the submesh
        if (length(inters)>0){
            prevLm <- t(specFull$vb[1:3, inters, drop=FALSE])
        }
    }

    # transfers vertex colors (if any)
    specFull2 <- colTransf(specFull, specFull2)

    # center the vertices of the submesh
    Trans2 <- apply(specFull2$vb[1:3, ], 1, mean)
    specFull2$vb[1:3,] <- specFull2$vb[1:3, ] - Trans2
    verts2 <- t(specFull2$vb[1:3,])
    norms2 <- specFull2$normals
    it2 <- matrix(as.integer(specFull2$it), nrow=3)

    # plot
    param3d <- par3d()
    vb <- specFull2$vb[1:3,] + Trans2
    if (grDev$zoomOptions$zoomPtsDraw){
        # projection of zoom extent on the decimated mesh
        points3d(vb[1, ], vb[2, ], vb[3, ], col = grDev$zoomOptions$zoomPtsCol)
    }
    if (grDev$winOptions$winNb == 1){
        next3d()
    } else {
        d2 <- open3d()
        par3d(windowRect = grDev$winOptions$winSize[2, ])
    }

    # plot zoomed mesh
    if (grDev$meshOptions$meshVertCol){
        Col <- NULL
        ColPts <- specFull2$material$color[match(1:ncol(specFull2$vb), specFull2$it)]
    } else{
        Col <- ColPts <- grDev$meshOptions$meshColor[2]
    }
    if (grDev$meshOptions$meshShade[2]){
        shade3d(specFull2, col = Col, alpha = grDev$meshOptions$meshAlpha[2], specular="black")
    }
    if (grDev$meshOptions$meshWire[2]){
        wire3d(specFull2, col = Col, alpha = grDev$meshOptions$meshAlpha[2])
    }
    if (grDev$meshOptions$meshPoints[2]){
        points3d(t(specFull2$vb[1:3, ]), col = ColPts, alpha = grDev$meshOptions$meshAlpha[2])
    }

    # draws eventually the intersections visibles on the submesh
    if (!is.null(idxPlanes)){
        cpt <- 0
        for (i in idxPlanes){
            cpt <- cpt + 1
            # Only intersection points closed to the intersection planes
            ddi <- sqrt(colSums((t(orthoplanes$vInter[[i]]) - Pt)^2))
            keep <- ddi < maxRad
            if (sum(keep, na.rm=TRUE)> 0){
                inter <- t(t(orthoplanes$vInter[[i]][keep, ]) - Trans2)
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
    tmp <- diff(apply(specFull2$vb[1:3, ], 1, range))
    grDev$spradius[, 2] <- grDev$spheresOptions$spheresRad[, 2] * mean(tmp)

    # plot other visible landmarks (if any) on the submesh
    if (!is.null(prevLm) & grDev$zoomOptions$zoomSeePrevLm){
        prevLm <- prevLm - matrix(Trans2, nrow(prevLm), 3, byrow = TRUE)
        spheres3d(prevLm, alpha = grDev$spheresOptions$spheresAlpha[2, 2],
                  radius = grDev$spradius[2, 2],
                  col = grDev$spheresOptions$spheresColor[2, 2])
        text3d(prevLm[, 1], prevLm[, 2], prevLm[, 3], texts = as.character(inters_idx[which(inters_idx > 0)]),
               cex = grDev$labelOptions$labelCex[2, 2],
               col = grDev$labelOptions$labelColor[2, 2],
               adj = grDev$labadj[2, 2,])
    }

    # Add the point on the zoomed mesh (if asked)
    if (!is.null(IdxPts) & grDev$zoomOptions$zoomSeeLm) {
        res2 <- SelectPoints3d(verts2, it2, norms2, modify, A, IdxPts, grDev, whichMesh = 2, #specFull2
                               prePlaced = Pt - Trans2)
    } else {
        res2 <- SelectPoints3d(verts2, it2, norms2, modify, A, IdxPts, grDev, whichMesh = 2)
    }
    IdxVert <- IdxVert[res2$the_idx]
    # redefine the landmark coordinates from the selected vertex coordinates (in some cases, not stricltly equivalent
    # due to rounding from c++ functions used in Rvcg:::vcgIsolated)
    res2$coords <- matrix(specFull$vb[1:3,IdxVert], 1, 3)

    if (grDev$winOptions$winNb > 1){
        # Adjust the orientation of the decimated mesh to correspond to the one of zoomed mesh
        par3d(dev = grDev$dev, userMatrix = par3d(dev = d2)$userMatrix)
        grDev$winOptions$winSize[2, ] <- par3d()$windowRect
        rgl.close()
    } else {
        tmp <- rgl.ids()
        delFromSubscene3d(ids = tmp[, 1])
        subS <- rgl.ids(type="subscene", subscene = 0)
        useSubscene3d(subS[2, 1])
    }

    # pop the projection of zoom extent on the decimated mesh (taking care to don't delete possible points for mesh
    # plotting)
    if (grDev$zoomOptions$zoomPtsDraw) {
        Ids <- rgl.ids()
        idxIdsPts <- which(Ids[, 2] == "points")
        li <- length(idxIdsPts)
        if (li > 1) {
            idxIdsPts <- idxIdsPts[li]
            rgl.pop("shapes", Ids[idxIdsPts, 1])
        } else {
            if (li == 1 & !grDev$meshOptions$meshPoints[1]){
                rgl.pop("shapes", Ids[idxIdsPts, 1])
            }
        }
    }

    return(list(coords = res2$coords, sp = res2$sp, tx = res2$tx, grDev = grDev, IdxVert = IdxVert))
}


