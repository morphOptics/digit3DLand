#' @title PlacePt
#' @description xxxx
#'
#' @details XXX
#' @param x Screen x-coordinate of the clicked point
#' @param y Screen y-coordinate of the clicked point
#' @param verts XXX
#' @param norms XXX
#' @param mesh 3d mesh
#' @param start xxxx
#'
#' @return list with the visible vertices and their index
#'
PlacePt <- function(x, y, verts, norms, mesh, start){
    # click point ----------------------------
    temp <- rgl.user2window(x = verts[, 1], y = verts[, 2],
                            z = verts[, 3], projection = start$projection)
    # conversion window coordinates in 2D screen coordinates
    X <- temp[, 1] * start$viewport[3]
    Y <- (1 - temp[, 2]) * start$viewport[4]

    # 3D window coordinates of cliked point
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
    if (length(selecMesh$vb) == 0 || !is.matrix(selecMesh$it) || length(int$vb)==0) {
        # clic in empty zone or only 1 face or empty submesh -> undefined visibles
        visibles <- idx <-NULL
    } else {
        isol <- vcgIsolated(selecMesh, split=TRUE)
        vd <- lapply(isol, distMin, int)
        vd <- matrix(unlist(vd), length(isol), 2, byrow=TRUE)
        idx <- which.min(vd[, 1])
        visibles <- cbind(isol[[idx]]$vb[1,], isol[[idx]]$vb[2,], isol[[idx]]$vb[3,])
        idx <- vd[idx, 2]
    }
    return(list(visibles=visibles, idx=idx))
}

distMin <- function(x, y){
    tmp <- sqrt(apply(sweep(x$vb, 1, y$vb)^2, 2, sum)) # To Remi why it was colMeans?
    return(c(min(tmp), which.min(tmp)))
}

#############################################################
rgl.select2<-function (button = c("left", "middle", "right"), verts, norms, mesh,
                       modify=FALSE, A=NULL, IdxPts=NULL, grDev) {

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
        # plot point/label
        Sp <<- spheres3d(visibles[idx, ], alpha = 0.5, radius=grDev$spradius)
        Tx <<- text3d(visibles[idx, ], texts = as.character(IdxPts),
                      cex=grDev$tcex, adj = rep(grDev$spradius, 2))
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
        Sp <<- spheres3d(visibles[idx, ], alpha=0.5, radius=grDev$spradius)
        Tx <<- text3d(visibles[idx, ], texts = as.character(IdxPts),
                    cex=grDev$tcex, adj = rep(grDev$spradius, 2))
    }
    # Finalizing ----------------------------
    End <- function(x,y){ }

    # Codes based on rgl.select()
    button <- match.arg(button)
    newhandler <- par3d("mouseMode")
    newhandler[button] <- "selecting"
    oldhandler <- par3d(mouseMode = newhandler)
    on.exit(par3d(mouseMode = oldhandler))

    # Modification of the user action when the user use right click
    # before it was a zoom, now track the surface by magnetism
    # (use of sub-functions Begin, Update and End)
    rMul <- rgl.setMouseCallbacks(2, Begin, Update, End)

    # Execution of the code is waiting until the user press ESC
    dev <- rgl.cur()
    while (dev == rgl.cur()) {
        if (!is.null(idx)) {
            result <- rgl:::rgl.selectstate()
            # if ESC -> get out
            if (result$state >= rgl:::msDONE)
                break
        }
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


#' @title SelectPoints3d
#' @description Function selects a point on the surface of a 3d mesh
#' @details A 3d point is selected using OSX:CMD+mouse, Win:?, Linux:? and
#' may be optionally moved by not releasing CMD. Final selection is done with ESC
#' @param mesh A 3d mesh as opened by (see \code{\link[Rvcg]{vcgPlyRead}})
#' @param modify logical allowing to modify continuously the point by tracking the surface
#' @param A Optional matrix that store all landmarks
#' @param IdxPts XXX
#' @param grDev Graphical parameters (see \code{\link[digit3DLand]{DigitFixed}})
#'
#' @return xxxx vvvvv

SelectPoints3d<-function (mesh, modify=FALSE, A=NULL, IdxPts=NULL, grDev) {
    # Do the transpose only once
    verts <- t(mesh$vb[1:3, ])
    norms <- t(mesh$normals[1:3, ])
    StopPts <- 0
    # Stop when landmark is validated by ESC or when the window is closed
    while (StopPts==0) {
        temp <- rgl.select2(button="right", verts = verts, norms = norms, mesh = mesh,
                            modify = modify, A = A, IdxPts = IdxPts, grDev = grDev)
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

#' @title SetPtZoom
#' @description Zoom around the selected point on the decimated mesh
#' @details Function zoomed on the selected point by opening a new 3d scene with
#' the full resolution mesh within some distance of the point
#' @param specFull Full resolution mesh
#' @param Pt Numeric vector with the xyz-coordinates of the clicked point
#' @param IdxPts xxx
#' @param orthoplanes xxx
#' @param percDist proportion of the maximu distance between the point and all vertices
#' @param modify if the point is modifiable
#' @param A Matrix of coordinates of landmarks
#' @param grDev Some graphical parameters
#'
#' @return xxx xxxx
#'
SetPtZoom <- function(specFull, Pt, IdxPts=NULL, orthoplanes,
                      percDist=0.15, modify=FALSE, A=NULL, grDev) {

    if (missing(grDev)) stop("grDev missing without default value. See 'DigitFixed' ")
    # Conserved only vertices at some distances (eg 15%) maximum of the cliked point
    dd <- sqrt(apply(sweep(specFull$vb[1:3,], 1, Pt)^2, 2, sum))
    keep <- dd < (percDist * max(dd))
    specFull2 <- subset(specFull, subset = keep)

    # if the submesh contains several isolated meshes, we keep the closest to the clicked point
    tmp <- vcgIsolated(specFull2, split = TRUE)
    vd <- lapply(tmp, distMin, list(vb = matrix(c(Pt, 1), 4, 1)))
    vd <- matrix(unlist(vd), length(tmp), 2, byrow=TRUE)
    specFull2 <- tmp[[which.min(vd[, 1])]]

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
            ddi <- sqrt(apply((t(orthoplanes$vInter[[i]]) - Pt)^2, 2, sum))
            keep <- ddi < (percDist * max(ddi))
            if (sum(keep)> 0){
                inter <- sweep(orthoplanes$vInter[[i]][keep, ], 2, Trans2)
                lines3d(inter, col="red", lwd=2)
            }
        }
    }

    # Adjust the orientation of the zoomed mesh to correspond to the one of the decim mesh
    rgl.viewpoint(userMatrix = param3d$userMatrix)

    if (is.null(grDev$spradius)) {
        tmp <- diff(apply(specFull2$vb[1:3,], 1, range))
        grDev$spradius <- (1/50)*min(tmp)
    }
    # Add the point on the zoomed mesh
      # if (!is.null(IdxPts)) {
      #     DrawSpheres(Pt - Trans2)
      #   #spheres3d(Pt - Trans2, alpha=0.10, color = "lightskyblue2", radius=8*grDev$spradius)
      # }

    res2 <- SelectPoints3d(specFull2, modify, A, IdxPts, grDev)
    res2$coords <- matrix(res2$coords + Trans2, 1, 3)
    # Adjust the orientation of the decimated mesh to correspond to the one of zoomed mesh
    par3d(dev=grDev$dev, userMatrix = par3d(dev=d2)$userMatrix)

    rgl.close()
    return(list(coords = res2$coords, sp = res2$sp, tx = res2$tx))
}

DrawSpheres <- function(Pt){
    alpha <- matrix(seq(0, 2*pi, by=pi/8), ncol=1)
    phi <- seq(-pi/2, pi/2, by=pi/50)
    r <- grDev$spradius * 5

    circl <- function(alpha, phi, r, Pt){
        x <- r * cos(alpha) * cos(phi) + Pt[1]
        y <- r * sin(alpha) * cos(phi)
        z <- r * sin(phi)
        lines3d(x, z + Pt[2], y + Pt[3], col="lightskyblue2")
        lines3d(x, y + Pt[2], z + Pt[3], col="lightskyblue2")
    }
    apply(alpha, 1, circl, phi, r, Pt)
}


