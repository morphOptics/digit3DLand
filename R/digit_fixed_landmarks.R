#' @title digitFixed
#' @description XXXX
#' @details cccccc
#' @param specFull full resolution mesh3d object
#' @param decim 0-1 number setting the amount of reduction relative to existing face number to decimate the full resolution mesh. Used for multiresolution. (see \code{\link[Rvcg]{vcgQEdecim}})
#' @param fixed number of landmarks to digitalize
#' @param index digitalization sequence of landmarks
#' @param templateFile template of 3D coordinates. Order of landmarks must be the
#' same than the one specified in index
#' @param idxPtsTemplate Indices of landmarks used to fit the template
#' @param orthoplane logical whether or not orthoplane are draw
#' @param percDist percentage of distance around the landmark use for zooming
#' @param grDev list of graphical parameters with ptSize the point size, windowRect xxxx, spradius and tcex
#' @param ... additional graphical parameters to be passed to plot (color, alpha, col)
#' @export
#' @return XXXX
#'
DigitFixed <- function (specFull, decim = 0.25, fixed, index = 1:fixed,
                        templateFile = NULL, idxPtsTemplate,
                        drawPCplane = FALSE, percDist = 0.15,
                        grDev = list(windowRect = rbind(c(0,50,830,904), c(840,50,1672,904)),
                                     ptSize = 1, spradius = NULL, tcex=2, nbWin=2), ...) {

    if (!(any(class(specFull) == "mesh3d")))
        stop("specFull must have class \"mesh3d\".")
    spec.name <- deparse(substitute(specFull))

    # check which setting of drawPCplane is called, and set idxPlanes consequently
    if(is.logical(drawPCplane)){
        if (drawPCplane){
            idxPlanes<-1:3
        }else{
            idxPlanes<-NULL
        }
    }else{
        V<-c("pc2-pc3","pc1-pc3","pc1-pc2")
        idxPlanes<-which(is.element(V,tolower(drawPCplane)))
    }

    # Correction if mesh has non-manifold faces (ie faces made of non-manifold edges, ie edges shared by more than 2 faces)
    # Correction needed for the ordering of the intersection points among mesh and planes
    specFull<-vcgClean(specFull,sel=2)

    if (decim != 1) {
        cat("\nMesh decimation for multiresolution view ----------\n")
        specDecim <- vcgQEdecim(specFull, percent = decim)
        specDecim <- vcgUpdateNormals(specDecim, silent = TRUE)
        # Correction if mesh has non-manifold faces
        specDecim<-vcgClean(specDecim,sel=2)
        cat("---------------------------------------------------\n")
    } else specDecim <- specFull

    if (is.null(specFull$material)) {
        specFull$material <- specDecim$material <- "gray"
    } else {
        specDecim$material <- "gray" # We may want to fix this by transfering vertex color etc
    }
    if (missing(fixed)) {
        stop("missing number of landmarks to digitalize")
    } else {
        if (!is.numeric(fixed) || length(fixed) > 1)
            stop("fixed must be a single number")
    }

    # Use or not of a template
    if (is.null(templateFile)){
        idxPtsTemplate <- 1:fixed
    } else {
        if (missing(idxPtsTemplate)) {
            idxPtsTemplate <- 1:4
            warning("idxPtsTemplate was missing.
                    First 4 landmarks will be used to align the template")
        } else {
            if (length(idxPtsTemplate) < 4)
                stop("idxPtsTemplate must contain at least 4 landmarks")
        }
        p1 <- length(idxPtsTemplate)
        template <- list()
        template$M <- templateFile
    }
    # Define default values for graphics interactivity
    grDev$vSp <- grDev$vTx <- Sp <- Tx <- rep(NA, fixed)
    if (is.null(grDev$windowRect))
        grDev$windowRect[1, ] <- c(0,50,830,904)
    if (length(grDev$windowRect) == 4)
        grDev$windowRect <- rbind(grDev$windowRect, c(840,50,1672,904))
    if (is.null(grDev$spradius)) {
        tmp <- diff(apply(specDecim$vb[1:3,], 1, range))
        grDev$spradius <- (1/50) * min(tmp)
    }

    # Centering of the meshes on the centroid of the decimated one
    tmp <- scale(t(specDecim$vb[-4, ]), scale = FALSE)
    specDecim$vb[-4, ] <- t(tmp)
    specFull$vb[-4, ] <- sweep(specFull$vb[-4, ], 1, attr(tmp, which="scaled:center"))

    # plot decimated mesh
    d1 <- Clear3d()
    par3d(windowRect = grDev$windowRect[1, ])
    if (grDev$nbWin==1){
        d1 <- currentSubscene3d()
        layout3d(t(c(1,2)), sharedMouse = TRUE)
        next3d()
    }
    shade3d(specDecim)
    grDev$dev <- rgl.cur()

    # Rotate the scene
    R <- rotMajorAxes(specDecim$vb[1:3, ])
    par3d(userMatrix = R)

    # plot of orthogonal planes: they are initialized as major axes of the mesh
    orthoplanes <- list(vInter=NULL, vPlanes = NULL)
    if (length(idxPlanes)>0){
        orthoplanes <- DrawOrthoplanes(specDecim,idxPlanes)
        if (ncol(specDecim$vb)!=ncol(specFull$vb)){
            # computation of intersections among full mesh and fixed planes
            orthoplanes <- DrawOrthoplanes(specFull,idxPlanes,planes=orthoplanes$vPlanes,interactive=FALSE,is.plot=FALSE)
        }
    }


    # Landmark selection - A is the individual configuration matrix [k x 3]
    A <- Adeci <- matrix(NA, fixed, 3, dimnames = list(1:fixed, c("x","y","z")))
    attr(A, which = "spec.name") <- spec.name

    Idx <- setdiff(1:fixed, index[idxPtsTemplate])
    for (i in 1:fixed){
        if (i <= length(idxPtsTemplate)){
            # Place 1st points require to adjust the template if it exist otherwise take all pts
            idx_pts <- index[idxPtsTemplate[i]]
            res <- SelectPoints3d(mesh = specDecim, A = A, IdxPts = idx_pts, grDev = grDev)

            Pt <- res$coords
            grDev$vSp[idx_pts] <- res$sp
            grDev$vTx[idx_pts] <- res$tx
        } else {
            # Selection of remaining landmarks (if any)
            idx_pts <- index[Idx[i-length(idxPtsTemplate)]]
            #Pt <- B[idx_pts, , drop = FALSE]
            Pt <- B[idx_pts, ]
        }
        # zoom on full resolution mesh around the selected landmark
        res2 <- SetPtZoom(specFull=specFull, Pt = Pt, IdxPts = idx_pts,
                          orthoplanes = orthoplanes, idxPlanes=idxPlanes, percDist = percDist, grDev=grDev)
        grDev<-res2$grDev
        # landmark coordinate on the full resolution mesh
        A[idx_pts, ] <- res2$coords
        # Projection of landmarks on decimated mesh for graphics
        Adeci[idx_pts,] <- project(res2$coords, specDecim, trans = TRUE)

        # Graphics
        grDev <- plot.landmark(Adeci[idx_pts, ], d1, idx_pts, grDev, exist = TRUE, ...)

        if(!is.null(templateFile) & i==length(idxPtsTemplate)){
            # all reference points of the template are placed
            # => impute missing landmarks
            B <- imputeCoords(A, template = template$M) #idx = idxPtsTemplate
            ptsB <- project(B, specDecim)
            B <- project(B, specFull, trans = TRUE)

            # plot points/labels of B not placed before
            vv <- index[Idx]
            for (ii in 1:length(vv)){
                grDev <- plot.landmark(t(ptsB$vb[1:3, vv[ii]]), d1, vv[ii], grDev,
                                       exist = FALSE, color = "blue", col = "red")
            }
            rgl.viewpoint(userMatrix = R)
        }
    }
    # Now, wait if the user want changed any landmark. Stop when the graphics is closed
    Stop <- 0
    grDev$dev <- d1
    while((Stop==0)){
        # clicks point on the decimated mesh
        res <- SelectPoints3d(mesh = specDecim, modify = TRUE, A = Adeci, grDev = grDev)
        if (res$isClosed) break

        idx_pts <- res$Idx
        # zoom on full resolution mesh
        res2 <- SetPtZoom(specFull=specFull, Pt = res$coords, IdxPts = idx_pts,
                          orthoplanes = orthoplanes, idxPlanes=idxPlanes, percDist = percDist, grDev =  grDev)
        grDev<-res2$grDev
        # landmark coordinate on the full resolution mesh
        A[idx_pts, ] <- res2$coords
        # Projection of the landmark on the decimated mesh for graphics
        Adeci[idx_pts, ] <- project(res2$coords, specDecim, trans = TRUE)

        # Graphics
        grDev$vSp[idx_pts] <- res$sp
        grDev$vTx[idx_pts] <- res$tx
        grDev <- plot.landmark(Adeci[idx_pts, ], d1, idx_pts, grDev, exist = TRUE, ...)
    }
    return(A)
}
