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
DigitFixed <- function (specFull, decim=0.5, fixed, idxFixed = 1:fixed, templateCoord = NULL, idxTemplate = NULL,
                        GrOpt=setDigitFixedOptions()) {

    if (!(any(class(specFull) == "mesh3d")))
        stop("specFull must have class \"mesh3d\".")
    spec.name <- deparse(substitute(specFull))

    # check which setting of GrOpt$PCplanesDraw is called, and set idxPlanes consequently
    if(is.logical(GrOpt$PCplanesOptions$PCplanesDraw)){
        if (GrOpt$PCplanesOptions$PCplanesDraw){
            idxPlanes<-1:3
        }else{
            idxPlanes<-NULL
        }
    }else{
        V<-c("pc2-pc3","pc1-pc3","pc1-pc2")
        idxPlanes<-which(is.element(V,tolower(GrOpt$PCplanesOptions$PCplanesDraw)))
    }

    # Correction if mesh has non-manifold faces (ie faces made of non-manifold edges, ie edges shared by more than 2 faces)
    # Correction needed for the ordering of the intersection points among mesh and planes
    specFull<-vcgClean(specFull,sel=2)

    # decimation
    if (decim != 1) {
        cat("\nMesh decimation for multiresolution view ----------\n")
        specDecim <- vcgQEdecim(specFull, percent = decim)
        specDecim <- vcgUpdateNormals(specDecim, silent = TRUE)
        # Correction if mesh has non-manifold faces
        specDecim<-vcgClean(specDecim,sel=2)
        cat("---------------------------------------------------\n")
    } else {
        specDecim <- specFull
    }
    # material<-"gray" ??? au lieu de material$color ???
    # if (is.null(specFull$material)) {
    #     specFull$material <- specDecim$material <- "gray"
    # } else {
    #     specDecim$material <- "gray" # We may want to fix this by transfering vertex color etc
    # }

    if (missing(fixed)) {
        stop("missing number of landmarks to digitalize")
    } else {
        if (!is.numeric(fixed) || length(fixed) > 1)
            stop("fixed must be a single number")
    }

    # Use or not of a template
    if (is.null(templateCoord)){
        idxTemplate <- 1:fixed
    } else {
        if (missing(idxTemplate)) {
            idxTemplate <- 1:4
            warning("idxTemplate was missing.
                    First 4 landmarks will be used to align the template")
        } else {
            if (length(idxTemplate) < 4)
                stop("idxTemplate must contain at least 4 landmarks")
        }
        p1 <- length(idxTemplate)
        template <- list()
        template$M <- templateCoord
    }

    # Define default values for graphics interactivity
    grDev<-GrOpt
    grDev$vSp <- grDev$vTx <- Sp <- Tx <- rep(NA, fixed)
    grDev$spradius <- GrOpt$spheresOptions$spheresRad
    tmp <- diff(apply(specDecim$vb[1:3,], 1, range))
    grDev$spradius[,1] <- GrOpt$spheresOptions$spheresRad[,1] * min(tmp)
    grDev$labadj <- GrOpt$labelOptions$labelAdj * min(tmp)

    # Centering of the meshes on the centroid of the decimated one
    tmp <- scale(t(specDecim$vb[-4, ]), scale = FALSE)
    specDecim$vb[-4, ] <- t(tmp)
    specFull$vb[-4, ] <- sweep(specFull$vb[-4, ], 1, attr(tmp, which="scaled:center"))

    # plot decimated mesh
    d1 <- Clear3d()
    par3d(windowRect = grDev$winOptions$winSize[1, ])
    if (grDev$winOptions$winNb==1){
        d1 <- currentSubscene3d()
        layout3d(t(c(1,2)), sharedMouse = grDev$winOptions$winSynchro)
        next3d()
    }
    if (grDev$meshOptions$meshShade[1]){
        shade3d(specDecim, col=grDev$meshOptions$meshColor[1], alpha=grDev$meshOptions$meshAlpha[1])
    }
    if (grDev$meshOptions$meshWire[1]){
        wire3d(specDecim, col=grDev$meshOptions$meshColor[1], alpha=grDev$meshOptions$meshAlpha[1])
    }
    if (grDev$meshOptions$meshPoints[1]){
        points3d(t(specDecim$vb[1:3,]), col=grDev$meshOptions$meshColor[1], alpha=grDev$meshOptions$meshAlpha[1])
    }
    grDev$dev <- rgl.cur()

    # Rotate the scene
    R <- rotMajorAxes(specDecim$vb[1:3, ])
    par3d(userMatrix = R)

    # plot of orthogonal planes: they are initialized as major axes of the mesh
    orthoplanes <- list(vInter=NULL, vPlanes = NULL)
    if (length(idxPlanes)>0){
        orthoplanes <- DrawOrthoplanes(mesh=specDecim, idxPlanes=idxPlanes, grDev=grDev)
        if (ncol(specDecim$vb)!=ncol(specFull$vb)){
            # computation of intersections among full mesh and fixed planes
            orthoplanes <- DrawOrthoplanes(mesh=specFull,idxPlanes=idxPlanes,planes=orthoplanes$vPlanes,
                                           interactive=FALSE,is.plot=FALSE, grDev=grDev)
        }
    }

    # Landmark selection - A is the individual configuration matrix [k x 3]
    A <- Adeci <- matrix(NA, fixed, 3, dimnames = list(1:fixed, c("x","y","z")))
    attr(A, which = "spec.name") <- spec.name

    Idx <- setdiff(1:fixed, idxFixed[idxTemplate])
    for (i in 1:fixed){
        if (i <= length(idxTemplate)){
            # Place 1st points require to adjust the template if it exist otherwise take all pts
            idx_pts <- idxFixed[idxTemplate[i]]
            res <- SelectPoints3d(mesh = specDecim, A = A, IdxPts = idx_pts, grDev = grDev, whichMesh = 1)

            Pt <- res$coords
            grDev$vSp[idx_pts] <- res$sp
            grDev$vTx[idx_pts] <- res$tx
        } else {
            # Selection of remaining landmarks (if any)
            idx_pts <- idxFixed[Idx[i-length(idxTemplate)]]
            #Pt <- B[idx_pts, , drop = FALSE]
            Pt <- B[idx_pts, ]
        }
        # zoom on full resolution mesh around the selected landmark
        Pt2 <- c(project(t(Pt), specFull, trans = TRUE)) # added
        res2 <- SetPtZoom(specFull=specFull, Pt = Pt2, IdxPts = idx_pts,
                          orthoplanes = orthoplanes, idxPlanes=idxPlanes, grDev=grDev)
        grDev<-res2$grDev
        # landmark coordinate on the full resolution mesh
        A[idx_pts, ] <- res2$coords
        # Projection of landmarks on decimated mesh for graphics
        # Adeci[idx_pts,] <- project(res2$coords, specDecim, trans = TRUE)
        Adeci[idx_pts,] <- specDecim$vb[1:3,which.min(colSums(abs(specDecim$vb[1:3,]-c(res2$coords))))]

        # Graphics
        grDev <- plot.landmark(Adeci[idx_pts, ], d1, idx_pts, grDev, exist = TRUE)

        if(!is.null(templateCoord) & i==length(idxTemplate)){
            # all reference points of the template are placed
            # => impute missing landmarks
            B <- imputeCoords(A, template = template$M) #idx = idxTemplate
            ptsB <- project(B, specDecim)
            B <- project(B, specFull, trans = TRUE)

            # plot points/labels of B not placed before
            vv <- idxFixed[Idx]
            for (ii in 1:length(vv)){
                grDev <- plot.landmark(t(ptsB$vb[1:3, vv[ii]]), d1, vv[ii], grDev, exist = FALSE)
            }
            rgl.viewpoint(userMatrix = R)
        }
    }
    # Now, wait if the user want changed any landmark. Stop when the graphics is closed
    Stop <- 0
    grDev$dev <- d1
    while((Stop==0)){
        # clicks point on the decimated mesh
        res <- SelectPoints3d(mesh = specDecim, modify = TRUE, A = Adeci, grDev = grDev, whichMesh = 1)
        if (res$isClosed) break

        idx_pts <- res$Idx
        # zoom on full resolution mesh
        Pt2 <- c(project(t(res$coords), specFull, trans = TRUE)) # added
        res2 <- SetPtZoom(specFull=specFull, Pt = Pt2, IdxPts = idx_pts,
                          orthoplanes = orthoplanes, idxPlanes=idxPlanes, grDev =  grDev)
        grDev<-res2$grDev
        # landmark coordinate on the full resolution mesh
        A[idx_pts, ] <- res2$coords
        # Projection of the landmark on the decimated mesh for graphics
        # Adeci[idx_pts, ] <- project(res2$coords, specDecim, trans = TRUE)
        Adeci[idx_pts,] <- specDecim$vb[1:3,which.min(colSums(abs(specDecim$vb[1:3,]-c(res2$coords))))]

        # Graphics
        grDev$vSp[idx_pts] <- res$sp
        grDev$vTx[idx_pts] <- res$tx
        grDev <- plot.landmark(Adeci[idx_pts, ], d1, idx_pts, grDev, exist = TRUE)
    }
    return(A)
}
