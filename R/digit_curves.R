#' @title Mesh Digitization
#' @description Generic function for mesh digitization. It invokes 2 particular methods depending on the class of the
#'              1st argument \code{M}. So far method is implemented only for the mesh3d object.
#' @details For details, user should refer to the method \code{\link{digitCurve.mesh3d}} (to digitize a single
#'          \code{mesh3d} object). \cr
#'              \cr
#'              \strong{WARNING}: For Mac users, \code{digitCurves} is currently not compatible with the RStudio interface.
#'                                Please, use the basic R interface instead. A call from RStudio will cause an error
#'                                and an exit from the function... \cr
#' @param M Either a mesh3d object (in this case user should refer to \code{\link{digitCurves.mesh3d}}).
#' @param ... Additional arguments (all are not optional!) needed for curve digitization.
#'
#' @return A numerical matrix (\code{\link{digitCurves.mesh3d}}) with attributes.
#' @seealso \code{\link{digitCurves.mesh3d}}.
#' @export
#'
digitCurves <- function(M, ...){
    UseMethod("digitCurves", M)
}
#' Digitizing 3d curves
#'
#' @param specFull mesh3d object
#' @param coords matrix of 3d coordinates
#' @param curves
#' @param npts
#' @param nsLds
#' @param checkMesh
#' @param smoothMesh
#' @param GrOpt
#' @param verbose
#' @param spec.name
#' @param ... additional smoothing parameters (see \code{\link{vcgSmooth}} function)
#'
#' @return a matrix of landmark and semilandmark 3d coordinates plus semiLds and curve attributes
#' @export
#'
#' @examples
digitCurves.mesh3d <- function(specFull, coords, curves, npts = 1000, nsLds,
                               checkMesh = TRUE, smoothMesh = NULL, GrOpt = setGraphicOptions(),
                               verbose = c(TRUE, TRUE),
                               spec.name = NULL, ...){

    # check OS and R GUI to avoid graphic incompatibilities related to mac os and Rstudio
    tmp <- checkOsGui(GrOpt$winOptions$winNb, GrOpt$winOptions$winSynchro)
    GrOpt$winOptions$winNb <- tmp[[1]]
    GrOpt$winOptions$winSynchro <- tmp[[2]]

    # check mesh
    if (!(any(class(specFull) == "mesh3d")))
        stop("specFull must have class \"mesh3d\".")
    if (is.null(spec.name)){
        spec.name <- deparse(substitute(specFull))
    }

    # check verbose
    verbose <- checkLogical(verbose, c(1,2))

    # Correction if mesh has non-manifold faces (ie faces made of non-manifold edges, ie edges shared by more than 2 faces)
    # Correction needed for the ordering of the intersection points among mesh and planes
    if (checkMesh) {
        if (verbose[1]){
            cat("\nChecking mesh: starts...")
            if (verbose[2]){
                cat("\n")
            }
        }
        specFull <- vcgUpdateNormals(specFull, silent = !verbose[2])
        specFull <- vcgClean(specFull, sel=2, silent = !verbose[2])

        if (verbose[1]){
            if (!verbose[2]){
                cat("\r")
            }
            cat("Checking mesh: done!    \n")
        }
    }
    if (verbose[1]){
        cat("\nInitializations for digitMesh.mesh3d: in progress...")
    }

    # check which setting of GrOpt$PCplanesDraw is called
    # and set idxPlanes consequently
    if (is.logical(GrOpt$PCplanesOptions$PCplanesDraw)) {
        if (GrOpt$PCplanesOptions$PCplanesDraw) {
            idxPlanes <- 1:3
        } else {
            idxPlanes <- NULL
        }
    } else {
        V <- c("pc2-pc3", "pc1-pc3", "pc1-pc2")
        idxPlanes <- which(is.element(V, tolower(GrOpt$PCplanesOptions$PCplanesDraw)))
    }

    # smoothing surface
    if (!is.null(smoothMesh)){
        specFull <- vcgSmooth(specFull, type = smoothMesh, ...)
    }
    # Define default values for graphics interactivity
    grDev <- GrOpt
    # check if the mesh is actually colored
    if (!GrOpt$meshOptions$meshVertCol)
        specFull$material$color <- specDecim$material$color <- NULL
    if (grDev$meshOptions$meshVertCol & is.null(specFull$material$color)){
        grDev$meshOptions$meshVertCol <- FALSE
    }

    # curves
    if(is.vector(curves))
        curves <- matrix(curves, nrow = 1)

    # number of landmarks + ctrl_pts
    fixed <- nrow(coords)
    n_ctrl <- sum(curves[, 2])

    grDev$vSp <- grDev$vTx <- Sp <- Tx <- rep(NA, fixed+n_ctrl)
    grDev$spradius <- GrOpt$spheresOptions$spheresRad

    tmp <- diff(apply(specFull$vb[1:3,], 1, range))
    grDev$spradius[, 1] <- GrOpt$spheresOptions$spheresRad[, 1] * mean(tmp)
    grDev$labadj <- GrOpt$labelOptions$labelAdj * mean(tmp)

    # defined values for curves
    grDev$vSp_sLds <- rep(NA, sum(nsLds))
    grDev$vLn <- rep(NA, nrow(curves))

    # Centering of the mesh and coordinates
    tmp <- scale(t(specFull$vb[-4, ]), scale = FALSE)
    Trans1 <- attr(tmp, which = "scaled:center")
    coords <- sweep(coords, 2, Trans1)
    specFull$vb[-4, ] <- t(tmp)

    if (verbose[1]){
        cat("\rInitializations for digitMesh.mesh3d: done!         \n")
        cat("\nPlotting decimated mesh: in progress...")
    }

    # plot decimated mesh
    d1 <- Clear3d()
    par3d(windowRect = grDev$winOptions$winSize[1, ])
    if (grDev$winOptions$winNb==1){
        d1 <- currentSubscene3d()
        layout3d(t(c(1,2)), sharedMouse = grDev$winOptions$winSynchro)
        next3d()
    }
    if (grDev$meshOptions$meshShade[1]){
        shade3d(specFull, col = grDev$meshOptions$meshColor[1],
                alpha=grDev$meshOptions$meshAlpha[1], specular="black")
    }
    if (grDev$meshOptions$meshWire[1]){
        wire3d(specFull, col=grDev$meshOptions$meshColor[1],
               alpha=grDev$meshOptions$meshAlpha[1])
    }
    if (grDev$meshOptions$meshPoints[1]){
        points3d(t(specFull$vb[1:3,]), col=grDev$meshOptions$meshColor[1],
                 alpha=grDev$meshOptions$meshAlpha[1])
    }
    grDev$dev <- rgl.cur()

    # Rotate the scene
    R <- rotMajorAxes(specFull$vb[1:3, ])
    par3d(userMatrix = R)

    if (verbose[1]) cat("\rPlotting decimated mesh: done!         \n")

    # plot of orthogonal planes: they are initialized as major axes of the mesh
    orthoplanes <- list(vInter=NULL, vPlanes = NULL)
    if (length(idxPlanes) > 0){
        if (verbose[1]) cat("\nPlotting mesh/plane intersections: starts...\n")

        orthoplanes <- DrawOrthoplanes(mesh=specFull, idxPlanes=idxPlanes,
                                                     grDev=grDev, verbose=verbose)
        if (verbose[1]) cat("\n Plotting mesh/plane intersections: done!\n")
    }

    # Landmark selection - A is the individual configuration matrix [(k+ctrl_pts) x 3]
    A <- rbind(coords, matrix(NA, n_ctrl, 3))
    attr(A, which = "spec.name") <- spec.name

    # loop over landmarks to get their Ids in the scene
    for (idx_pts in 1:fixed){
        # project to get the ids
        grDev <- plot.landmark(A[idx_pts, ], d1 = d1, idx_pts = idx_pts, grDev = grDev, exist = FALSE)
        # Graphics
        grDev <- plot.landmark(A[idx_pts, ], d1 = d1, idx_pts= idx_pts, grDev = grDev, exist = TRUE)
    }
    # set idx_pts @ max
    idx_pts <- fixed
    ncurves <- nrow(curves)
    ctrl_pts <- list()

    if (length(nsLds) == 1 & ncurves > 1) nsLds <- rep(nsLds, ncurves)
    if (length(nsLds) != ncurves) stop("length of semilandmarks and number of curves don't match")

    curve_pts <- smLds <- list()
    idx_sLds <- list(1:nsLds[1])
    if (ncurves > 1){
        idx_sLds <- c(idx_sLds, lapply(2:ncurves, function(idx)
            (cumsum(nsLds)[idx-1]+1):cumsum(nsLds)[idx]))
    }

    if (verbose[1]){
        cat("\nLoop for semilandmark digitization: starts...")
        cat("\nLeft click to rotate, scroll wheel to zoom, (for mac users: cmd +) right click to position a landmark.\n")
    }

    verts <- t(specFull$vb[-4,])

    for (ii in 1:ncurves){
        if (verbose[1]) cat("\nCurve digitization: click ", curves[ii, 2], "pts please on curve
                        defined between landmarks:", curves[ii, c(1, 3)], "\n")
        for (jj in 1:curves[ii, 2]){ # ?do a while condition to make the curve adaptive ?how get out
            idx_pts <- idx_pts + 1
            if (jj == 1) idx0 <- idx_pts
            # select ctrl points on the full resolution mesh
            #res <- SelectPoints3d(mesh = specFull, A = A, IdxPts = idx_pts, grDev = grDev, whichMesh = 1)
            res <- SelectPoints3d(verts = verts, norms=specFull$normals, it=specFull$it, A = A, IdxPts = idx_pts, grDev = grDev, whichMesh = 1)
            A[idx_pts, ] <- res$coords

            # cubic spline
            ctrl_pts[[ii]] <- c(curves[ii, 1], idx0:idx_pts, curves[ii, 3])
            curve_pts[[ii]] <- splineCurve(A[ctrl_pts[[ii]], ], npts, mesh = specFull)

            # Resampled curve
            smLds[[ii]] <- curve_pts[[ii]][(seq(1, npts, length = nsLds[ii] + 2))[-c(1, nsLds[ii] + 2)], ]

            # Graphics
            grDev$vSp[idx_pts] <- res$sp
            grDev$vTx[idx_pts] <- res$tx
            grDev <- plot.landmark(A[idx_pts, ], d1=d1, idx_pts=idx_pts,
                                                 grDev=grDev, exist = TRUE)
            grDev <- plot.curve(curve_pts[[ii]], smLds[[ii]], d1 = d1, idx_curv = ii,
                                idx_sLds = idx_sLds[[ii]], grDev = grDev, exist = TRUE)

        }

    }


    if (verbose[1]){
        cat("\n\nLoop for landmark digitization: ends.\n\n")
        cat("Now, you can:")
        cat("\n - validate your digitization by closing the graphic device\n")
        cat("\n - or modify some landmarks if necessary (just click a on a landmark to modify and redigitize it,\n")
        cat("\n   once all landmraks are correct, just close the graphic device).\n")
    }

    # Now, wait if the user want changed any landmark. Stop when the graphics is closed
    Stop <- 0
    grDev$dev <- d1
    while((Stop==0)){
        # clicks point on the decimated mesh
        #res <- SelectPoints3d(mesh = specFull, modify = TRUE, A = A, grDev = grDev, whichMesh = 1)
        res <- SelectPoints3d(verts = verts, norms=specFull$normals, it=specFull$it, modify = TRUE, A = A, grDev = grDev, whichMesh = 1)
        if (res$isClosed){
            if (verbose[1]){
                cat("\n")
            }
            break
        }

        idx_pts <- res$Idx
        A[idx_pts, ] <- res$coords
        if (verbose[1]){
            cat("\n")
            txt<-paste0("Landmark to modify: ",idx_pts)
            cat(txt)
        }


        # Graphics
        grDev$vSp[idx_pts] <- res$sp
        grDev$vTx[idx_pts] <- res$tx
        grDev <- plot.landmark(A[idx_pts, ], d1, idx_pts, grDev, exist = TRUE)
        # is ctrl points or fixed landmarks ?
        if (idx_pts > fixed){
            # recompute curve and semilandmarks
            # get curve on which idx_pts lie
            idxcrv <- which(unlist(lapply(ctrl_pts, function(x) any(idx_pts %in% x))))
            pts <- A[ctrl_pts[[idxcrv]], ]
            curve_pts[[idxcrv]] <- splineCurve(pts, npts, mesh = specFull)
            smLds[[idxcrv]] <- curve_pts[[idxcrv]][(seq(1, npts, by = nsLds[idxcrv] + 2))[-c(1, nsLds[idxcrv] + 2)], ]

            # Graphics
            grDev <- plot.curve(curve_pts[[idxcrv]], smLds[[idxcrv]], d1 = d1, idx_curv = idxcrv,
                                idx_sLds = idx_sLds[[idxcrv]], grDev = grDev, exist = TRUE)

        }
        if (verbose[1]){
            cat("\r")
            txt<-paste0("Landmark ",idx_pts," has been modified.")
            cat(txt)
        }

    }

    # restore semilandmarks, curves and save
    smLds <- do.call("rbind", smLds)
    # ctrl <- A[-fixed, ]
    # curve_pts <- do.call("rbind", curve_pts)
    ctrl <- lapply(ctrl_pts, function(idx) A[idx, ] + matrix(Trans1, nrow = length(idx), ncol = 3, byrow = TRUE))
    A <- rbind(A[1:fixed, ], smLds)
    A <- A + matrix(Trans1, nrow = nrow(A), ncol = 3, byrow = TRUE)
    #ctrl <- ctrl + matrix(Trans1, nrow = nrow(ctrl), ncol = 3, byrow = TRUE)
    curve_pts <- lapply(curve_pts, function(x) x + matrix(Trans1, nrow = nrow(x), ncol = 3, byrow = TRUE))

    if (verbose[1])
        cat("\nMesh digitization is ended!\n")
    # add semilandmark indices as attributes
    attr(A, which = "semiLds") <- c(1:nrow(A))[-c(1:fixed)]
    # add curves as attributes
    attr(A, which = "curve") <- list(ctrl_pts = ctrl, curve = curve_pts)
    return(A)
}
# Chord distance for spline curves
cumchord <- function(M){
    cumsum(sqrt(apply((M - rbind(M[1, ], M[-nrow(M),]))^2, 1, sum)))
}
# cubic spline
spline3d <- function(pts, npts = 100, proj = FALSE, mesh){
    z <- cumchord(pts)
    zz <- seq(0, max(z), length.out = npts)
    sdata <- data.frame(
        x = splinefun(z, pts[, 1])(zz),
        y = splinefun(z, pts[, 2])(zz),
        z = splinefun(z, pts[, 3])(zz))
    if (proj) {
        sdata <- projRead(as.matrix(sdata), mesh)
        sdata <- vert2points(sdata)
    }
    sdata
}
#
splineCurve <- function(pts, npts, mesh){
    # need to sort the ctrl_pts, so clicks don't have to be ordered
    crv <- sortCurve(pts, start = 1)$xsorted
    # 3d spline reprojected on mesh
    crv <- spline3d(crv, npts = npts, proj = TRUE, mesh = mesh)
}

plot.curve <- function(curve, slandmark, d1, idx_curv, idx_sLds, grDev, exist = FALSE,...){
    if (ncol(slandmark) != 3)
        stop("landmark should be a xyz point")

    # Graphical parameters
    alpha <- grDev$spheresOptions$spheresAlpha[2, 1]
    color <- grDev$spheresOptions$spheresColor[2, 1]
    rad <- grDev$spradius[2, 1]
    col <- grDev$labelOptions$labelColor[2, 1]
    cex <- grDev$labelOptions$labelCex[2, 1]
    adj <- grDev$labelOptions$labelAdj[2, 1,]

    # plot
    if (grDev$winOptions$winNb == 2){
        rgl.set(d1)
    }else{
        subS<-rgl.ids(type = "subscene", subscene = 0)
        useSubscene3d(subS[2, 1])
    }

    if (exist & !is.na(grDev$vSp_sLds[idx_sLds[1]])) {
        rgl.pop("shapes", grDev$vLn[idx_curv])
    }
    grDev$vLn[idx_curv] <- lines3d(curve, col = "blue", lwd = 2)   # after

    if (exist & !is.na(grDev$vSp_sLds[idx_sLds[1]])) {
        for (ii in idx_sLds){
            rgl.pop("shapes", grDev$vSp_sLds[ii])
        }
    }
    for (ii in 1:length(idx_sLds)){
        grDev$vSp_sLds[idx_sLds[ii]] <- spheres3d(slandmark[ii, ], col = 'cyan', radius = grDev$spradius[1]/3)
    }
    return(grDev)
}
