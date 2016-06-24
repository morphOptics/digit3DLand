#' @title digitFixed
#' @description XXXX
#' @details cccccc
#' @param specFull full resolution mesh3d object
#' @param specDecim decimated mesh3d object
#' @param decim 0-1 number setting the amount of reduction relative to existing face number to decimate the full resolution mesh. Used for multiresolution and only if specDecim is NULL. (see \code{\link[Rvcg]{vcgQEdecim}})
#' @param fixed number of landmarks to digitalize
#' @param index Digitalization sequence of landmarks
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
                        orthoplane = TRUE, percDist = 0.15,
                        grDev = list(windowRect = rbind(c(0,50,830,904), c(840,50,1672,904)),
                                     ptSize = 1, spradius = NULL, tcex=2), ...) {

    if (class(specFull) != "mesh3d")
        stop("specFull must have class \"mesh3d\".")

    if (is.null(specFull$material))
        specFull$material <- "gray"

    if (decim != 1) {
        cat("\nMesh decimation for multiresolution view ----------\n")
        specDecim <- vcgQEdecim(specFull, percent = decim)
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

    # Centering of the meshesb on the centroid of the decimated one
    tmp <- scale(t(specDecim$vb[-4, ]), scale = FALSE)
    specDecim$vb[-4, ] <- t(tmp)
    specFull$vb[-4, ] <- sweep(specFull$vb[-4, ], 1, attr(tmp, which="scaled:center"))

    # plot decimated mesh
    d1 <- Clear3d()
    par3d(windowRect = grDev$windowRect[1, ])
    shade3d(specDecim)
    grDev$dev <- rgl.cur()

    # plot of orthogonal planes: they are initialized as major axes of the mesh
    orthoplanes <- list(vInter=NULL, vPlanes = NULL)
    if (orthoplane) orthoplanes <- DrawOrthoplanes(specDecim)

    # Landmark selection - A is the individual configuration matrix [k x 3]
    A <- Adeci <- matrix(NA, fixed, 3, dimnames = list(1:fixed, c("x","y","z")))
    attr(A, which = "spec.name") <- deparse(substitute(specFull))

    Idx <- setdiff(1:fixed, index[idxPtsTemplate])
    for (i in 1:fixed){
        if (i <= length(idxPtsTemplate)){
            # Place 1st points require to adjust the template if it exist otherwise take all pts
            idx_pts <- index[idxPtsTemplate[i]]
            res <- SelectPoints3d(mesh = specDecim, A = A, IdxPts = idx_pts, grDev = grDev)
#??? why A is passed ?? no sense
            # zoom on full resolution mesh around the selected landmark
            res2 <- SetPtZoom(specFull=specFull, Pt = res$coords, IdxPts = idx_pts,
                              orthoplanes = orthoplanes, percDist = percDist, grDev=grDev)
            # landmark coordinate on the full resolution mesh
            A[idx_pts, ] <- res2$coords
            # Projection of landmarks on decimated mesh for graphics
            Adeci[idx_pts,] <- project(res2$coords, specDecim, trans = TRUE)

            # Graphics
            Sp[idx_pts] <- res$sp
            Tx[idx_pts] <- res$tx
            grDev <- plot.landmark(Adeci[idx_pts,], d1, Sp, Tx, idx_pts, grDev, ...)

            if(!is.null(templateFile) & i==length(idxPtsTemplate)){
                # all reference points of the template are placed => fits the template to the mesh
                # defines 3 matrices of configurations:
                # - configA : points placed on the mesh
                # - configB : points within the templateFile
                # - configC : points within the templateFile and in configA
                configA <- A[!is.na(A[,1]), , drop = FALSE]
                configB <- template$M
                configC <- configB[idxPtsTemplate, ]

                # transation & scaling of configA and configC
                transA <- apply(configA, 2, mean)
                AA <- sweep(configA, 2, transA)
                scaleA <- 1 / sqrt(sum(AA^2))
                AA <- AA * scaleA

                transB <- apply(configC, 2, mean)
                BB <- sweep(configC, 2, transB)
                scaleB <- 1/sqrt(sum(BB^2))
                BB <- BB * scaleB

                # rotation of configC on configA
                AB <- crossprod(AA, BB)
                sv <- svd(AB)
                sig <- sign(det(AB))
                sv$u[,3] <- sig * sv$u[,3]
                rot <- tcrossprod(sv$v, sv$u)

                # apply translation, rotation and scaling compute from configC to configB
                BB <- sweep(configB, 2, transB)
                BB <- BB * scaleB
                BB <- BB %*% rot

                # scaling and translation in order configB find a size and position comparable to configA
                BB <- BB / scaleA
                B <- sweep(BB, 2, transA, FUN = "+")

                # Copy of already placed landmarks in A into B
                B[!is.na(A[,1]), ] <- A[!is.na(A[,1]), ]
                ptsB <- project(B, specDecim)
                B <- project(B, specFull, trans = TRUE)

                # plot points/labels of B not placed before
                vv <- index[Idx]
                for (ii in 1:length(vv)){
                    grDev <- plot.landmark(t(ptsB$vb[1:3, vv[ii]], d1, Sp=NULL, Tx=NULL,
                                             vv[ii], grDev, color = "blue", col = "red")
                }
            }
        } else {
            # Selection of remaining landmarks (if any)
            idx_pts <- index[Idx[i-length(idxPtsTemplate)]]
            # zoom on full resolution mesh around the selected landmark
            res2 <- SetPtZoom(specFull=specFull, Pt = B[idx_pts, , drop = FALSE], IdxPts = idx_pts,
                              orthoplanes = orthoplanes, percDist = percDist, grDev=grDev)
            # landmark coordinate on full resolution mesh
            A[idx_pts, ] <- res2$coords
            # Projection of the landmark on the decimated mesh for graphics
            Adeci[idx_pts, ] <- project(res2$coords, specDecim, trans = TRUE)
            # plot
            grDev <- plot.landmark(Adeci[idx_pts, ], d1, grDev$vSp, grDev$vTx, idx_pts, grDev, ...)
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
        # zoom sur le mesh complet autour du point plac? + pla?age du point plus pr?cis
        res2 <- SetPtZoom(specFull=specFull, Pt = res$coords, IdxPts = idx_pts,
                          orthoplanes = orthoplanes, percDist = percDist, grDev =  grDev)
        # landmark coordinate on the full resolution mesh
        A[idx_pts, ] <- res2$coords
        # Projection of the landmark on the decimated mesh for graphics
        Adeci[idx_pts, ] <- project(res2$coords, specDecim, trans = TRUE)
        # Graphics
        Sp[idx_pts] <- res$sp
        Tx[idx_pts] <- res$tx
        grDev <- plot.landmark(Adeci[idx_pts, ], d1, Sp, Tx, idx_pts, grDev, ...)
    }
    return(A)
}







