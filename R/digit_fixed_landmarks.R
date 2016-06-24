#' @title digitFixed
#' @description XXXX
#' @details cccccc
#' @param specFull full resolution mesh3d object
#' @param specDecim decimated mesh3d object
#' @param decim 0-1 number setting the amount of reduction relative to existing face number to decimate the full resolution mesh. Used for multiresolution and only if specDecim is NULL. See vcg:::vcgQEdecim
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
DigitFixed <- function (specFull, specDecim = NULL, decim = 0.5, fixed, index = 1:fixed,
                        templateFile = NULL, idxPtsTemplate,
                        orthoplane = TRUE, percDist = 0.15,
                        grDev = list(windowRect = rbind(c(0,50,830,904), c(840,50,1672,904)),
                                     ptSize = 1, spradius = NULL, tcex=2), ...) {

    if (class(specFull) != "mesh3d")
        stop("specFull must have class \"mesh3d\".")

    if (is.null(specFull$material))
        specFull$material <- "gray"

    if (is.null(specDecim)) {
        cat("\nMesh decimation for multiresolution view ----------\n")
        specDecim <- vcgQEdecim(specFull, percent = decim)
        cat("---------------------------------------------------\n")
    }

    if (is.null(specDecim$material))
        specDecim$material <- "gray"

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

    # Centering of the mesh (not anymore an option, simplify back tracking of translations)
    tmp <- scale(t(specDecim$vb[-4, ]), scale = FALSE)
    specDecim$vb[-4, ] <- t(tmp)
    specFull$vb[-4, ] <- sweep(specFull$vb[-4, ], 1, attr(tmp, which="scaled:center"))

    # plot decimated mesh
    d1 <- Clear3d()
    par3d(windowRect = grDev$windowRect[1, ])
    # Don't need to plot all vertices for the decimated mesh (important is the full resolution)
    # ids1 <- plot3d(specDecim$vb[1, ], specDecim$vb[2, ], specDecim$vb[3, ],
    #                size = grDev$ptSize, aspect = FALSE,
    #                axes = F, box = F, xlab="", ylab="", zlab="")
    shade3d(specDecim)
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
            dd <- sqrt(colSums((specFull$vb[1:3,] - Pt)^2))

            # zoom on full resolution mesh around the selected landmark
            res2 <- SetPtZoom(dd, specFull=specFull, Pt=Pt, IdxPts = idx_pts,
                              orthoplanes = orthoplanes, percDist = percDist, grDev=grDev)
            # A coordinates
            A[idx_pts, ] <- res2$coords + res2$Trans2
            # Projection of landmarks on decimated mesh for graphics
            Adeci[idx_pts,] <- project(A[idx_pts, ,drop=FALSE], specDecim, trans = TRUE)
            # plot
            grDev <- plot.landmark(Adeci[idx_pts,], d1, Sp, Tx, idx_pts, grDev, ...)

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
                B[!is.na(A[,1]), ] <- A[!is.na(A[,1]), ]
                ptsB <- project(B, specDecim)
                # for (ii in 1:nrow(B)){ #DEBUG
                #    spheres3d(B[ii, ], color = "orange", alpha=0.5,radius=4*grDev$spradius)
                #}
                #sweep(B, 2, Trans - res2$Trans2)
                B <- project(B, specFull, trans = TRUE)
                # for (ii in 1:nrow(B)){ #DEBUG
                #    spheres3d(B[ii, ], color = "red", alpha=0.5,radius=4*grDev$spradius)
                #}
                # plot points/labels of B not placed before
                vv <- index[Idx]
                for (ii in 1:length(vv)){
                    #grDev <- plot.landmark(t(ptsB$vb[1:3, vv[ii]], d1, Sp, Tx, vv[ii], grDev, ...)
                    grDev$vSp[vv[ii]] <- spheres3d(t(ptsB$vb[1:3, vv[ii]]),
                                                 color = "blue", alpha=0.5, radius = grDev$spradius)
                    grDev$vTx[vv[ii]] <- text3d(t(ptsB$vb[1:3, vv[ii]]),
                                              texts=as.character(vv[ii]), col="red",
                                              cex = grDev$tcex, adj = rep(grDev$spradius, 2))
                }
            }
        } else {
            # Selection of remaining landmarks (if any)
            idx_pts <- index[Idx[i-length(idxPtsTemplate)]]
            # distances full resolution mesh to automatic template landmark
            Pt <- B[idx_pts, , drop = FALSE] #sweep(, 2, Trans, "+")
            dd <- sqrt(colSums((sweep(specFull$vb[1:3,], 1, Pt))^2))

            # zoom on full resolution mesh around the selected landmark
            res2 <- SetPtZoom(dd, specFull=specFull, Pt=Pt, IdxPts = idx_pts,
                              orthoplanes = orthoplanes, percDist = percDist, grDev=grDev)
            # A coordinates
            A[idx_pts, ] <- res2$coords + res2$Trans2
            # Projection of landmarks on decimated mesh for graphics
            Adeci[idx_pts,] <- project(t(A[idx_pts, ]), specDecim)$vb[1:3]
            # plot
            grDev <- plot.landmark(Adeci[idx_pts,], d1, grDev$vSp, grDev$vTx, idx_pts, grDev, ...)
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
        dd <- sqrt(colSums((specFull$vb[1:3,] - Pt)^2))

        # zoom sur le mesh complet autour du point plac? + pla?age du point plus pr?cis
        res2 <- SetPtZoom(dd, specFull=specFull, Pt = Pt, IdxPts=idx_pts,
                          orthoplanes = orthoplanes, percDist = percDist, grDev =  grDev)
        # A coordinates
        A[idx_pts,] <- res2$coords + res2$Trans2
        # Projection of landmarks on decimated mesh for graphics
        Adeci[idx_pts,] <- project(t(A[idx_pts,]), specDecim, sign=FALSE)$vb[1:3]
        # plot
        grDev <- plot.landmark(Adeci[idx_pts,], d1, Sp, Tx, idx_pts, grDev, ...)

    }
    return(A)
}
