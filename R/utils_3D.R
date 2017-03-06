# utils_3D.R
# copy from rgl mouseTrackball() - see rgl mouseCallbacks.R
# utilities used by selectPlanes()
xprod <- function(a, b){
    c(a[2]*b[3] - a[3]*b[2],
      a[3]*b[1] - a[1]*b[3],
      a[1]*b[2] - a[2]*b[1])
}

vlen <- function(a) sqrt(sum(a^2))

angle <- function(a,b) {
    dot <- sum(a*b)
    acos(dot/vlen(a)/vlen(b))
}
#
Clear3d <- function(type = c("shapes", "bboxdeco", "material"), defaults = getr3dDefaults(), subscene = 0) {
    d <- .check3d()
    rgl.clear(type, subscene = subscene)
    return(d)
}
# copied and simplified from Morpho::projRead() for convenience
# project points onto the closest point on a mesh
project <- function (lm, mesh, sign = TRUE, trans = FALSE) {
    data <- vcgClost(lm, mesh, smoothNormals = FALSE, sign = sign, borderchk = FALSE)
    if (trans) data <- t(data$vb[1:3, ])
    return(data)
}

rotMajorAxes <- function(mat) {
    if (dim(mat)[1] == 3) {
        u <- svd(mat)$u
    } else {
        u <- svd(mat)$v
    }
    R <- tcrossprod(u)
    if (det(R) < 0) {
        u[, 3] <- -1 * u[, 3]
        R <- tcrossprod(u)
    }
    tmp<-diag(rep(1,4))
    tmp[1:3,1:3]<-R
    return(tmp)
}

#' imputeCoords
#' @description Function imputes position of landmarks given a full template
#' and a subset of these landmarks
#'
#' @param landmark matrix with NA or position of landmarks on the specimen
#' @param landmark coordinates on the template
#'
#' @return A matrix of imputed landmark coordinates
#' @export
#'
imputeCoords <- function(A, template) {
    # defines 3 matrices of configurations:
    # - configA : points placed on the mesh
    # - configB : points within the templateFile
    # - configC : points within the templateFile and in configA
    configA <- A[!is.na(A[, 1]), , drop = FALSE]
    configB <- template
    configC <- configB[!is.na(A[, 1]), , drop=FALSE] #idx

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
    sv$u[, 3] <- sig * sv$u[,3]
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
    return(B)
}
#' @title plot.landmark
#' @description Function takes a landmark and plots it onto the surface of a 3D mesh
#' @param landmark a 3D landmark
#' @param d1 a rgl scene
#' @param idx_pts index of the landmark
#' @param grDev graphical parameters spradius, ctex, and vSp, vTx for object ID
#' in the rgl scene (as returned by \code{\link[Rvcg]{spheres3d}} or
#' \code{\link[Rvcg]{text3d}}) and used by \code{\link[Rvcg]{rgl.pop}} to modify the landmark
#' if it already exist in the scene
#' @param exist logical whether or not the landmark has existing labels in the scene
#' @param ... optional graphical arguments (color and alpha for 3dspheres, col for 3dtext)
plot.landmark <- function(landmark, d1, idx_pts, grDev, exist = FALSE,...){
    if (length(landmark) != 3)
        stop("landmark should be a xyz point")

    # Graphical parameters - default
    alpha <- 0.5
    color <- "cyan"
    col <- "magenta"
    # Modify them in they are in the optional args
    argin <- list(...)
    if (length(argin)) {
        if ("alpha" %in% names(argin)) alpha <- argin$alpha
        if ("color" %in% names(argin)) color <- argin$color
        if ("col" %in% names(argin)) col <- argin$col
    }
    # plot
    rgl.set(d1)
    if (exist) {
        rgl.pop("shapes", grDev$vSp[idx_pts])
        rgl.pop("shapes", grDev$vTx[idx_pts])
    }
    grDev$vSp[idx_pts] <- spheres3d(landmark, color = color, alpha = alpha, radius = grDev$spradius)
    grDev$vTx[idx_pts] <- text3d(landmark, texts = as.character(idx_pts), col = col,
                                 cex = grDev$tcex, adj = rep(grDev$spradius, 2))

    return(grDev)
}
#' @title subset.mesh3d
#' @description Extracts a submesh
#' @details Function build a submesh from the mesh3d object according to some kept vertices
#' @param mesh a mesh3d object
#' @param subset expression indicating vertices to keep
#' @param select expression indicating fields to select.
#' @return Return of mesh3d object
#' @export
subset.mesh3d <- function(mesh, subset, select=c("vb", "normals", "it", "material")) {

    if (missing(subset)) {
        stop("'subset' not defined and without default value")
    } else {
        if (!is.logical(subset))
            stop("'subset' must be logical")
    }

    subMesh <- list()
    if ("vb" %in% select)
        subMesh$vb <- mesh$vb[, subset]

    if (any(c("norm", "normals") %in% select))
        subMesh$normals <- mesh$normals[, subset]

    if (any(c("face", "faces", "it") %in% select)) {
        idx_subset <- which(subset)
        idxV <- is.element(mesh$it, idx_subset)
        idxV <- matrix(idxV, nrow=3, ncol=dim(mesh$it)[2])
        idx <- (colSums(idxV) == 3)
        subMesh$it <- matrix(match(mesh$it[, idx], idx_subset), nrow=3, ncol=sum(idx))
    }

    if (any(c("mat", "material") %in% select)) {
        if (!is.null(mesh$material)) {
            subMesh$material <- mesh$material
            if (is.list(mesh$material) && any(c("face", "faces", "it") %in% select)) {
                subMesh$material$color <- mesh$material$color[, idx] #per vertex color attribute
            }
        } else subMesh$material <- "gray"
    }
    class(subMesh) <- "mesh3d"
    return(subMesh)
}



