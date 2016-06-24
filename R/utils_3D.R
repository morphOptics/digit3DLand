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
#' @title plot of one landmark in a 3D scene
#' @description Function takes a landmark and plots it onto the surface of a 3D mesh
#' @param landmark a 3D landmark
#' @param d1 a rgl scene
#' @param Sp xxx
#' @param Tx xxxx
#' @param idx_pts index of the landmark
#' @param grDev xxxx
plot.landmark <- function(landmark, d1, Sp, Tx, idx_pts, grDev, ...){
    if (length(landmark) != 3)
        stop("landmark should be a xyz point")

    # Graphical parameters - default
    alpha <- 0.5
    color <- "green"
    col <- "red"
    # Modify them in they are in the optional args
    argin <- list(...)
    if (length(argin)) {
        if (!("alpha" %in% names(argin))) alpha <- 0.5 #argin$alpha
        if (!("color" %in% names(argin))) color <- "green" #argin$color
        if (!("col" %in% names(argin))) col <- "red" #argin$col
    }
    # plot
    rgl.set(d1)
    rgl.pop("shapes", Sp[idx_pts])
    rgl.pop("shapes", Tx[idx_pts])
    grDev$vSp[idx_pts] <- spheres3d(landmark, color = color, alpha = alpha, radius = grDev$spradius)
    grDev$vTx[idx_pts] <- text3d(landmark, texts = as.character(idx_pts), col = col,
                                 cex = grDev$tcex, adj = rep(grDev$spradius, 2))

    return(grDev)
}
#' @title submesh
#' @description Extracts a submesh
#' @details Function build a submesh from the mesh3d object according to some kept vertices
#' @param mesh a mesh3d object
#' @param subset expression indicating vertices to keep
#' @param select expression indicating fields to select.
#' @return Return of mesh3d object
#' @export
subset.mesh <- function(mesh, subset, select=NULL) {
    if (class(mesh) != "mesh3d") stop("mesh should be a 'mesh3d' object")
    if (missing(subset))
        stop("'subset' not defined and without default value")
    else if (!is.logical(subset))
        stop("'subset' must be logical")

    subMesh <- list()
    idx_subset <- which(subset)
    # extraction of vertices and normals
    subMesh$vb <- mesh$vb[, subset]
    if (is.null(select)) {
        subMesh$normals <- mesh$normals[, subset]
        # extraction of faces
        idxV <- is.element(mesh$it, idx_subset)
        idxV <- matrix(idxV, nrow=dim(mesh$it)[1], ncol=dim(mesh$it)[2])
        idx <- (colSums(idxV) == 3)
        M <- mesh$it[, idx]
        subMesh$it <- matrix(match(M, idx_subset), nrow=dim(M)[1], ncol=dim(M)[2])
        # Optional extraction of "material"
        if (!is.null(mesh$material)) {
            subMesh$material <- mesh$material
            if (is.list(mesh$material)){ #per vertex color attribute
                subMesh$material$color <- mesh$material$color[, idx]
            }
        }
    }
    class(subMesh) <- "mesh3d"
    return(subMesh)
}



