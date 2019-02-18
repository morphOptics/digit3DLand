# utils_3D ----
# copy from rgl mouseTrackball() - see rgl mouseCallbacks.R
# utilities used by selectPlanes()
xprod <- function(a, b){
    c(a[2]*b[3] - a[3]*b[2],
      a[3]*b[1] - a[1]*b[3],
      a[1]*b[2] - a[2]*b[1])
}

vlen <- function(a) sqrt(sum(a^2))

angle <- function(a, b) {
    dot <- sum(a*b)
    acos(dot/vlen(a)/vlen(b))
}

# Cleaning 3D scene ----
Clear3d <- function(type = c("shapes", "bboxdeco", "material"),
                    defaults = getr3dDefaults(), subscene = 0) {
    d <- .check3d()
    rgl.clear(type, subscene = subscene)
    return(d)
}

# Utilities for rotating, projecting, imputing and plotting 3D coordinates ----
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
    tmp <- diag(rep(1,4))
    tmp[1:3, 1:3] <- R
    return(tmp)
}

###########################
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
    B[!is.na(A[, 1]), ] <- A[!is.na(A[, 1]), ]
    return(B)
}

###########################
plot.landmark <- function(landmark, d1, idx_pts, grDev, exist = FALSE,...){
    if (length(landmark) != 3)
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

    if (exist) {
        rgl.pop("shapes", grDev$vSp[idx_pts])
        rgl.pop("shapes", grDev$vTx[idx_pts])
    }
    grDev$vSp[idx_pts] <- spheres3d(landmark, color = color, alpha = alpha, radius = rad)
    grDev$vTx[idx_pts] <- text3d(landmark, texts = as.character(idx_pts), col = col,
                                 cex = cex, adj = adj)

    return(grDev)
}

# Intersection between mesh and planes ----
meshPlaneIntersect2 <- function (mesh, v1, v2 = NULL, v3 = NULL, normal = NULL) {

    # modified function from Morpho:::meshPlaneIntersect
    # As for Morpho:::meshPlaneIntersect, this function computes the intersection points between
    # a mesh and a plane, but additionnally to return the intersection point coordinates (out), it
    # also returns infos on edges being intersected (edgesTot)

    # unchanged lines from meshPlaneIntersect: determine face numbers intersected by the plane
    pointcloud <- vert2points(mesh)
    updown <- cutSpace(pointcloud, v1 = v1, v2 = v2, v3 = v3, normal = normal)
    upface <- getFaces(mesh, which(updown))
    downface <- getFaces(mesh, which(!updown))
    nit <- 1:ncol(mesh$it)
    facesInter <- which(as.logical((nit %in% upface) * nit %in% downface))

    # unchanged lines from meshPlaneIntersect: define the mesh of intersection
    mesh$it <- mesh$it[, facesInter]
    mesh <- rmUnrefVertex(mesh, silent = TRUE)

    # Modified line from meshPlaneIntersect. It gives:
    # - the edges composing the mesh of intersection (2 1st columns of edgesTot)
    # - at each face the edge belongs (3rd column)
    # - if the edge is at the border of the mesh of intersection (4th column)
    # With the 2nd argument set to FALSE in vcgGetEdge() (contrary to the setting to TRUE
    # in the original function meshPlaneIntersect), all occurences of the edges are reported
    # when faces are browsed, meaning that for manifold mesh, edges can be appeared once (when
    # they belongs to mesh borders) or twice (in the other case).
    # This duplicated info will be needed later to compute how intersection points should be linked.
    # At this step, the edges at the intersection mesh borders correspond to the edges at the
    # borders in the whole mesh, but also to all edges delimiting that submesh. But by definition,
    # this second class of edges won't bear any intersection points (contrary to the first one which could).
    edgesTot <- as.matrix(vcgGetEdge(mesh, FALSE))

    # contrary to the "edges" in meshPlaneIntersect, the following "edges" is bigger (each non-border
    # edge being duplicated)
    edges <- edgesTot[, 1:2]

    # unchanged lines from meshPlaneIntersect: extract (duplicated) vertex coordinates from the mesh of intersection
    pointcloud <- vert2points(mesh)

    # Modified line from meshPlaneIntersect, calling edgePlaneIntersect2 rather than edgePlaneIntersect,
    # and extracting the edge indexes being intersected (idx_edges), in addition to the intersection coordinates (out),
    # both being duplicated
    tmp <- edgePlaneIntersect2(pointcloud, edges, v1 = v1, v2 = v2, v3 = v3, normal = normal)
    out <- tmp[, 1:3]
    idx_edges <- tmp[, 4]

    # additional line from MeshPlaneIntersect. Extracts the intersected edges from edgesTot.
    edgesTot <- edgesTot[idx_edges, ]


    return(list(out = out, edgesTot = edgesTot))
}

edgePlaneIntersect2 <- function (pointcloud, edges, v1, v2 = NULL, v3 = NULL, normal = NULL) {

    # modified function from Morpho:::edgePlaneIntersect
    # The only change corresonds to the call of the Rccp function"edgePlane2" rather than "edgePlane"

    if (!is.null(normal)) {
        tangent <- tangentPlane(normal)
        v2 <- v1 + tangent$z
        v3 <- v1 + tangent$y
    }
    e1 <- v2 - v1
    e2 <- v3 - v1
    e1 <- e1/sqrt(sum(e1^2))
    e2 <- e2/sqrt(sum(e2^2))
    normal <- crossProduct(e1, e2)
    normal <- normal/sqrt(sum(normal^2))
    e2a <- crossProduct(e1, normal)
    e2a <- e2a/sqrt(sum(e2a^2))
    Ep <- cbind(e1, e2a)
    edges <- edges - 1
    pointcloud0 <- sweep(pointcloud, 2, v1)
    orthopro <- t(Ep %*% t(Ep) %*% t(pointcloud0))
    diff <- orthopro - pointcloud0

    # modified line: calls "edgePlane2"
    out <- .Call("edgePlane2", pointcloud, diff, edges)

    return(out)
}
