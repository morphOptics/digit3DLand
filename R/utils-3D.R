# utils_3D.R
Clear3d <- function(type = c("shapes", "bboxdeco", "material"), defaults = getr3dDefaults(), subscene = 0) {
    d <- .check3d()
    rgl.clear(type, subscene = subscene)
    return(d)
}
# copy and simplify from Morpho::projRead() for convenience
# project points onto the closest point on a mesh
project <- function (lm, mesh, sign = TRUE) {
    data <- vcgClost(lm, mesh, smoothNormals = FALSE, sign = sign, borderchk = FALSE)
    return(data)
}

