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
# copy and simplify from Morpho::projRead() for convenience
# project points onto the closest point on a mesh
project <- function (lm, mesh, sign = TRUE) {
    data <- vcgClost(lm, mesh, smoothNormals = FALSE, sign = sign, borderchk = FALSE)
    return(data)
}

plotLandmarks <- function(d1, pts, Sp, Tx, idx_pts, grDev){
    rgl.set(d1)
    rgl.pop("shapes", Sp[idx_pts])
    rgl.pop("shapes", Tx[idx_pts])
    grDev$vSp[idx_pts]<-spheres3d(pts,color = "green",alpha=0.5,radius=grDev$spradius)
    grDev$vTx[idx_pts]<-text3d(pts,texts=as.character(idx_pts),col="red",cex=grDev$tcex)

    return(grDev)
}

#############################################################
submesh<-function(spec,keep){

    # Fonction qui construit un sous mesh ? partir de l'objet mesh3d "mesh" selon les points ? conserver contenus dans le vecteur logique keep
    # Retourne un objet mesh3d

    spec2<-list()

    # extraction des vertices
    spec2$vb<-spec$vb[,keep]

    # extraction des edges
    idxV<-is.element(spec$it,which(keep))
    idxV<-matrix(idxV,dim(spec$it)[1],dim(spec$it)[2],byrow=FALSE)
    M<-spec$it[,which(colSums(idxV)==3)]
    spec2$it<-matrix(match(M,which(keep)),dim(M)[1],dim(M)[2])

    # extraction des normales
    spec2$normals<-spec$normals[,keep]

    # extraction du "material"
    spec2$material<-spec$material
    if (length(spec$material)>0){
        spec2$material$color<-spec$material$color[,which(colSums(idxV)==3)]
    }

    # d?finition de la class de spec2
    class(spec2) <- "mesh3d"

    return(spec2)
}



