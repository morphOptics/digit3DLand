# utils.R
#############################################################
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
