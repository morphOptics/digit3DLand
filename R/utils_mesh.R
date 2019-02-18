# Subsetting and color transfert ----
#' Color transfert
#' @title Transfers Colors between Meshes
#' @description This function transfers vertex colors (if any) from one mesh to a one.
#' @details This function transfers vertex colors (if any) from the original mesh to the decimated one following
#'          2 steps: \cr
#'          - converts of vertex colors from the first mesh \code{M1} to face colors by averaging the colors from
#'          the 3 vertices of each face, \cr
#'          - estimates the correspondence among vertices from the second mesh \code{M2} to the first mesh \code{M1}
#'          by means of the \code{\link[Rvcg]{vcgClostKD}} function).
#' @param M1 A \code{mesh3d} object to transfer vertex colors from.
#' @param M2 A \code{mesh3d} object to transfer vertex colors to.
#'
#' @return A decimated \code{mesh3d} object.
#' @export
#'
#' @examples
#'
colTransf <- function(M1, M2){

    if (class(M1) != "mesh3d" | class(M2) != "mesh3d"){
        stop("This function only accepts mesh3d objects...")
    }

    Col <- M1$material$color

    if (!is.null(Col)){
        # material$color should be present...
        if (is.matrix(Col)){
            #... be stored as a matrix...
            nCol <- ncol(Col)
            if (length(unique(Col[, sample(1:nCol, min(nCol, 500))])) > 1){
                # ...and contain different color values (test through a random sample of 500 faces)

                # converting vertex colors for M1 in face colors (=> Mcol)
                Acol <- array(col2rgb(Col), c(3, 3, nCol)) # converts hewadecimal colors to rgb values
                mcol <- (Acol[, 1, ] + Acol[, 2, ] + Acol[, 3, ]) / 3 # averaging vertex colors for each face to get a single value per face
                mcol <- round(mcol)
                mcol <- rgb(mcol[1, ], mcol[2, ], mcol[3, ], maxColorValue = 255) # converts back to hexadecimal values
                Mcol <- matrix(mcol, nrow = 3, ncol = nCol, byrow = TRUE)

                # computing correspondence among vertices of M2 and faces of M1
                kd <- vcgClostKD(M2, M1)

                # transferring face colors from M (Mcol) to decMesh
                nc <- ncol(M2$it)
                idxFaces <- matrix(kd$faceptr[M2$it], nrow = 3, ncol = nc) # indexes of faces from M1 matching with the vertices of M2
                M2$material$color <- matrix(Mcol[1, idxFaces], nrow = 3, ncol = nc)

            }
        }
    }

    return(M2)
}
# Subset Mesh
subset.mesh3d <- function(mesh, subset, select = c("vb", "normals", "it", "material")) {

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
        idxV <- matrix(idxV, nrow = 3, ncol = dim(mesh$it)[2])
        idx <- (colSums(idxV) == 3)
        subMesh$it <- matrix(match(mesh$it[, idx], idx_subset), nrow = 3, ncol = sum(idx))
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

# Mesh decimation utilities ----
#' Mesh Decimation
#' @description Generic function for mesh decimation. It invokes 2 particular methods depending on the class of the
#'              1st argument \code{M}.
#' @details For details, user should refer to the 2 methods \code{\link{decimMesh.mesh3d}} (to decimate a single
#'          \code{mesh3d} object) and \code{\link{decimMesh.character}} (to decimate a set of mesh files).
#' @param M Either a \code{mesh3d} object (in this case user should refer to \code{\link{decimMesh.mesh3d}}) or a
#'          character value (in this case user should refer to \code{\link{decimMesh.character}})
#' @param ... Additional arguments for mesh decimation.
#'
#' @return Either \code{mesh3d} object (\code{\link{decimMesh.mesh3d}}), or a list of \code{mesh3d} objects
#'         (\code{\link{decimMesh.character}})
#' @seealso \code{\link{decimMesh.mesh3d}} and \code{\link{decimMesh.character}}.
#' @export
#'
decimMesh<-function(M, ...){
    UseMethod("decimMesh", M)
}

#' @title Decimates a Mesh
#' @description Wrapper for the \code{\link[Rvcg]{vcgQEdecim}} function making additionnal mesh cleaning needed for
#'              \code{\link{digitMesh}}. Additionally to the \code{\link[Rvcg]{vcgQEdecim}} function, this function
#'              transfers vertex colors (if any) from the original mesh to the decimated one.
#' @details This function decimates a given mesh \code{M} with \code{\link[Rvcg]{vcgQEdecim}}, then updates mesh
#'          normals, removes non manifold faces in the decimated mesh and transfers vertex colors (if any) from the
#'          original mesh to the decimated one (done firstly by converting vertex colors from the original mesh to
#'          face colors, and then by computing the correspondence among vertices from the decimated mesh to the given
#'          one by means of the \code{\link[Rvcg]{vcgClostKD}} function).
#' @param M A \code{mesh3d} object.
#' @param tarface An integer numerical value for the targetted number of faces for decimation (if superior or equal
#'                to the number of faces in \code{M}, then no decimation is performed).
#' @param silent A logical value indicating if console output should be issued once the calls of
#'               \code{\link[Rvcg]{vcgQEdecim}}, \code{\link[Rvcg]{vcgUpdateNormals}} and
#'               \code{\link[Rvcg]{vcgClean}} are done.
#' @param ... Optional arguments for \code{\link[Rvcg]{vcgQEdecim}}.
#'
#' @return A decimated \code{mesh3d} object.
#' @seealso \code{\link{decimMesh.character}}.
#' @export
#'
#' @examples
#'
decimMesh.mesh3d<-function(M, tarface=15000, silent = FALSE, ...){
    # decimation
    decMesh <- vcgQEdecim(M, tarface, silent = silent, ...)
    # update normals
    decMesh <- vcgUpdateNormals(decMesh, silent = silent)
    # Correction if mesh has non-manifold faces
    decMesh<-vcgClean(decMesh, sel = 2, silent = silent)

    # color transfer (in case of non uniform color for M)
    decMesh <- colTransf(M, decMesh)
}

#' @title Decimates Several Meshes
#' @description Wrapper for the function \code{\link{decimMesh.mesh3d}} to decimate a set of mesh files contained
#'              within a directory.
#' @details This function is a wrapper for \code{\link{decimMesh.mesh3d}} (itself a wrapper for
#'              \code{\link[Rvcg]{vcgQEdecim}}), calling it to treat a list of full mesh files, which will be
#'              decimated the one after the other. It looks in the directory \code{sdir} for either a mesh filename
#'              specified by the character value stored in \code{M}, or a subdirectory whose name stored in \code{M}
#'              that contains the mesh files to decimate. By default, the level of decimation is set by the
#'              \code{target} argument, but other settings are possible (see \code{\link[Rvcg]{vcgQEdecim}}).
#'              Decimated mesh are then saved either in the same directory than full meshes, by adding
#'              a suffix in decimated mesh filenames (to distinguish them from full mesh filenames), or in a separate
#'              subfolder (contained in the directory for full meshes). Both ply and stl mesh files are supported
#'              (but decimated meshes will be uniquely saved as ply).
#' @usage
#' \method{decimMesh}{character}(M, tarface=15000, sdir=getwd(), patt=".ply", deci.suffix=NULL,
#'           deci.dir="DecimMesh", verbose = TRUE, \dots)
#' @param M A character value indicating either a mesh filename to decimated within \code{sdir}, or a directory name
#'          within \code{sdir} containing the subdirectory M with the mesh files to decimate.
#' @param tarface An integer numerical value for the targetted number of faces for decimation (if superior or equal
#'                to the number of faces in \code{M}, then no decimation is performed).
#' @param sdir A character value indicating the directoy name containg eiher the mesh filename \code{M} to decimate,
#'             or the subdirectory name with the mesh files to decimate.
#' @param patt A character value within \{\code{".ply"},\code{".stl"}\} indicating the kind of mesh file to open.
#' @param deci.suffix A character value indicating the suffix added at the end of the decimated mesh filenames and
#'                    before the \code{".ply"} exension. Set to \code{NULL} by default.
#' @param deci.dir A character value indicating the name of the subdirectory in \code{M} where the decimated mesh will
#'                 be saved. By default, \code{deci.dir} is set to \code{"DecimMesh"}, meaning that a subdirectory
#'                 named \code{"DecimMesh"} will be automatically created within \code{M} (if it doesn't exist yet).
#' @param verbose Possible settings are: \cr
#'                - a logical value: in this case this value should be recycled in a 2 length vector indicating
#'                  for 2 levels of verbose if comments should be printed or not on screen as the computations are
#'                  processed. The firs level corresponds to comments specific to the functions from the
#'                  \code{digit3DLand} library, and the second one to comments specific to the functions from the
#'                  \code{Rvcg} library. \cr
#'                - a 2-length logical vector standing for the 2 possible levels of verbose.
#'                a logical vector
#' @param ... Optional arguments for \code{\link[Rvcg]{vcgQEdecim}}).
#'
#' @return A list of decimated \code{mesh3d} objects.
#' @seealso \code{\link{decimMesh.mesh3d}}.
#' @export
#'
#' @examples
#'
#' ## Not run:
#'
#' # 1st example: decimating a mesh file
#' # A basic call consists in giving the filename of the mesh to decimate (contained in the working directory):
#' L <- decim("mesh2decimate.ply")
#'
#' # 2nd example: decimating mesh files in a given directory
#' # If the working directory contains a subdirectory named "fold" with the full resolution mesh files to decimate, a basic call of the function is:
#' L <- decim("fold")
#'
#' # 3rd example: decimating mesh files following a targetted face number percentage (see the help of vcgQEdecim for details):
#' L <- decim("fold", percent = 0.3)
#'
#' # 4th example: decimating stl files (contined in the "fold4stl" folder):
#' L <- decim("fold4stl", patt=".stl")
#'
#' # 5th example: saving decimating meshes in the same folder as the full meshes by sepcifying their suffix:
#' L <- decim("fold", deci.suffix="suffix4decimatedMesh", deci.dir=NULL)
#'
#' ## End(Not run)
#'
decimMesh.character<-function(M, tarface = 15000, sdir = getwd(), patt = ".ply",
                              deci.suffix = NULL, deci.dir = "DecimMesh",
                              verbose = c(TRUE, TRUE), ...) {
    curdir<-getwd()
    # check verbose
    verbose<-checkLogical(verbose, c(1, 2))

    if (verbose[1])
        cat("\nChecking arguments for decimMesh: in progress...")

    setDecimOptions(tarface = tarface, patt = patt)
    FiOpt <- setFileOptions(M = M, sdir = sdir, patt = patt,
                            deci.suffix = deci.suffix, deci.dir = deci.dir)

    full.files <- FiOpt$full.files
    full.dir <- FiOpt$full.dir
    deci.dir <- FiOpt$deci.dir

    ID <- unlist(strsplit(full.files, patt)) # id for individual

    if (verbose[1])
        cat("\rChecking arguments for decimMesh: done!         \n")

    Ldeci <- list()
    attr(Ldeci, "full.dir") <- full.dir
    attr(Ldeci, "deci.dir") <- deci.dir
    deci.files <- rep(NA, length(full.files))
    for (i in 1:length(full.files)) {
        if (verbose[1]){
            header <- paste0("---------- Mesh to decimate: ", full.files[i], " ----------")
            cat("\n", header, "\n")
            if (verbose[2]) cat("\n")
        }

        setwd(full.dir)
        full <- vcgImport(full.files[i], updateNormals = TRUE,
                          readcolor = TRUE, clean = TRUE, silent = !verbose[2])
        if (!is.null(tarface) & tarface > dim(full$it)[2]){
            tarface <- dim(full$it)[2]
        }
        deci <- decimMesh(full, tarface, silent = !verbose[2], ...)

        if (verbose[1])
            cat("\nExporting decimated mesh: in progress...")

        setwd(deci.dir)
        if (!is.null(deci.suffix)){
            deci.files[i] <- paste0(ID[i], deci.suffix)
        } else {
            deci.files[i] <- ID[i]
        }
        while(file.exists(paste0(deci.files[i], patt))) {
            deci.files[i] <- paste0(deci.files[i], "_deci")
        }

        attr(deci,"name") <- ID[i]
        attr(deci,"full.file") <- full.files[i]
        attr(deci,"deci.file") <- paste0(deci.files[i], patt)
        Ldeci[[i]] <- deci

        vcgPlyWrite(deci, filename = deci.files[i])

        if (verbose[1]){
            cat("\rExporting decimated mesh: done!         \n")
            cat(rep("-", nchar(header)), sep="")
            cat("\n")
        }

    }

    setwd(curdir)
    return(Ldeci)
}
