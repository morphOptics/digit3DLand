#' @title Mesh Digitization
#' @description Generic function for mesh digitization. It invokes 2 particular methods depending on the class of the
#'              1st argument \code{M}.
#' @details For details, user should refer to the 2 methods \code{\link{digitMesh.mesh3d}} (to digitize a single
#'          \code{mesh3d} object) and \code{\link{digitMesh.character}} (to digitize a set of mesh files). \cr
#'              \cr
#'              \strong{WARNING}: For Mac users, \code{digitMesh} is currently not compatible with the RStudio interface.
#'                                Please, use the basic R interface instead. A call from RStudio will cause an error
#'                                and an exit from the function... \cr
#' @param M Either a mesh3d object (in this case user should refer to \code{\link{digitMesh.mesh3d}}) or a character
#'          value (in this case user should refer to \code{\link{digitMesh.character}}).
#' @param ... Additional arguments (all are not optional!) needed for mesh digitization.
#'
#' @return Either a numerical matrix (\code{\link{digitMesh.mesh3d}}) or a numerical array
#'         (\code{\link{digitMesh.character}}).
#' @seealso \code{\link{digitMesh.mesh3d}} and \code{\link{digitMesh.character}}.
#' @export
#'
digitMesh <- function(M, ...){
    UseMethod("digitMesh", M)
}

#' @title Digitizes a Mesh
#' @description Interactive digitization of a mesh3d object.
#' @details This function allows to interactively digitize \emph{p}=\code{fixed} landmarks on the surface of a mesh,
#'          using two versions of this mesh: \cr
#'          - a decimated version (\code{specDecim}) to grossly positionned in a first time a given landmark, \cr
#'          - and then the full resolution version of the mesh (\code{specFull}) to finely positionned this landmark on
#'            a zoomed area around the previously positionned landmark \cr
#'          User should interact first on the decimated mesh, and then, landmark location will be validated on the
#'          full resoltution mesh.
#'
#'          The rationale to adopt a two steps approach to position landmarks through the use of a
#'          decimated mesh was motivated by several causes: \cr
#'          1/ The user interactions with the mesh (as the rotation, the zoom, and also the landmark location) can
#'             become time-consuming for heavy meshes, making not fluid at all the mesh digitizing. \cr
#'          2/ The automatic display of the zoomed zone allowed to question user on its landmark positionning if it
#'             has been made with a low magnfication. \cr
#'          3/ Because mesh translation is not possible with this function, the zoomed zone allows to digitize more
#'             easily landmarks distant from the mesh (and rotation) center.
#'
#'          \strong{a. Basic process for landmark digitization}
#'
#'          Mouse and keyboard interactions for landmark digitization: \cr
#'          - left click: mesh rotation \cr
#'          - scroll wheel: zoom and dezoom \cr
#'          - right click for Windows users, and cmd + right click for Mac users: digitize a landmark \cr
#'          - Escap key press: validate the landmark positionning (and then pass to the next landmark digitizing) \cr
#'
#'          In its basic version, the function process divides into 2 steps: \cr
#'          1/ Landmark digitizing \cr
#'          2/ Once all the landmarks are digitized, user can modify (or not) any of them as many time as wanted. \cr
#'
#'          During step 1/, user should first grossly position a landmark on the decimated mesh (step 1a/). While the
#'          landmark is not validated by the Escape key, the user is free to re-position as many time as wanted the
#'          landmark (rotation and zoom can be changed as well). Once the Escape key is pressed, the zoomed full
#'          resolution mesh appears (step 1b/), to finely positioned the landmark. As before, the user can change the
#'          landmark positioning, and will definitely validate it by the Escape key to digitize then the next landmark.
#'
#'          Once all the landmarks are positionned, the step 2/ allow user to modify (if necessary) any landmark.
#'          Just click on (or near) the landmark to modify and then, the process is the same as the one described
#'          for step 1/ (1a/ & 1b/). Once the user estimates that all landmarks are correctly positionned, The
#'          digitization of the mesh is validated by closing the graphic device (red cross). If no landmark needs
#'          modification, the device can be closed as soon as the step 1/ is done.
#'
#'          \strong{WARNING}: be carreful that the window should be closed once all landmark are validated, meaning that
#'          step 1b/ should be finished, and that no one landmark modification during step 2/ is in progress. To
#'          ensure of this, the zoomed mesh shouldn't be visible anymore, and all landmarks on the decimated mesh
#'          should have the same color (blue by default). Otherwise, the landmark coordinates won't be exported.
#'
#'          \strong{b. First refinement: using a template configuration}
#'
#'          A configuration matrix of landmark coordinates (\code{templateCoord}) can be used as a template to fasten
#'          the digitization process. In this case, user should first digitize few landmarks (set with
#'          \code{idxTemplate}) on the treated mesh, and once they are positionned, the function computes first rigid
#'          transformations (translation, scaling, rotation) to fit the template first landmarks onto the first
#'          digitized landmarks from the digitized mesh, then apply those transformations on the full template
#'          configuration, and finally project the template landmark onto the digitzed mesh.
#'
#'          Those projected landmarks from template on the mesh are expected (at least for small shape variablity cases,
#'          and with a judicious choice for the template configuration) to be positionned near the actual landmark to
#'          digitize. Consequently, the function uses this information of approximative position for the remainging
#'          landmarks as an assesment of the zone where the final landmarks should be positionned. Concretely,
#'          during the digitization process, the template projection of remaining landmark allows to automatically
#'          process the step 1a/ without need of user interaction, and the zoomed full mesh is automatically assessed
#'          for each landmark, and user needs only to finely position them on during step 1b/.
#'
#'          \strong{Note 1}: It can occur that the assessed zoomed mesh doesn't correspond to the actual zone where the
#'                           landmark should be positioned. If so, and because step 1b/ will process automatically each
#'                           landmark the ones after the others, you can simply incorrectly positioned this landmark,
#'                           and modify it later during step 2/.
#'
#'          \strong{Note 2}: The first few digitized landmarks used to fit the template on the mesh are decisive to
#'                           obtain a good positioning of the remainging landmarks. We advise user to choose at least 4
#'                           landmarks, sufficiently spaced each one with other and describing the whole mesh in its 3
#'                           dimensions (so avoiding to choose landmarks positionned in a single plane).
#'
#'          \strong{Note 3}: In its current version, the function doesn't allow to modify these first landmarks before
#'                           to fit the tempalte. So, if some of those landmarks are uncorrectly positionned, it can
#'                           highly impact the projection of the remaining landmarks and then the assessment of the
#'                           zoomed zones. Thus, we can only encourage user to be carreful on the positioning of those
#'                           landmarks to avoid to have to modifiy most of the landmarks during the step 2/...
#'
#'          \strong{c. Second refinement: using mesh/plane intersection as a guideline to digitize landmark}
#'
#'          An optional preliminary step allow user to interactively rotate plane(s) intersecting the mesh, and
#'          to keep a record of this intersection during the landmark digitizing. It could be of interest when
#'          some landmarks to digitize are located along a symmetry axis. To process this step, user should set
#'          options of the GrOpt argument (through the function \code{\link{setGraphicOptions}}, see associated help
#'          for details).
#'
#'          Mouse and keyboard interactions for plane rotation: \cr
#'          - left click: mesh & plane(s) rotation as a single block \cr
#'          - scroll wheel: zoom and dezoom \cr
#'          - right click for Windows users, and cmd + right click for Mac users: rotation of the mesh while the
#'            plane(s) stay fixed. \cr
#'          - Escap key press: validate the plane(s) positioning (and then pass to the landmark digitization)
#'
#'          With this option, the function process divides into 3 steps: \cr
#'          0/ Plane(s) rotation \cr
#'          1/ Landmark digitizing \cr
#'          2/ Once all the landmarks are digitized, user can modify (or not) any of them as many times as necessary.
#'
#'          During step 0/, user is free to rotate as many time as wanted the plane(s) relative to the mesh. The
#'          rotation is validated by pressing the Escap key. Then, steps 1/ and 2/ will be processed as described
#'          before. Note that the proposed planes are limited to the planes made by major axes of the mesh, so that
#'          from 1 to 3 planes can be drawn, they are orthogonal each one to the other, and they are centred on the
#'          mesh centroid.
#' @usage
#' \method{digitMesh}{mesh3d}(specFull, specDecim, fixed, idxFixed = 1:fixed, templateCoord = NULL,
#'           idxTemplate = NULL, GrOpt = setGraphicOptions(), verbose = TRUE)
#' @param specFull Full resolution mesh3d object.
#' @param specDecim Decimated resolution mesh3d object, as obtained through \code{\link{decimMesh.mesh3d}} for example.
#'           If missing, the mesh will be decimated to the \code{tarface} target value.
#' @param fixed Number of landmarks to digitize.
#' @param idxFixed Numeric vector with \code{fixed} positive integers specifing the landmark ordering in which the
#'                 landmarks will be digitized. \cr
#'                 Default: \code{1:fixed}, meaning that landmarks will be digitized following the ordering of their
#'                          numbers.
#' @param templateCoord Numeric matrix with three columns (x,y,z) and at least four lines indicating the
#'                      coordinates of at least four 3D points needed to fit the template configuration on the mesh. \cr
#'                      \strong{Warning}: The landmarks in the template configurations should be sorted following their
#'                                        numbers, even if the landmarks used to fit it on the mesh aren't stored in
#'                                        the first lines of this matrix. In such a case, it should be specified with
#'                                        \code{idxTemplate}. \cr
#'                      Default:  \code{NULL} => no template will be used.
#' @param idxTemplate Numeric vector with positive integers indicating the numbers of landmarks of the template used to
#'                    fit it on the mesh, sorted in the order with which they will be digitized on the mesh. For
#'                    example, if the landmarks used are the landmarks numbered 10, 12, 17 and 23, and digitized in the
#'                    following order on the mesh: 12, 10, 23, 17, \code{idxTemplate} should be set to:
#'                    \code{c(12, 10, 23, 17)}. \cr
#'                    Default: \code{NULL} => no template will be used, but corrected to \code{1:4} if
#'                             \code{templateCoord} is provided, but not \code{idxTemplate}.
#' @param GrOpt List defining options for graphic rendering. See \code{\link{setGraphicOptions}} for details.
#' @param verbose Possible settings are: \cr
#'                - a logical value: in this case this value should be recycled in a 2 length vector indicating
#'                  for 2 levels of verbose if comments should be printed or not on screen as the computations are
#'                  processed. The firs level corresponds to comments specific to the functions from the
#'                  \code{digit3DLand} library, and the second one to comments specific to the functions from the
#'                  \code{Rvcg} library. \cr
#'                - a 2-length logical vector standing for the 2 possible levels of verbose.
#' @param tarface Number of target faces to decimated the mesh, used if specDecim is missing.
#' @param spec.name Attribute for the returned A array indicating the specimen name. Possible settings are: \cr
#'                  - NULL (default): array name is set as the mesh3d object name \cr
#'                  - a character value given the array name.
#' @return A numeric matrix with \code{fixed} lines and 3 columns containing the 3D coordinates of the digitized
#'         landmarks.
#' @seealso \code{\link{digitMesh.character}}.
#' @export
#'
digitMesh.mesh3d <- function (specFull, specDecim = NULL, fixed, idxFixed = 1:fixed,
                              templateCoord = NULL, idxTemplate = NULL,
                              GrOpt = setGraphicOptions(), verbose = c(TRUE, TRUE),
                              tarface = 15000, spec.name=NULL) {

    # check OS and R GUI to avoid graphic incompatibilities related to mac os and Rstudio
    tmp <- checkOsGui(GrOpt$winOptions$winNb, GrOpt$winOptions$winSynchro)
    GrOpt$winOptions$winNb <- tmp[[1]]
    GrOpt$winOptions$winSynchro <- tmp[[2]]

    # check mesh
    if (!(any(class(specFull) == "mesh3d"))) {
        stop("specFull must have class \"mesh3d\".")
    }
    if (is.null(spec.name)){
        spec.name <- deparse(substitute(specFull))
    }

    # check verbose
    verbose <- checkLogical(verbose, c(1, 2))

    # Correction if mesh has non-manifold faces (ie faces made of non-manifold edges,
    # ie edges shared by more than 2 faces)
    # Correction needed for the ordering of the intersection points among mesh and planes
    if (verbose[1]){
        cat("\nChecking full & decimated meshes: starts...")
        if (verbose[2]){
            cat("\n")
        }
    }
    specFull <- vcgUpdateNormals(specFull, silent = !verbose[2])
    specFull <- vcgClean(specFull, sel = 2, silent = !verbose[2])

    # checking for the second argument (either specDecim or Nb landmark)
    test1 <- test2 <- FALSE
    if (is.null(specDecim)){
      test1 <- TRUE
    } else if (!(any(class(specDecim) == "mesh3d"))){
      test2 <- TRUE
    }
    if ( test1 | test2 ){
        if (test2){
            # we check later on if fixed is set correctly as a single number
            fixed <- specDecim
        }
        warning(paste("specDecim was missing.
                      Decimates the full mesh to tarface = ", tarface), immediate. = TRUE)
        specDecim <- decimMesh(specFull, tarface = tarface, silent = FALSE)
    }

    specDecim <- vcgUpdateNormals(specDecim, silent = !verbose[2])
    specDecim <- vcgClean(specDecim, sel = 2, silent = !verbose[2])

    if (verbose[1]){
        if (!verbose[2]) cat("\r")
        cat("Checking full & decimated meshes: done!    \n")
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

    if (missing(fixed)) {
        stop("missing number of landmarks to digitalize")
    } else {
        if (!is.numeric(fixed) || length(fixed) > 1)
            stop("fixed must be a single number")
    }

    # Use or not of a template
    if (is.null(templateCoord)){
        idxTemplate <- 1:fixed
    } else {
        if (missing(idxTemplate)) {
            idxTemplate <- 1:4
            warning("idxTemplate was missing.
                    First 4 landmarks will be used to align the template")
        } else {
            if (length(idxTemplate) < 4) {
                stop("idxTemplate must contain at least 4 landmarks")
            }
        }
        p1 <- length(idxTemplate)
        template <- list()
        template$M <- templateCoord
    }

    # Define default values for graphics interactivity
    grDev <- GrOpt
    # check if the mesh is actually colored
    if (!GrOpt$meshOptions$meshVertCol)
        specFull$material$color <- specDecim$material$color <- NULL
    if (grDev$meshOptions$meshVertCol & is.null(specFull$material$color)){
        grDev$meshOptions$meshVertCol <- FALSE
    }
    grDev$vSp <- grDev$vTx <- Sp <- Tx <- rep(NA, fixed)
    grDev$spradius <- GrOpt$spheresOptions$spheresRad
    tmp <- diff(apply(specDecim$vb[1:3, ], 1, range))
    grDev$spradius[, 1] <- GrOpt$spheresOptions$spheresRad[, 1] * mean(tmp)
    grDev$labadj <- GrOpt$labelOptions$labelAdj * mean(tmp)

    # Centering of the meshes on the centroid of the decimated one
    tmp <- scale(t(specDecim$vb[-4, ]), scale = FALSE)
    Trans1 <- attr(tmp, which = "scaled:center")
    specDecim$vb[-4, ] <- t(tmp)
    specFull$vb[-4, ] <- sweep(specFull$vb[-4, ], 1, Trans1)

    if (verbose[1]){
        cat("\rInitializations for digitMesh.mesh3d: done!         \n\n")
        cat("Plotting decimated mesh: in progress...")
    }

    # plot decimated mesh
    d1 <- Clear3d()
    par3d(windowRect = grDev$winOptions$winSize[1, ])
    if (grDev$winOptions$winNb == 1){
        d1 <- currentSubscene3d()
        layout3d(t(c(1, 2)), sharedMouse = grDev$winOptions$winSynchro)
        next3d()
    }
    if (grDev$meshOptions$meshShade[1]){
        shade3d(specDecim, col = grDev$meshOptions$meshColor[1],
                alpha = grDev$meshOptions$meshAlpha[1])
    }
    if (grDev$meshOptions$meshWire[1]){
        wire3d(specDecim, col = grDev$meshOptions$meshColor[1],
               alpha = grDev$meshOptions$meshAlpha[1])
    }
    if (grDev$meshOptions$meshPoints[1]){
        points3d(t(specDecim$vb[1:3,]), col = grDev$meshOptions$meshColor[1],
                 alpha = grDev$meshOptions$meshAlpha[1])
    }
    grDev$dev <- rgl.cur()

    # Rotate the scene
    R <- rotMajorAxes(specDecim$vb[1:3, ])
    par3d(userMatrix = R)
    if (verbose[1]) cat("\rPlotting decimated mesh: done!         \n")

    # plot of orthogonal planes: they are initialized as major axes of the mesh
    orthoplanes <- list(vInter = NULL, vPlanes = NULL)
    if (length(idxPlanes) > 0){
        if (verbose[1]) cat("\nPlotting mesh/plane intersections: starts...\n")

        orthoplanes <- DrawOrthoplanes(mesh = specDecim, idxPlanes = idxPlanes,
                                       grDev = grDev, verbose = verbose)
        if (ncol(specDecim$vb) != ncol(specFull$vb)) {
            # computation of intersections among full mesh and fixed planes
            orthoplanes <- DrawOrthoplanes(mesh = specFull, idxPlanes = idxPlanes,
                                           planes = orthoplanes$vPlanes,
                                           interactive = FALSE, is.plot = FALSE,
                                           grDev = grDev, verbose = verbose)
        }
        if (verbose[1]) cat("\n Plotting mesh/plane intersections: done!\n")
    }

    # Landmark selection - A is the individual configuration matrix [k x 3]
    A <- Adeci <- matrix(NA, fixed, 3, dimnames = list(1:fixed, c("x", "y", "z")))
    attr(A, which = "spec.name") <- spec.name

    if (verbose[1]) {
        cat("\nLoop for landmark digitization: starts...\n")
        cat("Left click to rotate, scroll wheel to zoom, (for mac users: cmd +) right click to position a landmark.\n")
    }
    Idx <- setdiff(1:fixed, idxFixed[idxTemplate])
    for (i in 1:fixed) {
        if (i <= length(idxTemplate)) {
            # Place 1st points require to adjust the template if it exist otherwise take all pts
            idx_pts <- idxFixed[idxTemplate[i]]
            if (verbose[1])
                cat(paste0("\nPlease digitize landmark: ", idx_pts))

            res <- SelectPoints3d(mesh = specDecim, A = A, IdxPts = idx_pts, grDev = grDev, whichMesh = 1)

            Pt <- res$coords
            grDev$vSp[idx_pts] <- res$sp
            grDev$vTx[idx_pts] <- res$tx
        } else {
            # Selection of remaining landmarks (if any)
            idx_pts <- idxFixed[Idx[i - length(idxTemplate)]]
            if (verbose[1])
                cat(paste0("\nPlease digitize landmark: ", idx_pts))

            #Pt <- B[idx_pts,, drop = FALSE]
            Pt <- B[idx_pts, ]
        }
        # zoom on full resolution mesh around the selected landmark
        Pt2 <- c(project(t(Pt), specFull, trans = TRUE)) # added
        res2 <- SetPtZoom(specFull = specFull, Pt = Pt2, IdxPts = idx_pts,
                          orthoplanes = orthoplanes, idxPlanes = idxPlanes, grDev = grDev)
        grDev <- res2$grDev
        # landmark coordinate on the full resolution mesh
        A[idx_pts, ] <- res2$coords
        # Projection of landmarks on decimated mesh for graphics
        # Adeci[idx_pts,] <- project(res2$coords, specDecim, trans = TRUE)
        Adeci[idx_pts, ] <- specDecim$vb[1:3, which.min(colSums(abs(specDecim$vb[1:3, ] - c(res2$coords))))]

        # Graphics
        grDev <- plot.landmark(Adeci[idx_pts, ], d1, idx_pts, grDev, exist = TRUE)

        if (verbose[1])
            cat(paste0("\nLandmark ", idx_pts, " has been digitized."))

        if (!is.null(templateCoord) & i == length(idxTemplate)){

            # all reference points of the template are placed
            # => impute missing landmarks
            B <- imputeCoords(A, template = template$M) #idx = idxTemplate
            ptsB <- project(B, specDecim)
            B <- project(B, specFull, trans = TRUE)

            # plot points/labels of B not placed before
            vv <- idxFixed[Idx]
            if (length(vv) > 0){
                for (ii in 1:length(vv)){
                    grDev <- plot.landmark(t(ptsB$vb[1:3, vv[ii]]), d1, vv[ii],
                                           grDev, exist = FALSE)
                }
                #rgl.viewpoint(userMatrix = R)
            }
        }
    }

    if (verbose[1]){
        cat("\n\nLoop for landmark digitization: ends.\n\n")
        cat("Now, you can:\n")
        cat(" - validate your digitization by closing the graphic device\n")
        cat(" - or modify some landmarks if necessary (just click a on a landmark to modify and redigitize it, \n")
        cat("   once all landmraks are correct, just close the graphic device).\n")
    }

    # Now, wait if the user want changed any landmark. Stop when the graphics is closed
    Stop <- 0
    grDev$dev <- d1
    while (Stop == 0){
        # clicks point on the decimated mesh
        res <- SelectPoints3d(mesh = specDecim, modify = TRUE, A = Adeci, grDev = grDev, whichMesh = 1)
        if (res$isClosed) {
            if (verbose[1]) cat("\n")
            break
        }

        idx_pts <- res$Idx

        if (verbose[1])
            cat(paste0("\nLandmark to modify: ", idx_pts))

        # zoom on full resolution mesh
        Pt2 <- c(project(t(res$coords), specFull, trans = TRUE)) # added
        res2 <- SetPtZoom(specFull = specFull, Pt = Pt2, IdxPts = idx_pts,
                          orthoplanes = orthoplanes, idxPlanes = idxPlanes, grDev =  grDev)
        grDev <- res2$grDev
        # landmark coordinate on the full resolution mesh
        A[idx_pts, ] <- res2$coords
        # Projection of the landmark on the decimated mesh for graphics
        # Adeci[idx_pts, ] <- project(res2$coords, specDecim, trans = TRUE)
        idx <- which.min(colSums(abs(specDecim$vb[1:3, ] - c(res2$coords))))
        Adeci[idx_pts, ] <- specDecim$vb[1:3, idx]

        # Graphics
        grDev$vSp[idx_pts] <- res$sp
        grDev$vTx[idx_pts] <- res$tx
        grDev <- plot.landmark(Adeci[idx_pts, ], d1, idx_pts, grDev, exist = TRUE)

        if (verbose[1])
            cat(paste0("\rLandmark ", idx_pts, " has been modified."))

    }
    # restore position
    A <- A + matrix(Trans1, fixed, 3, byrow = TRUE)

    if (verbose[1])
        cat("\nMesh digitization is ended!\n")

    return(A)
}

#' @title Digitizes several Meshes
#' @description Interactive digitization of a mesh or a set of meshes from mesh files (either ply or stl).
#' @details This function is a wrapper for \code{\link{digitMesh.mesh3d}}, calling it to treat a list of mesh files,
#'          which will be digitized one after the other. Options for mesh file opening and landmark coordinate saving
#'          are settable with the \code{FiOpt} argument. A preliminary decimation step can be processed to create the
#'          decimated version of the meshes if needed. The settings for this decimation step are set through the
#'          \code{DeOpt} argument. Moreover, to fasten the digitizing process, the use of a template configuration
#'          is possible whose options can be set through the \code{TeOpt} argument. Basically, with default arguments,
#'          the function will sequentially open the mesh files in a given directory, decimate them, allow user to
#'          digitize landmark (see \code{\link{digitMesh.mesh3d}} for details), store landmark coordinates in an
#'          array, and export them into a tps file. The first digitized mesh will serve as template for the
#'          digitization of the following meshes. After each mesh digitization a message displayed in the console asks
#'          the user to digitize or not the next mesh file.
#' @usage
#' \method{digitMesh}{character}(sdir, fixed, idxFixed = 1:fixed, GrOpt=setGraphicOptions(),
#'           FiOpt=setFileOptions(sdir), DeOpt=setDecimOptions(), TeOpt=setTemplOptions(fixed),
#'           verbose = TRUE, \dots)
#' @param sdir A character value indicating either a mesh filename to decimate stored within the working directory,
#'             or a directory name within the working directory containing the subdirectory\code{M} with the mesh
#'             files to decimate.
#' @param fixed Number of landmarks to digitize.
#' @param idxFixed Numeric vector with \code{fixed} integers specifing the landmark labelling in which the landmarks
#'                 will be digitized.
#' @param GrOpt List defining options for graphic rendering. See \code{\link{setGraphicOptions}} for details.
#' @param FiOpt List defining options for file opening and saving. See \code{\link{setFileOptions}} for details.
#' @param DeOpt List defining options for mesh decimation. See \code{\link{setDecimOptions}} for details.
#' @param TeOpt List defining options for template definition and use. See \code{\link{setTemplOptions}} for details.
#' @param verbose Possible settings are: \cr
#'                - a logical value: in this case this value should be recycled in a 2 length vector indicating
#'                  for 2 levels of verbose if comments should be printed or not on screen as the computations are
#'                  processed. The firs level corresponds to comments specific to the functions from the
#'                  \code{digit3DLand} library, and the second one to comments specific to the functions from the
#'                  \code{Rvcg} library. \cr
#'                - a 2-length logical vector standing for the 2 possible levels of verbose.
#' @param ... Optional arguments used for mesh decimation. See \code{\link[Rvcg]{vcgQEdecim}} for details.
#' @return An array with \code{fixed} lines, 3 columns and \emph{n} slices (one for each treated mesh) containing the
#'         3D coordinates of the digitized landmarks.
#' @seealso \code{\link{digitMesh.mesh3d}}.
#' @export
#'
#' @examples
#'
#' ## Not run:
#' # Below some possible uses of the digitMesh.character() function (examples are not exhausive, see in
#' # particular the helps of setGraphicOptions(), setFileOptions(), setDecimOptions() and setFileOptions() for details).
#'
#' # 1st example: digitizing a mesh file
#' # A basic call consists in giving the filename of the mesh to digitize (contained in the working directory):
#' A <- digitMesh("mesh2digitize.ply", 10)
#'
#' # 2nd example: digitizing mesh files in a given directory
#' # If the working directory contains a subdirectory named "fold" with the full resolution mesh files to digitize, a basic call of the function is:
#' A <- digitMesh("fold", 10)
#'
#' # 3rd example: pursuing a previous digitiztion session
#' # If during a first digitization session, not all the mesh files were processed, and that landmark
#' # coordinates were saved into a tps file, it is possible to follow the mesh digitization with untreated
#' # mesh by specifying the already treated mesh files with the tps file (named "TPS_FileName.tps"):
#' sdir <- "fold"
#' FiOpt <- setFileOptions(sdir, saveTPS="TPS_FileName.tps", append=TRUE)
#' A <- digitMesh(sdir, 10, FiOpt=FiOpt)
#'
#' # 4th example: digitizing mesh files already decimated, contained in the subfolder "DM" itself within the
#' # folder full resolution meshes (decimated mesh filemanes distinguish also from full mesh filenames by the suffix "_D"):
#' sdir <- "fold"
#' deci.dir<- "DM"
#' FiOpt <- setFileOptions(sdir, deci.suffix="_D", deci.dir="DM")
#' DeOpt <- setDecimOptions(makeDecimation=FALSE)
#' A <- digitMesh(sdir, 10, FiOpt=FiOpt, DeOpt=DeOpt)
#'
#' # 5th example: using "percent" (percentage of face number reduction, see Rvg:::vcgQEdeci) in place of
#' # "tarface" (targetted number of faces) for decimation:
#' A <- digitMesh("fold", 10, percent=0.2)
#'
#' # 6th example: processing all mesh decimation before mesh digitization
#' DeOpt <- setDecimOptions(sequential=FALSE)
#' A <- digitMesh("fold", 10, DeOpt=DeOpt)
#'
#' # 7th example: defining a given individual for the template configuration
#' TeOpt<-setTemplOptions(10, template="mesh4tempalte.ply")
#' A <- digitMesh("fold", 10, TeOpt=TeOpt)
#'
#' # 8th example: using a coordinate matrix M (from an already digitized mesh) as template
#' TeOpt<-setTemplOptions(10, template=M)
#' A <- digitMesh("fold", 10, TeOpt=TeOpt)
#'
#' # 9th example: processing stl files (contained in the sudirectory "fold4stl")
#' FiOpt <- setFileOptions("stl", patt=".stl")
#' A <- digitMesh("fold4stl", 10, FiOpt=FiOpt)
#'
#' # 10th example: adjusting and drawing 1st major plane before mesh digitization
#' GrOpt <- setGraphicOptions(PCplanesDraw="pc1-pc2")
#' A <- digitMesh("fold4stl", 10, GrOpt=GrOpt)
#' Numerous graphical options are settable: see the help of setGraphicOptions() for details.
#'
#' ## End(Not run)
digitMesh.character <- function(sdir, fixed, idxFixed = 1:fixed,
                                GrOpt = setGraphicOptions(),
                                FiOpt = setFileOptions(sdir),
                                DeOpt = setDecimOptions(),
                                TeOpt = setTemplOptions(fixed),
                                verbose = c(TRUE, TRUE), ...) {

    curdir <- getwd()

    # check verbose
    verbose <- checkLogical(verbose, c(1, 2))

    if (verbose[1])
        cat("\nChecking arguments for digitMesh: in progress...")

    # extract decimation options
    makeDecimation <- DeOpt$makeDecimation
    sequential <- DeOpt$sequential
    tarface <- DeOpt$tarface

    # extract file options
    sdir <- FiOpt$sdir
    patt <- FiOpt$patt
    deci.suffix <- FiOpt$deci.suffix
    deci.dir <- FiOpt$deci.dir
    full.dir <- FiOpt$full.dir
    full.files <- FiOpt$full.files
    saveTPS <- FiOpt$saveTPS
    append <- FiOpt$append
    overwrite <- FiOpt$overwrite

    if (verbose[1])
        cat("\rChecking arguments for digitMesh: done!         \n")

    # get back template coordinates (if any)
    # if user go back to a partially processed directory
    ongoing <- FALSE
    if (is.character(saveTPS)) {
        if (append & !overwrite) {
            if (verbose[1])
                cat("\nRecovering previous digitized data: in progress...")

            ongoing <- TRUE
            coord <- read.tps(saveTPS, sdir = FiOpt$full.dir, quiet = !verbose[2])
            tmp <- dimnames(coord)[[3]]
            done.files.full <- paste0(tmp, patt)

            if (is.character(TeOpt$template)) {
                # template is a particular mesh already digitized during the previous session
                TeOpt$template <- coord[,, which(done.files.full == TeOpt$template)]
            }
            if (is.logical(TeOpt$template)) {
                if (TeOpt$template) {
                    # template is the 1st mesh already digitized during the previous session
                    TeOpt$template <- coord[,, 1]
                }
            }
            # full.files: remaining mesh to digitize
            full.files <- setdiff(full.files, done.files.full)
            if (!makeDecimation | (makeDecimation & !sequential)){
                deci.files <- paste0(gsub(patt, "", full.files), deci.suffix, patt)
            }

            if (verbose[1]){
                if (verbose[2]){
                    cat("Recovering previous digitized data: done!\n")
                } else {
                    cat("\rRecovering previous digitized data: done!         \n")
                }
            }
        }
    }

    if (verbose[1])
        cat("\nChecking template options: in progress...")

    # Now that full.files is defined,
    # checks if the supplied template filenemame (if any) is within those files
    TeOpt <- setTemplOptions(fixed, template = TeOpt$template,
                             idxTemplate = TeOpt$idxTemplate, full.files = full.files)
    # extract template options
    template <- TeOpt$template
    idxTemplate <- TeOpt$idxTemplate
    makeTempl <- TeOpt$makeTempl

    if (verbose[1])
        cat("\rChecking template options: done!         \n")

    if (makeDecimation & !sequential){
        # all meshes are decimated before to be digitized...
        if (!ongoing) {
            if (verbose[1])
                cat("\nDecimation of all meshes: starts...")

            #... but only if the folder is browsed for the first time
            if (identical(full.dir, sdir)){
                ar1 <- full.files
            } else {
                ar1 <- strsplit(full.dir, paste0(sdir, "/"))[[1]][2]
            }
            Ldeci <- decimMesh(ar1, tarface = tarface, sdir = sdir,
                             patt = patt, deci.suffix = deci.suffix,
                             deci.dir = strsplit(deci.dir, paste0(full.dir, "/"))[[1]][2],
                             verbose = verbose, ...)

            if (verbose[1])
                cat("\nDecimation of all meshes: ends!")
        }
    }

    if (!makeDecimation | (makeDecimation & !sequential)){
        # either decimation is not needed, or it is processed in a single pass
        # before digitization
        if (verbose[1])
            cat("\nChecking for filename correspondence among full and decimated meshes: in progress...")

        if (deci.dir == full.dir){
            # decimated and full meshes are stored in the same folder,
            # differenciating by a suffix in the filenames for decimated meshes
            setwd(deci.dir)
            ply.files <- list.files(pattern = patt, ignore.case = TRUE)
            deci.files <- list.files(pattern = paste0(deci.suffix, patt), ignore.case = TRUE)
            if (!ongoing){
                # folder browsed for the first time
                full.files <- setdiff(ply.files, deci.files)
            }
        } else {
            # decimated and full mesh files are in 2 separate subfolders within sdir
            setwd(full.dir)

            # extract identifiers for full meshes
            ID1 <- unlist(strsplit(full.files, patt))

            setwd(deci.dir)
            if (!ongoing) {
                # folder browsed for the first time
                deci.files <- list.files(pattern = patt, ignore.case = TRUE)
            }
            # extract identifiers for decimated meshes
            ID2 <- unlist(strsplit(deci.files, paste0(deci.suffix, patt)))
        }

        # check the correspondence among identifiers for full and decimated meshes
        ID1 <- unlist(strsplit(full.files, patt))
        ID2 <- unlist(strsplit(deci.files, paste0(deci.suffix, patt)))
        ID <- checkID(ID1, ID2)

        if (verbose[1])
            cat("\rChecking for filename correspondence among full and decimated meshes: done!         \n")

    } else {
        # sequential decimation: will be processed separately for each mesh to digitize
        setwd(full.dir)
        if (verbose[1])
            cat("\nExtracting mesh ID: in progress...")

        # extract identifiers for meshes
        ID <- unlist(strsplit(full.files, patt))
        if (verbose[1])
            cat("\rExtracting mesh ID: done!         \n")
    }

    # check if at least one mesh file was found
    n <- length(ID)
    if (n < 1){
        stop("No mesh files found...")
    }

    # if a template is used, and not being the first individual, find which file will be used, and define the order
    # for mesh digitizing (begining)
    idxMesh <- 1:n
    if (makeTempl & is.character(template)) {
        if (verbose[1])
            cat("\nExtracting individual for template: in progress...")

        idx_tpl <- which(full.files == template)
        idxMesh <- c(idx_tpl, sort(setdiff(1:n, idx_tpl)))
        if (verbose[1])
            cat("\rExtracting individual for template: done!         \n")
    }

    # loop for mesh digitizing (and possibly mesh decimation if sequential=TRUE)
    if (verbose[1])
        cat("\nLoop to digitize all meshes: starts...\n")

    A <- array(NA, c(fixed, 3, n))
    cpt <- 0
    Vspec.name <- rep(0, n)
    interrupt <- FALSE
    for (i in idxMesh){
        cpt <- cpt + 1

        # import full mesh
        setwd(full.dir)
        ff <- full.files[i]
        if (verbose[1]){
            header <- paste0("********** Mesh to digitize: ", ff, " **********")
            cat("\n", rep("*", nchar(header)), sep="")
            cat("\n", header)
            cat("\n", rep("*", nchar(header)), sep="")
            cat("\n\n")
            cat("Full resolution mesh opening: starts...")
            if (verbose[2])
                cat("\n\n")
        }
        full <- vcgImport(ff, updateNormals = TRUE, readcolor = TRUE,
                          clean = TRUE, silent = !verbose[2])
        if (verbose[1]){
            if (verbose[2]) {
                cat("\nFull resolution mesh opening: done!\n")
            } else {
                cat("\rFull resolution mesh opening: done!    \n")
            }
        }

        # create or import decimated mesh
        if (makeDecimation & sequential){
            # sequential decimation: decimate full mesh
            if (verbose[1])
                cat("\nFull resolution mesh decimation: starts...\n")

            deci <- decimMesh(full.files[i], tarface = tarface, sdir = full.dir,
                              patt = patt, deci.suffix = deci.suffix,
                              deci.dir = strsplit(deci.dir, paste0(full.dir, "/"))[[1]][2],
                              verbose = verbose, ...)
            deci <- deci[[1]]
            if (verbose[1])
                cat("\nFull resolution mesh decimation: done!\n")

        } else {
            # import decimated mesh
            if (verbose[1])
                cat("\nImporting decimated mesh: in progress...\n")

            setwd(deci.dir)
            deci <- vcgImport(deci.files[i], updateNormals = TRUE, readcolor = TRUE,
                              clean = TRUE, silent = !verbose[2])
            if (verbose[1])
                cat("\nImporting decimated mesh: done!         \n")
        }
        setwd(curdir)

        # use or not template
        if (makeTempl & is.character(template)) {
            # use template (template being a filename)
            if (cpt == 1){
                # the first indivual is the template
                tmpA <- digitMesh(full, deci, fixed = fixed, idxFixed = idxFixed,
                                  GrOpt = GrOpt, verbose = verbose,
                                  spec.name = gsub(".stl", "", gsub(".ply", "", ff)))
                # we store its coordinates for use with the next meshes to digitize
                tpl <- tmpA
            } else {
                # the other ones use this template
                tmpA <- digitMesh(full, deci, fixed = fixed, idxFixed = idxFixed,
                                templateCoord = tpl, idxTemplate = idxTemplate,
                                GrOpt = GrOpt, verbose = verbose,
                                spec.name = gsub(".stl", "", gsub(".ply", "", ff)))
            }
        } else {
            if (is.matrix(template)) {
                # use template (template being a matrix)
                tmpA <- digitMesh(full, deci, fixed = fixed, idxFixed = idxFixed,
                                  templateCoord = template, idxTemplate = idxTemplate,
                                  GrOpt = GrOpt, verbose = verbose,
                                  spec.name = gsub(".stl", "", gsub(".ply", "", ff)))
            } else {
                # don't use template
                tmpA <- digitMesh(full, deci, fixed = fixed, idxFixed = idxFixed,
                                  GrOpt = GrOpt, verbose = verbose,
                                  spec.name = gsub(".stl", "", gsub(".ply", "", ff)))
            }
        }

        # get spec.name
        A[,, cpt] <- tmpA
        Vspec.name[cpt] <- attr(tmpA, "spec.name")

        # saving coordinates in a TPS file
        if (is.character(saveTPS)){
            save.tps(A[,, i], ID[i], saveTPS, sdir = full.dir,
                     app = ifelse(cpt == 1, append, TRUE),
                     over.write = ifelse(cpt == 1, overwrite, FALSE))
        }

        # ask user if next mesh (if any) should be digitize
        cat("\n")
        if (cpt < length(idxMesh)) {
            ans <- readline(prompt = "Digitize next mesh ? Type y (for yes) or n (for no): ")
            if (ans == "n") {
                interrupt <- TRUE
                A <- A[,, 1:cpt, drop = FALSE]
                Vspec.name <- Vspec.name[1:cpt]
                dimnames(A) <- list(NULL, NULL, Vspec.name)
                if (verbose[1]) {
                    cat("\nLoop to digitize all meshes: stops before the end ...\n")
                    cat("\n=> The following files will remain to be digitized:")
                    rem.files <- full.files[idxMesh[(cpt + 1):n]]
                    cat(paste0("\n", " - ", rem.files, "\n\n"))
                }
                break
            }
        } else {
            dimnames(A) <- list(NULL, NULL, Vspec.name)
            if (verbose[1])
                cat("Last file reached: digitization loop is ending...\n\n")
        }
    }

    setwd(curdir)

    return(A)
}
