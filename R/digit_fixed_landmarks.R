###########################
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
#'          rotation is validated by pressing the Escap key. Then , steps 1/ and 2/ will be processed as described
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
digitMesh.mesh3d <- function (specFull, specDecim, fixed, idxFixed = 1:fixed, templateCoord = NULL, idxTemplate = NULL,
                              GrOpt=setGraphicOptions(), verbose=c(TRUE,TRUE), tarface = 15000, spec.name=NULL) {


    # check OS and R GUI
    os <- Sys.info()[1]
    gui <- .Platform$GUI
    if (os == "Darwin"){
        if (gui == "RStudio"){
            # not supported
            stop("The function is not supported with the RStudio interface in Mac OS
                 Please, use the R gui instead.")
        } else{
            # in Mac OS only 2 separate windows are supported
            # => forcing graphic options...
            if (GrOpt$winOptions$winNb == 1){
                warning('with mac OS, multiple interactive subscenes are not supported.
                        winWb option was set to 2')
                GrOpt <- setGraphicOptions(winNb = 2, winSynchro = FALSE)
            }
        }
    }

    # check mesh
    if (!(any(class(specFull) == "mesh3d")))
        stop("specFull must have class \"mesh3d\".")
    if (is.null(spec.name)){
        spec.name <- deparse(substitute(specFull))
    }

    # check verbose
    verbose<-checkLogical(verbose,c(1,2))

    # Correction if mesh has non-manifold faces (ie faces made of non-manifold edges, ie edges shared by more than 2 faces)
    # Correction needed for the ordering of the intersection points among mesh and planes
    if (verbose[1]){
        cat("\n")
        cat("Checking full & decimated meshes: starts...")
        if (verbose[2]){
            cat("\n")
        }
    }
    specFull <- vcgUpdateNormals(specFull, silent = !verbose[2])
    specFull <- vcgClean(specFull, sel=2, silent=!verbose[2])

    if (missing(specDecim)) {
        warning('specDecim was missing with no default. Run decimMesh')
        specDecim <- decimMesh(M, tarface = tarface, silent = FALSE)
    }

    specDecim <- vcgUpdateNormals(specDecim, silent = !verbose[2])
    specDecim <- vcgClean(specDecim, sel=2, silent=!verbose[2])

    if (verbose[1]){
        if (!verbose[2]){
            cat("\r")
        }
        cat("Checking full & decimated meshes: done!    ")
        cat("\n")
        cat("\n")
        cat("Initializations for digitMesh.mesh3d: in progress...")
    }


    # check which setting of GrOpt$PCplanesDraw is called, and set idxPlanes consequently
    if(is.logical(GrOpt$PCplanesOptions$PCplanesDraw)){
        if (GrOpt$PCplanesOptions$PCplanesDraw){
            idxPlanes<-1:3
        }else{
            idxPlanes<-NULL
        }
    }else{
        V<-c("pc2-pc3","pc1-pc3","pc1-pc2")
        idxPlanes<-which(is.element(V,tolower(GrOpt$PCplanesOptions$PCplanesDraw)))
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
            if (length(idxTemplate) < 4)
                stop("idxTemplate must contain at least 4 landmarks")
        }
        p1 <- length(idxTemplate)
        template <- list()
        template$M <- templateCoord
    }

    # Define default values for graphics interactivity
    grDev<-GrOpt
    grDev$vSp <- grDev$vTx <- Sp <- Tx <- rep(NA, fixed)
    grDev$spradius <- GrOpt$spheresOptions$spheresRad
    tmp <- diff(apply(specDecim$vb[1:3,], 1, range))
    grDev$spradius[,1] <- GrOpt$spheresOptions$spheresRad[,1] * mean(tmp)
    grDev$labadj <- GrOpt$labelOptions$labelAdj * mean(tmp)

    # Centering of the meshes on the centroid of the decimated one
    tmp <- scale(t(specDecim$vb[-4, ]), scale = FALSE)
    Trans1 <- attr(tmp, which="scaled:center")
    specDecim$vb[-4, ] <- t(tmp)
    specFull$vb[-4, ] <- sweep(specFull$vb[-4, ], 1, Trans1)

    if (verbose[1]){
        cat("\r")
        cat("Initializations for digitMesh.mesh3d: done!         ")
        cat("\n")
        cat("\n")
        cat("Plotting decimated mesh: in progress...")
    }

    # plot decimated mesh
    d1 <- Clear3d()
    par3d(windowRect = grDev$winOptions$winSize[1, ])
    if (grDev$winOptions$winNb==1){
        d1 <- currentSubscene3d()
        layout3d(t(c(1,2)), sharedMouse = grDev$winOptions$winSynchro)
        next3d()
    }
    if (grDev$meshOptions$meshShade[1]){
        shade3d(specDecim, col=grDev$meshOptions$meshColor[1], alpha=grDev$meshOptions$meshAlpha[1])
    }
    if (grDev$meshOptions$meshWire[1]){
        wire3d(specDecim, col=grDev$meshOptions$meshColor[1], alpha=grDev$meshOptions$meshAlpha[1])
    }
    if (grDev$meshOptions$meshPoints[1]){
        points3d(t(specDecim$vb[1:3,]), col=grDev$meshOptions$meshColor[1], alpha=grDev$meshOptions$meshAlpha[1])
    }
    grDev$dev <- rgl.cur()

    # Rotate the scene
    R <- rotMajorAxes(specDecim$vb[1:3, ])
    par3d(userMatrix = R)

    if (verbose[1]){
        cat("\r")
        cat("Plotting decimated mesh: done!         ")
        cat("\n")
    }

    # plot of orthogonal planes: they are initialized as major axes of the mesh
    orthoplanes <- list(vInter=NULL, vPlanes = NULL)
    if (length(idxPlanes)>0){
        if (verbose[1]){
            cat("\n")
            cat("Plotting mesh/plane intersections: starts...")
            cat("\n")
        }
        orthoplanes <- DrawOrthoplanes(mesh=specDecim, idxPlanes=idxPlanes, grDev=grDev, verbose=verbose)
        if (ncol(specDecim$vb)!=ncol(specFull$vb)){
            # computation of intersections among full mesh and fixed planes
            orthoplanes <- DrawOrthoplanes(mesh=specFull,idxPlanes=idxPlanes,planes=orthoplanes$vPlanes,
                                           interactive=FALSE,is.plot=FALSE, grDev=grDev, verbose=verbose)
        }
        if (verbose[1]){
            cat("\n")
            cat("Plotting mesh/plane intersections: done!")
            cat("\n")
        }
    }

    # Landmark selection - A is the individual configuration matrix [k x 3]
    A <- Adeci <- matrix(NA, fixed, 3, dimnames = list(1:fixed, c("x","y","z")))
    attr(A, which = "spec.name") <- spec.name

    if (verbose[1]){
        cat("\n")
        cat("Loop for landmark digitization: starts...")
        cat("\n")
        cat("Left click to rotate, scroll wheel to zoom, (for mac users: cmd +) right click to position a landmark.")
        cat("\n")
    }
    Idx <- setdiff(1:fixed, idxFixed[idxTemplate])
    for (i in 1:fixed){
        if (i <= length(idxTemplate)){
            # Place 1st points require to adjust the template if it exist otherwise take all pts
            idx_pts <- idxFixed[idxTemplate[i]]
            if (verbose[1]){
                cat("\n")
                txt<-paste0("Please digitize landmark: ",idx_pts)
                cat(txt)
            }

            res <- SelectPoints3d(mesh = specDecim, A = A, IdxPts = idx_pts, grDev = grDev, whichMesh = 1)

            Pt <- res$coords
            grDev$vSp[idx_pts] <- res$sp
            grDev$vTx[idx_pts] <- res$tx
        } else {
            # Selection of remaining landmarks (if any)
            idx_pts <- idxFixed[Idx[i-length(idxTemplate)]]
            if (verbose[1]){
                cat("\n")
                txt<-paste0("Please digitize landmark: ",idx_pts)
                cat(txt)
            }
            #Pt <- B[idx_pts, , drop = FALSE]
            Pt <- B[idx_pts, ]
        }
        # zoom on full resolution mesh around the selected landmark
        Pt2 <- c(project(t(Pt), specFull, trans = TRUE)) # added
        res2 <- SetPtZoom(specFull=specFull, Pt = Pt2, IdxPts = idx_pts,
                          orthoplanes = orthoplanes, idxPlanes=idxPlanes, grDev=grDev)
        grDev<-res2$grDev
        # landmark coordinate on the full resolution mesh
        A[idx_pts, ] <- res2$coords
        # Projection of landmarks on decimated mesh for graphics
        # Adeci[idx_pts,] <- project(res2$coords, specDecim, trans = TRUE)
        Adeci[idx_pts,] <- specDecim$vb[1:3,which.min(colSums(abs(specDecim$vb[1:3,]-c(res2$coords))))]

        # Graphics
        grDev <- plot.landmark(Adeci[idx_pts, ], d1, idx_pts, grDev, exist = TRUE)

        if (verbose[1]){
            cat("\r")
            txt<-paste0("Please digitize landmark: ",idx_pts)
            txt<-paste0("Landmark ",idx_pts," has been digitized.")
            cat(txt)
        }

        if(!is.null(templateCoord) & i==length(idxTemplate)){

            # all reference points of the template are placed
            # => impute missing landmarks
            B <- imputeCoords(A, template = template$M) #idx = idxTemplate
            ptsB <- project(B, specDecim)
            B <- project(B, specFull, trans = TRUE)

            # plot points/labels of B not placed before
            vv <- idxFixed[Idx]
            if (length(vv)>0){
                for (ii in 1:length(vv)){
                    grDev <- plot.landmark(t(ptsB$vb[1:3, vv[ii]]), d1, vv[ii], grDev, exist = FALSE)
                }
                #rgl.viewpoint(userMatrix = R)
            }
        }
    }

    if (verbose[1]){
        cat("\n")
        cat("\n")
        cat("Loop for landmark digitization: ends.")
        cat("\n")
        cat("\n")
        cat("Now, you can:")
        cat("\n")
        cat(" - validate your digitization by closing the graphic device")
        cat("\n")
        cat(" - or modify some landmarks if necessary (just click a on a landmark to modify and redigitize it,")
        cat("\n")
        cat("   once all landmraks are correct, just close the graphic device).")
        cat("\n")
    }

    # Now, wait if the user want changed any landmark. Stop when the graphics is closed
    Stop <- 0
    grDev$dev <- d1
    while((Stop==0)){
        # clicks point on the decimated mesh
        res <- SelectPoints3d(mesh = specDecim, modify = TRUE, A = Adeci, grDev = grDev, whichMesh = 1)
        if (res$isClosed){
            if (verbose[1]){
                cat("\n")
            }
            break
        }

        idx_pts <- res$Idx

        if (verbose[1]){
            cat("\n")
            txt<-paste0("Landmark to modify: ",idx_pts)
            cat(txt)
        }

        # zoom on full resolution mesh
        Pt2 <- c(project(t(res$coords), specFull, trans = TRUE)) # added
        res2 <- SetPtZoom(specFull=specFull, Pt = Pt2, IdxPts = idx_pts,
                          orthoplanes = orthoplanes, idxPlanes=idxPlanes, grDev =  grDev)
        grDev<-res2$grDev
        # landmark coordinate on the full resolution mesh
        A[idx_pts, ] <- res2$coords
        # Projection of the landmark on the decimated mesh for graphics
        # Adeci[idx_pts, ] <- project(res2$coords, specDecim, trans = TRUE)
        Adeci[idx_pts,] <- specDecim$vb[1:3,which.min(colSums(abs(specDecim$vb[1:3,]-c(res2$coords))))]

        # Graphics
        grDev$vSp[idx_pts] <- res$sp
        grDev$vTx[idx_pts] <- res$tx
        grDev <- plot.landmark(Adeci[idx_pts, ], d1, idx_pts, grDev, exist = TRUE)

        if (verbose[1]){
            cat("\r")
            txt<-paste0("Landmark ",idx_pts," has been modified.")
            cat(txt)
        }

    }

    # restore position
    A <- A + matrix(Trans1, fixed, 3, byrow=TRUE)

    if (verbose[1]){
        cat("\n")
        cat("Mesh digitization is ended!")
        cat("\n")
    }

    return(A)
}
