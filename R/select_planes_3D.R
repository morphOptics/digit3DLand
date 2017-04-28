DrawOrthoplanes <- function(mesh,idxPlanes=1:3, planes=NULL, interactive=TRUE, is.plot=TRUE) {

    # Two posssible usages:
    # 1/ DrawOrthoplanes(mesh)
    # Interactive plot and rotation of intersection among mesh and planes
    #
    # 2/ DrawOrthoplanes(mesh, planes, interactive=FALSE, is.plot=FALSE)
    # In case of already computed intersections among a decimated mesh and planes, and to compute
    # intersections among full mesh and planes. Needs infos on previously determined planes, and no
    # interactive plotting is allowed.
    # NOTE: in this usage, intersections among planes and full mesh involve planes determined on decimated mesh,
    # ie their centers are the centroid of the decimated mesh, and they correspond to major planes computed on
    # the decimated mesh.

    if (class(mesh) != "mesh3d") stop("mesh should be an object \"mesh3d\".")

    if (is.null(planes)){
        # centering of vertex coordinates
        pts <- t(mesh$vb[1:3, ])
        A <- apply(pts, 2, mean)
        pts <- sweep(pts, 2, A)

        # Planes
        sv <- svd(crossprod(pts))
    }

    # Compute and draw the intersection mesh/planes
    ptsPlanes <- array(NA, c(3,3,3))
    vInter<- list()
    for (i in idxPlanes){

        # plane PC2-PC3 when i=1, PC1-PC3 when i=2, PC1-PC2 when i=3
        ci <- setdiff(1:3, i)

        # planes/mesh intersection
        if (is.null(planes)){
            tmp <- meshPlaneIntersect2(mesh, A, A+sv$u[,ci[1]], A+sv$u[,ci[2]])
        }else{
            tmp <- meshPlaneIntersect2(mesh, planes[1,,i], planes[2,,i], planes[3,,i])
        }

        inter<- tmp[[1]]
        edgesTot<-tmp[[2]]
        edges_in_cplx<-edgesTot[,1]+1i*edgesTot[,2] # expression of edges of intersected faces in a complex form for facility
        faces_in<-edgesTot[,3] # intersected faces
        is_bd<-as.logical(edgesTot[,4]) # indicates if edges of the mesh of intersection are or not border edges

        # Sort the intersection points.
        # As returned by meshPlaneIntersect2, the intersection points aren't sorted in a way that a use of lines3d
        # will give proper links among those points.
        # The Morpho:::sortCurve function could be used in this aim, but it stays an approximation in our case on how points
        # must be linked in that sense that the sorting is based only on a proximity among points computed from their coordinates.
        # The connectivity of the mesh faces on which intersection points lay on are not taking into account in this sorting.
        # The sortCurveMesh.cpp function uses this info of face connectivity to determine the point sorting.
        out <- .Call("sortCurveMesh", edges_in_cplx,faces_in,is_bd)

        # Convert "out" into a list where each element correpond to a submesh with at each time a vector indicating the order
        # of points. A NA value will be append at the end of each vector to allow the plot of lines 3d in a single pass (NA
        # value not being linked by lines3d)
        Lve<-list()
        for (j in 1:max(out[[2]])){
            Lve[[j]]<-c(which(out[[2]]==j)[order(out[[1]][out[[2]]==j])],NA)
        }
        inter<-inter[unlist(Lve),]

        if (is.plot){
            # plot of the intersection
            lines3d(inter, col = "red", lwd=2)
        }

        # stockage
        vInter[[i]]<-inter
        if (is.null(planes)){
            ptsPlanes[, , i] <- rbind(A, A+sv$u[,ci[1]], A+sv$u[,ci[2]])
        }else{
            ptsPlanes[, , i] <- planes[,,i]
        }

        # for each plane, contains 3 points in the plane (needed just after)
    }

    if (is.plot){
        # Draw orthogonal planes via planes3d()
        # coefficients a,b,c,d of the plane equation can be determined by a point of the plane
        # and a normal vector to this plane
        # https://fr.wikipedia.org/wiki/Plan_%28math%C3%A9matiques%29#D.C3.A9finition_par_un_vecteur_normal_et_un_point

        for (i in idxPlanes) {
            # a normal vector to a plane can be computed through the vector cross product of 1 orthogonal vectors
            # contained in this plane
            A <- ptsPlanes[1, , i]
            v1<-ptsPlanes[2,,i]-A
            v2<-ptsPlanes[3,,i]-A
            n<-xprod(v1,v2)
            d <- -t(n) %*% A
            planes3d(n[1], n[2], n[3], d, alpha = 0.7, col="cyan")
        }
    }

    if (interactive){
        # User interaction: manual rotation of the mesh (the orthogonal planes being fixed)
        # until the wanted orientation
        return(RotateMeshPlane3d(mesh, planes = ptsPlanes, vInter, idxPlanes))
    }else{
        return(list(vInter=vInter))
    }

}

#############################################################
RotateMeshPlane3d <- function(mesh, planes, vInter, idxPlanes) {
    Stop <- 0
    # stop when the rotation is validated by ESC
    while (Stop==0) {
        temp <- selectPlanes(button="right", mesh=mesh, planes=planes, vInter=vInter, idxPlanes=idxPlanes)
        if (temp$isDone){
            # because of rgl.... need to close twice the window
            if (temp$isClosed){
                rgl.close()
            }
            break
        }
    }
    return(temp)
}

#############################################################
selectPlanes<-function (button = c("left", "middle", "right"), mesh, planes, vInter, idxPlanes) {

    # Re-use of the rgl:::mouseTrackball (from binary) normally used to manually rotate the mesh, and disorted here
    # so that depending on this manual rotation (through the mouse right button), plane(s) seem to stay fix relative
    # to this rotation (which will affect the mesh). In true, and even for the initial rgl:::mouseTrackball function,
    # the mesh stays fix in that its absolute coordinates don't change (axes rotate with the mesh), whereas the new
    # plane coordinates have to be updated at each user interaction.

    # initializations: re-used from rgl:::mouseTrackball
    width <- height <- rotBase <- NULL
    userMatrix <- list()
    cur <- rgl.cur()

    # our initializations
    is.complete<-TRUE
    uplanes<-NULL
    vPlanes<-planes

    # useful function from rgl:::mouseTrackball
    screenToVector <- function(x, y) {
        # Seems to compute a vector between 2 successive mouse positions
        # Called twice in the process, and the 2 vectors allow then to compute an angle and a 3D rotation axis.
        radius <- max(width, height)/2
        centre <- c(width, height)/2
        pt <- (c(x, y) - centre)/radius
        len <- vlen(pt)

        if (len > 1.e-6) pt <- pt/len

        maxlen <- sqrt(2)
        angle <- (maxlen - len)/maxlen*pi/2
        z <- sin(angle)
        len <- sqrt(1 - z^2)
        pt <- pt * len
        return (c(pt, z))
    }

    # function called when the right mouse button starts to be pressed
    trackballBegin <- function(x, y) {

        if(is.complete){
            # test needed when computations for plane intersections are long: in this case, a first held right click can start
            # the intersection computations in trackballUpdate(), but before the end of those computations, we need to prevent
            # another run of computation with a second held right click. Else, a first plane could be computed depending to a
            # first rotatio, and a second one could depend from another rotation. This test prevent to begin new computations
            # before the actual ones are finished.

            # unchanged lines from rgl:::mouseTrackball()
            vp <- par3d("viewport")
            width <<- vp[3]
            height <<- vp[4]
            cur <<- rgl.cur()
            for (i in dev) {
                if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
                else userMatrix[[i]] <<- par3d("userMatrix")
            }
            rgl.set(cur, TRUE)
            rotBase <<- screenToVector(x, height - y)

            # added initializations for plane(s)
            if (is.null(uplanes)){
                # 1st rotation done: "uplanes" is initialized to initial "planes" (array with 3 points lying in the plane along 2 orthogonal vectors lying alos in the plane)
                uplanes<<-planes
                # "vInter" initialization which will contain necessary data for plane plottings
                vInter<<-list()
            }else{
                # nth rotation done : "planes" is initialized to "uplanes" (the last done rotation)
                planes<<-uplanes
            }

        }

    }

    # function called when the press of the right mouse button is held
    trackballUpdate <- function(x,y) {

        if (is.complete){

            is.complete<<-FALSE

            # unchanged lines from rgl:::mouseTrackball()
            rotCurrent <- screenToVector(x, height - y)
            angle <- angle(rotBase, rotCurrent) # rotation angle computation...
            axis <- xprod(rotBase, rotCurrent)  # ...& 3D rotation axis

            mouseMatrix <- rotationMatrix(angle, axis[1], axis[2], axis[3]) # corresponding rotation matrix
            for (i in dev) {
                if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
                else par3d(userMatrix = mouseMatrix %*% userMatrix[[i]]) # actually, the user point of view is rotated (userMatrix) whith the previous rotation matrix
                # => absolute mesh coordinates stay unchanged
            }
            rgl.set(cur, TRUE)

            # added lines

            # Because "mouseMatrix" affects "userMatrix", it seems that the axis coordinates aren't expressed in mesh absolute coordinates, but were already
            # partially processed with the 3D coordinate treatment (cf help for par3d()).
            # More precisely, we should be at the step where coordinates ae multiplied by "userMatrix".
            # To sum up the different steps:
            # - mesh absolute coordinates: vector (x,y,z)
            # - v = (x,y,z,1)
            # - v'=v*scale : scale depending from axis scaling: don't care about this here...
            # - v''= userMatrix*v*scale
            # so:
            # axis = userMatrix%*%axis_abs_coord
            # inv(userMatrix)%*%axis = axis_abs_coord
            # Express the axis rotation into absolute coordinates is needed to rotates the plane(s)

            axis<-c(axis,1)
            axis<-solve(par3d("userMatrix"))%*%axis
            axis<-axis[1:3]

            # computation of rotation matrix needed to rotate absolute coordinates (only the rotation axis needs this transformation, rotation angle kepping the same)
            uMat <- rotationMatrix(angle, axis[1], axis[2], axis[3])

            # delete all previous lines and planes (if any) for plot updating
            Ids<-rgl.ids()
            rgl.pop("shapes",Ids[Ids[,2]=="linestrip",1])
            rgl.pop("shapes",Ids[Ids[,2]=="planes",1])

            # loop to update and plot the mesh/plane intersections
            for (i in idxPlanes){

                # plane PC2-PC3 when i=1, PC1-PC3 when i=2, PC1-PC2 when i=3
                ci<-setdiff(1:3,i)

                # extract points coordinates from given plane
                A<-p1<-planes[1,,i]
                p2<-planes[2,,i]
                p3<-planes[3,,i]

                # compute plane vectors
                p2<-p2-p1
                p3<-p3-p1

                # uMat rotation
                p2<-c(t(p2)%*%uMat[1:3,1:3])+p1
                p3<-c(t(p3)%*%uMat[1:3,1:3])+p1

                # storage for later use in a next function call
                uplanes[,,i]<<-rbind(p1,p2,p3)

                # planes/mesh intersection
                tmp <- meshPlaneIntersect2(mesh, p1,p2,p3)
                inter<- tmp[[1]]
                edgesTot<-tmp[[2]]
                edges_in_cplx<-edgesTot[,1]+1i*edgesTot[,2] # expression of edges of intersected faces in a complex form for facility
                faces_in<-edgesTot[,3] # intersected faces
                is_bd<-as.logical(edgesTot[,4]) # indicates if edges of the mesh of intersection are or not border edges

                out <- .Call("sortCurveMesh", edges_in_cplx,faces_in,is_bd)

                # Convert "out" into a list where each element correpond to a submesh with at each time a vector indicating the order
                # of points. A NA value will be append at the end of each vector to allow the plot of lines 3d in a single pass (NA
                # value not being linked by lines3d)
                Lve<-list()
                for (j in 1:max(out[[2]])){
                    Lve[[j]]<-c(which(out[[2]]==j)[order(out[[1]][out[[2]]==j])],NA)
                }
                inter<-inter[unlist(Lve),]

                # plot of the intersection
                lines3d(inter, col = "red", lwd=2)

                # storage
                vInter[[i]]<<-inter
            }
            vPlanes<<-uplanes

            # loop for plane plotting
            for (i in idxPlanes) {
                # a normal vector to a plane can be computed through the vector cross product of 1 orthogonal vectors
                # contained in this plane
                A <- uplanes[1, , i]
                v1<-uplanes[2,,i]-A
                v2<-uplanes[3,,i]-A
                n<-xprod(v1,v2)

                # plane parameter
                d <- -t(n) %*% A

                # plane plotting
                planes3d(n[1], n[2], n[3], d, alpha = 0.7, col="cyan")
            }

        }
        is.complete<<-TRUE
    }

    # unchanged lines from rgl:::rgl.select()
    button <- match.arg(button)
    newhandler <- par3d("mouseMode")
    newhandler[button] <- "selecting"
    oldhandler <- par3d(mouseMode = newhandler)
    on.exit(par3d(mouseMode = oldhandler))

    # modification of mouse action when user right clicks: up to now right click zooms, and now it allows dissociated
    # rotation of mesh and planes.
    rMul<-rgl.setMouseCallbacks(2, begin = trackballBegin, update = trackballUpdate, end = NULL)

    # code execution is now stopped, with this new definition of right click action, until user presses escap to
    # validate plane positions
    dev<-rgl.cur()

    while (dev==rgl.cur()){
        result <- rgl:::rgl.selectstate()
        # if ESC -> get-out
        if (result$state >= rgl:::msDONE) break
    }

    # if the window is closed (shouldn't occur): we exit this function
    if (dev!=rgl.cur())
    {
        if (modify){
            isDone<-TRUE
            return(list(isDone=isDone,isClosed=TRUE))
        }
    }

    # Otherwise, the mouse action for right click is reinitialized to zoom
    rgl.setMouseCallbacks(2)

    # unchanged line from rgl:::rgl.select :
    rgl:::rgl.setselectstate("none")

    # Exports
    isDone <- TRUE
    if (result$state == rgl:::msDONE) isDone <- FALSE

    Ids<-rgl.ids()
    rgl.pop("shapes",Ids[Ids[,2]=="planes",1])

    return(list(isDone=isDone,isClosed=FALSE,vInter=vInter,vPlanes=vPlanes))
}

#############################################################

