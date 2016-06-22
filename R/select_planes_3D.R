DrawOrthoplanes <- function(mesh) {
    if (class(mesh) != "mesh3d") stop("mesh should be an object \"mesh3d\".")
    # centering of vertex coordinates
    pts <- t(mesh$vb[1:3, ])
    A <- apply(pts, 2, mean)
    pts <- sweep(pts, 2, A)
    # Planes
    sv <- svd(crossprod(pts))

    # Compute and draw the intersection mesh/planes
    ptsPlanes <- array(NA, c(3,3,3))
    Pl <- pl <- list()
    for (i in 1:3){
        # plane PC2-PC3 when i=1, PC1-PC3 when i=2, PC1-PC2 when i=3
        ci <- setdiff(1:3, i)
        # intersection planes/mesh
        inter <- meshPlaneIntersect(mesh, A, A+sv$u[,ci[1]], A+sv$u[,ci[2]])
        # Centering
        int <- sweep(inter, 2, A)
        # Rotation in order that the intersection pts lied in the plan
        int2 <- (int %*% svd(crossprod(int))$u)[, 1:2]
        # Order according to the angle of each vector (needed before to used lines3d())
        idx <- order(Arg(int2[,1] + 1i * int2[,2]))
        inter <- inter[idx, ]
        # plot of the intersection
        Pl[[i]] <- lines3d(inter, col = "red", lwd=2)
        # stockage
        #vInter[[i]] <- inter
        ptsPlanes[, , i] <- rbind(A, A+sv$u[,ci[1]], A+sv$u[,ci[2]])
        # for each plane, contains 3 points in the plane (needed just after)
    }

    # Draw orthogonal planes via planes3d()
    # coefficients a,b,c,d of the plane equation can be determined by a point of the plane
    # and a normal vector to this plane
    # https://fr.wikipedia.org/wiki/Plan_%28math%C3%A9matiques%29#D.C3.A9finition_par_un_vecteur_normal_et_un_point
    A <- ptsPlanes[1, , 1]
    for (i in 1:3) {
        if (i == 1) {
            # vect norm to plane 2-3
            n <- ptsPlanes[2, , 3] - ptsPlanes[1, , 3] #Nico ??? why not A and 3
        } else {
            # vect norm to plan 1-2 and plan 1-3
            n <- ptsPlanes[i, , 1] - A
        }
        d <- -t(n) %*% A
        pl[[i]] <- planes3d(n[1], n[2], n[3], d, alpha = 0.7, col="cyan")
    }
    # User interaction: manual rotation of the mesh (the orthogonal planes being fixed)
    # until the wanted orientation
    return(RotateMeshPlane3d(mesh, planes = ptsPlanes, Pl=Pl, pl=pl))
}

RotateMeshPlane3d <- function(mesh, planes, Pl, pl) {
    Stop <- 0
    # stop when the rotation is validated by ESC
    while (Stop==0) {
        temp <- selectPlanes(button="right", mesh=mesh, planes=planes, Pl=Pl, pl=pl)
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
selectPlanes<-function (button = c("left", "middle", "right"), mesh, planes,
                        Pl=NULL, pl=NULL, dev = rgl.cur(), uplanes=NULL) {

    # R?cup?ration de la fonction mouseTrackball() de rgl (via les binary) utilis?e pour effectuer la rotation manuel du mesh
    # compl?t?e ici pour qu'ne fonction de cette rotation (via le bouton droit de la souris), les plans orthogonaux restent "fixes" par rapport ? cette rotation.
    # En r?alit?, c'est bien le mesh qui reste fixe en coordonn?es absolues (car avec le mesh ,les axes tournent), et ce sont des nouveaux plans qui sont calcul?s ? chaque fois.

    # initialisations r?cup?r?es de mouseTrackball()
    width <- height <- rotBase <- NULL
    userMatrix <- list()
    cur <- rgl.cur()

    vInter<-NULL
    vPlanes<-NULL

    # fonctions utilitaires r?cup?r?es de mouseTrackball()
    screenToVector <- function(x, y) {
        # doit permettre ? partir d'un point 2D ?cran de calculer un vecteur
        # entre 2 positions successives de la souris, fonction appel?e 2 fois, et les deux vecteurs permettent de d?terminer un angle et un axe de rotation 3D
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

    # fonction de d?but de clic (droit)
    trackballBegin <- function(x, y) {

        # lignes inchang?es de mouseTrackball()
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

        # initialisations rajout?es pour les plans orthogonaux
        if (is.null(uplanes)){
            # premi?re rotation effectu? : on initialise uplanes (updated planes) ? planes initial (array contenant pour chaque plan des points contenus dans le plan le long de 2 vecteurs orthogonaux)
            uplanes<<-planes
            # vInter et vPlanes contiendront les donn?es n?cessaires pour dessiner les plans et intersections
            vInter<<-list()
            vPlanes<<-matrix(NA,3,4)
        }else{
            # ni?me rotation effectu? : planes est initialis? ? uplanes (derni?re rotation effectu?e)
            planes<<-uplanes
        }
    }

    # fonction de maintien du clic (droit)
    trackballUpdate <- function(x,y) {

        # lignes inchang?es de mouseTrackball()
        rotCurrent <- screenToVector(x, height - y)
        angle <- angle(rotBase, rotCurrent) # d?termination de l'angle...
        axis <- xprod(rotBase, rotCurrent)  # ...et de l'axe de rotation

        mouseMatrix <- rotationMatrix(angle, axis[1], axis[2], axis[3]) # calcul de la matrice de rotation correspondante
        for (i in dev) {
            if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
            else par3d(userMatrix = mouseMatrix %*% userMatrix[[i]]) # en fait ce qui tourne : c'est le point de vue utilisateur (userMatrix) par la matrice de rotation pr?c?dente
            # => les coordonn?es r?elles du mesh restent inchang?es car les axes tournent en m?me temps
        }
        rgl.set(cur, TRUE)

        # lignes rajout?es

        # mouseMatrix affectant la userMatrix, je suppose donc que les coordonn?es d'axe ne sont pas en coordonn?es absolues du mesh, mais ont d?j? subit une partie du traitement des coordonn?es 3D (cf aide de par3d())
        # Plus pr?cis?ment, on doit en ?tre ? l'?tape o? les coordonn?es sont multipli?es par userMatrix
        # Si on reprend les ?tapes :
        # - coordonn?es absolues mesh : vector (x,y,z)
        # - v = (x,y,z,1)
        # - v'=v*scale : scale d?pendant de scaling des axes, pas d'importance ici...
        # - v''= userMatrix*v*scale
        # donc :
        # axis = userMatrix%*%axis_abs_coord
        # inv(userMatrix)%*%axis = axis_abs_coord
        axis<-c(axis,1)
        axis<-solve(par3d("userMatrix"))%*%axis
        axis<-axis[1:3]
        # d?termination de la matrice de rotation ? appliquer aux coordonn?es absolues (axis modifi?, angle restant le m?me)
        uMat <- rotationMatrix(angle, axis[1], axis[2], axis[3])

        # boucle pour d?terminer et dessiner les nouvelles intersections
        for (i in 1:3){

            # plan consid?r? : plan PC2-PC3 si i=1, PC1-PC3 si i=2, PC1-PC2 si i=3
            ci<-setdiff(1:3,i)

            # extraction des valeurs du plan consid?r?
            A<-p1<-planes[1,,i]
            p2<-planes[2,,i]
            p3<-planes[3,,i]

            # calcul des vecteurs du plan
            p2<-p2-p1
            p3<-p3-p1

            # rotation par uMat
            p2<-c(t(p2)%*%uMat[1:3,1:3])+p1
            p3<-c(t(p3)%*%uMat[1:3,1:3])+p1

            # stcokage pour utilisation ult?rieure d'un nouvel appel de la fonction
            uplanes[,,i]<<-rbind(p1,p2,p3)

            # intersection plan/mesh
            inter<- meshPlaneIntersect(mesh,p1,p2,p3)

            # les points fournis par meshPlane Intersect ne sont pas tri?s (pas g?nant quand on utilise points3d(), mais ?a l'est avec lines3d()...)
            # => on trie les points en fonction de l'angle qu'ils forment dans le plan de l'intersection
            # centrage des points
            int<-inter-matrix(A,nrow=dim(inter)[1],ncol=3,byrow=TRUE)
            # rotation pour que les points d'intersection soient dans un plan
            int2<-(int%*%svd(t(int)%*%int)$u)[,1:2]
            # d?termination du rang de l'angle form? par chaque vecteur point d'intersection (via les complexes)
            idx<-order(Arg(int2[,1]+1i*int2[,2]))
            # tri des points
            inter<-inter[idx,]

            # effacement des plans et intersections pr?c?dents
            rgl.pop("shapes",Pl[[i]])
            rgl.pop("shapes",pl[[i]])

            # plot des nouvelles intersections et stockage
            Pl[[i]]<<-lines3d(inter,col="red",lwd=2)
            vInter[[i]]<<-inter

        }

        # boucle pour dessiner les plans
        for (i in 1:3){
            A <- uplanes[1,,1]
            if (i==1){
                # vect norm ? plan 2-3
                n <- uplanes[2,,3] - uplanes[1,,3]
            } else {
                n <- uplanes[i,,1] - A
            }

            # calcul des param?tres du plan
            d<- -t(n)%*%A
            # plot des nouveaux plans et stockage
            pl[[i]]<<-planes3d(n[1],n[2],n[3],d,alpha=0.7,col="cyan")
            vPlanes[i,]<<-c(n[1:3],d)
        }
    }

    # lignes suivantes : reprises de rgl.select()
    button <- match.arg(button)
    newhandler <- par3d("mouseMode")
    newhandler[button] <- "selecting"
    oldhandler <- par3d(mouseMode = newhandler)
    on.exit(par3d(mouseMode = oldhandler))

    # modification de l'action d'interface quand l'utilisateur clique sur le bouton droit de la souris :
    # avant zoom, maintenant pla?age d'un point par magn?tisme (utilisation des sous-fonctions Begin, Update et End)
    rMul<-rgl.setMouseCallbacks(2, begin = trackballBegin, update = trackballUpdate, end = NULL)

    # il faut maintenant suspendre l'?xecution du code et maintenir cette d?finition du clic droit tant que l'utilisateur n'a pas appuy? sur la touche ?chappe
    dev<-rgl.cur()

    while (dev==rgl.cur()){
        result <- rgl:::rgl.selectstate()
        # if ESC -> get-out
        if (result$state >= rgl:::msDONE) break
    }

    # si la fen?tre a ?t? ferm?e (?a ne devrait pas arriver) : on sort de cette fonction
    if (dev!=rgl.cur())
    {
        if (modify){
            isDone<-TRUE
            return(list(isDone=isDone,isClosed=TRUE))
        }
    }

    # Sinon
    # on red?finit l'action du clic droit par l'action par d?faut comme ?tant le zoom
    rgl.setMouseCallbacks(2)
    # ligne de rgl.select :
    rgl:::rgl.setselectstate("none")

    # Exports
    isDone <- TRUE
    if (result$state == rgl:::msDONE) isDone <- FALSE
    return(list(isDone=isDone,isClosed=FALSE,vInter=vInter,vPlanes=vPlanes))
}

#############################################################

