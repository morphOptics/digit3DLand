###############################
checkColor<-function(M,message=as.character(deparse(substitute(M)))){
    t1<-is.element(c(M),colors())
    if (sum(t1)<length(t1)){
        stop(paste0(message," should be chosen within the values from colors()..."))
    }
}
###############################
checkAlpha<-function(M,message=as.character(deparse(substitute(M)))){
    if (!is.numeric(M) | any(is.na(M))){
        stop(paste0(message," should contain numeric values..."))
    }
    if (min(M)<0 | max(M)>1){
        stop(paste0(message," should contain values in [0,1]..."))
    }
}

###############################
checkLength<-function(M,len,message=as.character(deparse(substitute(M))),repV=TRUE){
    if (!is.element(length(M),len)){
        stop(paste0(message," should have length within {",paste(len,collapse=","),"}..."))
    }
    if (repV & length(M)>1 & length(M)<max(len)){
        stop(paste0(message," should have length within {",paste(range(len),collapse=","),"}..."))
    }
    if (repV & length(M)==1){
        M<-rep(M,max(len))
    }
    return(M)
}

###############################
checkLogical<-function(M,len,message=as.character(deparse(substitute(M)))){
    if(!is.logical(M) | any(is.na(M))){
        stop(paste0(message," should contain logical values..."))
    }
    M<-checkLength(M,len,message=message)
    return(M)
}

###############################
checkMat<-function(M,message=as.character(deparse(substitute(M)))){
    #?! (see Algebraic notation (chess) for meaning...)

    if (!is.element(class(matrix(1,2,2)),c("matrix","numeric"))){
        stop(paste0(message," should be either a matrix or vector..."))
    }
    if (length(M)>4 | length(M)==3){
        stop(paste0("Wrong value number specified for ",message,"..."))
    }
    if (length(M)==1){
        M<-matrix(M,2,2)
    }
    if (length(M)==2){
        M<-matrix(c(M),2,2)
    }
    return(M)
}



###############################
setDigitFixedOptions<-function( winNb=1, winSize= rbind(c(0,50,830,904), c(840,50,1672,904)), winSynchro=TRUE,
                                meshColor=rep("gray",2), meshAlpha=rep(1,2), meshShade=rep(TRUE,2),
                                meshPoints=c(FALSE,TRUE), meshWire=rep(FALSE,2),
                                PCplanesDraw=TRUE, PCplanesColor= "cyan", PCplanesAlpha=0.7,
                                intersectLines=TRUE, intersectPoints=FALSE, intersectColor="red",
                                spheresRad=1/50, spheresColor=matrix(c("black","blue"),2,2), spheresAlpha=1,
                                labelCex=2, labelColor="magenta", labelAdj=1/50,
                                zoomPercDist=0.15, zoomPtsDraw=TRUE, zoomPtsCol="red", zoomSeeLm=FALSE
                                )
{

    # Allow to set graphical options for DigitFixed.
    # The simplest use with no argument returns a list with default values for all settable parameters:
    # opt<-setDigitFixedOptions()
    #
    # Some or all parameters can be set individually (the ones not given will be set to default). Ex:
    # opt<-setDigitFixedOptions(winNb=2, winSynchro= FALSE, PCplanesDraw=c("pc2-pc3"))
    # opt<-setDigitFixedOptions(meshColor=c("gray","orange"), meshPoints=FALSE, zoomSeeLm=TRUE)
    # etc...
    #
    # Below possible settings for the different parameters:
    #
    # > window options
    #   * winNb: number of devices, possible settings:
    #     - 1 for a single device subdivided into 2 parts (one for the decimated mesh, the other for the zoomed mesh)
    #     - 2 for two separate devices
    #     Default: 1.
    #   * winSize: size & location of the device(s), possible settings:
    #     - a vector with 4 positive values indicating the left, top, right and bottom (in pixels, see windowRect
    #       parameter in rgl:::par3d()) for the device (when winNb=1)
    #     - a 2*4 matrix of positive values, each line indicating as before the left, top, right and bottom
    #       (in pixels) for each device (when winNb=2)
    #     Default: rbind(c(0,50,830,904), c(840,50,1672,904).
    #   * winSynchro: logical value indicating if user interaction (zoom, rotation) applied on a mesh should be
    #     synchronously applied on the second one (decimated or detailled mesh). Only Works for winNb=1.
    #     Default: TRUE
    #
    # > mesh options
    #   * meshColor: a character vector of length 1 or 2 indicating the color(s) for mesh plotting. Note: for the
    #     meshShade and the meshWire options, but not for the mesPoints option (see details after), this color
    #     won't overwrite the color stored in the material$color (if any) from the mesh object for plotting
    #     This color will be used only if the material$color is not given. Values for meshColor should be taken
    #     from colors().
    #     If two colors are given, the 2 meshes (the decimated one and the zoomed one) will be plotted with those
    #     two different colors, otherwise, the same color will be used for both ones.
    #     Default: rep("gray",2)
    #   * meshAlpha: a character vector of length 1 or 2 indicating the alpha value(s) for mesh plotting. It will
    #     overwrite the alpha value stored in material$alpha (if any). Values for meshAlpha should be taken within
    #     [0,1].
    #     If two alphas are given, the 2 meshes (the decimated one and the zoomed one) will be plotted with those
    #     two different alphas, otherwise, the same alpha will be used for both ones.
    #     Default: rep(1,2)
    #   * meshShade, meshWire, meshPoints: logical vector of length 1 or 2 indicating if the meshes should be plotted
    #     using this kind of representation. meshShade will use rgl:::shade3d, meshWire rgl:wire3d, meshPoints
    #     rgl:::points3d. The two values stand for the decimated and the zoomed mesh. If only one value is given,
    #     it will be recycled for the second mesh.
    #     Only one or any combinations of those kind of mesh representations are possible. If for one or both mesh(es)
    #     all values are set to FALSE, the meshShade representation will be used.
    #     Default: rep(TRUE,2) for MeshShade, rep(FALSE,2) for meshWire, c(FALSE,TRUE) for meshPoints
    #
    # > PC planes options
    #   * PCplanesDraw: logical or character vector of length 1, 2 or 3 indicating if major planes (from mesh principal
    #     components) as well as their intersections with the mesh should be plotted, and interactively set by user
    #     before the landmark placement step. Possible settings:
    #     - TRUE (or FALSE): in this case all (or no one of) the 3 major planes will be (or won't be) plotted.
    #     - any combination of 1 or 2 values within c("pc1-pc2","pc1-pc3","pc2-pc3") indicating which particular
    #       plane(s) should be plotted.
    #     Default: TRUE.
    #   * PCplanesColor, PCplanesAlpha : a character/numeric vector taking values within colors()/[1,2] indicating
    #     with which color/alpha each plane should be plotted. If only one value is given for more than one plane,
    #     color/alpha value are recycled.
    #     Default: "cyan", 0.7.
    #
    # > intersection line options
    #   * intersectLines, intersectPoints: a logical vector indicating for each intersection plane if the plotting
    #     of the intersection should use this kind of representation (using rgl:::lines3d/rgl:::points3d). If only
    #     one value is given for more than one plane, the logical values is recycled. Both representation can be used
    #     for the same plane. If for one or more plane, all values are set to FALSE, the intersectLines reprensentation
    #     will be used.
    #     Default: TRUE, FALSE.
    #     Note: for big meshes, the intersection plotting can be fasten by setting intersectLines to FALSE.
    #   * intersectColor: a character vector taking values within colors() indicating with which color each intersection
    #     should be plotted.
    #     Default: "red".
    #
    # > landmark sphere options
    #   * spheresRad, spheresColor, spheresAlpha: numerical/character/numerical vector or matrix taking values within
    #     [0,1]/colors()/[0,1] indicating with which radius/color/alpha the spheres for landmarks should be plotted.
    #     Possible settings:
    #     - a 2*2 matrix, the 1st line corresponding to the setting before user validation and the 2nd one after,
    #       the 1st column correspondng to the setting for the decimated mesh and the second one for the zoomed mesh
    #     - a vector of length 2 corresponding to the setting before user validation and the 2nd one after (those values
    #       will be recycled for the zoomed mesh)
    #     - a unique value for the setting before and after the user validation, and for the decimated and the zoomed mesh
    #     Default: 1/50, matrix(c("black","blue"),2,2), 1.
    #
    # > landmark label options
    #   * labelCex, labelColor: numerical/character vector or matrix taking values within [0,1]/colors() indicating
    #     with which size/color the landmark labels should be plotted. For details on possible settings: see e.g.
    #     spheresRad.
    #     Default: 2, "magenta".
    #   * labelAdj: numerical vector or matrix or array taking values within [0,1] indicating how to adjust the label
    #     location relative to the landmark sphere. Possible settings:
    #     - a 2*2*2 array, the 1st line corresponding to the setting before user validation and the 2nd one after,
    #       the 1st column correspondng to the setting for the decimated mesh and the second one for the zoomed mesh,
    #       the first slide corresponding to the horizontal adjustment end the second one to the vertical adjustment
    #     - a 2*2 matrix corresonding to the slide for horizontal adjustment which will be recycled for vertical adjustment
    #     - a vector of length 2 corresponding to the setting before user validation and the 2nd one after (those values
    #       will be recycled for the zoomed mesh and the vertical adjustment)
    #     - a unique value for the setting before and after the user validation (recyled for the zoomed mesh and the
    #       vertical adjustment)
    #     Default: 1/50.
    #
    # > zoom options
    #   * zoomPercDist: numerical value within [0,1] specifying the extent of the zoomed mesh. This extent is computed as
    #     the maximal distance between the clicked point and the mesh points multiplied by zoomPercDist.
    #     Default: 0.15.
    #     Note: for big meshes, the higher this value will be, the more slow down the computations will be (notably during
    #           user interaction through manual zoom and rotation with the zoomed mesh)
    #   * zoomPtsDraw: a logical value indicating if the exent of the zoomed mesh should be shown on the decimated mesh.
    #     This extent will be represented as a 3D point cloud.
    #     Default: TRUE.
    #   * zoomPtsCol: a character value taking values within colors() indicating with which color the 3D point cloud (if
    #     zoomPtsDraw is set to TRUE) should be plotted.
    #     Default: "red".
    #   * zoomSeeLm: a logical value indicating if the landmark placed by user on the decimated mesh should be visible
    #     on the zoomed mesh. It will fasten the process of manual digitizing enabling the direct validation (without
    #     any manual change) of the placed landmark, but at the risk of a more or less important approximation on the
    #     landmark positioning depending on the degree of decimation used for the decimated mesh.
    #     Default: FALSE.

    warn<-options()$warn
    options(warn=1)

    # Window options
    if (!is.numeric(winNb) | length(winNb)!=1 | !is.element(winNb[1],c(1,2))){
        stop("winNb should be a scalar taking values in {1,2}...")
    }
    if (!is.element(length(winSize),c(4,8))){
        message<-paste("With ",winNb," windows, winSize should contain ",4*winNb," integer values...",sep="")
        stop(message)
    }
    if (min(winSize)<0 | !is.numeric(winSize) | any(is.na(winSize))){
        stop("winSize should contain positive integer values...")
    }
    if (!all(floor(winSize) == winSize, na.rm = TRUE)){
        stop("winSize should contain positive integer values...")
    }
    if (winNb==2 & length(winSize)==4){
        tmp<-winSize
        winSize<-rbind(tmp,tmp+c(tmp[3],0,tmp[3],0))
        warning(paste0("The size and position for the second window wasn't specified... Set to: c(",toString(winSize[2,]),")"))
    }
    winSize<-matrix(c(t(winSize))[1:(4*winNb)],winNb,4,byrow=TRUE)
    if (min(winSize[,3]-winSize[,1])<=0 | min(winSize[,4]-winSize[,2])<=0){
        stop("Wrong specification for winSize: see the help for windowRect parameter in rgl:::par3d for details")
    }
    if (!is.logical(winSynchro) | is.na(winSynchro)){
        stop("winSynchro should be a logical value...")
    }
    if (winNb==2 & winSynchro){
        winSynchro<-FALSE
        warning("With 2 separate windows, winSynchro=TRUE is not supported... Set to FALSE")
    }
    winOptions<-list(winNb=winNb, winSize=winSize, winSynchro=winSynchro)

    # mesh options
    checkColor(meshColor)
    meshColor<-checkLength(meshColor,c(1,2))

    checkAlpha(meshAlpha)
    meshAlpha<-checkLength(meshAlpha,c(1,2))

    meshShade<-checkLogical(meshShade,c(1,2))
    meshPoints<-checkLogical(meshPoints,c(1,2))
    meshWire<-checkLogical(meshWire,c(1,2))
    Sm<-meshShade+meshPoints+meshWire
    if (min(Sm)==0){
        meshShade<-rep(TRUE,2)
        warning("Some mesh plot modes are missing! meshShade will be set to TRUE for both windows...")
    }
    meshOptions<-list(meshColor=meshColor,meshAlpha=meshAlpha,meshShade=meshShade,meshPoints=meshPoints,meshWire=meshWire)

    # PC planes options
    stopMessage<-"PCplanesDraw should be a logical value, or a character vector taking value within {\"pc1-pc2\",\"pc1-pc3\",\"pc2-pc3\"}..."
    if (any(is.na(PCplanesDraw))){
        stop(stopMessage)
    }
    if (!is.character(PCplanesDraw) & !is.logical(PCplanesDraw)){
        stop(stopMessage)
    }
    PCplanesDraw<-checkLength(PCplanesDraw,1:3,repV=FALSE)
    if (is.logical(PCplanesDraw) & length(PCplanesDraw)==1){
        nbPlanes<-0
        if (PCplanesDraw){
            nbPlanes<-3
        }
    }else{
        if (is.logical(PCplanesDraw)){
            stop(stopMessage)
        }
        V<-c("pc2-pc3","pc1-pc3","pc1-pc2")
        t1<-is.element(PCplanesDraw,V)
        nbPlanes<-sum(t1)
        if (sum(t1)!=length(t1)){
            stop(stopMessage)
        }
    }
    if (nbPlanes>0){
        checkColor(PCplanesColor)
        PCplanesColor<-checkLength(PCplanesColor,1:nbPlanes)
        checkAlpha(PCplanesAlpha)
        PCplanesAlpha<-checkLength(PCplanesAlpha,1:nbPlanes)
    }else{
        PCplanesColor<-PCplanesAlpha<-NULL
    }
    PCplanesOptions<-list(PCplanesDraw=PCplanesDraw,PCplanesColor=PCplanesColor,PCplanesAlpha=PCplanesAlpha)

    # intersect options
    if (nbPlanes>0){
        intersectLines<-checkLogical(intersectLines,1:nbPlanes)
        intersectPoints<-checkLogical(intersectPoints,1:nbPlanes)
        if (sum(intersectLines)==0 & sum(intersectPoints)==0){
            intersectLines<-rep(TRUE,nbPlanes)
            warning("To be visible, mesh/plane intersections should plotted at least by lines: instersectLines was set to TRUE...")
        }
        checkColor(intersectColor)
        intersectColor<-checkLength(intersectColor,1:nbPlanes)

    }else{
        intersectLines<-intersectPoints<-FALSE
        intersectColor<-NULL
    }
    intersectOptions<-list(intersectLines=intersectLines,intersectPoints=intersectPoints,intersectColor=intersectColor)

    # sphere options
    if (!is.numeric(spheresRad) | any(spheresRad < 0)){
        stop("spheresRad should contained positive numerical values...")
    }
    if (max(spheresRad)>1){
        warning("spheresRad should be a fractional number < 1 to allow readable plots")
    }
    spheresRad<-checkMat(spheresRad)
    checkColor(spheresColor)
    spheresColor<-checkMat(spheresColor)
    checkAlpha(spheresAlpha)
    spheresAlpha<-checkMat(spheresAlpha)
    spheresOptions<-list(spheresRad=spheresRad,spheresColor=spheresColor,spheresAlpha=spheresAlpha)

    # label options
    if (any(is.na(labelCex))){
        stop("labelCex should contained numerical values...")
    }
    if (!is.numeric(labelCex) | min(labelCex)<0){
        stop("labelCex should contained positive numerical values...")
    }
    labelCex<-checkMat(labelCex)
    checkColor(labelColor)
    labelColor<-checkMat(labelColor)
    if (!is.numeric(labelAdj) | any(is.na(labelAdj))){
        stop("labelAdj should be numeric...")
    }
    if (is.element(length(labelAdj),c(1,2,4))){
        labelAdj<-array(rep(c(labelAdj),8/length(labelAdj)),c(2,2,2))
    }else if(!identical(dim(labelAdj), as.integer(c(2,2,2)) )){
        stop("labelAdj should be a 2*2*2 array or a vector/matriw with 1, 2 or 4 values...")
    }
    if (max(labelAdj)>1){
        warning("labelAdj should contain fractional values to allow readable plots")
    }
    labelOptions<-list(labelCex=labelCex,labelColor=labelColor,labelAdj=labelAdj)

    # zoom options
    if (length(zoomPercDist)!=1 | !is.numeric(zoomPercDist) | min(zoomPercDist)<0 | max(zoomPercDist)>1){
        stop("zoomPercDist should be a numerical value between 0 and 1...")
    }
    if (length(zoomPtsDraw)!=1 | !is.logical(zoomPtsDraw) | any(is.na(zoomPtsDraw))){
        stop("zoomPtsDraw should be a logical value...")
    }
    if (zoomPtsDraw){
        if (length(zoomPtsCol)!=1 | any(!is.element(zoomPtsCol,colors()))){
            stop("zoomPtsCol should be a unique value chosen within the values from colors()...")
        }
    }else{
        zoomPtsCol<-NULL
    }
    if (length(zoomSeeLm)!=1 | !is.logical(zoomSeeLm) | any(is.na(zoomSeeLm))){
        stop("zoomSeeLm should be a logical value...")
    }
    zoomOptions<-list(zoomPercDist=zoomPercDist,zoomPtsDraw=zoomPtsDraw,zoomPtsCol=zoomPtsCol,zoomSeeLm=zoomSeeLm)

    options(warn=warn)

    return(list(winOptions=winOptions, meshOptions=meshOptions, PCplanesOptions=PCplanesOptions,
                intersectOptions=intersectOptions, spheresOptions=spheresOptions, labelOptions=labelOptions,
                zoomOptions=zoomOptions))
}
