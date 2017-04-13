#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

///////////////////////////////////////////////

IntegerVector which(LogicalVector vLog);

IntegerVector which(LogicalVector vLog){

    // equivalent to the R which function: convert a logical vector into TRUE indices

    int S = sum(vLog);

    IntegerVector idx(S);
    int cpt = 0;
    for (int i=0; i<vLog.length(); i++){
        if (vLog[i]){
            idx[cpt]=i;
            cpt++;
        }
    }

    return (idx);
}

///////////////////////////////////////////////

// TO DO: the 3 following functions checkVal, checkVal2 & checkValCplx are striclty identical in their rationale,
// ie testing in a vector the equality of a given value (Vector==given_value in R) and returning a logical vector,
// but differ only in the kind of input value:
// - integer vector & a case of an integer vector in checkVal
// - integer vector & a case of an integer in checkVal2
// - complex vector & a case of a complex vector in checkValCplx
// A fusion of those 3 functions would be good...

LogicalVector checkVal(IntegerVector vec, IntegerVector val);

LogicalVector checkVal(IntegerVector vec, IntegerVector val){

    int n=vec.length();
    LogicalVector Vout(n);

    for (int i=0; i<n; i++){
        if (vec[i]==val[0]){
            Vout[i]=TRUE;
        }
    }

    return(Vout);

}

///////////////////////////////////////////////

LogicalVector checkVal2(IntegerVector vec, int val);

LogicalVector checkVal2(IntegerVector vec, int val){

    int n=vec.length();
    LogicalVector Vout(n);

    for (int i=0; i<n; i++){
        if (vec[i]==val){
            Vout[i]=TRUE;
        }
    }

    return(Vout);

}

///////////////////////////////////////////////

LogicalVector checkValCplx(ComplexVector vec, ComplexVector val);

LogicalVector checkValCplx(ComplexVector vec, ComplexVector val){

    int n=vec.length();
    LogicalVector Vout(n);

    for (int i=0; i<n; i++){
        if (vec[i]==val[0]){
            Vout[i]=TRUE;
        }

    }

    return(Vout);

}

///////////////////////////////////////////////

RcppExport SEXP sortCurveMesh(SEXP edges_in_cplx_, SEXP faces_in_, SEXP is_bd_) {

    // From:
    // - a set of edges expressed in a complex form, all being intersected by a common plane (edges_in_cplx)
    // - the set of faces to which those edges belong (faces_in)
    // - the information for each edge if it is or not an edge border (is_bd)
    // this function determines in which order the edges should be sorted according to the connectivity of
    // the faces they belong
    // The algorithm works in this way:
    // 1) if some border edges exist:
    //    a. start from a border edge (combining info from edges_in_cplx & is_bd)
    //    b. determine to which face it belong (info from faces_in)
    //    c. determine the other edge defining this face (combining info from faces_in & edges_in_cplx)
    //    d. repeat steps a, b & c until another border edge is reached
    //    e. if all edges were browsed during this search, the stops, otherwise, repeat step 1) or 2) until all edges are browsed
    // 2) otherwise:
    //    a'. start from any non-border edge
    //    b'. determine to which face it belongs (info from faces_in)
    //    c'. determine the other edge defining this face (combining info from faces_in & edges_in_cplx)
    //    d'. repeat steps a', b' & c' until the initial non-border edge is reached
    //    e'. if all edges were browsed during this search, then stops, otherwise, repeat 2) until all edges are browsed
    //
    // !!!WARNING: This algorithm won't work with non-manifold mesh geometry, in particular with non-manifold edges, ie
    // edges being shared by more than 2 faces. In such a case, the execution of the R function calling sortCurveMesh will
    // generate a crash of RStudio
    //
    //Outputs: a list containing 2 vectors:
    // - index_submesh: for each edge, gives to which isolated submesh (if any) the edge belongs
    // - index_vertex: for each edge of a given submesh (if any), gives the order with which edges should be browsed

    // input declaration
    ComplexVector edges_in_cplx(edges_in_cplx_);
    IntegerVector faces_in(faces_in_);
    LogicalVector is_bd(is_bd_);

    // numbers of border edges
    IntegerVector num_edge_bd = which(is_bd);

    // nb of edges
    int N = is_bd.length();

    // declaration & initialization of some intermediate and result objets
    List ret;
    IntegerVector vOut1(N);
    IntegerVector vOut2(N);
    IntegerVector rest = seq(0,N-1);
    int nOut=0;
    bool Stop=FALSE;
    IntegerVector start_edge=0;
    LogicalVector cV;

    // main loop over all the nOut submeshes
    while(Stop==FALSE){

        // counter for nb of submeshes
        nOut++;

        // Check if some border edges exist...
        if (num_edge_bd.length()>0){
            // ...if so, the edge browsing will start by one of them
            start_edge=num_edge_bd[0];
            num_edge_bd[0]=-1; // this case set to -1 will be deleted in the following
        }else{
            // ...else, the edge browsing will start by any of the remaining edge
            start_edge=rest[which(rest>=0)[0]];
        }

        // declaration & initialization
        IntegerVector ve(N);
        ve=ve-1;
        bool Stop2=FALSE;
        IntegerVector start_edge_init=start_edge;
        ve[0]=start_edge[0]; // TO DO: convert start_edge into int to avoid [0]...
        int cpt=0;

        // second loop over all the edges in the considered submesh
        while(Stop2==FALSE){

            // counter for nb of edges
            cpt++;

            // update the "rest" vector recording which edges were browsed
            rest[start_edge[0]]=-1;

            // find where the face containing the edge "start_edge" appear in "faces_in"
            cV= checkVal(faces_in,faces_in[start_edge]);
            IntegerVector idxF = which(cV);

            // find the line number where this face appears (excluding the line with "start_edge")
            cV=!checkVal(idxF,start_edge);
            IntegerVector Tmp=which(cV);

            // update "start_edge" as the other edge of the considered face
            cV=!checkVal(idxF,start_edge);
            start_edge=idxF[which(cV)];

            if (start_edge.length()>1){
                // shouldn't occur...
                return wrap("Fail1");
            }else{

                // update the "rest" vector recording which edges were browsed
                rest[start_edge[0]]=-1;

                // store the current edge in "ve"
                ve[cpt]=start_edge[0];

                if (edges_in_cplx[start_edge[0]]==edges_in_cplx[start_edge_init[0]]){
                    // The current edge is the same as the initial edge: a complete circuit was performed.
                    // The second while loop can be stopped
                    Stop2=TRUE;
                }else{

                    // Test if the current edge is a border edge...
                    LogicalVector Cond=is_bd[start_edge];

                    if (Cond[0]){
                        // ... if so, the second while loop can be stopped
                        Stop2=TRUE;

                        // and "num_edge_bd" is updated with this current border edge
                        cV=checkVal(num_edge_bd,start_edge);
                        num_edge_bd[which(cV)]=-1;

                    }else{
                        // ... otherwise, the loop goes on

                        // find where the current edge "start_edge" appear in "edges_in_cplx"
                        cV=checkValCplx(edges_in_cplx,edges_in_cplx[start_edge]);
                        IntegerVector idxE = which(cV);

                        // find the line number where this edge appears (excluding the line with "start_edge")
                        cV=!checkVal(idxE,start_edge);
                        IntegerVector Tmp=which(cV);

                        // update "start_edge" as the other occurence of this edge (actually, start_edge stay the
                        // same, but at this line, start_edge will be assosicated to another face)
                        start_edge=idxE[which(!checkVal(idxE,start_edge))];
                    }
                }
            }
        }

        // cases of "ve" set to -1 are removed
        cV=!checkVal2(ve,-1);
        ve=ve[which(cV)];

        // store results of order of edges and to which submesh thy belongs on "vOut1" and "vOut2"
        IntegerVector vs=seq(1,ve.length());
        vOut1[ve]=vs;
        vOut2[ve]=nOut;

        // cases of "num_edge_bd" set to -1 are removed
        cV=!checkVal2(num_edge_bd,-1);
        num_edge_bd=num_edge_bd[which(cV)];

        // Test if all edges were browsed...
        if (max(rest)==-1){
            //... if so, the main while loop can be stopped
            Stop=TRUE;
        }
    }

    // return list of results
    ret["index_vertex"] = vOut1;
    ret["index_submesh"] = vOut2;

    return wrap(ret);
}

