#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

///////////////////////////////////////////////

double angcalcArma(colvec a, colvec b); // called within doozer.h file of Morpho src

double angcalcArma(vec a, vec b) { // unchanged function from doozer.cpp file of Morpho src
  try {
    double alen = norm(a,2);
    double blen = norm(b,2);
    if (alen > 0)
      a = a/alen;
    if (blen > 0)
      b = b/blen;
    vec diffvec = a-b;
    double angle = acos((dot(diffvec,diffvec)-2)/-2);
    return angle;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

///////////////////////////////////////////////

RcppExport SEXP edgePlane2(SEXP vb_, SEXP diff_, SEXP edges_) {

    // modified from edgePlane.cpp of Morpho src
    // this version adds the export of the edge numbers being intersected

    try {
        IntegerMatrix edges(edges_);
        NumericMatrix vb(vb_);
        NumericMatrix diff(diff_);
        unsigned int nedges = edges.nrow();
        mat out(nedges,3); out.zeros();
        std::vector<unsigned int> test;
        for (unsigned int i = 0; i < nedges; i++) {
            int i1 = edges(i,0);
            int i2 = edges(i,1);
            vec vb2 =  vb(i2,_);
            vec diff0 = diff(i2,_);
            double ancath = sqrt(dot(diff0,diff0));

            vec resvec = vb(i1,_)-vb(i2,_);
            double angle = angcalcArma(diff0,resvec);
            double lres = sqrt(dot(resvec,resvec));
            resvec = resvec/lres;
            double hypoth = ancath/cos(angle);
            if (hypoth <= lres && hypoth >= 0) {
                out.row(i) = conv_to<rowvec>::from(vb2+hypoth*resvec);
                test.push_back(i);
            }
        }
        uvec myinds = conv_to<uvec>::from(test);
        vec myidxs =  conv_to<vec>::from(myinds); // added from edgePlane.cpp
        out = out.rows(myinds);
        //a 4th column containg intersected edge numbers is added to "out":
        out.insert_cols(3, myidxs+1.0); // added from edgePlane.cpp

        return wrap(out);
    } catch (std::exception& e) {
        ::Rf_error( e.what());
    } catch (...) {
        ::Rf_error("unknown exception");
    }
}
