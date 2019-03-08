/** @file DEstimator_IFace.h
 * @author Mark J. Olah (mjo\@cs.unm DOT edu)
 * @author Peter Relich (physx.grad\@gmail.com)
 * @date 2014-2019
 * @brief The class declaration and inline and templated functions for DEstimator_IFace.
 */

#ifndef DESTIMATOR_DESTIMATOR_IFACE
#define DESTIMATOR_DESTIMATOR_IFACE

#include "MexIFace/MexIFace.h"
#include "DEstimator/DEstimator.h"

template<class FloatT>
class DEstimator_IFace : public mexiface::MexIFace,
                         public mexiface::MexIFaceHandler<destimator::Destimator<FloatT>>
{
public:
    DEstimator_IFace();
private:
    using mexiface::MexIFaceHandler<destimator::DEstimator<FloatT>>::obj;

    //Constructor
    void objConstruct() override;

    //Exposed method calls
    void objLLH();
    void objLLHdim();
    void objStaticLLH(const std::string &method);
};

template<class FloatT>
DEstimator_IFace<FloatT>::DEstimator_IFace()
{
    methodmap["LLH"] = std::bind(&DEstimator_IFace::objLLH, this);
    methodmap["LLHdim"] = std::bind(&DEstimator_IFace::objLLHdim, this);
    staticmethodmap["LLH_laplace1D"] = std::bind(&DEstimator_IFace::objStaticLLH, this, "laplace");
    staticmethodmap["LLH_recursive1D"] = std::bind(&DEstimator_IFace::objStaticLLH, this, "recursive");
    staticmethodmap["LLH_markov1D"] = std::bind(&DEstimator_IFace::objStaticLLH, this, "markov");
}

template<class FloatT>
void DEstimator_IFace<FloatT>::objConstruct()
{
    // [in] Obs: N x Ndim double.  Vector of observed positions
    // [in] T: N x 1 double.  Vector of observation times
    // [in] SE: N x Ndim double.  Vector of observation standard errors
    // [in] exposureT: scalar double.  The duration over which an observation is made
    checkNumArgs(1,4);
    auto Obs = getMat<FloatT>();
    auto T = getVec<FloatT>();
    auto SE = getMat<FloatT>();
    auto exposureT = getAsFloat<<FloatT>();
    this->outputHandle(new DEstimator<FloatT>(Obs, T, SE, exposureT));
}

template<class FloatT>
void DEstimator_IFace<FloatT>::objLLH()
{
    // [in] D: double vector of diffusion constants to estimate LLH for
    // [out] llh: double vector of log-likelihood for each given D value.
    checkNumArgs(1,1);
    auto D = getVec<FloatT>();
    auto llh = makeOutputArray<FloatT>(D.n_elem);
    if(D.n_elem==1) {
        llh(0) = obj->LLH(D(0));
    } else {
        obj->LLH(D, llh); //parallelized
    }
}

template<class FloatT>
void DEstimator_IFace<FloatT>::objLLHdim()
{
    // [in] D: double vector of diffusion constants to estimate LLH for
    // [in] dim: scalar integer giving index of dimension to estimate LLH for.  dim is zero-indexed:  0<=dim<Ndim;
    // [out] llh: double vector of log-likelihood for each given D value only considering dimension dim
    checkNumArgs(1,2);
    auto D = getVec<FloatT>();
    auto dim = getAsUnsigned<IdxT>();
    auto llh = makeOutputArray<FloatT>(D.n_elem);
    if(D.n_elem==1) {
        llh(0) = obj->LLHdim(D(0),dim);
    } else {
        obj->LLHdim(D, dim, llh); //parallelized
    }
}

template<class FloatT>
void DEstimator_IFace<FloatT>::objStaticLLH(const std::string &method)
{
    // This method will be exposed as static methods LLH_laplace1D, LLH_recursive1D, LLH_markov1D,
    // Each uses the same arguments but a different underlying algorithm
    // [in] D: double vector of diffusion constants to estimate LLH for
    // [in] Obs: N x 1 double.  Vector of observed 1D positions
    // [in] T: N x 1 double.  Vector of observation times
    // [in] SE: N x 1 double.  Vector of 1D observation standard errors
    // [in] exposureT: scalar double.  The duration over which an observation is made
    // [out] llh: double vector of log-likelihood for each given D value.
    checkNumArgs(1,5);
    auto D = getVec<FloatT>();
    auto Obs = getVec<FloatT>();
    auto T = getVec<FloatT>();
    auto SE = getVec<FloatT>();
    auto exposureT = getAsFloat<FloatT>();
    auto llh = makeOuputVer<FloatT>(D.n_elem);
    if(!method.compare("laplace")){
        DEstimator<FloatT>::LLH_laplace1D(D, Obs, T, SE, exposureT, llh);
    } else if(!method.compare("recursive")){
        DEstimator<FloatT>::LLH_recursive1D(D, Obs, T, SE, exposureT, llh);
    } else if(!method.compare("markov")){
        DEstimator<FloatT>::LLH_markov1D(D, Obs, T, SE, exposureT, llh);
    }
}

#endif /* DESTIMATOR_DESTIMATOR_IFACE */
