/** @file DEstimator_IFace.cpp
 * @author Mark J. Olah (mjo\@cs.unm DOT edu)
 * @date 2014-2019
 * @brief The class declaration and inline and templated functions for DEstimator_IFace.
 */
#include "DEstimator_IFace.h"

DEstimator_IFace<double> iface; /!< Global iface object provides a iface.mexFunction

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, const mxArray *rhs[])
{
    iface.mexFunction(nlhs, lhs, nrhs, rhs);
}
