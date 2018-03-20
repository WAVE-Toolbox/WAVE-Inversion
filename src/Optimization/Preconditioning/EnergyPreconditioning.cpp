#include "EnergyPreconditioning.hpp"
// #include <cmath> // for asinh()

/*! \brief Initialize private members of the object
 *
 *
 \param dist Distribution of all DenseVectors with model size
 \param config Configuration
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::init(scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config)
{
    approxHessian.allocate(dist);
    approxHessian.assign(0.0);
    
    wavefieldSeismogramType.allocate(dist);
    wavefieldSeismogramType.assign(0.0);
    
    saveApproxHessian = config.get<bool>("saveApproxHessian");  
    approxHessianName = config.get<std::string>("approxHessianName");  
    epsilonHessian = config.get<ValueType>("epsilonHessian");
 
    
    /* Calculate approximation of receiver Green's function (migration weight K3 after Plessix and Mulder, 2004)
       Attention: only valid, if z_r = 0 (depth) for all receivers! */
    
//     approxRecGreensFct.allocate(dist);
//     approxRecGreensFct.assign(0.0);
    
    /* Allocate approximation of the receiver Green's function: 
       It (x_r,min and x_r,max) depends on the shot if the shots use different reicevers
       A loop over the model is necessary:
       for ...
          index2coordinate(1Dcoordinate, NX,Ny,NZ);
          ( asinh() - asinh() )                        */
    
}

/*! \brief Calculate the approximation of the diagonal of the inverse of the Hessian for one shot
 *
 *
 \param wavefield Wavefield of one time step
 \param DT temporal sampling intervall
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::intSquaredWavefields(KITGPI::Wavefields::Wavefields<ValueType> &wavefield, ValueType DT)
{
    
    /* Function should use the wavefield of the seismogram type which is used  to calculate the misfit! */
    
//     receivers.getSeismogramTypes();
    
    wavefieldSeismogramType = wavefield.getRefP();
    wavefieldSeismogramType *= wavefieldSeismogramType;
    wavefieldSeismogramType *= DT;
    approxHessian += wavefieldSeismogramType;
    
}

/*! \brief Apply the approximation of the diagonal of the inverse of the Hessian for one shot
 *
 *
 \param gradientPerShot 
 \param shotNumber 
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, IndexType shotNumber)
{
//     approxHessian *= approxRecGreensFct;
//     sqrt(approxHessian) missing because of |u_i| (see old IFOS)??
    
    /* Stabilize Hessian for inversion (of diagonal matrix) and normalize Hessian */
    approxHessian += epsilonHessian*approxHessian.maxNorm(); 
    approxHessian *= 1 / approxHessian.maxNorm(); 
    
    if(saveApproxHessian==1){
        approxHessian.writeToFile(approxHessianName + ".shot_" + std::to_string(shotNumber+1) + ".mtx");}
    
    approxHessian.invert();
    gradientPerShot *= approxHessian;    // overload operator /= in gradient-class    
    
}

template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::resetApproxHessian()
{
    approxHessian.assign(0.0);
    
}

template class KITGPI::Preconditioning::EnergyPreconditioning<double>;
template class KITGPI::Preconditioning::EnergyPreconditioning<float>;
