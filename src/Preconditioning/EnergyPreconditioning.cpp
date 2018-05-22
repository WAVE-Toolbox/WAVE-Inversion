#include "EnergyPreconditioning.hpp"

/*! \brief Initialize private members of the object
 *
 *
 \param dist Distribution of all DenseVectors with model size
 \param config Configuration
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::init(scai::dmemo::DistributionPtr dist, KITGPI::Configuration::Configuration config)
{
    approxHessian.setSameValue(dist, 0 );

    wavefieldVX.setSameValue(dist, 0 );
    wavefieldVY.setSameValue(dist, 0 );
    wavefieldVZ.setSameValue(dist, 0 );
    
    saveApproxHessian = config.get<bool>("saveApproxHessian");  
    approxHessianName = config.get<std::string>("approxHessianName");  
    epsilonHessian = config.get<ValueType>("epsilonHessian");
    
    dimension = config.get<std::string>("dimension");
    equationType = config.get<std::string>("equationType");
    
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
            
    if(equationType!="SH"){
        wavefieldVX = wavefield.getRefVX();
        wavefieldVX *= wavefieldVX;
        wavefieldVX *= DT;
        approxHessian += wavefieldVX;}
        
    if(equationType!="SH"){
        wavefieldVY = wavefield.getRefVY();
        wavefieldVY *= wavefieldVY;
        wavefieldVY *= DT;
        approxHessian += wavefieldVY;}
    
    if(dimension=="3D" || equationType=="SH"){
        wavefieldVZ = wavefield.getRefVZ();
        wavefieldVZ *= wavefieldVZ;
        wavefieldVZ *= DT;
        approxHessian += wavefieldVZ;}
    
}

/*! \brief Apply the approximation of the diagonal of the inverse of the Hessian for one shot
 *
 *
 \param gradientPerShot 
 \param shotNumber 
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, scai::IndexType shotNumber)
{
//     sqrt(approxHessian) missing because of |u_i| (see old IFOS)??
    
    /* Stabilize Hessian for inversion (of diagonal matrix) and normalize Hessian */
    approxHessian += epsilonHessian*approxHessian.maxNorm(); 
    approxHessian *= 1 / approxHessian.maxNorm(); 
    
    if(saveApproxHessian==1){
        approxHessian.writeToFile(approxHessianName + ".shot_" + std::to_string(shotNumber+1) + ".mtx");}
    
    approxHessian = 1 / approxHessian;
    gradientPerShot *= approxHessian;    // overload operator /= in gradient-class    
    
}

template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::resetApproxHessian()
{
    approxHessian = 0;   // does not change size/distribution
}

template class KITGPI::Preconditioning::EnergyPreconditioning<double>;
template class KITGPI::Preconditioning::EnergyPreconditioning<float>;
