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
    approxHessian.setSameValue(dist, 0);

    wavefieldVX.setSameValue(dist, 0);
    wavefieldVY.setSameValue(dist, 0);
    wavefieldVZ.setSameValue(dist, 0);
    
    saveApproxHessian = config.get<bool>("saveApproxHessian");  
    approxHessianName = config.get<std::string>("approxHessianName");  
    epsilonHessian = config.get<ValueType>("epsilonHessian");
    
    dimension = config.get<std::string>("dimension");
    equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);
}

/*! \brief Calculate the approximation of the diagonal of the inverse of the Hessian for one shot
 *
 *
 \param wavefield Wavefield of one time step
 \param DT temporal sampling interval
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::intSquaredWavefields(KITGPI::Wavefields::Wavefields<ValueType> &wavefield, ValueType DT)
{
            
    if(equationType.compare("sh") != 0){
        wavefieldVX = wavefield.getRefVX();
        wavefieldVX *= wavefieldVX;
        wavefieldVX *= DT;
        approxHessian += wavefieldVX;        
    }
        
    if(equationType.compare("sh") != 0){
        wavefieldVY = wavefield.getRefVY();
        wavefieldVY *= wavefieldVY;
        wavefieldVY *= DT;
        approxHessian += wavefieldVY;        
    }
    
    if(dimension.compare("3d") == 0 || equationType.compare("sh") == 0){
        wavefieldVZ = wavefield.getRefVZ();
        wavefieldVZ *= wavefieldVZ;
        wavefieldVZ *= DT;
        approxHessian += wavefieldVZ;        
    }    
}

/*! \brief Apply the approximation of the diagonal of the inverse of the Hessian for one shot
 *
 *
 \param gradientPerShot 
 \param shotNumber 
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, scai::IndexType shotNumber, scai::IndexType fileFormat)
{
//     sqrt(approxHessian) missing because of |u_i| (see IFOS2D)?
    
    /* Stabilize Hessian for inversion (of diagonal matrix) and normalize Hessian */
    approxHessian += epsilonHessian*approxHessian.maxNorm(); 
    approxHessian *= 1 / approxHessian.maxNorm(); 
    
    if(saveApproxHessian==1){
        IO::writeVector(approxHessian, approxHessianName + ".shot_" + std::to_string(shotNumber), fileFormat);        
    }
        
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
