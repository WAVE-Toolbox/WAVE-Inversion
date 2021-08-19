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
    useEnergyPreconditioning = config.get<scai::IndexType>("useEnergyPreconditioning");  
    saveApproxHessian = config.get<bool>("saveApproxHessian");  
    approxHessianName = config.get<std::string>("approxHessianName");  
    epsilonHessian = config.get<ValueType>("epsilonHessian");
    
    dimension = config.get<std::string>("dimension");
    equationType = config.get<std::string>("equationType");
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);
    std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);
    isSeismic = Common::checkEquationType<ValueType>(equationType);
    
    approxHessian.setSameValue(dist, 0);

    if (isSeismic) {
        wavefieldVX.setSameValue(dist, 0);
        wavefieldVY.setSameValue(dist, 0);
        wavefieldVZ.setSameValue(dist, 0);
    } else {
        wavefieldEX.setSameValue(dist, 0);
        wavefieldEY.setSameValue(dist, 0);
        wavefieldEZ.setSameValue(dist, 0);       
    }
}

/*! \brief Calculate the approximation of the diagonal of the inverse of the Hessian for one shot
 *
 *
 \param wavefield Wavefield of one time step
 \param DT temporal sampling interval
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::intSquaredWavefields(KITGPI::Wavefields::Wavefields<ValueType> &wavefield, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldBack, ValueType DT)
{
    if (useEnergyPreconditioning != 0 && isSeismic) {               
        if(equationType.compare("sh") != 0){
            wavefieldVX = wavefield.getRefVX();
            wavefieldVX *= wavefieldVX;
            wavefieldVX *= DT;
            approxHessian += wavefieldVX;  
            if (useEnergyPreconditioning == 2) {
                wavefieldVX = wavefieldBack.getRefVX();
                wavefieldVX *= wavefield.getRefVX().maxNorm() / wavefieldVX.maxNorm();
                wavefieldVX *= wavefieldVX;
                wavefieldVX *= DT;
                approxHessian += wavefieldVX; 
            }      
        }
            
        if(equationType.compare("sh") != 0){
            wavefieldVY = wavefield.getRefVY();
            wavefieldVY *= wavefieldVY;
            wavefieldVY *= DT;
            approxHessian += wavefieldVY;   
            if (useEnergyPreconditioning == 2) {
                wavefieldVY = wavefieldBack.getRefVY();
                wavefieldVY *= wavefield.getRefVY().maxNorm() / wavefieldVY.maxNorm();
                wavefieldVY *= wavefieldVY;
                wavefieldVY *= DT;
                approxHessian += wavefieldVY;
            }             
        }
        
        if(dimension.compare("3d") == 0 || equationType.compare("sh") == 0){
            wavefieldVZ = wavefield.getRefVZ();
            wavefieldVZ *= wavefieldVZ;
            wavefieldVZ *= DT;
            approxHessian += wavefieldVZ;  
            if (useEnergyPreconditioning == 2) {
                wavefieldVZ = wavefieldBack.getRefVZ();
                wavefieldVZ *= wavefield.getRefVZ().maxNorm() / wavefieldVZ.maxNorm();
                wavefieldVZ *= wavefieldVZ;
                wavefieldVZ *= DT;
                approxHessian += wavefieldVZ; 
            }                   
        }    
    } else if (useEnergyPreconditioning != 0 && !isSeismic) {            
        if(equationType.compare("tmem") != 0 && equationType.compare("viscotmem") != 0){
            wavefieldEX = wavefield.getRefEX();
            wavefieldEX *= wavefieldEX;
            wavefieldEX *= DT;
            approxHessian += wavefieldEX;  
            if (useEnergyPreconditioning == 2) {
                wavefieldEX = wavefieldBack.getRefEX();
                wavefieldEX *= wavefield.getRefEX().maxNorm() / wavefieldEX.maxNorm();
                wavefieldEX *= wavefieldEX;
                wavefieldEX *= DT;
                approxHessian += wavefieldEX; 
            }
        }
            
        if(equationType.compare("tmem") != 0 && equationType.compare("viscotmem") != 0){
            wavefieldEY = wavefield.getRefEY();
            wavefieldEY *= wavefieldEY;
            wavefieldEY *= DT;
            approxHessian += wavefieldEY;
            if (useEnergyPreconditioning == 2) {
                wavefieldEY = wavefieldBack.getRefEY();
                wavefieldEY *= wavefield.getRefEY().maxNorm() / wavefieldEY.maxNorm();
                wavefieldEY *= wavefieldEY;
                wavefieldEY *= DT;
                approxHessian += wavefieldEY;
            }        
        }
        
        if(dimension.compare("3d") == 0 || equationType.compare("tmem") == 0 || equationType.compare("viscotmem") == 0){
            wavefieldEZ = wavefield.getRefEZ();
            wavefieldEZ *= wavefieldEZ;
            wavefieldEZ *= DT;
            approxHessian += wavefieldEZ;    
            if (useEnergyPreconditioning == 2) {
                wavefieldEZ = wavefieldBack.getRefEZ();
                wavefieldEZ *= wavefield.getRefEZ().maxNorm() / wavefieldEZ.maxNorm();
                wavefieldEZ *= wavefieldEZ;
                wavefieldEZ *= DT;
                approxHessian += wavefieldEZ; 
            }            
        } 
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
    if (useEnergyPreconditioning != 0) {  
    //     sqrt(approxHessian) missing because of |u_i| (see IFOS2D)?
        
        /* Stabilize Hessian for inversion (of diagonal matrix) and normalize Hessian */
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm(); 
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, approxHessianName + ".shot_" + std::to_string(shotNumber), fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        gradientPerShot *= approxHessian;    // overload operator /= in gradient-class    
    }
}

template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::resetApproxHessian()
{
    approxHessian = 0;   // does not change size/distribution
}

template class KITGPI::Preconditioning::EnergyPreconditioning<double>;
template class KITGPI::Preconditioning::EnergyPreconditioning<float>;
