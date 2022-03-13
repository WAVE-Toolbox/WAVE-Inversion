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
    if (useEnergyPreconditioning != 0) {
        saveApproxHessian = config.get<bool>("saveApproxHessian");  
        approxHessianName = config.get<std::string>("approxHessianName");  
        epsilonHessian = config.get<ValueType>("epsilonHessian");
        
        dimension = config.get<std::string>("dimension");
        equationType = config.get<std::string>("equationType");
        std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);
        std::transform(equationType.begin(), equationType.end(), equationType.begin(), ::tolower);
        isSeismic = Common::checkEquationType<ValueType>(equationType);
        
        approxHessian.setSameValue(dist, 0);
        if (useEnergyPreconditioning == 2 || useEnergyPreconditioning == 4)
            approxHessianAdjoint.setSameValue(dist, 0);

        wavefieldX.setSameValue(dist, 0);
        wavefieldY.setSameValue(dist, 0);
        wavefieldZ.setSameValue(dist, 0);
    }
}

/*! \brief Calculate the approximation of the diagonal of the inverse of the Hessian for one shot
 *
 *
 \param wavefield Wavefield of one time step
 \param DT temporal sampling interval
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::intSquaredWavefields(KITGPI::Wavefields::Wavefields<ValueType> &wavefield, ValueType DT, bool isAdjoint)
{
    if (isSeismic && ((!isAdjoint && (useEnergyPreconditioning == 1 || useEnergyPreconditioning == 2 || useEnergyPreconditioning == 4)) || (isAdjoint && (useEnergyPreconditioning == 2 || useEnergyPreconditioning == 4)))) {               
        if(equationType.compare("sh") != 0 && equationType.compare("viscosh") != 0){
            wavefieldX = wavefield.getRefVX();
            wavefieldX *= wavefieldX;
            wavefieldX *= DT;
            if (!isAdjoint) {
                approxHessian += wavefieldX;  
            } else {
                approxHessianAdjoint += wavefieldX; 
            }
        }
            
        if(equationType.compare("sh") != 0 && equationType.compare("viscosh") != 0){
            wavefieldY = wavefield.getRefVY();
            wavefieldY *= wavefieldY;
            wavefieldY *= DT;
            if (!isAdjoint) {
                approxHessian += wavefieldY;  
            } else {
                approxHessianAdjoint += wavefieldY; 
            }      
        }
        
        if(dimension.compare("3d") == 0 || equationType.compare("sh") == 0 || equationType.compare("viscosh") == 0){
            wavefieldZ = wavefield.getRefVZ();
            wavefieldZ *= wavefieldZ;
            wavefieldZ *= DT;
            if (!isAdjoint) {
                approxHessian += wavefieldZ;  
            } else {
                approxHessianAdjoint += wavefieldZ; 
            }
        }    
    } else if (!isSeismic && ((!isAdjoint && (useEnergyPreconditioning == 1 || useEnergyPreconditioning == 2 || useEnergyPreconditioning == 4)) || (isAdjoint && (useEnergyPreconditioning == 2 || useEnergyPreconditioning == 4)))) {   
        if(equationType.compare("tmem") != 0 && equationType.compare("viscotmem") != 0){
            wavefieldX = wavefield.getRefEX();
            wavefieldX *= wavefieldX;
            wavefieldX *= DT;
            if (!isAdjoint) {
                approxHessian += wavefieldX;  
            } else {
                approxHessianAdjoint += wavefieldX; 
            }
        }
            
        if(equationType.compare("tmem") != 0 && equationType.compare("viscotmem") != 0){
            wavefieldY = wavefield.getRefEY();
            wavefieldY *= wavefieldY;
            wavefieldY *= DT;
            if (!isAdjoint) {
                approxHessian += wavefieldY;  
            } else {
                approxHessianAdjoint += wavefieldY; 
            }
        }
        
        if(dimension.compare("3d") == 0 || equationType.compare("tmem") == 0 || equationType.compare("viscotmem") == 0){
            wavefieldZ = wavefield.getRefEZ();
            wavefieldZ *= wavefieldZ;
            wavefieldZ *= DT;
            if (!isAdjoint) {
                approxHessian += wavefieldZ;  
            } else {
                approxHessianAdjoint += wavefieldZ; 
            }          
        } 
    }
}

/*! \brief Apply the approximation of the diagonal of the inverse of the Hessian for one shot
 *
 \param gradientPerShot 
 \param shotNumber 
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::apply(KITGPI::Gradient::Gradient<ValueType> &gradientPerShot, scai::IndexType shotNumber, scai::IndexType fileFormat)
{
    if (useEnergyPreconditioning != 0 && useEnergyPreconditioning != 3) {  
    //     sqrt(approxHessian) missing because of |u_i| (see IFOS2D)?
    
        /* Stabilize Hessian for inversion (of diagonal matrix) and normalize Hessian */
        if (useEnergyPreconditioning == 2 || useEnergyPreconditioning == 4) {
            approxHessian *= approxHessianAdjoint;
            approxHessian = scai::lama::sqrt(approxHessian);
        }
        approxHessian += epsilonHessian*approxHessian.maxNorm(); 
        approxHessian *= 1 / approxHessian.maxNorm();   
        
        if(saveApproxHessian){
            IO::writeVector(approxHessian, approxHessianName + ".shot_" + std::to_string(shotNumber), fileFormat);        
        }
            
        approxHessian = 1 / approxHessian;
        gradientPerShot *= approxHessian;    // overload operator /= in gradient-class 
        
        if (useEnergyPreconditioning == 4) {
            gradientPerShot.applyEnergyPreconditioning(epsilonHessian, saveApproxHessian, approxHessianName + ".shot_" + std::to_string(shotNumber), fileFormat);
        }
    } else if (useEnergyPreconditioning == 3) {
        gradientPerShot.applyEnergyPreconditioning(epsilonHessian, saveApproxHessian, approxHessianName + ".shot_" + std::to_string(shotNumber), fileFormat);
    }
}

/*! \brief Reset approxHessian
 *
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::resetApproxHessian()
{
    if (useEnergyPreconditioning != 0) {
        approxHessian.setSameValue(wavefieldX.getDistributionPtr(), 0); // does not change size/distribution
        if (useEnergyPreconditioning == 2 || useEnergyPreconditioning == 4)
            approxHessianAdjoint.setSameValue(wavefieldX.getDistributionPtr(), 0);
    }
}

/*! \brief Transform approxHessian
 *
 */
template <typename ValueType>
void KITGPI::Preconditioning::EnergyPreconditioning<ValueType>::applyTransform(scai::lama::Matrix<ValueType> const &lhs)
{
    if (useEnergyPreconditioning != 0) {
        approxHessian = lhs * approxHessian;   // change size/distribution
        if (useEnergyPreconditioning == 2 || useEnergyPreconditioning == 4)
            approxHessianAdjoint = lhs * approxHessianAdjoint; 
    }
}

template class KITGPI::Preconditioning::EnergyPreconditioning<double>;
template class KITGPI::Preconditioning::EnergyPreconditioning<float>;
