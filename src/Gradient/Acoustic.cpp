#include "Acoustic.hpp"
#include <scai/lama/io/FileIO.hpp>
using namespace scai;
using namespace KITGPI;

/*! \brief Constructor that is using the Configuration class
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType>::Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "acoustic";
    init(ctx, dist, 0.0, 0.0);
}

/*! \brief Initialisation that is using the configuration class
 *
 *  Generates a homogeneous gradient, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist, 0.0, 0.0);
}


/*! \brief Initialisation that is generating a homogeneous gradient
 *
 *  Generates a homogeneous gradient, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param velocityP_const velocity gradients given as scalar value
 \param rho_const Density gradient given as scalar value
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType rho_const)
{
    this->initParameterisation(velocityP, ctx, dist, velocityP_const);
    this->initParameterisation(density, ctx, dist, rho_const);
}


//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType>::Acoustic(const Acoustic &rhs)
{
    equationType = rhs.equationType;
    velocityP = rhs.velocityP;
    density = rhs.density;
}

/*! \brief Write gradient to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::write(std::string filename, IndexType partitionedOut, KITGPI::Workflow::Workflow<ValueType> const &workflow) const
{
    if (workflow.getInvertForVp() == 1) {
        std::string filenameP = filename + ".vp.mtx";
        this->writeParameterisation(velocityP, filenameP, partitionedOut);
    }

    if (workflow.getInvertForDensity() == 1) {
        std::string filenamedensity = filename + ".density.mtx";
        this->writeParameterisation(density, filenamedensity, partitionedOut);
    }
    
};

/*! \brief Get equationType (acoustic)
 */
template <typename ValueType>
std::string KITGPI::Gradient::Acoustic<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Acoustic<ValueType>::getVelocityS()
{
    COMMON_THROWEXCEPTION("The S-wave velocity is not defined in an acoustic simulation.")
    return (velocityS);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Acoustic<ValueType>::getTauP()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Acoustic<ValueType>::getTauS()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::Acoustic<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::Acoustic<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Overloading * Operation
 *
 \param rhs scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> KITGPI::Gradient::Acoustic<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Gradient::Acoustic<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief non-member function to multiply (scalar as left operand)
 *
 \param lhs scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> operator*(ValueType lhs, KITGPI::Gradient::Acoustic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> &KITGPI::Gradient::Acoustic<ValueType>::operator*=(ValueType const &rhs)
{
    density *= rhs;
    velocityP *= rhs;

    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> KITGPI::Gradient::Acoustic<ValueType>::operator+(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    KITGPI::Gradient::Acoustic<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> &KITGPI::Gradient::Acoustic<ValueType>::operator+=(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    density += rhs.density;
    velocityP += rhs.velocityP;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Gradient which is subtractet.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> KITGPI::Gradient::Acoustic<ValueType>::operator-(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    KITGPI::Gradient::Acoustic<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> &KITGPI::Gradient::Acoustic<ValueType>::operator-=(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    density -= rhs.density;
    velocityP -= rhs.velocityP;
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::Acoustic<ValueType> &KITGPI::Gradient::Acoustic<ValueType>::operator=(KITGPI::Gradient::Acoustic<ValueType> const &rhs)
{
    // why does rhs.density not work (density = protected)
    velocityP = rhs.velocityP;
    density = rhs.density;
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::assign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    density = rhs.getDensity();
    velocityP = rhs.getVelocityP();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{

    density -= rhs.getDensity();
    velocityP -= rhs.getVelocityP();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{

    density += rhs.getDensity();
    velocityP += rhs.getVelocityP();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::timesAssign(ValueType const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::timesAssign(scai::lama::Vector<ValueType> const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract gradient.
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    scai::lama::DenseVector<ValueType> temp;
    temp = lhs.getVelocityP() - rhs.getVelocityP();
    lhs.setVelocityP(temp);
    temp = lhs.getDensity() - rhs.getDensity();
    lhs.setDensity(temp);
};

/*! \brief function for scaling the gradients with the model parameter 
 *
 \param model Abstract model.
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForVp()) {
        velocityP *= 1 / velocityP.maxNorm() * model.getVelocityP().maxNorm();
    }
    if (workflow.getInvertForDensity()) {
        density *= 1 / density.maxNorm() * model.getDensity().maxNorm();
    }
}

/*! \brief Function for normalizing the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::normalize()
{
    if (this->getNormalizeGradient()) {
        ValueType gradientMax = velocityP.maxNorm();
        if (gradientMax != 0)
            velocityP *= 1 / gradientMax;
        gradientMax = density.maxNorm();
        if (gradientMax != 0)
            density *= 1 / gradientMax;
    }
}

/*! \brief function for calculating the acoustic gradients from the cross-correlation and the model parameter 
 *
 \param model Abstract model.
 \param correlatedWavefields Abstract xCorr.
 \param DT Temporal discretization 
 \param workflow 
 *
 \f{eqnarray*}
   \nabla_K E &=& -\mathrm{d}t \frac{1}{v_{\mathrm{p}}^4 \rho^2} \cdot X_{\lambda} \\ 
   \nabla_{v_{\mathrm{p}}} E &=& 2 v_{\mathrm{p}} \rho \cdot \nabla_K E \\
   \nabla_{\rho} E &=& v_{\mathrm{p}}^2 \cdot \nabla_K E - \mathrm{d}t \cdot X_{\rho}
 \f}
 *
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr2Dacoustic<ValueType>::update} for the cross-correlations \f$ (X_{\lambda},X_{\rho}) \f$ in 2D 
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr3Dacoustic<ValueType>::update} for the cross-correlations \f$ (X_{\lambda},X_{\rho}) \f$ in 3D
 *
 */
template <typename ValueType>
void KITGPI::Gradient::Acoustic<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //dt should be in cross correlation!
    //gradBulk = -dt * Padj * dPfw/dt / lambda^2
    scai::lama::DenseVector<ValueType> gradBulk;
    gradBulk = scai::lama::pow(model.getVelocityP(), 2);
    gradBulk *= model.getDensity();
    gradBulk = scai::lama::pow(gradBulk, -2);
    Common::replaceInvalid<ValueType>(gradBulk,0.0);

    gradBulk *= correlatedWavefields.getXcorrLambda();

    gradBulk *= -DT;

    scai::hmemo::ContextPtr ctx = gradBulk.getContextPtr();
    scai::dmemo::DistributionPtr dist = gradBulk.getDistributionPtr();

    if (workflow.getInvertForVp()) {
        //grad_vp = 2 * rho *vp * gradBulk
        velocityP = 2 * gradBulk;
        velocityP *= model.getDensity();
        velocityP *= model.getVelocityP();
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }

    if (workflow.getInvertForDensity()) {
	//grad_densityRaw
        //grad_density' = vp^2 * gradBulk + (-) (sum(i) (dt Vadj,i * dVfw,i/dt))
        density = scai::lama::pow(model.getVelocityP(), 2);
        density *= gradBulk;
        density -= DT * correlatedWavefields.getXcorrRho();
    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }
}

template class KITGPI::Gradient::Acoustic<double>;
template class KITGPI::Gradient::Acoustic<float>;
