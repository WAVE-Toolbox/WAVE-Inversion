#include "Elastic.hpp"
#include <scai/dmemo/SingleDistribution.hpp>

using namespace scai;

/*! \brief Constructor that is using the Configuration class
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType>::Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "elastic";
    init(ctx, dist, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation with zeros
 *
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist, 0.0, 0.0, 0.0);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho Density given as Scalar
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho)
{
    this->initParameterisation(velocityP, ctx, dist, velocityP_const);
    this->initParameterisation(velocityS, ctx, dist, velocityS_const);
    this->initParameterisation(density, ctx, dist, rho);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType>::Elastic(const Elastic &rhs)
{
    equationType = rhs.equationType;
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
}

/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 \param fileFormat format of output file
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::write(std::string filename, IndexType fileFormat, KITGPI::Workflow::Workflow<ValueType> const &workflow) const
{
    if (workflow.getInvertForVp() == 1) {
        std::string filenameP = filename + ".vp";
        this->writeParameterisation(velocityP, filenameP, fileFormat);
    }

    if (workflow.getInvertForVs() == 1) {
        std::string filenameS = filename + ".vs";
        this->writeParameterisation(velocityS, filenameS, fileFormat);
    }

    if (workflow.getInvertForDensity() == 1) {
        std::string filenamedensity = filename + ".density";
        this->writeParameterisation(density, filenamedensity, fileFormat);
    }
};

/*! \brief Get equationType (elastic)
 */
template <typename ValueType>
std::string KITGPI::Gradient::Elastic<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Elastic<ValueType>::getTauP()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Gradient::Elastic<ValueType>::getTauS()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Gradient::Elastic<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Gradient::Elastic<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> KITGPI::Gradient::Elastic<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Gradient::Elastic<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief Free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> operator*(ValueType lhs, KITGPI::Gradient::Elastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> &KITGPI::Gradient::Elastic<ValueType>::operator*=(ValueType const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
    velocityS *= rhs;

    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> KITGPI::Gradient::Elastic<ValueType>::operator+(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    KITGPI::Gradient::Elastic<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Gradient which is added.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> &KITGPI::Gradient::Elastic<ValueType>::operator+=(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    density += rhs.density;
    velocityP += rhs.velocityP;
    velocityS += rhs.velocityS;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> KITGPI::Gradient::Elastic<ValueType>::operator-(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    KITGPI::Gradient::Elastic<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Gradient which is subtracted.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> &KITGPI::Gradient::Elastic<ValueType>::operator-=(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    density -= rhs.density;
    velocityP -= rhs.velocityP;
    velocityS -= rhs.velocityS;

    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Gradient which is copied.
 */
template <typename ValueType>
KITGPI::Gradient::Elastic<ValueType> &KITGPI::Gradient::Elastic<ValueType>::operator=(KITGPI::Gradient::Elastic<ValueType> const &rhs)
{
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;

    return *this;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::assign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    density = rhs.getDensity();
    velocityP = rhs.getVelocityP();
    velocityS = rhs.getVelocityS();
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::minusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{

    density -= rhs.getDensity();
    velocityP -= rhs.getVelocityP();
    velocityS -= rhs.getVelocityS();
}

/*! \brief Function for overloading += Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtractet.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::plusAssign(KITGPI::Gradient::Gradient<ValueType> const &rhs)
{

    density += rhs.getDensity();
    velocityP += rhs.getVelocityP();
    velocityS += rhs.getVelocityS();
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::timesAssign(ValueType const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
    velocityS *= rhs;
}

/*! \brief Function for overloading *= Operation (called in base class)
 *
 \param rhs Abstract gradient which is subtracted.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::timesAssign(scai::lama::Vector<ValueType> const &rhs)
{
    density *= rhs;
    velocityP *= rhs;
    velocityS *= rhs;
}

/*! \brief Function for overloading -= Operation (called in base class)
 *
 \param lhs Abstract model.
 \param rhs Abstract gradient which is assigned.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> &lhs, KITGPI::Gradient::Gradient<ValueType> const &rhs)
{
    scai::lama::DenseVector<ValueType> temp;
    temp = lhs.getVelocityP() - rhs.getVelocityP();
    lhs.setVelocityP(temp);
    temp = lhs.getVelocityS() - rhs.getVelocityS();
    lhs.setVelocityS(temp);
    temp = lhs.getDensity() - rhs.getDensity();
    lhs.setDensity(temp);
};

/*! \brief Function for summing the gradients of all shot domains
 *
 \param commInterShot inter shot communication pointer
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::sumShotDomain(scai::dmemo::CommunicatorPtr commInterShot)
{
    /*reduction between shot domains.
    each shot domain may have a different distribution of (gradient) vectors.
      This happens if geographer is used (different result for dist on each shot domain even for homogenous architecture)
      or on heterogenous architecture. In this case even the number of processes on each domain can vary.
    Therfore it is necessary that only one process per shot domain communicates all data.
    */
    
    //get information from distributed vector
    auto size = velocityP.size();
    auto dist = velocityP.getDistributionPtr();
    auto comm = dist->getCommunicatorPtr();
    
    // create single distribution, only master process owns the complete vector (no distribution).
    
    int shotMaster=0;
    auto singleDist = std::make_shared<dmemo::SingleDistribution>( size, comm, shotMaster );
    
    //redistribute vector to master process
    // (this may cause memory issues for big models)
    velocityP.redistribute(singleDist);
    
    //reduce local array (size of local array is !=0 only for master process)
    commInterShot->sumArray(velocityP.getLocalValues());
    
    //redistribute vector to former partition
    velocityP.redistribute(dist);
    
    velocityS.redistribute(singleDist);
    commInterShot->sumArray(velocityS.getLocalValues());
    velocityS.redistribute(dist);
    
    density.redistribute(singleDist);
    commInterShot->sumArray(density.getLocalValues());
    density.redistribute(dist);
}

/*! \brief Smoothen gradient by gaussian window
 \param gradient gradient model
 \param modelCoordinates coordinate class object of the subset
 \param NX NX in model
 \param NY NY in model
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::smoothGradient(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType NX, scai::IndexType NY)
{
    auto savedVelocityP = velocityP;
    auto savedVelocityS = velocityS;
    auto savedDensity = density;
    
    for (IndexType y = 0; y < NY; y++) {
        for (IndexType x = 0; x < NX; x++) {
            if (x == 0) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x+3, y, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x+2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x+1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x+1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x+2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x+3, y, 0)]*0.0055;
            }
            else if (x == 1) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x+1, y, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x-1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x+1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x+2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x+3, y, 0)]*0.0055;
            }
            else if (x == 2) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x-1, y, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x-2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x-1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x+1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x+2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x+3, y, 0)]*0.0055;
            }
            else if (x == NX-3) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x-3, y, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x-2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x-1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x+1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x+2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x+1, y, 0)]*0.0055;
            }
            else if (x == NX-2) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x-3, y, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x-2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x-1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x+1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x-1, y, 0)]*0.0055;
                }
            else if (x == NX-1) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x-3, y, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x-2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x-1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x-1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x-2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x-3, y, 0)]*0.0055;
            }
            else {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x-3, y, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x-2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x-1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x+1, y, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x+2, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x+3, y, 0)]*0.0055;
            }
        }
    }
    for (IndexType x = 0; x < NX; x++) {
        for (IndexType y = 0; y < NY; y++) {
            if (y == 0) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x, y+3, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x, y+2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y+1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x, y+1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y+2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y+3, 0)]*0.0055;
            }
            else if (y == 1) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x, y+1, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y-1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x, y+1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y+2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y+3, 0)]*0.0055;
            }
            else if (y == 2) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x, y-1, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x, y-2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y-1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x, y+1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y+2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y+3, 0)]*0.0055;
            }
            else if (y == NY-3) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x, y-3, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x, y-2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y-1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x, y+1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y+2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y+1, 0)]*0.0055;
            }
            else if (y == NY-2) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x, y-3, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x, y-2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y-1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x, y+1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y-1, 0)]*0.0055;
                }
            else if (y == NY-1) {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x, y-3, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x, y-2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y-1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x, y-1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y-2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y-3, 0)]*0.0055;
            }
            else {
                velocityS[modelCoordinates.coordinate2index(x, y, 0)] = savedVelocityS[modelCoordinates.coordinate2index(x, y-3, 0)]*0.0055 + savedVelocityS[modelCoordinates.coordinate2index(x, y-2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y-1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y, 0)]*0.383 + savedVelocityS[modelCoordinates.coordinate2index(x, y+1, 0)]*0.242 + savedVelocityS[modelCoordinates.coordinate2index(x, y+2, 0)]*0.061 + savedVelocityS[modelCoordinates.coordinate2index(x, y+3, 0)]*0.0055;
            }
        }
    }
}

/*! \brief Function for scaling the gradients with the model parameter
 *
 \param model Abstract model.
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::scale(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    if (workflow.getInvertForVp()) {
        velocityP *= 1 / velocityP.maxNorm() * model.getVelocityP().maxNorm();
    }

    if (workflow.getInvertForVs()) {
        velocityS *= 1 / velocityS.maxNorm() * model.getVelocityS().maxNorm();
    }

    if (workflow.getInvertForDensity()) {
        density *= 1 / density.maxNorm() * model.getDensity().maxNorm();
    }
}

/*! \brief Function for normalizing the gradient
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::normalize()
{
    if (this->getNormalizeGradient()) {
        ValueType gradientMax = velocityP.maxNorm();
        if (gradientMax != 0)
            velocityP *= 1 / gradientMax;
        gradientMax = velocityS.maxNorm();
        if (gradientMax != 0)
            velocityS *= 1 / gradientMax;
        gradientMax = density.maxNorm();
        if (gradientMax != 0)
            density *= 1 / gradientMax;
    }
}

/*! \brief Function for calculating the elastic gradients from the cross correlation and the model parameter
 *
 \param model Abstract model.
 \param correlatedWavefields Abstract xCorr.
 \param DT Temporal discretization
 \param workflow
 *
 \f{eqnarray*}
   \nabla_{\lambda} E &=& - \mathrm{d}t \frac{1}{(N \lambda+2\mu)^2} \cdot X_{\lambda} \\
   \nabla_{\mu} E &=& \mathrm{d}t \left[ \frac{N \lambda^2 + 4 \mu \lambda}{2 \mu^2 (N \lambda+2\mu)^2} - \frac{1}{2 \mu^2} \right] \cdot X_{\mu,A} + \mathrm{d}t \frac{N \lambda^2 + 4 \mu \lambda}{2 \mu^2 (N \lambda+2\mu)^2} \cdot X_{\mu,B} - \mathrm{d}t  \frac{1}{\mu^2} \cdot X_{\mu,C} \\ \\
   \nabla_{v_{\mathrm{p}}} E &=& 2 \rho v_{\mathrm{p}} \cdot \nabla_{\lambda} E \\
   \nabla_{v_{\mathrm{s}}} E &=& 2 \rho v_{\mathrm{s}} \cdot \nabla_{\mu} E - 4 \rho v_{\mathrm{s}} \nabla_{\lambda} E \\
   \nabla_{\rho} E &=& ( v_{\mathrm{p}}^2-2v_{\mathrm{s}}^2 ) \cdot \nabla_{\lambda} E + v_{\mathrm{s}}^2 \cdot \nabla_{\mu} E - \mathrm{d}t \cdot X_{\rho}
 \f}
 *
 * with \f$ N \f$ as the number of dimensions, \f$ \mu = \rho v_{\mathrm{s}}^2 \f$ and \f$ \lambda = \rho v_{\mathrm{p}}^2 - 2\mu\f$.
 *
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr2Delastic<ValueType>::update} for the cross-correlations \f$ (X_{\lambda},X_{\rho},X_{\mu,A},X_{\mu,B},X_{\mu,C}) \f$ in 2D
 \sa{KITGPI::ZeroLagXcorr::ZeroLagXcorr3Delastic<ValueType>::update} for the cross-correlations \f$ (X_{\lambda},X_{\rho},X_{\mu,A},X_{\mu,B},X_{\mu,C}) \f$ in 3D
 */
template <typename ValueType>
void KITGPI::Gradient::Elastic<ValueType>::estimateParameter(KITGPI::ZeroLagXcorr::ZeroLagXcorr<ValueType> const &correlatedWavefields, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, ValueType DT, KITGPI::Workflow::Workflow<ValueType> const &workflow)
{
    //dt should be in cross correlation!

    scai::lama::DenseVector<ValueType> gradLambda;
    scai::lama::DenseVector<ValueType> gradMu;
    scai::lama::DenseVector<ValueType> temp;
    scai::lama::DenseVector<ValueType> temp2;
    scai::lama::DenseVector<ValueType> lambda;
    scai::lama::DenseVector<ValueType> mu;

    //dimension
    int N = correlatedWavefields.getNumDimension();

    mu = scai::lama::pow(model.getVelocityS(), 2);
    mu *= model.getDensity();

    lambda = scai::lama::pow(model.getVelocityP(), 2);
    lambda *= model.getDensity();
    lambda -= 2 * mu;

    // Lambda and Mu gradients
    //-1/(N*lambda+2mu)^2)
    gradLambda = N * lambda;
    gradLambda += 2 * mu;
    gradLambda = scai::lama::pow(gradLambda, -2);

    gradLambda *= correlatedWavefields.getXcorrLambda();
    gradLambda *= -DT;
    Common::replaceInvalid<ValueType>(gradLambda, 0.0);

    scai::hmemo::ContextPtr ctx = gradLambda.getContextPtr();
    scai::dmemo::DistributionPtr dist = gradLambda.getDistributionPtr();

    if ((workflow.getInvertForVs()) || (workflow.getInvertForDensity())) {
        //(N*lambda^2+4mu*lambda)/(2mu^2(N*lambda+2mu)^2)

        //temp2=>B
        temp2 = scai::lama::pow(lambda, 2);
        temp2 *= N;
        temp = 4 * mu;
        temp *= lambda;
        temp2 += temp;

        temp = scai::lama::pow(mu, 2);
        temp *= 2;
        temp2 /= temp;

        temp = N * lambda;
        temp += 2 * mu;
        temp = scai::lama::pow(temp, 2);
        temp2 /= temp;

        gradMu = temp2;
        gradMu *= correlatedWavefields.getXcorrMuB();

        //temp2=>A
        temp = scai::lama::pow(mu, -2);
        temp /= -2;
        temp2 += temp;

        temp2 *= correlatedWavefields.getXcorrMuA();
        gradMu += temp2;

        //temp2=>C
        temp2 = mu;
        temp2 *= mu;
        temp2 = 1 / temp2;
        temp2 *= correlatedWavefields.getXcorrMuC();
        gradMu -= temp2;

        gradMu *= DT;
        Common::replaceInvalid<ValueType>(gradMu, 0.0);
    }

    // vp, vs , rho gradients

    if (workflow.getInvertForVp()) {
        //grad_vp = 2*rho*vp*grad_lambda
        velocityP = 2 * gradLambda;
        velocityP *= model.getDensity();
        velocityP *= model.getVelocityP();
    } else {
        this->initParameterisation(velocityP, ctx, dist, 0.0);
    }

    if (workflow.getInvertForVs()) {

        velocityS = -4 * gradLambda;
        velocityS *= model.getDensity();
        velocityS *= model.getVelocityS();

        temp = 2 * gradMu;
        temp *= model.getDensity();
        temp *= model.getVelocityS();

        velocityS += temp;
    } else {
        this->initParameterisation(velocityS, ctx, dist, 0.0);
    }

    if (workflow.getInvertForDensity()) {

        density = scai::lama::pow(model.getVelocityP(), 2);
        temp = scai::lama::pow(model.getVelocityS(), 2);
        temp *= 2;
        density -= temp;
        density *= gradLambda;

        temp = scai::lama::pow(model.getVelocityS(), 2);
        temp *= gradMu;
        density += temp;

        temp = DT * correlatedWavefields.getXcorrRho();
        density -= temp;

    } else {
        this->initParameterisation(density, ctx, dist, 0.0);
    }
}

template class KITGPI::Gradient::Elastic<float>;
template class KITGPI::Gradient::Elastic<double>;
