
// #include <iosfwd>

#include <scai/lama.hpp>
#include <scai/common/Walltime.hpp>

#include <Configuration/Configuration.hpp>
#include <ForwardSolver/ForwardSolver.hpp>
#include <ForwardSolver/Derivatives/DerivativesFactory.hpp>
#include <Acquisition/Receivers.hpp>
#include <Acquisition/Sources.hpp>
#include <Modelparameter/ModelparameterFactory.hpp>
#include <Wavefields/WavefieldsFactory.hpp>


template <typename ValueType>
class GradientCalculation{
    
public:
    
    /* Default constructor and destructor */
    GradientCalculation(){};
    ~GradientCalculation(){};

    void allocate(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr no_dist_NT, scai::dmemo::CommunicatorPtr comm, scai::hmemo::ContextPtr ctx);
    /* Calculate gradients */
    void calc(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, IndexType iteration);
    
    scai::lama::DenseVector<ValueType> grad_bulk;  // make it private and write getter and setter functions
    scai::lama::DenseVector<ValueType> grad_rho;
    scai::lama::DenseVector<ValueType> grad_vp;
    
private:
    
    scai::lama::DenseVector<ValueType> waveconv_v;
    scai::lama::DenseVector<ValueType> waveconv_p;
    scai::lama::DenseVector<ValueType> tmp;
    scai::lama::DenseVector<ValueType> waveconv_v_sum;
    scai::lama::DenseVector<ValueType> waveconv_p_sum;
    
    /* WavefieldRecords (does it make sense here? -> is it used somewhere else or only in calc?, if yes declare a pointer instead of class and call constructor in function?) */ 
    scai::lama::DenseMatrix<ValueType> wavefieldrecordvx;
    scai::lama::DenseMatrix<ValueType> wavefieldrecordvy;
    scai::lama::DenseMatrix<ValueType> wavefieldrecordp;
    
};
    

