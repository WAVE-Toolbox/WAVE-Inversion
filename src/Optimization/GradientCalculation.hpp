
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

    /* Calculate gradients */
    void calc(KITGPI::ForwardSolver::ForwardSolver<ValueType> &solver, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, KITGPI::Acquisition::Receivers<ValueType> &receivers, KITGPI::Acquisition::Sources<ValueType> &sources, KITGPI::Modelparameter::Modelparameter<ValueType> const &model, KITGPI::Wavefields::Wavefields<ValueType> &wavefields, KITGPI::Configuration::Configuration config, IndexType iteration);
    
    scai::lama::DenseVector<ValueType> grad_bulk;
    scai::lama::DenseVector<ValueType> grad_rho;
    scai::lama::DenseVector<ValueType> grad_vp;
    
// private:
    
    /* WavefieldRecords (does it make sense here? -> is it used somewhere else or only in calc?, if yes declare a pointer instead of class and call constructor in function?) */ 
//     scai::lama::DenseMatrix<ValueType> wavefieldrecordvx;
//     scai::lama::DenseMatrix<ValueType> wavefieldrecordvy;
//     scai::lama::DenseMatrix<ValueType> wavefieldrecordp;
    
    /* Adjoint sources  (does it make sense here? -> is it used somewhere else or only in calc?) */ 
//     KITGPI::Acquisition::Receivers<ValueType> adjoint;
//     KITGPI::Acquisition::Seismogram<ValueType> truedata;
//     KITGPI::Acquisition::Seismogram<ValueType> synthetic;
    
};
    

