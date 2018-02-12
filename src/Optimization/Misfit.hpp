#pragma once

#include <scai/lama.hpp>
#include <vector>

/*! \brief The class Misfit holds the misfits for inversion
 *
 */
// Does it need to be a template class?
template <typename ValueType>
class Misfit
{

  public:
    /* Default constructor and destructor */
    Misfit(){};
    ~Misfit(){};

    scai::lama::Scalar getMisfitSum(int iteration);
    scai::lama::Scalar getMisfitShot(int iteration, int shotNumber);
    void add(scai::lama::DenseVector<ValueType> vector);

    void l1();
    void l2();

  private:
    std::vector<scai::lama::DenseVector<ValueType>> misfitShot;
};
