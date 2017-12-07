#pragma once

#include <vector>
#include <scai/lama.hpp>

// Does it need to be a template class?
template <typename ValueType>
class Misfit{
    
public:
    
    /* Default constructor and destructor */
    Misfit(){};
    ~Misfit(){};
    
    scai::lama::Scalar get(int element);
    void set(int element, scai::lama::Scalar value);
    void add(scai::lama::Scalar value);
       
    void l1();
    void l2();

private:

    std::vector<scai::lama::Scalar> misfit; 
    
};
