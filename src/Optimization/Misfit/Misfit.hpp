#pragma once

#include <scai/lama.hpp>
#include <vector>

namespace KITGPI
{
    
    //! \brief Misfit namespace
    namespace Misfit
    {
        //! \brief Abstract class for the misfit
        /*!
         * As this class is an abstract class, all constructors are protected.
         */
        
        template <typename ValueType>
        class Misfit
        {

        public:

            virtual void calc() = 0;
            scai::lama::Scalar getMisfitSum(int iteration);
            scai::lama::Scalar getMisfitShot(int iteration, int shotNumber);
            void add(scai::lama::DenseVector<ValueType> vector);

            //! \brief Misfit pointer
            typedef std::shared_ptr<Misfit<ValueType>> MisfitPtr;
            
        protected:
            
            /* Default constructor and destructor */
            Misfit(){};
            ~Misfit(){};
            
            std::vector<scai::lama::DenseVector<ValueType>> misfitShot;
            
        };
    }
}
