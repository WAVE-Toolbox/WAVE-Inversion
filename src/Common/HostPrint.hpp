#pragma once

#define MASTERGPI 0

/*!
 \def HOST_PRINT(comm,msg)
 Print message only on master.
 */
#define HOST_PRINT(comm, msg)         \
    {                                 \
        int myRank = comm->getRank(); \
        if (myRank == MASTERGPI) {    \
            std::cout << msg;         \
        }                             \
    }
