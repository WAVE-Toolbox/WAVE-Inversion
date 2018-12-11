#pragma once

#define MASTERGPI 0

extern bool verbose; // global variable

// all macros defined here are to support HOST_PRINT with an optional third parameter

/*!
 \def HOST_PRINT_STANDARD(comm, msg)
 Print message only on master. Message is printed independent of verbosity level.
 */
#define HOST_PRINT_STANDARD(comm, msg) \
    {                             \
    int myRank = comm->getRank(); \
    if (myRank == MASTERGPI)      \
        std::cout << msg;         \
    }                             \

/*!
 \def HOST_PRINT_VERBOSE(comm, msg, msgVerbose)
 Print message only on master. Verbosity level decides if msgVerbose is printed.
 */
#define HOST_PRINT_VERBOSE(comm, msg, msgVerbose)  \
    {                                       \
        int myRank = comm->getRank();       \
        if (myRank == MASTERGPI) {          \
            std::cout << msg;               \
            if (verbose) {                  \
                std::cout << msgVerbose; }  \
        }                                   \
    }                                       \

/*!
 \def GET_4TH_ARG(arg1, arg2, arg3, arg4, ...)
 Returns the fourth argument of the argument list.
 */
#define GET_4TH_ARG(arg1, arg2, arg3, arg4, ...) arg4

/*!
 \def HOST_PRINT_MACRO_CHOOSER(...)
 Chooses eigther HOST_PRINT_STANDARD or HOST_PRINT_VERBOSE depending on the number of arguments it receives.
 By calling GET_4TH_ARG with the list of the arguments HOST_PRINT_MACRO_CHOOSER receives and the two macros the 4th argument is always the correct one. 
 */
#define HOST_PRINT_MACRO_CHOOSER(...) \
    GET_4TH_ARG(__VA_ARGS__, HOST_PRINT_VERBOSE, \
                HOST_PRINT_STANDARD, )

/*!
 \def HOST_PRINT(comm,msg, (msgVerbose))
 Print message only on master. The third argument is optional. 
 If the global variable verbose is set to true, msgVerbose is printed in addition to msg.
 */
#define HOST_PRINT(...) HOST_PRINT_MACRO_CHOOSER(__VA_ARGS__)(__VA_ARGS__)