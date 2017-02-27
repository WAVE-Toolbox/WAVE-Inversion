# Contributing Guidelines

**All contributions to this project are warmly welcome!**

In the following some guidelines are given which ensure a smooth and hassle-free cooperation.

## Write and adjust the source code

The source code itself is written in object-oriented C++11. This means, all features of the C++11 standard can be used and the used compilers have to be compatible to the C++11 standard. However, to increase compatibility, features of newer standards, such as C++14 or C++17, should not be used.

To ensure that the code will stay easy to maintain and extend, all developers are requested to stick to the five **SOLID principles**. More information on this principles can be found eg. on [Wikipedia](https://en.wikipedia.org/wiki/SOLID_(object-oriented_design).

#### Code formatting

For consistent source code formatting a clang-format style file is provided.
Please config your integrated development environment (IDE) to format source code according to the style file given by the `.clang-format` file located within the `src/` folder.
All major IDEs (e.g. Kdevelop, XCode, emacs, vim) support this file type.
Alternatively, a manual source code formatting is possible by the shell script, which is located in `src/Scripts/formatSourceCode.sh`.

You may be required to install the clang-format tool on your machine.
Since clang-format is part of llvm project, more information are given on [llvm.org](llvm.org) or [https://clang.llvm.org/docs/ClangFormat.html](https://clang.llvm.org/docs/ClangFormat.html).

#### Code documentation

The `doxygen` tool is used for automatic generation of source code documentation.
Thus, it is required that the whole code is annotated in `doxygen` syntax.
A detailed and complete source code documentation is of particular importance for new developers, since it allows to access and understand the interfaces and classes quickly.


More information on the `doxygen` syntax can be found on [www.doxygen.org](www.doxygen.org).

#### Code testing

Before you commit your changes make sure all tests pass without errors.
Two kinds of tests are provided:
1. Unit tests (based on the [Google Test framework](https://github.com/google/googletest))
2. Full integration tests

In order to use the Google Test framework, you have to set the environment variable `GTEST_DIR` to the location of the compiled Google Test library (`libgtest.*` and `libgtest_main.*`). The integration test has no special requirements apart LAMA and MPI. 

The unit tests will verify the functionality of individual functions and classes, while the integration tests will run the whole forward code and verify that the obtained result is identical to a reference solution.
