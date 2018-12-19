#!/bin/bash

# Change to src/ folder
cd ../

# Format *.cpp files
find $(pwd) -name '*.cpp' -exec clang-format -style=file -i {} \;

# Format *.hpp files
find $(pwd) -name '*.hpp' -exec clang-format -style=file -i {} \;
