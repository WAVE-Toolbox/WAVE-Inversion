## Geophysical inversion software **WAVE-Inversion**

Master: [![build status](https://git.scc.kit.edu/WAVE/WAVE-Inversion/badges/master/build.svg)](https://git.scc.kit.edu/WAVE/WAVE-Inversion/commits/master) Develop: [![build status](https://git.scc.kit.edu/WAVE/WAVE-Inversion/badges/develop/build.svg)](https://git.scc.kit.edu/WAVE/WAVE-Inversion/commits/develop)

Geophysical inversion software developed by the [WAVE-Project](http://wave-toolbox.org).

Compilation:
`mkdir build && cd build`
`SIMULATION_DIR=<WAVE-Simulation build dir> cmake ../src/ -DCMAKE_INSTALL_PREFIX=./`
`make install -j 4`

Documentation:
`make doc`
`make pdf`


Contributing:
- See the contributing guidelines provided in `CONTRIBUTING.md`

