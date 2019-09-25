## Geophysical inversion software **IFOS**

Master: [![build status](https://git.scc.kit.edu/WAVE/IFOS/badges/master/build.svg)](https://git.scc.kit.edu/WAVE/IFOS/commits/master) Develop: [![build status](https://git.scc.kit.edu/WAVE/IFOS/badges/develop/build.svg)](https://git.scc.kit.edu/WAVE/IFOS/commits/develop)

Geophysical inversion software developed by the [WAVE-Project](http://wave-toolbox.org).

Compilation:
`mkdir build && cd build`
`SOFI_DIR=<SOFI build dir> cmake ../src/ -DCMAKE_INSTALL_PREFIX=./`
`make install -j 4`

Documentation:
`make doc`
`make pdf`


Contributing:
- See the contributing guidelines provided in `CONTRIBUTING.md`

Download:
```
git clone https://git.scc.kit.edu/WAVE/IFOS.git
git submodule init
git submodule update
```
