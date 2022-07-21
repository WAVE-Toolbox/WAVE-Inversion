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

Publication materials:
- Synthetic models and data used in Qin et al., (2022a) are available in https://bwsyncandshare.kit.edu/s/YJMgn8z7Siez4dL
- Synthetic models and data used in Qin et al., (2022b) are available in https://bwsyncandshare.kit.edu/s/bDsTJGPjzjcRQ4n
- Initial models and processed field data used in Qin et al., (2022c) are available in https://bwsyncandshare.kit.edu/s/ZyppjtcNcT5k2B8
- Initial models and processed field data used in Qin et al., (2022d) are available in https://bwsyncandshare.kit.edu/s/SRRnNfn8CYgtj5d

Reference:
- Qin, Tan, Thomas Bohlen, and Yudi Pan. "Indirect joint petrophysical inversion of synthetic shallow-seismic and multi-offset ground-penetrating radar data." Geophysical Journal International 229.3 (2022a): 1770-1784.
- Qin, Tan, Thomas Bohlen, and Niklas Allroggen. "Full-waveform inversion of ground-penetrating radar data in frequency-dependent media involving permittivity attenuation." Geophysical Journal International (2022b): In review.
- Qin, Tan, Thomas Bohlen, and Niklas Allroggen. "Fast full-waveform inversion of ground-penetrating radar long profile using data approximation." GEOPHYSICS (2022c): In preparation.
- Qin, Tan, Thomas Bohlen, Yudi Pan, and Niklas Allroggen. "Consistent Imaging of Near-Surface Targets Using Indirect Joint Petrophysical Inversion of Love-Wave and Multi-Offset Surface Ground-Penetrating Radar Field Data." Geophysical Research Letters (2022d): In preparation.
