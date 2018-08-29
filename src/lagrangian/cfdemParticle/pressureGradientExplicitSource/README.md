# Modified `pressureGradientExplicitSource`
- Code modified so that the pressure gradient can write to the global IO object registry
- The value is then read by ShirgaonkarIB to set the pressure gradient force
- Note that the surface integral of pressure is equal to the volume integral of the pressure gradient (Gauss's theorem)

- To compile:
  1. overwrite the 'pressureGradientExplicitSource' folder in the OpenFOAM 'OpenFOAM-2.3.x/src/fvOptions/sources/derived/'
  2. navigate to 'OpenFOAM-2.3.1/src/fvOptions' 
  3. run `wmake`
