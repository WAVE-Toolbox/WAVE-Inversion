from writeMtx import write_mtx

OUTPUT_FILENAME = 'sources.mtx'

## Required Parameters
X = [50] # Coordinates in X (Grid points)
Y = [50] # Coordinates in Y (Grid points) (depth)
Z = [50] # Coordinates in Z (Grid points)

SOURCE_TYPE = [1]  # Source Type (1=P, 2=vX, 3=vY, 4=vZ)
WAVELET_TYPE = [1] # Wavelet Type (1=Synthetic)

## Optional Parameters
WAVELET_SHAPE = [1] # Wavelet Shape (1=Ricker, 2=Sinw, 3=sin^3, 4=FGaussian, 5=Spike/Delta, 6=integral sin^3)
FC = [5]            # Center Frequency in Hz
AMP = [5]           # Amplitude
TShift = [0]        # Time shift in s

## Write to file

# Create Matrix
SOURCE_FILE=[X, Y, Z, SOURCE_TYPE, WAVELET_TYPE, WAVELET_SHAPE, FC, AMP, TShift,];

#Write mtx file
writeMatrix2mtx(OUTPUT_FILENAME,SOURCE_FILE)

