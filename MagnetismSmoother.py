import matplotlib.pyplot as plt
import numpy as np 
import sys
from functions import *
from astropy.io import fits 
from astropy.wcs import WCS
from scipy.ndimage import label
np.set_printoptions(threshold=np.inf)


PolFracFileExists = True            # Variable to grab data from pfrac file instead of calculating it
CreateSeparatePolFracFile = False   # Determines if polCalc_FITS makes a new file in the directory
IterateMultipleBoxCars = True       # Should the code run itself multiple times on different boxcars
BoxCarIncrementRange = [1, 8]       # Box car values to try if IterateMultipleBoxCars
incrementAmount = 1                 # Amount to increment boxcar by
BoxCarDictionary = {}               # Hash table with key box car and std value of residuals
BoxCarSize = 3                      # If = 3 then boxcar is 3x3, will matter if ~IterateMultipleBoxCars
QuiverMagnifier = 9                 # Will change the size of the lines in the plot uniformly

qFile = 'HH212M.flagging_brightlines.selfcal_avgch.StokesQ.finalmfs.image.fits'
uFile = 'HH212M.flagging_brightlines.selfcal_avgch.StokesU.finalmfs.image.fits'
iFile = 'HH212M.flagging_brightlines.selfcal_avgch.StokesI.finalmfs.image.fits'


def ApplyBoxcarAngleSmoothing(polarizationAngleData, mask, boxCarSize):
    filteredPolarizationAngleData = np.zeros_like(polarizationAngleData)  
    residualData = np.zeros_like(polarizationAngleData)      

    sourceWindowScale = SourceWindowCalculator()

    windowDimensionsSize = int(np.round(sourceWindowScale * boxCarSize))
    halfWindow = windowDimensionsSize // 2
    print("Dimension of Boxcar: ", windowDimensionsSize)

    # Go over the data but ignore the edges to avoid out of bounds
    for i in range(halfWindow, polarizationAngleData.shape[0] - halfWindow):
        for j in range(halfWindow, polarizationAngleData.shape[1] - halfWindow):
            
            window = polarizationAngleData[i-halfWindow:i+halfWindow+1, j-halfWindow:j+halfWindow+1]
            maskWindow = mask[i-halfWindow:i+halfWindow+1, j-halfWindow:j+halfWindow+1]
            smoothedValue, residual = ProcessWindow(window, windowDimensionsSize, maskWindow, polarizationAngleData[i, j])

            filteredPolarizationAngleData[i, j] = smoothedValue
            residualData[i, j] = residual

    return filteredPolarizationAngleData, residualData


def SourceWindowCalculator():  
    qHDU = fits.open(qFile)
    bMajor = qHDU[0].header.get('BMAJ')
    bMinor = qHDU[0].header.get('BMIN')
    
    return np.sqrt(bMajor * bMinor) * 3600 / 0.1


def ProcessWindow(window, windowSize, maskWindow, originalValue):
    validValues = window[~np.isnan(window) & maskWindow]
    
    if len(validValues) >= (windowSize * windowSize) // 2:  
        smoothedValue, residual = GetSmoothAngle(validValues, originalValue)
    else:
        smoothedValue = np.nan
        residual = np.nan
    
    return smoothedValue, residual


def GetSmoothAngle(validValues, originalValue):
    anglesRad = np.radians(validValues)

    xCoords = np.cos(anglesRad)
    yCoords = np.sin(anglesRad)

    meanX = np.mean(xCoords)
    meanY = np.mean(yCoords)

    smoothedAngleRad = np.arctan2(meanY, meanX)
    smoothedAngle = np.degrees(smoothedAngleRad)

    residual = min(
        abs(originalValue - smoothedAngle),
        180 - abs(originalValue - smoothedAngle)
    )

    if np.abs(residual) > 90:
        residual = np.nan

    return smoothedAngle, residual


def BigMain(boxCarSize):
    qHDU = fits.open(qFile)
    uHDU = fits.open(uFile)
    iHDU = fits.open(iFile)
    
    I_rms = 0.00012 #float(argv[4])
    Q_rms = 0.000069 #float(argv[5])
    U_rms = 0.000069 #float(argv[6])

    qData = qHDU[0].data[0, 0, :, :]  
    uData = uHDU[0].data[0, 0, :, :]  
    iData = iHDU[0].data[0, 0, :, :] 

    polarizationAngleData = 0.5 * np.degrees(np.arctan2(uData, qData)) 

    wcs = WCS(iHDU[0].header)
    if wcs.naxis > 2:
        wcs = wcs.slice([0, 0])


    polarizationFractionData = np.sqrt(qData**2 + uData**2) / iData
    polarizationFractionData[iData <= 0] = np.nan

    if PolFracFileExists: 
        sourceName = qFile.split(".")[0]
        pfracFileName = sourceName + ".pfrac.fits"
        
        polarizationFractionData = fits.open(pfracFileName)[0].data[0, 0, :, :]
        
    elif CreateSeparatePolFracFile:
        sourceName = qFile.split("_")[0]
        pfracFileName = sourceName + "_robust0.5_2sig"
        
        polCalc_FITS_file(iFile, I_rms, qFile, Q_rms, uFile, U_rms, sigma_QU=2, 
                          sigma_I=None, set_min_pfrac_rms=False, output_name=pfracFileName)
        
        print("New file pfrac file created.")
        sys.exit()
        
    else:
        polarizationFractionData = polCalc_FITS_return(iFile, I_rms, qFile, Q_rms, uFile, U_rms, sigma_QU=2, 
                                                       sigma_I=None, set_min_pfrac_rms=False, output_name=None)[1][0, 0, :, :]

    polarizedIntensity = iData * polarizationFractionData


    connectedComponentArray, numFeatures = label(polarizedIntensity > 0.05 * np.nanmax(polarizedIntensity))

    # Label connected components
    connectedComponentArray, numFeatures = label(polarizedIntensity > 0.05 * np.nanmax(polarizedIntensity))

    # Find the component with the maximum polarized intensity value
    maxIntensityComponent = None
    maxIntensity = -np.inf
    for component in range(1, numFeatures + 1):
        componentPixels = np.column_stack(np.where(connectedComponentArray == component))
        componentIntensity = polarizedIntensity[componentPixels[:, 0], componentPixels[:, 1]].max()
        if componentIntensity > maxIntensity:
            maxIntensity = componentIntensity
            maxIntensityComponent = component

    # Create mask for the component with the highest polarized intensity
    mask = (connectedComponentArray == maxIntensityComponent)

    filteredPolarizationAngleData, residualData = ApplyBoxcarAngleSmoothing(polarizationAngleData, mask, boxCarSize)

    largestConnectedComponentResiduals = residualData[mask]

    validResiduals = largestConnectedComponentResiduals[~np.isnan(largestConnectedComponentResiduals)]
    stdConnectedStructure = np.std(validResiduals)

    print(f"Standard Deviation of Residuals in Central LargestConnectedStructure: {stdConnectedStructure:.4f} degrees")
    BoxCarDictionary[np.round(boxCarSize, 1)] = np.round(stdConnectedStructure, 3)


    fig = plt.figure(figsize=(10, 8))
    ax = plt.subplot(projection=wcs)

    # Make sure to * 1000 in order to make it mJy / beam rather than Jy / beam
    im = ax.imshow(polarizedIntensity * 1000, cmap='gnuplot', origin='lower')
    plt.title('Residual Angular Deviation over Polarized Intensity Map (OMC-1)')
    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Polarized 870 μm Intensity (mJy/beam)')

    X, Y = np.meshgrid(np.arange(polarizationAngleData.shape[1]), np.arange(polarizationAngleData.shape[0]))

    XSlice = X[mask]
    YSlice = Y[mask]
    residualSlice = residualData[mask]
    filteredPolarizationAngleSlice = filteredPolarizationAngleData[mask]
    polarizationFractionSlice = polarizationFractionData[mask]
    theta_slice = polarizationAngleData[mask]

    # Plot the quiver (vector field) for residuals
    # Use cos and sin of residuals to get X and Y components of vectors
    # Multiply by polarizationFractionSlice to correlate the length with the strength of detection
    ax.quiver(XSlice, YSlice, 
              polarizationFractionSlice * QuiverMagnifier * np.cos(np.radians(residualSlice)), 
              polarizationFractionSlice * QuiverMagnifier * np.sin(np.radians(residualSlice)),
              color='red', scale=50, headlength=0, headaxislength=0, 
              label='Residual Angles')

    ax.quiver(XSlice, YSlice, 
              polarizationFractionSlice * QuiverMagnifier * np.cos(np.radians(filteredPolarizationAngleSlice)), 
              polarizationFractionSlice * QuiverMagnifier * np.sin(np.radians(filteredPolarizationAngleSlice)),
              color='green', scale=50, headlength=0, headaxislength=0, 
              label='Smoothed Angles')

    ax.quiver(XSlice, YSlice, 
              polarizationFractionSlice * QuiverMagnifier * np.cos(np.radians(theta_slice)), 
              polarizationFractionSlice * QuiverMagnifier * np.sin(np.radians(theta_slice)),
              color='pink', scale=50, headlength=0, headaxislength=0, 
              label='Normal Angles')

    plt.grid(False)
    plt.show()

    qHDU.close()
    uHDU.close()
    iHDU.close()


def main():
    if IterateMultipleBoxCars:
        i = BoxCarIncrementRange[0]
        while (i <= BoxCarIncrementRange[1]):
            BigMain(i)
            i += incrementAmount
        print("Box Car Values and Residuals: ") 
        print(BoxCarDictionary)
    else:
        BigMain(BoxCarSize)

main()
