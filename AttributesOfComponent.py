import numpy as np


distanceCm = 1.197e+21
pixelSizeInArcSeconds = 4
effectiveBeamSize = 14.1

def GetNumberValidPixels(component):
    count = 0

    for value in component:
        if value != np.NaN:
            count += 1
    return count


def PrintNumberValidPixels(component):
    print(f'Number of valid pixels: {GetNumberValidPixels(component)}')


def AreaCalculation(component):
    count = GetNumberValidPixels(component)
    
    print(f'Angular area in arcseconds^2: {pixelSizeInArcSeconds**2 * count}')
    print(f'Area in cm^2: {distanceCm**2 * count * (pixelSizeInArcSeconds * np.pi / 648000)**2}')


def VolumeCalculation(component):
    count = GetNumberValidPixels(component)
    area = distanceCm**2 * count * (pixelSizeInArcSeconds * np.pi / 648000)**2

    volume = (4/3) * (np.pi **(-1/2)) * (area ** (3/2))
    print(f'Volume in cm^3: {volume}')


def FluxCalculation(mask, iHDU):
    data = np.squeeze(iHDU[0].data)

    beam_area = 1.133 * effectiveBeamSize**2
    pixel_area = pixelSizeInArcSeconds**2
    total_intensity = np.sum(data[mask])
    print(f'Intensity: {total_intensity}')

    jy_pixel_per_beam = total_intensity * 725
    flux_density = jy_pixel_per_beam * (pixel_area / beam_area)

    print(f'Flux Density in Jy: {flux_density}')

