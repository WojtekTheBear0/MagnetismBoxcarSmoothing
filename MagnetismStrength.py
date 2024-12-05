import numpy as np

# Source and boxcar size key, std of residual angles value
sourceAngleData = {
    "HH212M": {1: 3.613, 2: 5.728, 3: 6.373, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HH270IRS": {1: 8.94, 2: 13.047, 3: 18.868, 4: 16.594, 5: 17.846, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-10": {1: 2.704, 2: np.nan, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-11": {1: 15.901, 2: 21.753, 3: 21.316, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-12E": {1: 3.335, 2: 9.151, 3: 11.978, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-12W": {1: 7.788, 2: 13.102, 3: 15.298, 4: 18.363, 5: 19.272, 6: 19.337, 7: 18.937, 8: 11.515},
    "HOPS-50": {1: 2.921, 2: 3.366, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-53": {1: 15.799, 2: 11.509, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-60": {1: 10.994, 2: 20.311, 3: 2.686, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-78": {1: 8.812, 2: 14.777, 3: 17.529, 4: 20.946, 5: 24.49, 6: 21.964, 7: np.nan, 8: np.nan},
    "HOPS-81": {1: 3.609, 2: np.nan, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-84": {1: 15.114, 2: 6.451, 3: 4.103, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-87": {1: 1.789, 2: 3.774, 3: 5.276, 4: 6.379, 5: 7.679, 6: 8.62, 7: 8.942, 8: 8.448},
    "HOPS-88": {1: 4.037, 2: 7.088, 3: 8.969, 4: 8.622, 5: 6.434, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-96": {1: 3.16, 2: 7.715, 3: 10.031, 4: 9.752, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-124": {1: 21.562, 2: 28.085, 3: 27.788, 4: 25.333, 5: 4.335, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-153": {1: 3.833, 2: 1.807, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-168": {1: 10.137, 2: 11.642, 3: 9.227, 4: 8.297, 5: 8.094, 6: 8.552, 7: 9.024, 8: np.nan},
    "HOPS-169": {1: 7.602, 2: 13.194, 3: 17.623, 4: 19.523, 5: 21.116, 6: 22.628, 7: 20.869, 8: np.nan},
    "HOPS-182": {1: 9.878, 2: 13.65, 3: 16.263, 4: 18.276, 5: 20.832, 6: 23.117, 7: 24.216, 8: 23.413},
    "HOPS-203N": {1: 18.404, 2: 24.837, 3: 29.695, 4: 21.5, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-203S": {1: 15.821, 2: 18.613, 3: 22.003, 4: 22.651, 5: 14.309, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-224": {1: 5.267, 2: 7.785, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-247": {1: 2.141, 2: 1.183, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-250": {1: 3.468, 2: np.nan, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-288": {1: 11.37, 2: 15.737, 3: 18.464, 4: 19.132, 5: 19.462, 6: 19.188, 7: np.nan, 8: np.nan},
    "HOPS-303": {1: 11.147, 2: 20.786, 3: 24.756, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-310": {1: 4.334, 2: 5.005, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-317S": {1: 14.23, 2: 22.051, 3: 30.568, 4: 20.872, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-325": {1: 17.134, 2: 25.927, 3: 27.89, 4: 25.185, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-340": {1: 3.238, 2: np.nan, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-341": {1: 14.385, 2: 24.708, 3: 30.254, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-354": {1: 3.876, 2: np.nan, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-358": {1: 5.351, 2: 6.365, 3: 1.856, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-359": {1: 12.899, 2: 19.989, 3: 23.804, 4: 27.225, 5: 30.058, 6: 31.652, 7: 31.642, 8: 30.532},
    "HOPS-361N": {1: 8.183, 2: 13.174, 3: 17.476, 4: 19.734, 5: 20.945, 6: 22.526, 7: 22.683, 8: 22.59},
    "HOPS-370": {1: 4.503, 2: 6.687, 3: 7.819, 4: 8.702, 5: 9.399, 6: 8.982, 7: 8.857, 8: 8.541},
    "HOPS-373EW": {1: 9.064, 2: 5.254, 3: 3.442, 4: 2.622, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-383": {1: 3.729, 2: 2.629, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-384": {1: 11.775, 2: 16.489, 3: 19.648, 4: 21.385, 5: 22.72, 6: 23.834, 7: 25.203, 8: 26.829},
    "HOPS-395": {1: 18.652, 2: 23.245, 3: 19.692, 4: 14.227, 5: 17.27, 6: 14.557, 7: np.nan, 8: np.nan},
    "HOPS-398": {1: 10.834, 2: 18.159, 3: 22.42, 4: 23.08, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-399": {1: 6.561, 2: 8.557, 3: 10.664, 4: 12.172, 5: 13.514, 6: 14.767, 7: 15.912, 8: 16.418},
    "HOPS-400": {1: 3.852, 2: 8.387, 3: 11.653, 4: 12.926, 5: 12.036, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-401": {1: 15.917, 2: 22.577, 3: 24.64, 4: 19.438, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-402": {1: 4.949, 2: 5.734, 3: np.nan, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-403": {1: 13.552, 2: 16.297, 3: 17.455, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-404": {1: 9.386, 2: 14.045, 3: 13.931, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-407": {1: 2.447, 2: 4.13, 3: 4.926, 4: 5.354, 5: 5.887, 6: 6.219, 7: np.nan, 8: np.nan},
    "HOPS-408": {1: 14.846, 2: 5.772, 3: 1.987, 4: np.nan, 5: np.nan, 6: np.nan, 7: np.nan, 8: np.nan},
    "HOPS-409": {1: 2.395, 2: 4.471, 3: 6.774, 4: 8.836, 5: 7.99, 6: np.nan, 7: np.nan, 8: np.nan},
    "OMC1N-4-5-E": {1: 10.888, 2: 14.414, 3: 18.295, 4: 20.987, 5: 22.438, 6: 24.149, 7: 26.095, 8: 27.878},
    "OMC1N-6-7": {1: 9.934, 2: 13.35, 3: 17.57, 4: 19.433, 5: 20.627, 6: 22.275, 7: 23.42, 8: 24.18},
    "OMC1N-8-N": {1: 4.759, 2: 6.999, 3: 9.968, 4: 11.534, 5: 13.095, 6: 14.461, 7: 15.808, 8: 17.2}
}

# Source key, with 2 value array key, 0 is H2 density and 1 is std of velocity
sourceDictionary = {
    "HH212M": [6.6E+06, 0.73], "HH270IRS": [7.2E+06, 0.77], "HOPS-10": [2.6E+06, 0.68], 
    "HOPS-11": [5.4E+06, 0.72], "HOPS-124": [9.9E+06, 1.51], "HOPS-12E": [4.3E+06, 0.78], 
    "HOPS-12W": [6.1E+06, 0.63], "HOPS-153": [3.3E+06, 0.92], "HOPS-164": [4.0E+06, 0.77], 
    "HOPS-168": [3.5E+06, 0.90], "HOPS-169": [9.7E+06, 0.63], "HOPS-182": [8.7E+06, 1.15], 
    "HOPS-203N": [3.8E+06, 0.86], "HOPS-203S": [3.4E+06, 0.71], "HOPS-224": [7.1E+06, 0.78], 
    "HOPS-247": [7.1E+06, 0.84], "HOPS-250": [3.3E+06, 0.92], "HOPS-288": [8.9E+06, 1.38], 
    "HOPS-303": [7.8E+06, 0.73], "HOPS-310": [6.8E+06, 0.95], "HOPS-317S": [4.0E+07, 0.64], 
    "HOPS-325": [4.6E+06, 0.65], "HOPS-340": [4.2E+06, 0.70], "HOPS-341": [4.8E+06, 0.63], 
    "HOPS-354": [3.8E+06, 0.90], "HOPS-358": [4.6E+06, 1.23], "HOPS-359": [1.3E+07, 0.65], 
    "HOPS-361N": [7.5E+06, 0.92], "HOPS-361S": [4.8E+06, 0.88], "HOPS-370": [4.3E+06, 1.22], 
    "HOPS-373EW": [None, None], "HOPS-383": [2.9E+06, 0.89], "HOPS-384": [5.8E+06, 1.02], 
    "HOPS-395": [1.0E+07, 0.70], "HOPS-398": [1.2E+07, 0.47], "HOPS-399": [2.6E+07, 0.70], 
    "HOPS-400": [2.1E+07, 0.82], "HOPS-401": [7.7E+06, 0.63], "HOPS-402": [1.1E+07, 0.66], 
    "HOPS-403": [1.6E+07, 1.12], "HOPS-404": [1.1E+07, 0.70], "HOPS-407": [1.2E+07, 0.64], 
    "HOPS-408": [5.6E+06, 0.94], "HOPS-409": [4.8E+06, 0.78], "HOPS-50": [4.8E+06, 1.01], 
    "HOPS-53": [3.2E+06, 0.74], "HOPS-60": [5.0E+06, 0.87], "HOPS-78": [9.9E+06, 0.77], 
    "HOPS-81": [2.4E+06, 0.74], "HOPS-84": [4.2E+06, 0.92], "HOPS-87": [3.0E+07, 0.70], 
    "HOPS-88": [7.1E+06, 0.76], "HOPS-96": [6.0E+06, 0.57], "OMC1N-4-5-E": [None, None], 
    "OMC1N-6-7": [None, None], "OMC1N-8-N": [None, None]
}


# Print of a single source with a single boxcar
currentSource = "HOPS-407"
currentBoxcarSize = 1

densityOfH2 = sourceDictionary.get(currentSource)[0]
stdVelocity = sourceDictionary.get(currentSource)[1]
velocityDispersion = stdVelocity * np.sqrt(8 * np.log(2)) * 1000
stdTheta = sourceAngleData.get(currentSource, {}).get(currentBoxcarSize, np.nan)
meanMolecularMass = 2.37
G = 6.6743 * (10**-11)

magneticFieldStrength = 9.3 * np.sqrt(densityOfH2) * (velocityDispersion/stdTheta) * meanMolecularMass * G

print("Magnetism Strength: {:.6e} Gauss".format(magneticFieldStrength))


# Print of the avergae of all boxcars across every available source
averageForBoxcars = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0}
tempSumStrength = 0

for i in range(1, 9):
    numSourcesWithSize = 0
    for key in sourceDictionary:
        if not np.isnan(sourceAngleData.get(key, {}).get(i, np.nan)) and sourceDictionary.get(key)[0] is not None:
            densityOfH2 = sourceDictionary.get(key)[0]
            stdVelocity = sourceDictionary.get(key)[1]
            velocityDispersion = stdVelocity * np.sqrt(8 * np.log(2)) * 1000
            stdTheta = sourceAngleData.get(key, {}).get(i, np.nan)
            meanMolecularMass = 2.37
            G = 6.6743 * (10**-11)

            magneticFieldStrength = 9.3 * np.sqrt(densityOfH2) * (velocityDispersion/stdTheta) * meanMolecularMass * G

            averageForBoxcars[i] += magneticFieldStrength 
            numSourcesWithSize += 1
    averageForBoxcars[i] /= numSourcesWithSize

print("\nBoxcar Size Averages: ", averageForBoxcars)
