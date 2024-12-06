# Explanation of Scripts

MagnetismSmoother.py - Script that autmotates the process the applying a boxcar smoothing algorithm onto sources. MagnetismSmoother can be adjusted to run from either within the same directory as its target data or from a directory above.

MagnetismStrength.py - Script that calculates the magentic strength of source. MagnetismStrength runs on its own.

functions.py - Script that does specific calculations useful for creating a polarization fraction file. This script is not meant to be called on its own and is used by MagnetismSmoother, it must be in the same directory as MagnetismSmoother in order to work. This file is slightly different from the original version I was given.
