import numpy as np
import subprocess
import os
import shlex
import sys
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.io.fits as pyfits
from scipy.optimize import minimize
from scipy.special import iv


def sh(s, log=None) :
    if log :
        sh('touch {0}'.format(log))
        f = open(log, 'a')
        # Deal with parentheses and quotes
        s = s.replace( '"', '' )
        s = s.replace( '(', '"(' )
        s = s.replace( ')', ')"' )
        f.write(s)
        f.write('\n')
        f.close()

    return subprocess.Popen(shlex.split(s)).wait()

def sh_pipe(s, log=None, output=False) :
    if log :
        sh('touch {0}'.format(log))
        f = open(log, 'a')
        # Deal with parentheses and quotes
        s = s.replace( '"', '' )
        s = s.replace( '(', '"(' )
        s = s.replace( ')', ')"' )
        f.write(s)
        f.write('\n')
        f.close()

    proc = subprocess.Popen(shlex.split(s), stdout=subprocess.PIPE, \
stderr=subprocess.PIPE)

    stdout,stderr = proc.communicate('stdin')

    # Print stdout and stderr
    if output :
        print(stdout)
        print(stderr)

    if not len(stderr) :
        return stdout.split('\n'),stderr
    else :
        return stdout.split('\n'),stderr.split('\n')

def comment(s, log) :
    sh('touch {0}'.format(log))
    f = open(log, 'a')
    f.write(s)
    f.write('\n')
    f.close()

def HMS_to_decimal(s) :
    # NOTE: THIS IS ONLY GOOD FOR RA!
    hms = list(map(float, s.split(':')))
    decimal = hms[0] + hms[1]/60 + hms[2]/3600 
    return decimal

def decimal_to_HMS(d) :
	# NOTE: THIS IS ONLY GOOD FOR RA!
	
	# Stay between 0 and 24
	if d < 0 :
		d = (24 - abs(d)%24) % 24
	elif d >= 24 :
		d %= 24
	    
	h = int(np.trunc(d))
	rem1 = (d - h)*60
	m = int(np.trunc(rem1))
	rem2 = (rem1 - m)*60
	s = int(np.trunc(round(rem2)))
	if s > 59 :
		s %= 60
	m += 1
	if m > 59 :
		m %= 60
		h += 1
	if h > 24 :
		h %= 24
	hms = '%02d:%02d:%02d'%(h,m,s)
	return hms
	
def isfloat(x) :
	try:
		float(x)
		return True
	except:
		return False

def iscomment(s):
    return s.startswith('#')

def text2ps(file, out_path='./', out_file=None) :
    prefix = file.split('.')[0]
    if not out_file :
        out_file = prefix + '.ps'
    sh('enscript -r {0} -o {1}{2}'.format(file, out_path, out_file))

def pdf2png(in_file, out_file=None) :
    if not out_file :
        out_file = in_file.strip('.pdf')+'.png'
    sh('gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 \
-sOutputFile={0} {1}'.format(out_file, in_file))
    # sh('mogrify -rotate 90 {0}'.format(out_file))

def vec_color(color_ps, old_color="1024 512 0", new_color="0 0 0") :
    f = open(color_ps, 'r')
    lines = f.readlines()
    f.close()
    # Replace RGB vector color with grayscale
    for i in range(len(lines)) :
        if "{0} K".format(old_color) in lines[i] :
            lines[i] = lines[i].replace("{0} K".format(old_color), "{0} K".format(new_color))
    f = open(color_ps, 'w')
    for line in lines :
        f.write(line)
    f.close()

def beam_color(color_ps, old_color="0 0 1024", new_color="0 0 0") :
    f = open(color_ps, 'r')
    lines = f.readlines()
    f.close()
    # Replace RGB beam color with grayscale
    j = 0
    for i in range(len(lines)) :
        if "{0} K".format(old_color) in lines[i] :
            if j :
                lines[i] = lines[i].replace("{0} K".format(old_color), "{0} K".format(new_color))
            j += 1
    f = open(color_ps, 'w')
    for line in lines :
        f.write(line)
    f.close()

def ps2eps(ps) :
    # This function will work, but will rotate
    # the resulting file...using the ps2eps version
    # installed in the OS is better
    name = ps.split('.')[0]
    sh('ps2epsi {0}.ps'.format(name))
    sh('mv {0}.epsi {0}.eps'.format(name))

def get_activeScopes(vis_file) :
    # Get the active antennas
    log_file_1 = 'listobs.log'   
    # Delete any previous temporary listobs log file
    sh('rm -f {0}'.format(log_file_1)) 
    sh('listobs vis={0} log={1}'.format(vis_file,log_file_1)) 
    listobs_log = open(log_file_1,'r')
    lines = listobs_log.readlines()
    listobs_log.close()

    active_scopes=[]
    for line in lines:
        if line.find(':') != -1 and line.find('Antenna') != -1:
            # Find active telescopes
            active_scopes.append(int(line.split(':')[0].split()[1]))
    print('\nThese antennas were used in the observation:')
    print(active_scopes)
    return active_scopes

def get_uvindex(vis_file) :
    log_file = 'uvindex.log'   
    # Delete any previous temporary listobs log file
    sh('rm -f {0}'.format(log_file)) 
    sh('uvindex vis={0} log={1}'.format(vis_file,log_file)) 
    return log_file

def get_listobs(vis_file) :
    log_file = 'listobs.log'   
    # Delete any previous temporary listobs log file
    sh('rm -f {0}'.format(log_file)) 
    sh('listobs vis={0} log={1}'.format(vis_file,log_file)) 
    return log_file

def get_obsDate(uvindex_log='uvindex.log') :
    f = open(uvindex_log)
    lines = f.readlines()
    f.close()

    months = [['JAN','01'],
              ['FEB','02'],
              ['MAR','03'],
              ['APR','04'],
              ['MAY','05'],
              ['JUN','06'],
              ['JUL','07'],
              ['AUG','08'],
              ['SEP','09'],
              ['OCT','10'],
              ['NOV','11'],
              ['DEC','12']]

    for i in range(len(lines)) :
        if len(lines[i].split()) > 1 :
            if lines[i].split()[1] == 'Total' :
                date = lines[i].split()[0]
                year = '20' + date[0:2]
                for m in months :
                    if date[2:5] == m[0] :
                        month = m[1]
                        break
                day = date[5:7]
    if len(date) :
        return '{0}{1}{2}'.format(year,month,day)
    else :
        raise Exception( 'ERROR: Date not found in uvindex.log' )

def get_RA_DEC ( vis_file, source ) :
    # Get RA, DEC of source
    lines,err = sh_pipe('uvlist vis={0} options=variables select="source({1})"'.format(vis_file, source))

    RA = []
    DEC = []
    for line in lines:
        for i in range(len(line.split())) :
            if line.split()[i] == 'dec' :
                DEC = line.split()[i+2]
            elif line.split()[i] == 'ra' :
                RA = line.split()[i+2]
    return RA,DEC

def get_bandInfo ( uvindex_log='uvindex.log' ) :
    '''
    Returns for both wide and narrow bands:
    -- Selection strings for windows 
    -- Frequencies for each band
    -- Number of channels in each band
    '''

    # Get information for wideband windows
    win_wb = ''
    chan_wb = ''
    freq_wb = ''
    velres_wb = ''

    win_nb = ''
    chan_nb = ''
    freq_nb = ''
    velres_nb = ''
    velres_wb = ''
    # Speed of light in km/s
    c = 2.99 * 10**5

    f = open(uvindex_log)
    lines = f.readlines()
    f.close()
    for i in range(len(lines)) :
        if len(lines[i].split()) > 0 :
            if lines[i].split()[0] == 'Frequency' and lines[i].split()[2] == '2' :
                break
            if lines[i].split()[0] == 'Spectral' :
                for j in range(1,17) :
                    # Calculate bandwidth to determine WB or NB
                    freq = float( lines[i+j].split()[1] )
                    increment = abs( float(lines[i+j].split()[2]) )
                    BW = int(lines[i+j].split()[0]) * increment
                    if BW > 0.48 :
                        win_wb += '{0},'.format(j)
                        chan_wb += lines[i+j].split()[0] + ','
                        freq_wb += lines[i+j].split()[1] + ','
                        velres_wb += str( round(c*(increment/freq), 2) ) + ','
                    elif 0 < BW < 0.48 :
                        win_nb += '{0},'.format(j)
                        chan_nb += lines[i+j].split()[0] + ','
                        freq_nb += lines[i+j].split()[1] + ','
                        velres_nb += str( round(c*(increment/freq), 2) ) + ','
    win_wb = win_wb.strip(',')
    win_LSB = win_wb[0:len(win_wb)/2].strip(',')
    win_USB = win_wb[len(win_wb)/2 : len(win_wb)].strip(',')
    win_nb = win_nb.strip(',')
    chan_wb = chan_wb.strip(',')
    chan_nb = chan_nb.strip(',')
    freq_wb = freq_wb.strip(',')
    freq_nb = freq_nb.strip(',')
    velres_wb = velres_wb.strip(',')
    velres_nb = velres_nb.strip(',')

    # Calculate number of wideband and narrowband windows
    num_wb = len(chan_wb.split(','))
    num_nb = len(chan_nb.split(','))

    # Calculate frequencies for grand average, LSB, and USB (wideband only)
    freq_ALL = np.mean( list(map(float, freq_wb.split(','))) )
    freq_LSB = np.mean( list(map(float, freq_wb.split(',')[0:len(freq_wb.split(','))/2])) )
    freq_USB = np.mean( list(map(float, freq_wb.split(',')[len(freq_wb.split(','))/2 : len(freq_wb.split(','))])) )

    return win_wb,chan_wb,freq_wb,win_nb,chan_nb,freq_nb,win_LSB,win_USB,freq_LSB,freq_USB,freq_ALL,num_wb,num_nb,velres_nb,velres_wb

def formatPar ( parString ) :
    '''
    formatPar ( parString )
 
    Returns parString as string, bool, float or int, 
    depending on what it is:
    
    -- bool is 'True', 'true', 'False' or 'false'
    -- int is all digits
    -- float has ONE '.' in it and the pieces on either side of 
       the '.' are digits
    -- Otherwise, it returns a string 
    '''

    if parString in ['True','true'] :
        return True
    elif parString in ['False','false'] :
        return False
    elif parString.isdigit() :
        return int(parString)
    elif len(parString.split('.')) == 2 : 
        s = parString.split('.')
        if s[0].isdigit() and s[1].isdigit() :
            return float(parString)
    return parString


def peakData ( map_name, region="arcsec,box(7)" ) :
    if not os.path.exists(map_name) :
        raise Exception('{0} not found.'.format(map_name))
    lines,err = sh_pipe('imstat in={0} region={1}'.format(map_name, region))

    for i in range(len(lines)) :
        if len(lines[i].split()) > 0 :
            if "Total" in lines[i] :
                min_flux = float( lines[i+1][61:71] )
                max_flux = float( lines[i+1][51:61] )
                fluxes = [ min_flux, max_flux ] 
                peak_index = list( np.abs(fluxes) ).index( max(np.abs(fluxes)) )
                peak = float( fluxes[peak_index] )
                print('\nPeak of %s = %.4f Jy/bm\n'%(map_name, peak))

    return peak


def noiseList ( map_name, region=None, moment=False ) :
    # box(xmin, ymin, xmax, ymax)
    regions = []
    vals = []    
    
    # Sample four different regions
#    r1 = 'arcsec,box(30,-25,45,25)(1)'
#    r2 = 'arcsec,box(-45,-25,-30,25)(1)'
#    r3 = 'arcsec,box(-25,30,25,45)(1)'
#    r4 = 'arcsec,box(-25,-45,25,-30)(1)'

    r1 = 'arcsec,box(15,-15,30,15)(1)'
    r2 = 'arcsec,box(-30,-15,-15,15)(1)'
    r3 = 'arcsec,box(-15,25,15,30)(1)'
    r4 = 'arcsec,box(-15,-30,15,-15)(1)'

    if not region :
        regions = [r1, r2, r3, r4]
    else :
        regions = [region]

    # Looks for noise value in cleaned maps
    print('\n')
    for r in regions :
        if not os.path.exists(map_name) :
            raise Exception('{0} not found.'.format(map_name))
        lines,err = sh_pipe('imlist options=stat in={0} region={1}'.format(map_name, r))
        for line in lines:
            if len(line.split()) > 0 :
                if 'K/Jy:' in line.split() and isfloat(line.split()[-1]) :
                    KJy = float(line.split()[-1])
                    print('Conversion factor = {0} K/Jy\n'.format(KJy))
                if isfloat(line.split()[0]) :
                    val = float(line.split()[6])
                    vals.append(val)
                    print('RMS of %s = %.5f Jy\n'%(map_name, val))
                    # if region : print '\n'

    if not moment :
        rms = np.median(vals)
    else :
        rms = np.min(vals)

    if not region :
        if not moment :
            print('\nMedian RMS of %s = %.5f Jy/bm\n'%(map_name, rms))
        else :
            print('\nRMS of %s in selected velocity range = %.5f Jy/bm\n'%(map_name, rms))

    return rms,KJy


def fits ( map_name, type ) :
    sh('fits in={0} op={1} out={0}.fits'.format(map_name, type))


def coords ( source, region, path='./' ) :
    # Fit Gaussian (size of synthesized beam) to dust peak
    lines,err = sh_pipe('imfit in={0}{1}.dust.I.cm region={2} object=point'.format(
            path, source, region))
    for line in lines :
        for i,val in enumerate(line.split()) :
            if val == "Right" :
                RA = line.split()[-1]
            if val == "Declination:" :
                DEC = line.split()[-1]

    return RA,DEC


def beam_size ( source, path='./' ) :
    lines,err = sh_pipe('imlist in={0}{1}.dust.I.cm'.format(path, source))

    for line in lines :
        for i,val in enumerate(line.split()) :
            if val == "bmin" :
                bmin = float( line.split()[i+2].split(':')[-1] )
            if val == "bmaj" :
                bmaj = float( line.split()[i+2].split(':')[-1] )
            if val == "bpa" :
                bpa = float( line.split()[i+2].split(':')[-1] )

    return bmin,bmaj,bpa


def avgAngle ( source, region, path='./', weight=None, type=None ) :

    results = []
    angles = []
    errors = []
    intensity = []
    pixel=0

    for file in ['{0}{1}.dust.pa.cm'.format(path,source),
                 '{0}{1}.dust.pa_err.cm'.format(path,source),
                 '{0}{1}.dust.I.cm'.format(path,source)] :

        logfile = file.split('.')[2] + "_" + source + ".log"
        out,err = sh_pipe('imlist in={0} options=data region={1} log={2}{3}'.format(file, region, path, logfile), output=True)
        f = open(path+logfile, "r")
        lines = f.readlines()
        f.close()
    
        for i in range(9, len(lines)-2) :
            for num in lines[i].split()[0:len(lines[i].split())-1] :
                if isfloat(num) :
                    if 'pa' in file.split('.') :
                        angles.append(float(num))
                        pixel += 1                        
                    elif 'pa_err' in file.split('.') :
                        errors.append(float(num))
                    elif 'I' in file.split('.') :
                        intensity.append(float(num))
                
                # If there is no pol detection, set a flag
                else :
                    if 'pa' in file.split('.') :
                        angles.append('X')
                    if 'pa_err' in file.split('.') :
                        errors.append('X')

    # Set intensity to zero where there is no pol detection
    for i in range(len(angles)) :
        if angles[i] == 'X' :
            intensity[i] = 0
            angles[i] = 0
            errors[i] = 0

    # WEIGHTED AVERAGE

    # Factor required for averaging angles with different ranges
    # Polars with 180 deg degeneracy require a factor of two
    if type == 360 :
        factor = 1.
    elif type == 180 :
        factor = 2.
    else : raise Exception('You must use either type=180 or type=360')

    # Average cosine and sine to avoid problems around +-90 deg
    c = cos( factor*radians(angles) )
    s = sin( factor*radians(angles) )
    c_err = cos( factor*radians(errors) )
    s_err = sin( factor*radians(errors) )

    # Weighted angles and errors by Stokes I intensity
    if weight.lower() == 'intensity' :
        c_avg = sum( array(c) * array(intensity) ) / sum( array(intensity) )
        s_avg = sum( array(s) * array(intensity) ) / sum( array(intensity) )

        c_avg_err = sum( array(c_err) * array(intensity) ) / sum( array(intensity) )
        s_avg_err = sum( array(s_err) * array(intensity) ) / sum( array(intensity) )

    elif weight.lower() == 'variance' :
        c_avg = sum( array(c) * array(errors)**2 ) / sum( array(errors)**2 )
        s_avg = sum( array(s) * array(errors)**2 ) / sum( array(errors)**2 )

        c_avg_err = sum( array(c_err) * array(errors)**2 ) / sum( array(errors)**2 )
        s_avg_err = sum( array(s_err) * array(errors)**2 ) / sum( array(errors)**2 )        

    else : raise Exception('You must supply a weighting scheme.')

    results.append( degrees( arctan2( s_avg, c_avg ) ) / factor )
    results.append( degrees( arctan2( s_avg_err, c_avg_err ) ) / factor )

    print("\nSource = {0}".format(source))
    print("Number of pixels with pol detection = {0}".format(pixel))
    print("Average B-field angle = %.1f +- %.1f deg\n\n"%(results[0]+90, results[1]))

    # Return B-field angle
    return results[0]+90, results[1]


def polDebias ( image_array, rms ) :
     
    print('\nDebiasing polarization image...\n')

    pol = image_array
     
    # print 'Image shape'
    # print pol.shape
     
    # Normalize by the rms
    pol /= rms
     
    # Array of debiased values
    real = []
     
    for obs in pol.ravel() : 
        if isnan(obs) :
            real.append(0)
        elif obs > 9 : 
            # Big values: use usual debiasing limit of P_c = sqrt(P^2 - rms^2)
            real.append( np.sqrt(obs**2 - 1) )
        else :
            # Debias using PDF of p_true given p_obs
            def PDF( p_true, p_obs=obs, rms=1 ) :
                bessel_param = p_obs*p_true/rms**2
                bessel_param [ bessel_param > 700 ] = 0
                # This is the PDF of p_true given p_obs
                P = p_obs/rms**2 * iv(0, bessel_param ) * np.exp(-(p_obs**2 + p_true**2)/(2*rms**2))
                # Print the negative of the function so it can be minimized
                return -P

            # Guess for minimization
            if obs < 3 :
                guess = 1
            else :
                guess = obs

            # Minimize the PDF; the most likely p_true value is at the minimum of the (negated) PDF
            P_min = minimize( PDF, guess )
            
            # Turn tiny negative values into 0
            if P_min.x[0] < 0 :
                real.append(0)
            else :
                real.append(P_min.x[0])
                     
    # Turn back to actual pol values
    real = np.array(real)*rms

    # Reshape to original array shape
    real = np.reshape(real, pol.shape)
     
    return real

    # Test values
    # print pol*rms
    # print nanmax(pol)
    # print real*rms
    # print real.max()


def polCalc_FITS_file ( I, I_rms, Q, Q_rms, U, U_rms, path='./', sigma_I=None, sigma_QU=None, output_name=None, set_min_pfrac_rms = True, min_pfrac_rms = 0.001 ) :
    # Prefix of files to write
    if not output_name :
        name = I.split('.')[0]+'.'+I.split('.')[1]
    else :
        name = output_name

    # Open FITS files
    # Get header for Stokes I
    h = pyfits.open (path+I)[0].header
    # Open I, Q, U files
    I = pyfits.open (path+I)[0].data 
    Q = pyfits.open (path+Q)[0].data 
    U = pyfits.open (path+U)[0].data 

    # Average Q,U rms
    rms_QU = np.average( [Q_rms, U_rms] )

    # Make masks
    PA_mask = np.empty_like(I)
    PA_mask[:] = np.NAN
    PA_rms_mask = np.empty_like(I)
    PA_rms_mask[:] = np.NAN
    pol_debias_mask = np.empty_like(I)
    pol_debias_mask[:] = np.NAN
    pfrac_mask = np.empty_like(I)
    pfrac_mask[:] = np.NAN
    pfrac_rms_mask = np.empty_like(I)
    pfrac_rms_mask[:] = np.NAN

    # Pol intensity (debiased) calculation
    pol = np.sqrt(Q**2 + U**2)
  
    # This tends toward plain old sigma when U_rms = Q_rms
    ##### NOTE: this is not correct for low-SNR measurements!
    # The correct way to get the RMS in the pol would be to take the area
    #  under the PDF that encompasses 68% of the probability.
    pol_rms = (1/pol) * np.hypot( (U_rms * U) , (Q_rms * Q) )

    # Debias pol intensity
    pol_debias = polDebias( pol, rms_QU )


    # Fractional polarization and rms 
    pfrac = pol_debias / I
    pfrac_rms = pfrac * np.sqrt( (pol_rms/pol_debias)**2 + (I_rms/I)**2 )
    
    # Angle calculation
    PA = 0.5 * np.degrees ( np.arctan2 (U,Q) )    
    PA_rms = 0.5 * degrees( sqrt( (Q*U_rms)**2 + (U*Q_rms)**2 ) / pol_debias**2 )

    ##### This is only true if you assume that sigmaQ = sigmaU
    # PA_rms = 0.5 * np.degrees( pol_rms / pol_debias )
 
    # Mask based on sigma_I only, sigma_QU only, or both
    if sigma_I and not sigma_QU :
        mask = I > sigma_I*I_rms
    elif sigma_QU and not sigma_I :
        mask = pol_debias > sigma_QU*pol_rms
    # Mask on BOTH sigma_QU and sigma_I :
    elif sigma_QU and sigma_I :
        mask = np.logical_and(I > sigma_I*I_rms, pol_debias > sigma_QU*pol_rms)
 
    # Consider minimum pfrac
    if set_min_pfrac_rms:
    	pfrac_rms[pfrac_rms < min_pfrac_rms] = min_pfrac_rms
    	mask = np.logical_and(mask,pfrac > sigma_QU*pfrac_rms)
     
    # Mask files to be written
    pol_debias_mask [ mask ] = pol_debias [ mask ]
    PA_mask [ mask ] = PA [ mask ]
    PA_rms_mask [ mask ] = PA_rms [ mask ]
    pfrac_mask [ mask ] = pfrac [ mask ]
    pfrac_rms_mask [ mask ] = pfrac_rms [ mask ]
     
    # Write pol intensity (debiased)
    os.system('rm -r {0}{1}.pol.fits'.format(path, name)) 
    pyfits.writeto('{0}{1}.pol.fits'.format(path, name), pol_debias_mask, h)
     
    # Write pfrac
    os.system('rm -r {0}{1}.pfrac.fits'.format(path, name)) 
    pyfits.writeto('{0}{1}.pfrac.fits'.format(path, name), pfrac_mask, h)
     
    # Write pfrac rms
    os.system('rm -r {0}{1}.pfrac_rms.fits'.format(path, name)) 
    pyfits.writeto('{0}{1}.pfrac_rms.fits'.format(path, name), pfrac_rms_mask, h)

    # Write PA
    os.system('rm -r {0}{1}.pa.fits'.format(path, name)) 
    pyfits.writeto('{0}{1}.pa.fits'.format(path, name), PA_mask, h)

    # Write PA rms
    os.system('rm -r {0}{1}.pa_rms.fits'.format(path, name)) 
    pyfits.writeto('{0}{1}.pa_rms.fits'.format(path, name), PA_rms_mask, h)


def polCalc_FITS_return ( I, I_rms, Q, Q_rms, U, U_rms, path='./', sigma_I=None, sigma_QU=None, output_name=None, set_min_pfrac_rms = True, min_pfrac_rms = 0.001 ) :
    # Prefix of files to write
    if not output_name :
        name = I.split('.')[0]+'.'+I.split('.')[1]
    else :
        name = output_name

    # Open FITS files
    # Get header for Stokes I
    h = pyfits.open (path+I)[0].header
    # Open I, Q, U files
    I = pyfits.open (path+I)[0].data 
    Q = pyfits.open (path+Q)[0].data 
    U = pyfits.open (path+U)[0].data 

    # Average Q,U rms
    rms_QU = np.average( [Q_rms, U_rms] )

    # Make masks
    PA_mask = np.empty_like(I)
    PA_mask[:] = np.NAN
    PA_rms_mask = np.empty_like(I)
    PA_rms_mask[:] = np.NAN
    pol_debias_mask = np.empty_like(I)
    pol_debias_mask[:] = np.NAN
    pfrac_mask = np.empty_like(I)
    pfrac_mask[:] = np.NAN
    pfrac_rms_mask = np.empty_like(I)
    pfrac_rms_mask[:] = np.NAN

    # Pol intensity (debiased) calculation
    pol = np.sqrt(Q**2 + U**2)
  
    # This tends toward plain old sigma when U_rms = Q_rms
    ##### NOTE: this is not correct for low-SNR measurements!
    # The correct way to get the RMS in the pol would be to take the area
    #  under the PDF that encompasses 68% of the probability.
    pol_rms = (1/pol) * np.hypot( (U_rms * U) , (Q_rms * Q) )

    # Debias pol intensity
    pol_debias = polDebias( pol, rms_QU )


    # Fractional polarization and rms 
    pfrac = pol_debias / I
    pfrac_rms = pfrac * np.sqrt( (pol_rms/pol_debias)**2 + (I_rms/I)**2 )
    
    # Angle calculation
    PA = 0.5 * np.degrees ( np.arctan2 (U,Q) )    
    PA_rms = 0.5 * degrees( sqrt( (Q*U_rms)**2 + (U*Q_rms)**2 ) / pol_debias**2 )

    ##### This is only true if you assume that sigmaQ = sigmaU
    # PA_rms = 0.5 * np.degrees( pol_rms / pol_debias )
 
    # Mask based on sigma_I only, sigma_QU only, or both
    if sigma_I and not sigma_QU :
        mask = I > sigma_I*I_rms
    elif sigma_QU and not sigma_I :
        mask = pol_debias > sigma_QU*pol_rms
    # Mask on BOTH sigma_QU and sigma_I :
    elif sigma_QU and sigma_I :
        mask = np.logical_and(I > sigma_I*I_rms, pol_debias > sigma_QU*pol_rms)
 
    # Consider minimum pfrac
    if set_min_pfrac_rms:
    	pfrac_rms[pfrac_rms < min_pfrac_rms] = min_pfrac_rms
    	mask = np.logical_and(mask,pfrac > sigma_QU*pfrac_rms)
     
    # Mask files to be written
    pol_debias_mask [ mask ] = pol_debias [ mask ]
    PA_mask [ mask ] = PA [ mask ]
    PA_rms_mask [ mask ] = PA_rms [ mask ]
    pfrac_mask [ mask ] = pfrac [ mask ]
    pfrac_rms_mask [ mask ] = pfrac_rms [ mask ]
     
    # Return data as a list of arrays
    return [pol_debias_mask, pfrac_mask, pfrac_rms_mask, PA_mask, PA_rms_mask]


def polCalc ( Q, Q_rms, U, U_rms, I, I_rms, FITS=False, rot90=False ) :
    Q = float(Q)
    Q_rms = float(Q_rms)
    U = float(U)
    U_rms = float(U_rms)
    I = float(I)
    I_rms = float(I_rms)
    
    # This is the amount of polarized radiation
    pol = np.hypot( Q, U )
    # This is the uncertainty in the amount of polarized radiation.  It's NOT a "percent error"; 
    # i.e. if there's a 10% error in the 5% polarization measurement, this program
    # will output (5+-0.5)%, NOT 5% +- 10%
    pol_rms = (1/pol) * np.hypot( (U_rms * U) , (Q_rms * Q) )

    # Note: this is ONLY true for high SNR!
    ##### Need to implement method from polDebias()
    pol_debias = np.sqrt( pol**2 - pol_rms**2 )
    print('\nNOTE: the debiased polarization np.sqrt( pol**2 - pol_rms**2 ) is only valid for high SNR!')

    # Angle calculation
    PA = 0.5 * np.degrees ( np.arctan2 (U,Q) )    
    PA_rms = 0.5 * degrees( sqrt( (Q*U_rms)**2 + (U*Q_rms)**2 ) / pol_debias**2 )

    ##### This is only true if you assume that sigmaQ = sigmaU
    # PA_rms = 0.5 * np.degrees( pol_rms / pol_debias )
    
    # Finally, calculate percentages
    pct = pol_debias / I * 100
    pct_rms = pct * np.sqrt( (pol_rms/pol)**2 + (I_rms/I)**2 ) * 100
    
    # In the end, if you get a three-sigma Q,U detection (i.e. hypot(Q,U)/rms_QU), you get a +- 10 deg in P.A.
    # Double the sigma (6-sigma), halve the uncertainty (+-5 deg)

    return pct,pct_rms,pa,pa_rms


def writeLeakSets ( leakTablePath='/o/chat/c_soft/scripts/ant_leaks/',
                    leakSetPath='/o/chat/c_data/leakages/') :
    """
    Purpose:
        Uses do_leakAnalyze to write out a different set of leakages for each day, 
        in a set of frequency ranges.

        It writes the leakage table into skeleton.mir, and then uses GPCOPY to create a table
        with an appropriate name.

        freq should be a range of frequencies
          freq=[226,230]
    """ 

    ##### NOTE: C7 currently uses D_terms for 239 GHz for any frequency >239 GHz.
    # Create README file
    f = open('{0}README.txt'.format(leakSetPath),'w')
    f.write('NOTE: C7 currently uses a 239 GHz leakage soltion for any frequency >239 GHz.\n')
    f.write('Last updated on 06/18/2012, by Chat Hull.\n')
    f.close()

    # Crucial date ranges:
    date_ranges = [ '00000000-20111025',
               '20111026-20111102',
               '20111103-20120117',
               '20120118-20120416',
               '20120417-current' ]

    obs_dates = [ 20111001,
                  20111101,
                  20120101,
                  20120201,
                  20120501 ]

    freq_ranges_old = [ '216,222',
                        '221,227',
                        '226,232' ]

    freq_ranges_new = [ '206,212',
                    '211,217',
                    '216,222',
                    '221,227',
                    '226,232',
                    '231,237',
                    '236,242',
                    '241,247',
                    '246,252',
                    '251,257',
                    '256,262',
                    '261,267',
                    '266,272']

    freqs_old = [ 219,
                  224,
                  229 ]

    freqs_new = [ 209,
                  214,
                  219,
                  224,
                  229,
                  234,
                  239,
                  244,
                  249,
                  254,
                  259,
                  264,
                  269 ]

    for i in range(len(date_ranges)) :
        if int( date_ranges[i].split('-')[0] ) >= 20120417 :
            ranges = freq_ranges_new
            freqs = freqs_new
        else :
            ranges = freq_ranges_old
            freqs = freqs_old
        for j in range(len(ranges)) :
            sh('delhd in={0}skeleton.mir/leakage'.format(leakTablePath))
            sh('do_leakAnalyze write=True interact=False verbose=False freq={0} obs_date={1} \
target={2}skeleton.mir write_leakset=True'.format(ranges[j], obs_dates[i], leakTablePath ))
            sh('gpcopy vis={0}skeleton.mir out={3}{1}.{2}.leak mode=create \
options=nocal,nopass'.format(leakTablePath, date_ranges[i], freqs[j], leakSetPath))
    os.system('tar -cvf leak_current.tar *.leak README.txt')
