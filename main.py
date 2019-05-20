import numpy as np
import matplotlib.pyplot as plt


def bracewell_buneman(xarray, length, log2length):
    '''
    bracewell-buneman bit reversal function
    inputs: xarray is array; length is array length; log2length=log2(length).
    output: bit reversed array xarray.
    '''
    muplus = int((log2length+1)/2)
    mvar = 1
    reverse = np.zeros(length, dtype = int)
    upper_range = muplus+1
    for _ in np.arange(1, upper_range):
        for kvar in np.arange(0, mvar):
            tvar = 2*reverse[kvar]
            reverse[kvar] = tvar
            reverse[kvar+mvar] = tvar+1
        mvar = mvar+mvar
    if (log2length & 0x01):
        mvar = mvar/2
    for qvar in np.arange(1, mvar):
        nprime = qvar-mvar
        rprimeprime = reverse[int(qvar)]*mvar
        for pvar in np.arange(0, reverse[int(qvar)]):
            nprime = nprime+mvar
            rprime = rprimeprime+reverse[pvar]
            temp = xarray[int(nprime)]
            xarray[int(nprime)] = xarray[int(rprime)]
            xarray[int(rprime)] = temp
    return xarray



def dif_fft0 (xarray, twiddle, log2length):
    '''
    radix-2 dif fft
    '''
    xarray = xarray.astype(np.complex_)
    b_p = 1
    nvar_p = xarray.size
    twiddle_step_size = 1
    for _ in range(0,  log2length):           # pass loop
        nvar_pp =  int(nvar_p/2)
        base_e = 0
        for _ in range(0,  b_p):       # block loop
            base_o = int(base_e + nvar_pp)
            for nvar in range(0, nvar_pp):   # butterfly loop

                evar =  xarray[int(base_e + nvar)]+xarray[int(base_o + nvar)]
                if nvar == 0:
                    ovar = xarray[int(base_e + nvar)]-xarray[int(base_o + nvar)]
                else:
                    twiddle_factor =  nvar*twiddle_step_size
                    ovar = (xarray[int(base_e + nvar)] \
                        -xarray[int(base_o + nvar)])*twiddle[twiddle_factor]
                xarray[int(base_e + nvar)] = evar
                xarray[int(base_o + nvar)] = ovar
            base_e = base_e+nvar_p
        b_p = b_p*2
        nvar_p = int(nvar_p/2)
        twiddle_step_size = 2*twiddle_step_size
    xarray = bracewell_buneman(xarray, xarray.size, log2length)
    return xarray


def test(time, yarray, samplefreq):
    '''
    Set up plot, call FFT function, plot result.
    Called  from testbench function.
    Inputs time:time vector, yarray: array, samplefreq: sampling rate.
    Outputs: none.
    '''
    plt.subplot(2, 1, 1)
    plt.title('Test of DIF FFT with 10 Hz Sine Input')
    plt.plot(time, yarray ,'k-')
    plt.xlabel('time')
    plt.ylabel('amplitude')
    plt.subplot(2, 1, 2)
    ylength = len(yarray)                       # length of the signal
    kvar = np.arange(ylength)
    tvar = ylength/samplefreq
    frq = kvar/tvar                        # two-sided frequency range
    freq = frq[list(range(int(ylength/2)))]           # one-sided frequency range
    twiddle = np.exp(-2j*np.pi*np.arange(0, 0.5, 1./ylength, dtype=np.complex_))
    y2array = abs(dif_fft0(yarray, twiddle, \
              int(np.log2(ylength)))/ylength)   # fft normalized magnitude
    y3array = y2array[list(range(ylength//2))]
    markerline, stemlines, baseline = plt.stem(freq, y3array, '--')
    plt.xlabel('freq (Hz)')
    plt.ylabel('|Y(freq)|')
    plt.ylim((0.0, 0.55))
    plt.setp(markerline, 'markerfacecolor', 'b')
    plt.setp(baseline, 'color', 'b', 'linewidth', 2)
    plt.show()
    return None


def testbench():
    '''
    no inputs, no outputs, calls test function.
    '''
    samplefreq = 512                                # sampling rate
    samplinginterval = 1.0/samplefreq               # sampling interval
    time = np.arange(0, 1, samplinginterval)        # time vector
    frequency = 100                                  # frequency of the signal
    yarray = 2 * np.sin(2 * np.pi * frequency * time + (4*np.pi/2))   # 10 hz sine wave signal
    test(time, yarray, samplefreq)                 # send sine to test
    return None


testbench()

def main():
    pass

if __name__ == '__main__':
    main()
