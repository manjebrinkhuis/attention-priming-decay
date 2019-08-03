from scipy.stats import gamma
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaln
from scipy import interp

def spm_Gpdf(x, h, l):
    return np.exp(h * np.log(l) + (h-1) * np.log(x) - (l*x) - gammaln(h))

def create_hrf(time_unit=3.0/99,
               response_delay=6.0,
               undershoot_delay=16.0,
               response_dispersion=1.0,
               undershoot_dispersion=1.0,
               resp_undershoot_ratio=6.0,
               onset=0.0,
               kernel_length=32.0,
               fmri_t=16.0):
    """ """

    dt = time_unit / fmri_t
    u = np.arange(kernel_length / dt + 1) - (onset / dt)
    
    with np.errstate(divide='ignore'):
        g1 = spm_Gpdf(u, 
                      response_delay / response_dispersion,
                      dt / response_dispersion)
    
        g2 = spm_Gpdf(u,
                      undershoot_delay / undershoot_dispersion,
                      dt / undershoot_dispersion) / resp_undershoot_ratio

    hrf = g1 - g2
    idx = np.arange(int(kernel_length/time_unit + 1)) * fmri_t
    hrf = hrf[idx.astype(int)]
    hrf = hrf / np.sum(hrf)
    
    return hrf 


def create_regressor(length, onsets, amplitudes,
                     duration=.2, time_unit=.01, 
                     hrf=None, TR=2.1):
    """"""
    
    x = np.arange(0, length, time_unit)
    reg = np.zeros(len(x))
    
    for i, onset in enumerate(onsets):
        reg[(x > onset) * (x <= onset+duration)] = amplitudes[i]

    if hrf is None:
        hrf = create_hrf(time_unit=time_unit)
        
    conv = np.convolve(reg, hrf, mode="full")[:len(x)]    
    
    xselect = np.arange(length/TR)*TR
    yselect = interp(xselect, x, conv)

    return xselect, yselect
    

if __name__ == "__main__":
    
    hrf = create_hrf(time_unit=0.01)
               
    plt.figure()
    
    pulse_lengths = [.2, 4, 8]
    for i in range(3):
        
        onsets = [np.random.rand()*240 for onset in range(100)]
        t1, reg1 = create_regressor(240, onsets, time_unit=.016)   
        
        plt.subplot(3, 1, i+1)
        plt.plot(t1, reg1)
        plt.ylim(-.1, 1.1)

        
        
    