import numpy as np

def merge_targetHDU(hdu1, hdu2):
    # re-index first HDU
    
    # append second HDU, but only lines that are different
    lines = [tuple(t1) for t1 in hdu1.data]
    i = 1 + len(old_id1)
    old_id2 = hdu2.TARGET_ID
    new_id2 = np.zeros_like(old_id2)
    for j, t2 in enumerate(hdu2.data):
        same = same_target(hdu1, t2)
        if any(same):
            k = np.argwhere(same)[0,0]
            new_id2[j] = old_id1[k] 
        else:
            new_id2[j] = i
            i += 1 
            x = np.copy(t2)
            x['TARGET_ID'] = new_id2[j]
            lines.append(t2.tolist())
    #
    merged = hdu1 + hdu2[lines]
    hdulist1 =  hdu1.get_referenceHDU()  
    if hdusl
    return hdu, map1, map2
