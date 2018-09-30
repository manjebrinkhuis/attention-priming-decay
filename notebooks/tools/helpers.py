import nibabel as nib

def leading_zeros(string, delimiter='_', n_digits=3, exclude=None):
    """ 
    Returns a new string with added leading zeros.
    Parse string and for every section (as is 
    indicated by delimiter) try to convert it to an
    integer and if possible, add zeros (n_digits). 
    It is possible to exclude part of the string 
    (e.g. the file extension) by passing a string 
    "exclude". Only the part of the string up to the 
    first occurance of that string is parsed.    
    """
    
    if exclude:
        fname = string[:string.find(exclude)]
        ext = exclude
    else:
        fname = string
        ext = ''
    
    fname_info = fname.split(delimiter)
    
    for i, info in enumerate(fname_info):                
        try: 
            # try converting string to int
            num = int(info)
                    
            # if succeeded convert int to string with n_digits
            nnum = '%0' + str(n_digits) + 'd'
            nnum %= num
            
            # change string entry
            fname_info[i] = nnum
        except:
            pass
        
    return '_'.join(fname_info) + ext    


def select_files(source, include=[], exclude=None, extensions=['.PAR']):
    """ 
    Returns a list of filenames in source. Select files by 
    identifiers. If include is true, filenames from source
    that include at least one identifier are selected. Else,
    names that include at least one identifier are not 
    returned (while others are).
    Only files with extension are included.
    """
    
    files = os.listdir(source)
    files = [f for f in files if (True in [ext in f for ext in extensions])]
    
    out = []
    for fname in files:
        if include:
            for inc in include:
                if inc.upper() in fname.upper():
                    out.append(fname)
                    break
        else:
            out.append(fname)
        
    to_remove = []
    for fname in out:
        if exclude:
            for exc in exclude:                
                if exc.upper() in fname.upper():                    
                    to_remove.append(fname)
    
    for rm in to_remove:
        out.remove(rm)
        
    return out


def parrec2nii(files,
               verbose=True,
               outdir='./nii',
               compressed=True,
               origin='fov',
               minv='parse',
               maxv='parse',
               store_header=False,
               scaling='dv',
               overwrite=False,
               cmd="parrec2nii"):
        
    """ 
    Wrapper to use parrec2nii from python
    and not terminal. Uses a subprocess. 
    """

    # Create folders
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    # Options
    opts = ['-v' if verbose else '',
            '-o', outdir,
            '-c' if compressed else '',
            '--origin', origin,
            '--minmax', minv, maxv,
            '--store-header' if store_header else '',
            '--scaling', scaling]
    
    # Remove empty flags from list.
    while '' in opts:
        opts.remove('')
        
    # Run subprocess and provide some output.
    print('Running parrec2nii...')
        
    for fname in files:
        print('-'*20)
        
        # create new filename to check if file exists
        par_file = fname.split('/')[-1]
        extension = '.nii.gz' if compressed else '.nii'
        nii_file = par_file[:par_file.find('.PAR')] + extension

        nii_lead0 = leading_zeros(nii_file, exclude=extension)
        files_out = os.listdir(outdir)
        if overwrite or nii_lead0 not in files_out:
            print('convert "' + fname + '" ...')
            print('to "' + outdir + nii_lead0 + '".')
            subprocess.call([cmd] + opts + [fname])
        else:
            print('File "' +nii_lead0+ '" exists in "' +outdir+ '".')
            print('Skipped.')
            
        # rename
        if os.path.exists(os.path.join(outdir, nii_file)):
            os.rename(os.path.join(outdir, nii_file),
                      os.path.join(outdir, nii_lead0))
        
        print('-'*20)
        
    print('Done. Nii-files are in "' + outdir + '".')
    
def convert_to_nifti(in_path, out_path):
    img = nib.load(in_path)
    new = nib.Nifti1Image(img, affine=img.affine, header=img.header.as_analyze_map())       
    nib.save(new, out_path)
    
def save_nii(data, aff, fname):
    # Create nifti image.
    img = nib.Nifti1Image(data, aff)
    # Save as nifti.
    nib.save(img, fname)
    
def mni2index(mni):
    x_, y_, z_ = mni
    x = (90 - x_) / (90+90) * 90
    y = (y_ + 126) / (126+90) * 108
    z = (z_ + 72) / (72+108) * 90
    return int(x), int(y), int(z)

def index2mni(index):
    x_, y_, z_ = index    
    x = (x_ * ((90+90) / 90) - 90) * -1
    y = (y_ * ((126+90) /108) - 126)
    z = (z_ * ((72+108) / 90) - 72)    
    return int(x), int(y), int(z)
    