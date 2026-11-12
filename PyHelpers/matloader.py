import numpy as np
import scipy.io as sio
from scipy.io.matlab import mat_struct

def check_keys(d):
    """Recursively convert mat_struct objects to nested dictionaries."""
    for key in d:
        if isinstance(d[key], mat_struct) or isinstance(d[key], sio.matlab._mio5_params.mat_struct):
            d[key] = to_dict(d[key])
    return d

def to_dict(matobj):
    d = {}
    for field in matobj._fieldnames:
        elem = getattr(matobj, field)
        if isinstance(elem, mat_struct) or isinstance(elem, sio.matlab._mio5_params.mat_struct):
            d[field] = to_dict(elem)
        elif isinstance(elem, str):
            d[field] = elem
        else:
            if isinstance(elem, int) or isinstance(elem, float):
                d[field] = elem
            else:
                d[field] = np.array(elem)
    return d

def loadmat(filename):
    """Load a MATLAB .mat file (v7 or earlier) into nested Python dicts."""
    md_dict = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    md_dict = check_keys(md_dict)
    return md_dict

def construct_model_key(md_dict):
    return f"ts{md_dict['time_step']}_ct{md_dict['calving_timescale']}_cd{md_dict['calving_crevasse_threshold']}_{md_dict['mesh_res']}"