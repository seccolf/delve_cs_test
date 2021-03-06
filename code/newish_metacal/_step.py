from __future__ import print_function, absolute_import
import logging
import os

import numpy as np
import joblib
import esutil as eu
import fitsio
from ngmix import ObsList, MultiBandObsList
from ngmix.gexceptions import GMixRangeError

from ngmix.medsreaders import MultiBandNGMixMEDS, NGMixMEDS
#from .metacal import MetacalFitter
from metacal import MetacalFitter
#from .ngmix_compat import NGMIX_V2
from ngmix_compat import NGMIX_V2
#from eastlake.step import Step
#from eastlake.utils import safe_mkdir

import pdb

logger = logging.getLogger(__name__)

# always and forever
MAGZP_REF = 30.0

CONFIG = {
    'metacal': {
        # check for an edge hit
        'bmask_flags': 2**30,

        'metacal_pars': {
            'psf': 'fitgauss',
            'types': ['noshear', '1p', '1m', '2p', '2m'],
        },

        'model': 'gauss',

        'max_pars': {
            'ntry': 2,
            'pars': {
                'method': 'lm',
                'lm_pars': {
                    'maxfev': 2000,
                    'xtol': 5.0e-5,
                    'ftol': 5.0e-5,
                }
            }
        },

        'priors': {
            'cen': {
                'type': 'normal2d',
                'sigma': 0.263
            },

            'g': {
                'type': 'ba',
                'sigma': 0.2
            },

            'T': {
                'type': 'two-sided-erf',
                'pars': [-1.0, 0.1, 1.0e+06, 1.0e+05]
            },

            'flux': {
                'type': 'two-sided-erf',
                'pars': [-100.0, 1.0, 1.0e+09, 1.0e+08]
            }
        },

        'psf': {
            'model': 'gauss',
            'ntry': 2,
            'lm_pars': {
                'maxfev': 2000,
                'ftol': 1.0e-5,
                'xtol': 1.0e-5
            }
        }
    },
}


def _run_metacal(meds_files, seed):
    """Run metacal on a tile.

    Parameters
    ----------
    meds_files : list of str
        A list of the meds files to run metacal on.
    seed : int
        The seed for the global RNG.
    """
    with NGMixMEDS(meds_files[0]) as m:
        cat = m.get_cat() #LS: gets all the relevant info as a recarray (eg cat['ra'] contains the arrays of all cutouts)
    logger.info(' meds files %s', meds_files)

    #LS: the stuff below seems to be setting up parallel execution and nothing else.
    #a 'chunk' is apparently a subset of the objects in the meds file
    n_cpus = joblib.externals.loky.cpu_count()
    n_chunks = max(n_cpus, 60)
    n_obj_per_chunk = int(cat.size / n_chunks)
    if n_obj_per_chunk * n_chunks < cat.size:
        n_obj_per_chunk += 1
    assert n_obj_per_chunk * n_chunks >= cat.size
    logger.info(
        ' running metacal for %d objects in %d chunks', cat.size, n_chunks)

    seeds = np.random.RandomState(seed=seed).randint(1, 2**30, size=n_chunks)

    jobs = []
    for chunk in range(n_chunks):
        start = chunk * n_obj_per_chunk
        end = min(start + n_obj_per_chunk, cat.size)
        jobs.append(joblib.delayed(_run_mcal_one_chunk)( ###LS: seems like the function being run is _run_mcal_one_chunk
            meds_files, start, end, seeds[chunk]))

    with joblib.Parallel(
        n_jobs=n_cpus, backend='multiprocessing',
        verbose=100, max_nbytes=None
    ) as p:
        outputs = p(jobs)

    assert not all([o is None for o in outputs]), (
        "All metacal fits failed!")

    output = eu.numpy_util.combine_arrlist( #LS: no idea what's going on here, but seems to be simply organizing the outputs
        [o for o in outputs if o is not None])
    logger.info(' %d of %d metacal fits worked!', output.size, cat.size)

    return output


def _run_mcal_one_chunk(meds_files, start, end, seed): #LS: this is where the magic happens I guess
    """Run metcal for `meds_files` only for objects from `start` to `end`.

    Note that `start` and `end` follow normal python indexing conventions so
    that the list of indices processed is `list(range(start, end))`.

    Parameters
    ----------
    meds_files : list of str
        A list of paths to the MEDS files.
    start : int
        The starting index of objects in the file on which to run metacal.
    end : int
        One plus the last index to process.
    seed : int
        The seed for the RNG.

    Returns
    -------
    output : np.ndarray
        The metacal outputs.
    """
    rng = np.random.RandomState(seed=seed)

    # seed the global RNG to try to make things reproducible
    np.random.seed(seed=rng.randint(low=1, high=2**30))

    output = None
    mfiles = []
    data = []
    try:
        # get the MEDS interface
        for m in meds_files:
            mfiles.append(NGMixMEDS(m)) #LS: puts all open meds files in a list
        mbmeds = MultiBandNGMixMEDS(mfiles) #LS: not sure what a multiband meds is
        cat = mfiles[0].get_cat()

        for ind in range(start, end):
            o = mbmeds.get_mbobs(ind)
            #LS: pdb
            pdb.set_trace()
            o = _strip_coadd(o) #LS: removes the first image of an object, since that's the coadd I think (and afaik we jointly fit shapes to the single epochs, not the coadd - why?)
            pdb.set_trace()
            o = _strip_zero_flux(o) #LS: removes images with zero total flux - when would that happen??
            pdb.set_trace()
            if not NGMIX_V2:
                # ngmix v1 worked in surface brightness, not flux
                o = _apply_pixel_scale(o) #LS: no idea

            skip_me = False
            for ol in o:
                if len(ol) == 0:
                    logger.debug(' not all bands have images - skipping!')
                    skip_me = True
            if skip_me:
                continue

            o.meta['id'] = ind
            o[0].meta['Tsky'] = 1
            o[0].meta['magzp_ref'] = MAGZP_REF
            o[0][0].meta['orig_col'] = cat['orig_col'][ind, 0]
            o[0][0].meta['orig_row'] = cat['orig_row'][ind, 0]

            nband = len(o)
            mcal = MetacalFitter(CONFIG, nband, rng)

            try:
                mcal.go([o])
                res = mcal.result
            except GMixRangeError as e:
                logger.debug(" metacal error: %s", str(e))
                res = None

            if res is not None:
                data.append(res)

        if len(data) > 0:
            output = eu.numpy_util.combine_arrlist(data)
    finally:
        for m in mfiles:
            m.close()

    return output


def _strip_coadd(mbobs):
    _mbobs = MultiBandObsList()
    _mbobs.update_meta_data(mbobs.meta)
    for ol in mbobs:
        _ol = ObsList()
        _ol.update_meta_data(ol.meta)
        for i in range(1, len(ol)):
            _ol.append(ol[i])
        _mbobs.append(_ol)
    return _mbobs


def _strip_zero_flux(mbobs):
    _mbobs = MultiBandObsList()
    _mbobs.update_meta_data(mbobs.meta)
    for ol in mbobs:
        _ol = ObsList()
        _ol.update_meta_data(ol.meta)
        for i in range(len(ol)):
            if np.sum(ol[i].image) > 0:
                _ol.append(ol[i])
        _mbobs.append(_ol)
    return _mbobs


def _apply_pixel_scale(mbobs):
    for ol in mbobs:
        for o in ol:
            scale = o.jacobian.get_scale()
            scale2 = scale * scale
            scale4 = scale2 * scale2
            o.image = o.image / scale2
            o.weight = o.weight * scale4
    return mbobs
