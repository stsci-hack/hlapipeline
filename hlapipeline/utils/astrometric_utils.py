"""Utilities to support creation of astrometrically accurate reference catalogs

The function, create_astrometric_catalog, allows the user to query an
astrometric catalog online to generate a catalog of astrometric sources that
should fall within the field-of-view of all the input images.

This module relies on the definition of an environment variable to specify
the URL of the astrometric catalog to use for generating this
reference catalog.

    ASTROMETRIC_CATALOG_URL  -- URL of web service that can be queried to
                                obtain listing of astrometric sources,
                                sky coordinates, and magnitudes.

"""
import os
from io import BytesIO

import csv
import requests
from lxml import etree
import inspect

import numpy as np

import stwcs
from stwcs.distortion import utils
from stwcs import wcsutil
from stsci.tools import fileutil as fu
from stsci.tools import parseinput

from astropy import units as u
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy.io import fits as pf
from astropy.io import ascii
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import photutils
from photutils import detect_sources, source_properties
from photutils import Background2D, MedianBackground
from scipy.spatial import distance_matrix
from scipy import ndimage

import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

import pysynphot as S

from hlapipeline.utils import bitmask

ASTROMETRIC_CAT_ENVVAR = "ASTROMETRIC_CATALOG_URL"
DEF_CAT_URL = 'http://gsss.stsci.edu/webservices'

if ASTROMETRIC_CAT_ENVVAR in os.environ:
    SERVICELOCATION = os.environ[ASTROMETRIC_CAT_ENVVAR]
else:
    SERVICELOCATION = DEF_CAT_URL

ISOLATION_LIMIT = 0.55

MODULE_PATH = os.path.dirname(inspect.getfile(inspect.currentframe()))
VEGASPEC = os.path.join(os.path.dirname(MODULE_PATH),
                        'data','alpha_lyr_stis_008.fits')

__all__ = ['create_astrometric_catalog', 'compute_radius', 'find_gsc_offset',
            'extract_sources', 'find_isolated_source', 'find_best_offset']

def create_astrometric_catalog(inputs, **pars):
    """Create an astrometric catalog that covers the inputs' field-of-view.

    Parameters
    ===========
    input : str
        Filenames of images to be aligned to astrometric catalog

    catalog : str, optional
        Name of catalog to extract astrometric positions for sources in the
        input images' field-of-view. Default: GSC241. Options available are
        documented on the catalog web page.

    output : str, optional
        Filename to give to the astrometric catalog read in from the master
        catalog web service.  If 'None', no file will be written out.
        Default: ref_cat.ecsv

    gaia_only : bool, optional
        Specify whether or not to only use sources from GAIA in output catalog
        Default: False

    note ::
        This function will point to astrometric catalog web service defined
        through the use of the ASTROMETRIC_CATALOG_URL environment variable.

    Returns
    =======
    ref_table : object
        Astropy Table object of the catalog

    """
    # interpret input parameters
    catalog = pars.get("catalog", 'GSC241')
    output = pars.get("output", 'ref_cat.ecsv')
    gaia_only = pars.get("gaia_only", False)
    table_format = pars.get("table_format", 'ascii.ecsv')

    inputs, _ = parseinput.parseinput(inputs)
    # start by creating a composite field-of-view for all inputs
    wcslist = []
    for img in inputs:
        nsci = fu.countExtn(img)
        for num in range(nsci):
            extname = '{}[sci,{}]'.format(img, num+1)
            wcslist.append(stwcs.wcsutil.HSTWCS(extname))

    # This default output WCS will have the same plate-scale and orientation
    # as the first chip in the list, which for WFPC2 data means the PC.
    # Fortunately, for alignment, this doesn't matter since no resampling of
    # data will be performed
    outwcs = utils.output_wcs(wcslist)
    radius = compute_radius(outwcs)
    ra, dec = outwcs.wcs.crval

    # perform query for this field-of-view
    ref_dict = get_catalog(ra, dec, sr=radius, catalog=catalog)
    colnames = ('ra','dec', 'mag', 'objID', 'GaiaID')
    col_types = ('f8', 'f8', 'f4', 'U25', 'U25')
    ref_table = Table(names = colnames, dtype=col_types)

    # extract just the columns we want...
    num_sources = 0
    for source in ref_dict:
        if 'GAIAsourceID' in source:
            g = source['GAIAsourceID']
            if gaia_only and g.strip() is '':
                continue
        else:
            g = -1  # indicator for no source ID extracted
        r = float(source['ra'])
        d = float(source['dec'])
        m = float(source['mag'])
        o = source['objID']
        num_sources += 1
        ref_table.add_row((r,d,m,o,g))

    # Write out table to a file, if specified
    if output:
        ref_table.write(output, format=table_format)
    print("Created catalog '{}' with {} sources".format(output, num_sources))

    return ref_table


def get_catalog(ra, dec, sr=0.1, fmt='CSV', catalog='GSC241'):
    """ Extract catalog from VO web service.

    Parameters
    ----------
    ra : float
        Right Ascension (RA) of center of field-of-view (in decimal degrees)

    dec : float
        Declination (Dec) of center of field-of-view (in decimal degrees)

    sr : float, optional
        Search radius (in decimal degrees) from field-of-view center to use
        for sources from catalog.  Default: 0.1 degrees

    fmt : str, optional
        Format of output catalog to be returned.  Options are determined by
        web-service, and currently include (Default: CSV):
            VOTABLE(default) | HTML | KML | CSV | TSV | JSON | TEXT

    catalog : str, optional
        Name of catalog to query, as defined by web-service.  Default: 'GSC241'

    Returns
    -------
    csv : obj
        CSV object of returned sources with all columns as provided by catalog

    """
    serviceType = 'vo/CatalogSearch.aspx'
    spec_str = 'RA={}&DEC={}&SR={}&FORMAT={}&CAT={}&MINDET=5'
    headers = {'Content-Type': 'text/csv'}

    spec = spec_str.format(ra, dec, sr, fmt, catalog)
    serviceUrl = '{}/{}?{}'.format(SERVICELOCATION, serviceType,spec)
    rawcat = requests.get(serviceUrl, headers=headers)
    r_contents = rawcat.content.decode()  # convert from bytes to a String
    rstr = r_contents.split('\r\n')
    # remove initial line describing the number of sources returned
    # CRITICAL to proper interpretation of CSV data
    del rstr[0]
    r_csv = csv.DictReader(rstr)

    return r_csv


def compute_radius(wcs):
    """Compute the radius from the center to the furthest edge of the WCS."""

    ra,dec = wcs.wcs.crval
    img_center = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    wcs_foot = wcs.calc_footprint()
    img_corners = SkyCoord(ra=wcs_foot[:,0]*u.degree,
                           dec=wcs_foot[:,1]*u.degree)
    radius = img_center.separation(img_corners).max().value

    return radius

def find_gsc_offset(image, input_catalog='GSC1', output_catalog='GAIA'):
    """Find the GSC to GAIA offset based on guide star coordinates

    Parameters
    ----------
    image : str
        filename of image to be processed

    Returns
    -------
    delta_ra,delta_dec : tuple of floats
        Offset in decimal degrees of image based on correction to guide star
        coordinates relative to GAIA
    """
    serviceType = "GSCConvert/GSCconvert.aspx"
    spec_str = "TRANSFORM={}-{}&IPPPSSOOT={}"

    if 'rootname' in pf.getheader(image):
        ippssoot = pf.getval(image, 'rootname').upper()
    else:
        ippssoot = fu.buildNewRootname(image).upper()

    spec = spec_str.format(input_catalog, output_catalog, ippssoot)
    serviceUrl = "{}/{}?{}".format(SERVICELOCATION, serviceType,spec)
    rawcat = requests.get(serviceUrl)
    if not rawcat.ok:
        print("Problem accessing service with:\n{{}".format(serviceUrl))
        raise ValueError

    delta_ra = delta_dec = None
    tree = BytesIO(rawcat.content)
    for _,element in etree.iterparse(tree):
        if element.tag == 'deltaRA':
            delta_ra = float(element.text)
        elif element.tag == 'deltaDEC':
            delta_dec = float(element.text)

    return delta_ra,delta_dec


def extract_sources(img, fwhm=3.0, threshold=None, source_box=7,
                    classify=True, output=None, plot=False, vmax=None):
    """Use photutils to find sources in image based on segmentation."""
    if threshold is None:
        bkg_estimator = MedianBackground()
        bkg = Background2D(img, (50, 50), filter_size=(3, 3),
                           bkg_estimator=bkg_estimator)
        threshold = bkg.background + (3. * bkg.background_rms)
    sigma = fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=source_box, y_size=source_box)
    kernel.normalize()
    segm = detect_sources(img, threshold, npixels=source_box,
                          filter_kernel=kernel)
    cat = source_properties(img, segm)
    print("Total Number of detected sources: {}".format(len(cat)))
    if classify:
        # Remove likely cosmic-rays based on central_moments classification
        goodsrcs = np.where(classify_sources(cat) == 1)[0].tolist()
        newcat = photutils.segmentation.properties.SourceCatalog([])
        for src in goodsrcs:
            newcat._data.append(cat[src])
    else:
        newcat = cat

    tbl = newcat.to_table()
    print("Final Number of selected sources: {}".format(len(newcat)))
    if output:
        tbl['xcentroid'].info.format = '.10f'  # optional format
        tbl['ycentroid'].info.format = '.10f'
        tbl['source_sum'].info.format = '.10f'
        tbl['cxy'].info.format = '.10f'
        tbl['cyy'].info.format = '.10f'
        if not output.endswith('.cat'):
            output += '.cat'
        tbl.write(output, format='ascii.no_header')
        print("Wrote source catalog: {}".format(output))
        
    if plot:
        norm = None
        if vmax is None:
            norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        ax1.imshow(img, origin='lower', cmap='Greys_r', norm=norm, vmax=vmax)
        ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345))

    return tbl, segm

def classify_sources(catalog, sources=None):
    """ Convert moments_central attribute for source catalog into star/cr flag

    This algorithm interprets the central_moments from the source_properties
    generated for the sources as more-likely a star or a cosmic-ray.  It is not
    intended or expected to be precise, merely a means of making a first cut at
    removing likely cosmic-rays or other artifacts.

    Parameters
    -----------
    catalog : object
        The photutils.SourceCatalog object for the image/chip

    sources : tuple
        Range of objects from catalog to process as a tuple of (min, max).
        Default: None which simply processes all sources.

    Returns
    -------
    srctype : ndarray
        An ndarray where a value of 1 indicates a likely valid, non-cosmic-ray
        source, and a value of 0 indicates a likely cosmic-ray.
    """
    moments = catalog.moments_central
    if sources is None:
        sources = (0,len(moments))
    num_sources = sources[1] - sources[0]
    srctype = np.zeros((num_sources,),np.int32)
    for src in range(sources[0],sources[1]):
        x,y = np.where(moments[src] == moments[src].max())
        if (x[0] > 1) and (y[0] > 1):
            srctype[src] = 1

    return srctype


def build_source_catalog(image, refwcs, **kwargs):
    dqname = kwargs.get('dqname','DQ')
    output = kwargs.get('output',None)
    # Build source catalog for entire image
    master_cat = None
    numSci = countExtn(image, extname='SCI')
    for chip in range(numSci):
        chip += 1
        # find sources in image
        if output:
            rootname = image[0].header['rootname']
            outroot = '{}_sci{}_src'.format(rootname, chip)
            kwargs['output'] = outroot
        imgarr = image['sci',chip].data
        photmode = image['sci',chip].header['photmode']
        # apply any DQ array, if available
        if image.index_of(dqname):
            dqarr = image[dqname,chip].data
            dqmask = bitmask.bitfield_to_boolean_mask(dqarr, good_mask_value=False)
            imgarr[dqmask] = 0. # zero-out all pixels flagged as bad
        seg_tab, segmap = extract_sources(imgarr, **kwargs)
        seg_tab_phot = compute_photometry(seg_tab,photmode)
        # Convert to sky coordinates
        chip_wcs = wcsutil.HSTWCS(image,ext=('sci',chip))
        seg_ra,seg_dec = chip_wcs.all_pix2world(seg_tab_phot['xcentroid'],seg_tab_phot['ycentroid'],1)
        # Convert sky positions to pixel positions in the reference WCS frame
        seg_xy_out = refwcs.all_world2pix(seg_ra,seg_dec,1)
        seg_tab_phot['xcentroid'] = seg_xy_out[0]
        seg_tab_phot['ycentroid'] = seg_xy_out[1]
        if master_cat is None:
            master_cat = seg_tab_phot
        else:
            master_cat = vstack([master_cat, seg_tab_phot])

    return master_cat

def compute_photometry(catalog, photmode):
    """ Compute magnitudes for sources from catalog based on observations photmode


    Returns
    --------
    phot_cat : object
        Astropy Table object of input source catalog with added column for
        VEGAMAG photometry (in magnitudes).
    """
    # Determine VEGAMAG zero-point using pysynphot for this photmode
    photmode = photmode.replace(' ',',')
    vega = S.FileSpectrum(VEGASPEC)
    bp = S.ObsBandpass(photmode)
    vegauvis = S.Observation(vega,bp)
    vegazpt = 2.5*np.log10(vegauvis.countrate())

    # Use zero-point to convert flux values from catalog into magnitudes
    source_phot = vegazpt - 2.5*np.log10(catalog['source_sum'])
    source_phot.name = 'vegamag'
    # Now add this new column to the catalog table
    catalog.add_column(source_phot)

    return catalog

def filter_catalog(catalog, **kwargs):
    """ Create a new catalog selected from input based on photometry

    Parameters
    ----------
    bright_limit : float
        Fraction of catalog based on brightness that should be retained.
        Default: 1.00 (full catalog)

    max_bright : int
        Maximum number of sources to keep regardless of `bright_limit`
        Default: 100

    min_bright : int
        Minimum number of sources to keep regardless of `bright_limit`
        Default: 20

    colname : string
        Name of column to use for selection/sorting. Default: 'vegamag'

    Returns
    --------
    new_catalog : Table
        New table which only has the sources that meet the selection criteria
    """
    # interpret input pars
    bright_limit = kwargs.get('bright_limit',1.00)
    max_bright = kwargs.get('max_bright',None)
    min_bright = kwargs.get('min_bright',20)
    colname = kwargs.get('colname','vegamag')

    # sort by magnitude
    phot_column = catalog[colname]
    num_sources = len(phot_column)
    sort_indx = np.argsort(phot_column)
    if max_bright is None:
        max_bright = num_sources

    # apply limits, insuring no more than full catalog gets selected
    limit_num = max(int(num_sources*bright_limit), min_bright)
    limit_num = min(max_bright, limit_num, num_sources)

    # Extract sources identified by selection
    new_catalog = catalog[sort_indx[:limit_num]]

    return new_catalog

def countExtn(fimg, extname='SCI'):
    """
    Return the number of 'extname' extensions, defaulting to counting the
    number of SCI extensions.
    """

    closefits = False
    if isinstance(fimg, str):
        fimg = pf.open(fimg)
        closefits = True

    n = 0
    for e in fimg:
        if 'extname' in e.header and e.header['extname'] == extname:
            n += 1

    if closefits:
        fimg.close()

    return n


def build_self_reference(filename,wcslist=None, clean_wcs=False):
    """ This function creates a reference, undistorted WCS that can be used to
    apply a correction to the WCS of the input file.

    PARAMETERS
    ----------
    filename : str
        Filename of image which will be corrected, and which will form the basis
        of the undistorted WCS

    wcslist : list, optional
        List of HSTWCS objects for all SCI extensions in input file.
        If None, will create using `wcsutil.HSTWCS` on each identified SCI
        extension identified in `filename`

    clean_wcs : bool
        Specify whether or not to return the WCS object without any distortion
        information, or any history of the original input image.  This converts
        the output from `utils.output_wcs()` into a pristine `HSTWCS` object.

    Returns
    --------
    customwcs : object
        HSTWCS object which contains the undistorted WCS representing the entire
        field-of-view for the input image

    Syntax
    -------
    This function can be used with the following syntax to apply a shift/rot/scale
    change to the same image:

    >>> import buildref
    >>> from drizzlepac import updatehdr
    >>> filename = "jce501erq_flc.fits"
    >>> wcslin = buildref.build_self_reference(filename)
    >>> updatehdr.updatewcs_with_shift(filename,wcslin,xsh=49.5694, ysh=19.2203, rot = 359.998, scale = 0.9999964)

    """
    if 'sipwcs' in filename:
        sciname = 'sipwcs'
    else:
        sciname = 'sci'
    if wcslist is None:
        numSci = countExtn(filename, extname=sciname.upper())
        wcslist = []
        for extnum in range(numSci):
            extname = (sciname, extnum+1)
            if sciname == 'sci':
                extwcs = wcsutil.HSTWCS(filename, ext=extname)
            else:
                # Working with HDRLET as input and do the best we can...
                extwcs = read_hlet_wcs(filename, ext=extname)
            wcslist.append(extwcs)
    wcslin = utils.output_wcs(wcslist)
    if clean_wcs:
        wcsbase = wcslin.wcs
        customwcs = build_hstwcs(wcsbase.crval[0],wcsbase.crval[1],wcsbase.crpix[0],wcsbase.crpix[1],wcslin._naxis1,wcslin._naxis2,wcslin.pscale,wcslin.orientat)
    else:
        customwcs = wcslin
    return customwcs

def read_hlet_wcs(filename, ext):
    """Insure HSTWCS includes all attributes of a full image WCS.

    For headerlets, the WCS does not contain information about the size of the
    image, as the image array is not present in the headerlet.
    """
    hstwcs = wcsutil.HSTWCS(filename, ext=ext)
    if hstwcs.naxis1 is None:
        hstwcs.naxis1 = int(hstwcs.wcs.crpix[0]*2.) # Assume crpix is center of chip
        hstwcs.naxis2 = int(hstwcs.wcs.crpix[1]*2.)

    return hstwcs

def build_hstwcs(crval1, crval2, crpix1, crpix2, naxis1, naxis2, pscale, orientat):
    """ Create an HSTWCS object for a default instrument without distortion
        based on user provided parameter values.

        .. note :: COPIED from drizzlepac.wcs_functions
    """
    wcsout = wcsutil.HSTWCS()
    wcsout.wcs.crval = np.array([crval1,crval2])
    wcsout.wcs.crpix = np.array([crpix1,crpix2])
    wcsout.naxis1 = naxis1
    wcsout.naxis2 = naxis2
    wcsout.wcs.cd = buildRotMatrix(orientat)*[-1,1]*pscale/3600.0
    # Synchronize updates with PyWCS/WCSLIB objects
    wcsout.wcs.set()
    wcsout.setPscale()
    wcsout.setOrient()
    wcsout.wcs.ctype = ['RA---TAN','DEC--TAN']

    return wcsout

def buildRotMatrix(theta):
    _theta = DEGTORAD(theta)
    _mrot = np.zeros(shape=(2,2), dtype=np.float64)
    _mrot[0] = (np.cos(_theta), np.sin(_theta))
    _mrot[1] = (-np.sin(_theta), np.cos(_theta))

    return _mrot

def DEGTORAD(deg):
    return (deg * np.pi / 180.)

def RADTODEG(rad):
    return (rad * 180. / np.pi)


def find_isolated_source(catalog, columns=None):
    """Find the source in the catalog which is most isolated from all others.

    PARAMETERS
    -----------
    catalog : obj
        list of source positions with (at a minimum) ra,dec or x,y positions
        Catalog could come from photutils.source_properties, for example.

    columns : list, optional
        List of column names for source positions as provided by the catalog.
        If specified (not None), convert input catalog using these column names.
        Default: None.

    Returns
    --------
    index : int
        Index of source from catalog list for source furthest from all others

    nearest_dist : float
        Distance from isolated source to it's nearest neighbor in the catalog

    .. note ::
        This function assumes that the catalog has already been limited to
        only sources which overlap the actual image area.

    """
    if columns:
        cat1 = np.column_stack((catalog[columns[0]], catalog[columns[1]]))
    else:
        cat1 = catalog
    cat_dist = distance_matrix(cat1, cat1)
    y_max = cat_dist.sum(axis=1)
    iso_indx = y_max.argsort()[-1]

    # get index to next nearest neighbor to this isolated source
    next_nearest = cat_dist[iso_indx].argsort()[1]
    # now remember the distance to that source from the isolated source
    iso_dist = cat_dist[iso_indx][next_nearest]

    return iso_indx, iso_dist

def find_central_isolated_source(catalog, wcs, columns=None):
    """Find the most central, isolated source in the catalog

    PARAMETERS
    -----------
    catalog : obj
        list of source positions with (at a minimum) ra,dec or x,y positions
        Catalog could come from photutils.source_properties, for example.

    wcs : obj
        WCS of the field for the catalog.

    columns : list, optional
        List of column names for source positions as provided by the catalog.
        If specified (not None), convert input catalog using these column names.
        Default: None.

    Returns
    --------
    index : int
        Index of source from catalog list for source furthest from all others

    nearest_dist : float
        Distance from isolated source to it's nearest neighbor in the catalog

    .. note ::
        This function assumes that the catalog has already been limited to
        only sources which overlap the actual image area.

    """
    refpix = wcs.wcs.crpix
    if columns:
        cat1 = np.column_stack((catalog[columns[0]], catalog[columns[1]]))
    else:
        cat1 = catalog
    # Compute distance from WCS (field-of-view) center
    central_distance = distance_matrix([refpix], cat1)[0]

    # Compute separations from every member in catalog
    cat_dist = distance_matrix(cat1, cat1)
    y_max = cat_dist.sum(axis=1)
    iso_indx = y_max.argsort()
    max_iso = y_max[iso_indx[0]]
    # Separation weighting: 0.0 - 1.0 for min-max separation
    iso_wht = 1.0/((max_iso*1.01 - y_max)/(max_iso-y_max[iso_indx[-1]]))

    # Determine final sorting based on central_distance weighted by separation
    central_distance *= iso_wht

    # Now, select based on weighted, central distance
    dist_indx = central_distance.argsort()
    iso_indx = dist_indx[0]
    # get index to next nearest neighbor to this isolated source
    next_nearest = dist_indx[1]
    # now remember the distance to that source from the isolated source
    iso_dist = cat_dist[iso_indx][next_nearest]

    return iso_indx, iso_dist

def find_common_isolated_source(ref_catalog, img_catalog, columns=None):
    """Find the most isolated source in the densest part of the source catalog

    PARAMETERS
    -----------
    ref_catalog, img_catalog : obj
        list of source positions with (at a minimum) ra,dec or x,y positions
        Catalog could come from photutils.source_properties, for example.

    columns : list, optional
        List of column names for source positions as provided by the catalog.
        If specified (not None), convert input catalog using these column names.
        Default: None.

    Returns
    --------
    index : int
        Index of source from catalog list for source furthest from all others

    nearest_dist : float
        Distance from isolated source to it's nearest neighbor in the catalog

    .. note ::
        This function assumes that the catalog has already been limited to
        only sources which overlap the actual image area.

    """
    if columns:
        cat1 = np.column_stack((ref_catalog[columns[0]], ref_catalog[columns[1]]))
    else:
        cat1 = ref_catalog
    # find mean position of sources in the image source catalog
    img_mean = (img_catalog[:,0].mean(),img_catalog[:,1].mean())

    # Compute distance from WCS (field-of-view) center
    central_distance = distance_matrix([img_mean], cat1)[0]

    # Compute separations from every member in catalog
    cat_dist = distance_matrix(cat1, cat1)
    y_max = cat_dist.sum(axis=1)
    iso_indx = y_max.argsort()
    max_iso = y_max[iso_indx[0]]
    # Separation weighting: 0.0 - 1.0 for min-max separation
    iso_wht = 1.0/((max_iso*1.01 - y_max)/(max_iso-y_max[iso_indx[-1]]))

    # Determine final sorting based on central_distance weighted by separation
    central_distance *= iso_wht

    # Now, select based on weighted, central distance
    dist_indx = central_distance.argsort()
    iso_indx = dist_indx[0]
    # get index to next nearest neighbor to this isolated source
    next_nearest = dist_indx[1]
    # now remember the distance to that source from the isolated source
    iso_dist = cat_dist[iso_indx][next_nearest]

    return iso_indx, iso_dist


def within_footprint(img, wcs, x, y):
    """Determine whether input x,y fall in the science area of the image.

    Parameters
    -----------
    img : ndarray
        ndarray of image where non-science areas are marked with value of NaN

    wcs : obj
        HSTWCS or WCS object with naxis terms defined

    x,y : arrays
        arrays of x,y positions for sources to be checked

    Returns
    -------
    x,y : arrays
        New arrays which have been trimmed of all sources that fall outside
        the science areas of the image

    """
    # start with limits of WCS shape
    if hasattr(wcs, 'naxis1'):
        naxis1 = wcs.naxis1
        naxis2 = wcs.naxis2
    elif hasattr(wcs, 'pixel_shape'):
        naxis1, naxis2 = wcs.pixel_shape
    else:
        naxis1 = wcs._naxis1
        naxis2 = wcs._naxis2
    maskx = np.bitwise_or(x<0, x>naxis1)
    masky = np.bitwise_or(y<0, y>naxis2)
    mask = ~np.bitwise_or(maskx,masky)
    x = x[mask]
    y = y[mask]

    # Now, confirm that these points fall within actual science area of WCS
    img_mask = create_image_footprint(img, wcs, border=1.0)
    inmask = np.where(img_mask[y.astype(np.int32),x.astype(np.int32)])[0]
    x = x[inmask]
    y = y[inmask]
    return x,y

def create_image_footprint(image, refwcs, border=0.):
    """ Create the footprint of the image in the reference WCS frame

    Parameters
    ----------
    image : HDUList or filename
        Image to extract sources for matching to
        the external astrometric catalog

    refwcs : object
        Reference WCS for coordinate frame of image

    border : float
        Buffer (in arcseconds) around edge of image to exclude astrometric
        sources. Default: 0.

    """
    # Interpret input image to generate initial source catalog and WCS
    if isinstance(image, str):
        image = pf.open(image)
    numSci = countExtn(image, extname='SCI')
    ref_x = refwcs._naxis1
    ref_y = refwcs._naxis2
    # convert border value into pixels
    border_pixels = int(border/refwcs.pscale)

    mask_arr = np.zeros((ref_y,ref_x),dtype=int)

    for chip in range(numSci):
        chip += 1
        # Build arrays of pixel positions for all edges of chip
        chip_y,chip_x = image['sci',chip].data.shape
        chipwcs = wcsutil.HSTWCS(image,ext=('sci',chip))
        xpix = np.arange(chip_x)+1
        ypix = np.arange(chip_y)+1
        edge_x = np.hstack([[1]*chip_y,xpix,[chip_x]*chip_y,xpix])
        edge_y = np.hstack([ypix,[1]*chip_x,ypix,[chip_y]*chip_x])
        edge_ra,edge_dec = chipwcs.all_pix2world(edge_x,edge_y,1)
        edge_x_out,edge_y_out = refwcs.all_world2pix(edge_ra,edge_dec,0)
        edge_x_out = np.clip(edge_x_out.astype(np.int32),0,ref_x-1)
        edge_y_out = np.clip(edge_y_out.astype(np.int32),0,ref_y-1)
        mask_arr[edge_y_out, edge_x_out] = 1

    # Fill in outline of each chip
    mask_arr = ndimage.binary_fill_holes(ndimage.binary_dilation(mask_arr,iterations=2))

    if border > 0.:
        mask_arr = ndimage.binary_erosion(mask_arr, iterations=border_pixels)

    return mask_arr

def compute_hist2d(img_coords, ref_coords, limit=None):
    deltas_full = distance_matrix(img_coords, ref_coords)


def find_isolated_offset(filename, reference,  refwcs = None, refnames=['ra', 'dec'],
                     match_tolerance=5., isolation_limit = ISOLATION_LIMIT,
                     min_match=10, classify=True):
    """Iteratively look for the best cross-match between the catalog and ref.

    Parameters
    ----------
        filename : HDUList or filename
            Image to extract sources for matching to
            the external astrometric catalog

        reference : str or object
            Reference catalog, either as a filename or ``astropy.Table``
            containing astrometrically accurate sky coordinates for astrometric
            standard sources

        refnames : list
            List of table column names for sky coordinates of astrometric
            standard sources from reference catalog

        match_tolerance : float
            Tolerance (in pixels) for recognizing that a source position matches
            an astrometric catalog position.  Larger values allow for lower
            accuracy source positions to be compared to astrometric catalog
            Default: 5 pixels

        isolation_limit : float
            Fractional value (0-1.0) of distance from most isolated source to
            next source to use as limit on overlap for source matching.
            Default: 0.55

        classify : bool
            Specify whether or not to use central_moments classification to
            ignore likely cosmic-rays/bad-pixels when generating the source
            catalog.  Default: True

    Returns
    -------
        best_offset : tuple
            Offset in input image pixels between image source positions and
            astrometric catalog positions that results in largest number of
            matches of astrometric sources with image sources
    """
    # Interpret input image to generate initial source catalog and WCS
    if isinstance(filename, str):
        image = pf.open(filename)
        rootname = filename.split("_")[0]
    else:
        image = filename
        rootname = image[0].header['rootname']

    # check to see whether reference catalog can be found
    if not os.path.exists(reference):
        print("Could not find input reference catalog: {}".format(reference))
        raise FileNotFoundError

    # Extract reference WCS from image
    if refwcs is None:
        refwcs = build_self_reference(image, clean_wcs=True)

    # read in reference catalog
    if isinstance(reference, str):
        refcat = ascii.read(reference)
    else:
        refcat = reference

    ref_ra = refcat[refnames[0]]
    ref_dec = refcat[refnames[1]]

    # Build source catalog for entire image
    iso_cat = build_source_catalog(image, refwcs, classify=classify)

    # Create selective catalog for determining offset
    #   - Limit to brightest 25% of sources (min: 10)
    #max_ref = max(int(0.5*len(refcat)),min_match)
    #iso_cat = filter_catalog(master_cat, bright_limit=0.25,
    #                         min_bright=min_match)

    # Retreive source XY positions in reference frame
    seg_xy = np.column_stack((iso_cat['xcentroid'], iso_cat['ycentroid']))
    seg_xy = seg_xy[~np.isnan(seg_xy[:, 0])]

    # Translate reference catalog positions into input image coordinate frame
    xref, yref = refwcs.all_world2pix(ref_ra, ref_dec, 1)

    # look for the most isolated reference source to serve as a
    # zero-point/anchor
    xref, yref = within_footprint(image, refwcs, xref, yref)
    ref_xy = np.column_stack((xref, yref))
    print("Working with {} astrometric sources for this field".format(len(ref_xy)))

    # write out astrometric reference catalog that was actually used
    ref_ra_img, ref_dec_img = refwcs.all_pix2world(xref, yref, 1)
    ref_tab = Table([ref_ra_img,ref_dec_img, xref, yref],names=['ra','dec', 'x', 'y'])
    ref_tab.write(reference.replace('.cat','_{}.cat'.format(rootname)),
                  format='ascii.fast_commented_header', overwrite=True)

    # Try full 2d hist method
    hist2d_offset = compute_hist2d(ref_xy, seg_xy)

    # Now, determine what catalog to use for computing the offset
    iso_xy = ref_xy
    img_xy = seg_xy

    # assumption: catalog is column_stack((xc,yc))
    #ref_iso, iso_dist = find_isolated_source(iso_xy)
    ref_iso, iso_dist = find_central_isolated_source(iso_xy, refwcs)
    print("ref_iso source at: {}".format(iso_xy[ref_iso]))

    # compute match limit based on distance to neighbor
    match_limit = iso_dist*isolation_limit
    print("ref_iso: {} with limit of {}".format(ref_iso, match_limit))

    # look for nearest matches to isolated reference source
    distxy = distance_matrix([iso_xy[ref_iso]], img_xy)[0]
    maskdist = (distxy < match_limit)

    # identify x,y positions of sources nearest isolated reference
    # close_matches = seg_xy[maskdist]

    # Now, start iterating through all combinations
    delta_refs = img_xy[maskdist] - iso_xy[ref_iso]
    num_matches = []
    for delta in delta_refs:
        mdist = distance_matrix(ref_xy+delta, img_xy)
        num_matches.append(len(np.where(mdist < match_tolerance)[0]))
    max_matches = max(num_matches)
    best_offset = delta_refs[num_matches.index(max_matches)]
    if max_matches < int(len(iso_xy)*0.1):
        best_offset = None
        print("No valid offset found for {}".format(rootname))
    print('best offset {} based on {} cross-matches'.format(best_offset, max_matches))
    return best_offset


def find_best_offset(image, wcs, reference, refnames=['ra', 'dec'],
                     match_tolerance=5., isolation_limit = ISOLATION_LIMIT):
    """Iteratively look for the best cross-match between the catalog and ref
        for a drizzled/resampled image.

    Parameters
    ----------
        image : ndarray
            Image (as ndarray) of image to extract sources for matching to
            the external astrometric catalog.  This image needs to be a
            resampled/undistorted mosaic of the calibrated input image.

        wcs : object
            WCS object which provides translation from sky coordinates to
            image coordinates for input image

        reference : str or object
            Reference catalog, either as a filename or ``astropy.Table``
            containing astrometrically accurate sky coordinates for astrometric
            standard sources

        refnames : list
            List of table column names for sky coordinates of astrometric
            standard sources from reference catalog

        match_tolerance : float
            Tolerance (in pixels) for recognizing that a source position matches
            an astrometric catalog position.  Larger values allow for lower
            accuracy source positions to be compared to astrometric catalog
            Default: 5 pixels

        isolation_limit : float
            Fractional value (0-1.0) of distance from most isolated source to
            next source to use as limit on overlap for source matching.
            Default: 0.55

    Returns
    -------
        best_offset : tuple
            Offset in input image pixels between image source positions and
            astrometric catalog positions that results in largest number of
            matches of astrometric sources with image sources
    """
    # read in reference catalog
    if isinstance(reference, str):
        refcat = ascii.read(reference)
    else:
        refcat = reference
    ref_ra = refcat[refnames[0]]
    ref_dec = refcat[refnames[1]]
    xref, yref = wcs.all_world2pix(ref_ra, ref_dec, 1)

    # look for the most isolated reference source to serve as a
    # zero-point/anchor
    xref, yref = within_footprint(image, wcs, xref, yref)
    xyarr = np.column_stack((xref, yref))
    # assumption: catalog is column_stack((xc,yc))
    ref_iso, iso_dist = find_isolated_source(xyarr)
    print("Found isolated source at position : {}".format(xyarr[ref_iso]))
    print("  with separation from neighbor of: {}".format(iso_dist))

    # compute match limit based on distance to neighbor
    match_limit = iso_dist*isolation_limit

    # find sources in image
    seg_cat, segmap = extract_sources(image)
    seg_tab = seg_cat.to_table()
    seg_xy = np.column_stack((seg_tab['xcentroid'], seg_tab['ycentroid']))
    seg_xy = seg_xy[~np.isnan(seg_xy[:, 0])]

    # look for nearest matches to isolated reference source
    distxy = distance_matrix(xyarr[ref_iso], seg_xy)
    maskdist = (distxy < match_limit)[0]
    # identify x,y positions of sources nearest isolated reference
    # close_matches = seg_xy[maskdist]

    # Now, start iterating through all combinations
    delta_refs = seg_xy[maskdist].value - xyarr[ref_iso]
    num_matches = []
    for delta in delta_refs:
        mdist = distance_matrix(xyarr+delta, seg_xy)
        num_matches.append(len(np.where(mdist < match_tolerance)[0]))
    max_matches = max(num_matches)
    best_offset = delta_refs[num_matches.index(max_matches)]
    if max_matches < int(len(xyarr)*0.1):
        best_offset = None
    return best_offset
