import sys
import traceback
import os
import pytest
from astropy.table import Table

import tweakwcs
from base_test import BaseHLATest
import hlapipeline.utils.astrometric_utils as amutils
from hlapipeline.alignimages import generate_source_catalogs
from ci_watson.artifactory_helpers import get_bigdata


class TestPipeline(BaseHLATest):
    """ Tests which validate whether pipeline code will meet requirements.

        Characeteristics of these tests:
          * A single astrometric catalog was obtained with both GAIA and non-GAIA
            (PanSTARRS?) sources for the entire combined field-of-view using the GSSS
            server.
          * This test only determines whether enough sources
            extracted from the images can match the objects included in
            the astrometric catalog.
          * Tests are included to cover most observation modes
            (ACS/SBC, WFC3/IR,...)
          * The WCS information for the input exposures do not get updated in this test.
          * No mosaic gets generated.

        Success Criteria:
          * Success criteria for source extraction tests
            * For observations with >3 astrometric sources, at least 3 sources
            * For observations with 3 or fewer sources, all astrometric sources
              were identified in the image

    """

    ref_loc = ['truth']

    @pytest.mark.parametrize("input_filenames, truth_file",
                                [('j8ep04lwq','j8ep04lwq_sky.cat')]
                            )
    def test_generate_catalog(self,input_filenames, truth_file):
        """ Verify whether sources from astrometric catalogs can be extracted from images.

        The following datasets are used in these tests:

            * ACS/SBC dataset with >10 sources (in image and in GAIADR2)

        """
        self.input_loc = 'catalog_tests'
        self.curdir = os.getcwd()
        truth_path = [self.input_repo, self.tree, self.input_loc, *self.ref_loc]

        if not isinstance(input_filenames, list):
            input_filenames = [input_filenames]

        try:
            # Make local copies of input files
            local_files = []
            for infile in input_filenames:
                downloaded_files = self.get_input_file(infile, docopy=True)
                local_files.extend(downloaded_files)

            reference_wcs = amutils.build_reference_wcs(local_files)
            input_catalog_dict = generate_source_catalogs([local_files[0]], reference_wcs)
            imcat = input_catalog_dict[local_files[0]]['catalog_table']
            imcat.rename_column('xcentroid', 'x')
            imcat.rename_column('ycentroid', 'y')
            
            # create FITS WCS corrector object
            wcs_corrector = tweakwcs.FITSWCS(reference_wcs)

            # get reference catalog as 'truth' files
            reference_catalog = get_bigdata(*truth_path, truth_file, docopy=True)
            reference_table = Table.read(reference_catalog, format='ascii.fast_commented_header')
            num_expected = len(reference_table)

            # Perform matching
            match = tweakwcs.TPMatch(searchrad=5, separation=0.1, tolerance=5, use2dhist=True)
            ridx, iidx = match(reference_table, imcat, wcs_corrector)
            nmatches = len(ridx)

        except Exception:
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
            sys.exit()

        assert (nmatches > 0.8*num_expected)
