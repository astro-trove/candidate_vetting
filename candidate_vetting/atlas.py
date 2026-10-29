from tom_dataproducts.processors.atlas_processor import AtlasProcessor
from candidate_vetting.public_catalogs.phot_catalogs import ATLAS_Forced_Phot

class ClippedStackedAtlasProcessor(AtlasProcessor):
    """
    Applies the sigma clipping and stacking recommended by David Young:
    https://github.com/thespacedoctor/plot-results-from-atlas-force-photometry-service/blob/main/plot_atlas_fp.py
    """

    def _process_photometry_from_plaintext(self, data_product):
        """
        Processes the photometric data from a plaintext file into a list of dicts. File is read using astropy as
        specified in the below documentation. The file is expected to be a multi-column delimited space delimited
        text file, as produced by the ATLAS forced photometry service at https://fallingstar-data.com/forcedphot
        See https://fallingstar-data.com/forcedphot/resultdesc/ for a description of the output format.

        The header looks like this:
        ###MJD   m   dm  uJy   duJy F err chi/N   RA  Dec   x   y  maj  min   phi  apfit mag5sig Sky   Obs

        :param data_product: ATLAS Photometric DataProduct which will be processed into a list of dicts
        :type data_product: DataProduct

        :returns: python list containing the photometric data from the DataProduct
        :rtype: list
        """
        with open(data_product.data.path) as f:
            filecontent = f.read()
        atlas_fp = ATLAS_Forced_Phot(name="ATLAS FP")
        data = atlas_fp._ATLAS_stack(filecontent)
        photometry = atlas_fp._to_magnitudes(data)
        return photometry
