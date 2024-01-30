from astropy.io import fits
import numpy as np

@


def test_alignement(self, method='correlation', correct_shift_solar_rotation=False):

        self.lonlims = None
        self.latlims = None
        self.shape = None
        self.reference_date = None
        self.function_to_apply = self._interpolate_on_large_data_grid


        self.method = method
        self.coordinate_frame = "helioprojective"
        f_large = fits.open(self.large_fov_known_pointing)
        f_small = fits.open(self.small_fov_to_correct)

        self.data_large = np.array(f_large[self.large_fov_window].data.copy(), dtype=np.float64)
        self.hdr_large = f_large[self.large_fov_window].header.copy()
        self._recenter_crpix_in_header(self.hdr_large)

        self.hdr_small = f_small[self.small_fov_window].header.copy()
        self._recenter_crpix_in_header(self.hdr_small)
        self.data_small = np.array(f_small[self.small_fov_window].data.copy(), dtype=np.float64)
        f_large.close()
        f_small.close()

        results = self._find_best_header_parameters()

        return results