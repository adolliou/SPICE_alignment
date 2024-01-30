import numpy as np
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
import scipy.interpolate
import cv2
import scipy.ndimage
import astropy.io.fits as Fits
from . import c_correlate
from . import rectify
from astropy.wcs import WCS
import astropy.units as u
from ..plot import plot
import warnings
from ..utils import Util
from astropy.time import Time


class Alignement:

    def __init__(self, large_fov_known_pointing: str, small_fov_to_correct: str, lag_crval1: np.array,
                 lag_crval2: np.array, lag_cdelta1: object, lag_cdelta2: object, lag_crota: object, lag_solar_r: object = None,
                 small_fov_value_min: object = None,
                 parallelism: object = False, use_tqdm: object = False,
                 small_fov_value_max: object = None, counts_cpu_max: object = 40, large_fov_window: object = -1, small_fov_window: object = -1,
                 path_save_figure: object = None, ) -> object:
        """

        @param large_fov_known_pointing: path to the reference file fits (most of the time an imager or a synthetic raster)
        @param small_fov_to_correct: path to the fits file to align. Only the header values will be changed.
        @param lag_crval1: (arcsec) array of header CRVAL1 lags.
        @param lag_crval2: (arcsec) array of header CRVAL2 lags.
        @param lag_cdelta1: (arcsec) array of header CDELT1 lags.
        @param lag_cdelta2: (arcsec) array of header CDELT2 lags.
        @param lag_crota: (deg) array of header CROTA lags. the PC1_1/2 matrixes will be updated accordingly.
        @param lag_solar_r: ([1/Rsun]) set 1.004 most of the time. Only needed if apply carrington transformation
        @param small_fov_value_min: min value (optional)
        @param small_fov_value_max: max value (optional)
        @param parallelism: set true to allow parallelism.
        @param use_tqdm: show advancement bar in terminal.
        @param counts_cpu_max: allow max number of cpu for the parallelism.
        @param large_fov_window: (str or int) HDULIST window for the reference file
        @param small_fov_window: (str or int) HDULIST window for the fits to align
        @param path_save_figure: folder where to save figs following the alignement (optional, will increase computational time)
        """
        self.large_fov_known_pointing = large_fov_known_pointing
        self.small_fov_to_correct = small_fov_to_correct
        self.lag_crval1 = lag_crval1
        self.lag_crval2 = lag_crval2
        self.lag_cdelta1 = lag_cdelta1
        self.lag_cdelta2 = lag_cdelta2

        self.lag_crota = lag_crota
        self.lag_solar_r = lag_solar_r
        self.lonlims = None
        self.latlims = None
        self.shape = None
        self.reference_date = None
        self.parallelism = parallelism
        self.small_fov_window = small_fov_window
        self.large_fov_window = large_fov_window

        self.crval1_ref = None
        self.crval2_ref = None
        self.crota_ref = None
        self.cdelta_ref = None
        self.data_large = None
        self.counts = counts_cpu_max
        self.data_small = None
        self.hdr_small = None
        self.hdr_large = None
        self.method = None
        self.rat_wave = {'171': '171', '193': '195', '211': '195', '131': '171', '304': '304', '335': '304',
                         '94': '171', '174': '171'}
        self.small_fov_value_min = small_fov_value_min
        self.small_fov_value_max = small_fov_value_max
        self.path_save_figure = path_save_figure
        self.use_tqdm = use_tqdm
        self.marker = False

    def _shift_header(self, hdr, **kwargs):
        if 'd_crval1' in kwargs.keys():
            if self.unit_lag == hdr["CUNIT1"]:
                hdr['CRVAL1'] = self.crval1_ref + kwargs["d_crval1"]
            else:
                hdr['CRVAL1'] = u.Quantity(self.crval1_ref, self.unit_lag).to(hdr["CUNIT1"]).value \
                                + u.Quantity(kwargs["d_crval1"], self.unit_lag).to(hdr["CUNIT1"]).value
        if 'd_crval2' in kwargs.keys():
            if self.unit_lag == hdr["CUNIT2"]:
                hdr['CRVAL2'] = self.crval2_ref + kwargs["d_crval2"]
            else:
                hdr['CRVAL2'] = u.Quantity(self.crval2_ref, self.unit_lag).to(hdr["CUNIT2"]).value \
                                + u.Quantity(kwargs["d_crval2"], self.unit_lag).to(hdr["CUNIT2"]).value
        if 'd_cdelta' in kwargs.keys():
            hdr['CDELT1'] = self.cdelta1_ref + kwargs["d_cdelta1"]
            hdr['CDELT2'] = self.cdelta2_ref + kwargs["d_cdelta2"]
        if 'd_crota' in kwargs.keys():
            if 'CROTA' in hdr:
                hdr['CROTA'] = self.crota_ref + kwargs["d_crota"]
                crot = hdr['CROTA']
            elif 'CROTA2' in hdr:
                hdr['CROTA2'] = self.crota_ref + kwargs["d_crota"]
                crot = hdr['CROTA2']
            else:
                if kwargs["d_crota"] != 0.0:
                    crot = np.rad2deg(np.arccos(hdr["PC1_1"]))
                    hdr["CROTA"] = crot
                    # raise NotImplementedError
            if ("PC1_1" in hdr) & (kwargs["d_crota"] != 0):
                theta = u.Quantity(crot, "deg").to("radian").value
                lam = hdr["CDELT2"] / hdr["CDELT1"]
                hdr["PC1_1"] = np.cos(theta)
                hdr["PC2_2"] = np.cos(theta)
                hdr["PC1_2"] = -lam * np.sin(theta)
                hdr["PC2_1"] = (1 / lam) * np.sin(theta)

    def _iteration_step_along_crval2(self, d_crval1, d_cdelta1, d_cdelta2, d_crota, d_solar_r, method: str, ):

        results = np.zeros(len(self.lag_crval2), dtype=np.float64)
        if self.use_tqdm:
            for ii, d_crval2 in enumerate(tqdm(self.lag_crval2, desc='crval1 = %.2f' % (d_crval1))):
                results[ii] = self._step(d_crval2=d_crval2, d_crval1=d_crval1,
                                         d_cdelta1=d_cdelta1, d_cdelta2=d_cdelta2, d_crota=d_crota,
                                         method=method, d_solar_r=d_solar_r,
                                         )
        else:
            for ii, d_crval2 in enumerate(self.lag_crval2):
                results[ii] = self._step(d_crval2=d_crval2, d_crval1=d_crval1,
                                         d_cdelta1=d_cdelta1, d_cdelta2=d_cdelta2, d_crota=d_crota,
                                         method=method, d_solar_r=d_solar_r,
                                         )

        return results

    def _step(self, d_crval2, d_crval1, d_cdelta1, d_cdelta2, d_crota, d_solar_r, method: str, ):
        # print(hdr_small['CRVAL1'])
        # print(self.crval1_ref)

        hdr_small_shft = self.hdr_small.copy()
        self._shift_header(hdr_small_shft, d_crval1=d_crval1, d_crval2=d_crval2,
                           d_cdelta1=d_cdelta1, d_cdelta2=d_cdelta2,
                           d_crota=d_crota)

        data_small = self.function_to_apply(d_solar_r=d_solar_r, data=self.data_small,
                                            hdr=hdr_small_shft)

        condition_1 = np.ones(len(data_small.ravel()), dtype='bool')
        condition_2 = np.ones(len(data_small.ravel()), dtype='bool')

        if self.small_fov_value_min is not None:
            condition_1 = np.array(data_small.ravel() > self.small_fov_value_min, dtype='bool')
        if self.small_fov_value_max is not None:
            condition_2 = np.array(data_small.ravel() < self.small_fov_value_max, dtype='bool')

        if method == 'correlation':

            lag = [0]
            is_nan = np.array((np.isnan(self.data_large.ravel(), dtype='bool')
                               | (np.isnan(data_small.ravel(), dtype='bool'))),
                              dtype='bool')
            return c_correlate.c_correlate(self.data_large.ravel()[(~is_nan) & (condition_1) & (condition_2)],
                                           data_small.ravel()[(~is_nan) & (condition_1) & (condition_2)],
                                           lags=lag)

        elif method == 'residus':
            norm = np.sqrt(self.data_large.ravel())
            diff = (self.data_large.ravel() - data_small.ravel()) / norm
            return np.std(diff[(condition_1) & (condition_2)])
        else:
            raise NotImplementedError

    def align_using_carrington(self, lonlims: list, latlims: list, shape: tuple, reference_date, method='correlation'):

        self.lonlims = lonlims
        self.latlims = latlims
        self.shape = shape
        self.reference_date = reference_date
        self.function_to_apply = self._carrington_transform
        self.method = method
        self.coordinate_frame = "carrington"

        f_large = Fits.open(self.large_fov_known_pointing)
        f_small = Fits.open(self.small_fov_to_correct)

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

    def align_using_helioprojective(self, method='correlation', correct_shift_solar_rotation=False):

        self.lonlims = None
        self.latlims = None
        self.shape = None
        self.reference_date = None
        self.function_to_apply = self._interpolate_on_large_data_grid


        self.method = method
        self.coordinate_frame = "helioprojective"
        f_large = Fits.open(self.large_fov_known_pointing)
        f_small = Fits.open(self.small_fov_to_correct)

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


    def _find_best_header_parameters(self):

        data_large = self.data_large
        naxis1, naxis2 = self._get_naxis(self.hdr_large)

        if (self.hdr_large["CRPIX1"] != (naxis1 + 1) / 2) or (self.hdr_large["CRPIX2"] != (naxis2 + 1) / 2):
            crpix1_old = self.hdr_large["CRPIX1"]
            crpix2_old = self.hdr_large["CRPIX2"]
            self._recenter_crpix_in_header(self.hdr_large)
            crpix1_new = self.hdr_large["CRPIX1"]
            crpix2_new = self.hdr_large["CRPIX2"]
            warnings.warn(f"in hdr_large : CRPIX1 or 2 not in center of FOV."
                          f"\nReplacing CRPIX1 from {crpix1_old} to {crpix1_new}"
                          f"\nReplacing CRPIX2 from {crpix2_old} to {crpix2_new}")

        naxis1, naxis2 = self._get_naxis(self.hdr_small)

        if (self.hdr_small["CRPIX1"] != (naxis1 + 1) / 2) or (self.hdr_small["CRPIX2"] != (naxis2 + 1) / 2):
            crpix1_old = self.hdr_small["CRPIX1"]
            crpix2_old = self.hdr_small["CRPIX2"]
            self._recenter_crpix_in_header(self.hdr_small)
            crpix1_new = self.hdr_small["CRPIX1"]
            crpix2_new = self.hdr_small["CRPIX2"]

            warnings.warn(f"in hdr_small : CRPIX1 or 2 not in center of FOV."
                          f"\nReplacing CRPIX1 from {crpix1_old} to {crpix1_new}"
                          f"\nReplacing CRPIX2 from {crpix2_old} to {crpix2_new}")

        self.crval1_ref = self.hdr_small['CRVAL1']
        self.crval2_ref = self.hdr_small['CRVAL2']
        self.use_crota = True
        if 'CROTA' in self.hdr_small:
            self.crota_ref = self.hdr_small['CROTA']
        elif 'CROTA2' in self.hdr_small:
            self.crota_ref = self.hdr_small['CROTA2']
        else:
            self.crota_ref = np.rad2deg(np.arccos(self.hdr_small['PC1_1']))
            self.hdr_small["CROTA"] = np.rad2deg(np.arccos(self.hdr_small['PC1_1']))
            # self.use_crota = False
        self.cdelta1_ref = self.hdr_small['CDELT1']
        self.cdelta2_ref = self.hdr_small['CDELT2']

        self.unit1 = self.hdr_small["CUNIT1"]
        self.unit2 = self.hdr_small["CUNIT2"]

        if (self.hdr_large["CRPIX1"] !=
            (self.hdr_large["NAXIS1"] + 1) / 2) or (self.hdr_large["CRPIX2"] !=
                                                    (self.hdr_large["NAXIS2"] + 1) / 2):
            raise ValueError("in hdr_large : CRPIX1 or 2 not in center of FOV")

        if "arcsec" in self.unit1:
            warnings.warn("Unit of headers in arcsec : must provide arcsec units for dcrval")
            self.unit_lag = "arcsec"

        elif "deg" in self.unit1:
            warnings.warn("Units of headers in deg: Modyfying inputs units to deg.")
            self.lag_crval1 = Util.ang2pipi(u.Quantity(self.lag_crval1, "arcsec")).to("deg").value
            self.lag_crval2 = Util.ang2pipi(u.Quantity(self.lag_crval2, "arcsec")).to("deg").value
            self.lag_cdelta1 = Util.ang2pipi(u.Quantity(self.lag_cdelta1, "arcsec")).to("deg").value
            self.lag_cdelta2 = Util.ang2pipi(u.Quantity(self.lag_cdelta2, "arcsec")).to("deg").value
            self.unit_lag = "deg"
        if self.lag_solar_r is None:
            self.lag_solar_r = np.array([1.004])
        results = np.empty((len(self.lag_crval1), len(self.lag_crval2), len(self.lag_cdelta1), len(self.lag_cdelta2),
                            len(self.lag_crota), len(self.lag_solar_r)), dtype=np.float64)
        if self.parallelism:
            # if (len(self.lag_cdelta1) > 1) or (len(self.lag_crota) > 1) or (len(self.lag_cdelta2) > 1):
            #     raise NotImplementedError
            # else:

            for kk, d_solar_r in enumerate(self.lag_solar_r):
                if self.coordinate_frame == "carrington":
                    self.data_large = self.function_to_apply(d_solar_r=d_solar_r, data=data_large,
                                                             hdr=self.hdr_large)
                elif self.coordinate_frame == "helioprojective":
                    self.data_large = self._create_submap_of_large_data(data_large=data_large)
                for ii, d_cdelta1 in enumerate(self.lag_cdelta1):
                    for ll, d_cdelta2 in enumerate(self.lag_cdelta2):
                        for jj, d_crota in enumerate(self.lag_crota):

                            count = self.counts
                            if self.counts > mp.cpu_count(): count = mp.cpu_count() - 10
                            pool = mp.Pool(count)
                            results[:, :, ii, ll, jj, kk] = pool.map(partial(self._iteration_step_along_crval2,
                                                                             d_cdelta1=d_cdelta1,
                                                                             d_cdelta2=d_cdelta2,
                                                                             d_crota=d_crota,
                                                                             method=self.method,
                                                                             d_solar_r=d_solar_r, ),
                                                                     self.lag_crval1)
                            pool.close()

        else:
            for hh, d_solar_r in enumerate(self.lag_solar_r):
                if self.coordinate_frame == "carrington":
                    self.data_large = self.function_to_apply(d_solar_r=d_solar_r, data=data_large,
                                                             hdr=self.hdr_large)
                elif self.coordinate_frame == "helioprojective":
                    self.data_large = self._create_submap_of_large_data(data_large=data_large)

                for ii, d_crval1 in enumerate(self.lag_crval1):
                    for jj, d_crval2 in enumerate(tqdm(self.lag_crval2)):
                        for kk, d_cdelta1 in enumerate(self.lag_cdelta1):
                            for mm, d_cdelta2 in enumerate(self.lag_cdelta2):
                                for ll, d_crota in enumerate(self.lag_crota):
                                    results[ii, jj, kk, mm, ll, hh] = self._step(d_crval2=d_crval2, d_crval1=d_crval1,
                                                                                 d_cdelta1=d_cdelta1,
                                                                                 d_cdelta2=d_cdelta2,
                                                                                 d_crota=d_crota,
                                                                                 method=self.method, d_solar_r=d_solar_r
                                                                                 )

        return results

    def _recenter_crpix_in_header(self, hdr):
        w = WCS(hdr)
        naxis1, naxis2 = self._get_naxis(hdr)
        x_mid = (naxis1 - 1) / 2
        y_mid = (naxis2 - 1) / 2
        lon_mid, lat_mid = w.pixel_to_world(np.array([x_mid]), np.array([y_mid]))
        lon_mid = lon_mid[0].to(hdr["CUNIT1"]).value
        lat_mid = lat_mid[0].to(hdr["CUNIT2"]).value
        hdr["CRVAL1"] = lon_mid
        hdr["CRVAL2"] = lat_mid
        hdr["CRPIX1"] = (naxis1 + 1) / 2
        hdr["CRPIX2"] = (naxis2 + 1) / 2

        return hdr

    def _carrington_transform(self, d_solar_r, data, hdr):

        spherical = rectify.CarringtonTransform(hdr, radius_correction=d_solar_r,
                                                reference_date=self.reference_date,
                                                rate_wave=self.rat_wave[
                                                    '%i' % (self.hdr_large['WAVELNTH'])])
        spherizer = rectify.Rectifier(spherical)
        image = spherizer(data, self.shape, self.lonlims, self.latlims, opencv=False, order=1, fill=-32762)
        image = np.where(image == -32762, np.nan, image)
        if Fits.HeaderDiff(hdr, self.hdr_large).identical:
            if self.path_save_figure is not None:
                plot.PlotFunctions.plot_fov(image, show=False,
                                            path_save='%s/image_large.pdf' % (self.path_save_figure))
                spherical = rectify.CarringtonTransform(self.hdr_small, radius_correction=d_solar_r,
                                                        reference_date=self.reference_date,
                                                        rate_wave=self.rat_wave[
                                                            '%i' % (self.hdr_large['WAVELNTH'])])
                spherizer = rectify.Rectifier(spherical)

                image_small = spherizer(self.data_small, self.shape, self.lonlims, self.latlims, opencv=False,
                                        order=1, fill=-32762)
                image_small = np.where(image_small == -32762, np.nan, image_small)

                plot.PlotFunctions.plot_fov(image_small, show=False,
                                            path_save='%s/image_small.pdf' % (self.path_save_figure))

        return image

    def _create_submap_of_large_data(self, data_large):
        if self.path_save_figure is not None:
            plot.PlotFunctions.simple_plot(self.hdr_large, self.data_large, show=False,
                                           path_save='%s/large_fov_before_cut.pdf' % (self.path_save_figure))

        hdr_cut = self.hdr_small.copy()

        longitude_cut, latitude_cut, dsun_obs_cut = Util.extract_EUI_coordinates(hdr_cut)
        w_xy_large = WCS(self.hdr_large.copy())
        x_cut, y_cut = w_xy_large.world_to_pixel(longitude_cut, latitude_cut)
        image_large_cut = Util.CommonUtil.interpol2d(np.array(data_large, dtype=np.float64), x=x_cut, y=y_cut, order=1, fill=-32768)
        image_large_cut[image_large_cut == -32768] = np.nan
        self.hdr_large = hdr_cut.copy()

        w_xy_small = WCS(self.hdr_small.copy())
        x_cut, y_cut = w_xy_small.world_to_pixel(longitude_cut, latitude_cut)
        image_small_cut = Util.CommonUtil.interpol2d(np.array(self.data_small.copy(), dtype=np.float64), x=x_cut, y=y_cut,
                                     order=1, fill=-32768)
        image_small_cut[image_small_cut == -32768] = np.nan

        self.data_small = image_small_cut
        self.hdr_small = hdr_cut.copy()
        levels = [0.15 * np.nanmax(self.data_small)]
        if self.path_save_figure is not None:
            date_small = self.hdr_small["DATE-AVG"]
            date_small = date_small.replace(":", "_")
            plot.PlotFunctions.simple_plot(self.hdr_large, image_large_cut, show=False,
                                           path_save='%s/large_fov_%s.pdf' % (self.path_save_figure, date_small))
            plot.PlotFunctions.simple_plot(self.hdr_small, self.data_small, show=False,
                                           path_save='%s/small_fov_%s.pdf' % (self.path_save_figure, date_small))
            plot.PlotFunctions.contour_plot(self.hdr_large, image_large_cut, self.hdr_small, self.data_small,
                                            show=False, path_save='%s/compare_plot_%s.pdf' % (self.path_save_figure,
                                                                                              date_small),
                                            levels=levels)
        self.step_figure = False
        return np.array(image_large_cut)

    def _interpolate_on_large_data_grid(self, d_solar_r, data, hdr):

        w_xy_small = WCS(hdr)
        longitude_large, latitude_large, dsun_obs_large = Util.extract_EUI_coordinates(self.hdr_large)
        x_large, y_large = w_xy_small.world_to_pixel(longitude_large, latitude_large)
        image_small_shft = Util.CommonUtil.interpol2d(np.array(data, dtype=np.float64), x=x_large, y=y_large, order=1, fill=-32768)
        image_small_shft = np.where(image_small_shft == -32768, np.nan, image_small_shft)


        return image_small_shft

    @staticmethod
    def _get_naxis(hdr):
        if "ZNAXIS1" in hdr:
            naxis1 = hdr["ZNAXIS1"]
            naxis2 = hdr["ZNAXIS2"]
        else:
            naxis1 = hdr["NAXIS1"]
            naxis2 = hdr["NAXIS2"]
        return naxis1, naxis2


class Util:
    @staticmethod
    def ang2pipi(ang):
        """ put angle between ]-180, +180] deg """
        pi = u.Quantity(180, 'deg')
        return - ((- ang + pi) % (2 * pi) - pi)

    @staticmethod
    def extract_EUI_coordinates(hdr):
        w = WCS(hdr)
        idx_lon = np.where(np.array(w.wcs.ctype, dtype="str") == "HPLN-TAN")[0][0]
        idx_lat = np.where(np.array(w.wcs.ctype, dtype="str") == "HPLT-TAN")[0][0]
        x, y = np.meshgrid(np.arange(w.pixel_shape[idx_lon]), np.arange(w.pixel_shape[idx_lat]), )  # t dépend de x,
        # should reproject on a new coordinate grid first : suppose slits at the same time :
        longitude, latitude = w.pixel_to_world(x, y)
        if "DSUN_OBS" in w.to_header():
            dsun_obs_large = w.to_header()["DSUN_OBS"]
        else:
            dsun_obs_large = None
        return Util.ang2pipi(longitude), Util.ang2pipi(latitude), dsun_obs_large