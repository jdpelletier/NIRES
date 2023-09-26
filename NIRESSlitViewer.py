import os, time, sys, threading, math
from os import listdir
from os.path import abspath, isfile, join
from pathlib import Path
import datetime
import csv
import copy

import numpy as np
from astropy.io import fits
from astropy import wcs
import astropy.units as u
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.modeling import models, fitting
import PIL.Image as PILimage

from ginga import Bindings
from ginga.misc import log
from ginga.qtw.QtHelp import QtGui, QtCore
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
from ginga.util import iqcalc
from ginga.util.loader import load_data
from ginga.AstroImage import AstroImage

import ktl

class ScannerSignals(QtCore.QObject):
    load = QtCore.Signal(object)

class Scanner(QtCore.QRunnable):
    '''
    Scanner thread
    '''
    def __init__(self, fn, *args, **kwargs):
        super(Scanner, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = ScannerSignals()
        self.kwargs['file_callback'] = self.signals.load

        # Add the callback to our kwargs
    @QtCore.Slot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        self.fn(*self.args, **self.kwargs)

class UpdateControlWindowSignals(QtCore.QObject):
    load = QtCore.Signal()

class UpdateControlWindow(QtCore.QRunnable):
    '''
    Control Window thread
    '''
    def __init__(self, fn, *args, **kwargs):
        super(UpdateControlWindow, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = UpdateControlWindowSignals()
        self.kwargs['file_callback'] = self.signals.load

        # Add the callback to our kwargs
    @QtCore.Slot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        self.fn(*self.args, **self.kwargs)


class FitsViewer(QtGui.QMainWindow):

    def __init__(self, logger):
        super(FitsViewer, self).__init__()
        self.logger = logger

        self.cachedFiles = None
        #KTL stuff
        #Cache KTL keywords
        self.slit_filename = ktl.cache('nids', 'FILENAME')
        # self.trickxpos = ktl.cache('tds', 'TRKRO1X')
        # self.trickypos = ktl.cache('tds', 'TRKRO1Y')
        # self.trickxsize = ktl.cache('tds', 'TRKRO1SX')
        # self.trickysize = ktl.cache('tds', 'TRKRO1SY')
        # self.progress = ktl.cache('tds', 'progress')
        # self.roipixels = ktl.cache('ao', 'TRKRO1PX')
        # self.roipixels.monitor()
        # self.tkenrup = ktl.cache('ao', 'tkenrup')
        # self.tkcrxs = ktl.cache('ao','tkcrxs')
        # self.tkcrys = ktl.cache('ao','tkcrys')
        # self.tkcrevxp = ktl.cache('ao','tkcrevxp')
        # self.tkcrevyp = ktl.cache('ao','tkcrevyp')
        # self.tkcrevxo = ktl.cache('ao','tkcrevxo')
        # self.tkcrevyo = ktl.cache('ao','tkcrevyo')
        # self.tkcxim = ktl.cache('ao','tkcxim')
        # self.tkcyim = ktl.cache('ao','tkcyim')

        self.rawfile = ''
        self.mode = ''

        self.img = AstroImage()

        self.threadpool = QtCore.QThreadPool()

        self.iqcalc = iqcalc.IQCalc(self.logger)

        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params('zscale')
        fi.enable_autozoom('on')
        # fi.set_callback('drag-drop', self.drop_file)
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        self.fitsimage = fi

        # enable some user interaction
        self.bd = fi.get_bindings()
        self.bd.enable_all(True)
        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setObjectName("vbox")
        viewer_hbox = QtGui.QHBoxLayout()
        viewer_hbox.setObjectName("viewer_hbox")
        w = fi.get_widget()
        w.setMinimumSize(QtCore.QSize(512, 700))
        viewer_hbox.addWidget(w)
        hw = QtGui.QWidget()
        hw.setLayout(viewer_hbox)
        vbox.addWidget(hw)
        readout_hbox = QtGui.QHBoxLayout()
        readout_hbox.setObjectName("readout_hbox")
        self.readout = QtGui.QLabel("X:                 Y:                    Value:")
        self.readout.setObjectName("readout")
        self.readout.setMinimumSize(QtCore.QSize(350, 0))
        readout_hbox.addWidget(self.readout)
        self.wcut = QtGui.QComboBox()
        for name in fi.get_autocut_methods():
            self.wcut.addItem(name)
        self.wcut.currentIndexChanged.connect(self.cut_change)
        readout_hbox.addWidget(self.wcut)
        self.wcolor = QtGui.QComboBox()
        for name in fi.get_color_algorithms():
            self.wcolor.addItem(name)
        self.wcolor.currentIndexChanged.connect(self.color_change)
        readout_hbox.addWidget(self.wcolor)
        hw = QtGui.QWidget()
        hw.setLayout(readout_hbox)
        vbox.addWidget(hw)
        file_hbox = QtGui.QHBoxLayout()
        file_hbox.setObjectName("file_hbox")
        self.file_info = QtGui.QLabel("file: ")
        self.file_info.setObjectName("file_info")
        file_hbox.addWidget(self.file_info)
        self.box_readout = QtGui.QLabel("Amplitude:                  FWHM: ")
        self.box_readout.setMinimumSize(QtCore.QSize(350, 0))
        self.box_readout.setObjectName("box_readout")
        file_hbox.addWidget(self.box_readout)
        hw = QtGui.QWidget()
        hw.setLayout(file_hbox)
        vbox.addWidget(hw)
        buttons_hbox = QtGui.QHBoxLayout()
        buttons_hbox.setObjectName("buttons_hbox")
        buttons_vbox_left = QtGui.QVBoxLayout()
        buttons_vbox_left.setContentsMargins(QtCore.QMargins(0, 0, 10, 0))
        buttons_vbox_left.setObjectName("buttons_vbox_left")
        self.wcenter = QtGui.QPushButton("Center")
        self.wcenter.setObjectName("wcenter")
        # self.wcenter.clicked.connect(self.center)
        buttons_vbox_left.addWidget(self.wcenter)
        self.wrecenter = QtGui.QPushButton("Re-center")
        self.wrecenter.setObjectName("wrecenter")
        # self.wrecenter.clicked.connect(self.recenter)
        buttons_vbox_left.addWidget(self.wrecenter)
        hw = QtGui.QWidget()
        hw.setLayout(buttons_vbox_left)
        buttons_hbox.addWidget(hw)
        buttons_vbox_cent = QtGui.QVBoxLayout()
        buttons_vbox_cent.setObjectName("buttons_vbox_cent")
        self.wstartscan = QtGui.QPushButton("Start Scan")
        self.wstartscan.setObjectName("wstartscan")
        self.wcenter.clicked.connect(self.start_scan)
        buttons_vbox_left.addWidget(self.wstartscan)
        self.wstopscan = QtGui.QPushButton("Stop Scan")
        self.wstopscan.setObjectName("wstopscan")
        self.wstopscan.clicked.connect(self.stop_scan)
        buttons_vbox_left.addWidget(self.wstopscan)
        hw = QtGui.QWidget()
        hw.setLayout(buttons_vbox_cent)
        buttons_hbox.addWidget(hw)
        buttons_vbox_right = QtGui.QVBoxLayout()
        buttons_vbox_right.setObjectName("buttons_vbox_right")
        self.wopen = QtGui.QPushButton("Open File")
        self.wopen.setObjectName("wopen")
        self.wopen.clicked.connect(self.open_file)
        buttons_vbox_right.addWidget(self.wopen)
        self.wquit = QtGui.QPushButton("Quit")
        self.wquit.setObjectName("wquit")
        self.wquit.clicked.connect(self.quit)
        buttons_vbox_right.addWidget(self.wquit)
        hw = QtGui.QWidget()
        hw.setLayout(buttons_vbox_right)
        buttons_hbox.addWidget(hw)
        hw = QtGui.QWidget()
        hw.setLayout(buttons_hbox)
        vbox.addWidget(hw)

        vw = QtGui.QWidget()
        self.setCentralWidget(vw)
        vw.setLayout(vbox)

        fi.set_callback('cursor-changed', self.motion_cb)
        fi.add_callback('cursor-down', self.btndown)

        # self.recdc, self.compdc = self.add_canvas()
        self.recdc = self.add_canvas()
        self.picktag = "pick-box"

        self.start_updating()

    def start_updating(self):
        self.updating = True
        updater = UpdateControlWindow(self.update)
        updater.signals.load.connect(self.update_gui)
        self.threadpool.start(updater)

    def update(self, file_callback):
        while self.updating:
            file_callback.emit()
            time.sleep(1)

    def stop_updating(self):
        self.updating = False

    def update_gui(self):
        name = self.slit_filename.read()
        self.file_info.setText(f"Name: {name}")

    def add_canvas(self, tag=None):
        # add a canvas to the view
        my_canvas = self.fitsimage.get_canvas()
        RecCanvas = my_canvas.get_draw_class('rectangle')
        # CompCanvas = my_canvas.get_draw_class('compass')
        # return RecCanvas, CompCanvas
        return RecCanvas


    def cut_change(self):
        self.fitsimage.set_autocut_params(self.wcut.currentText())

    def color_change(self):
        self.fitsimage.set_color_algorithm(self.wcolor.currentText())

    def motion_cb(self, viewer, button, data_x, data_y):

        # Get the value under the data coordinates
        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            value = viewer.get_data(int(data_x + 0.5), int(data_y + 0.5))

        except Exception:
            value = None

        fits_x, fits_y = data_x, data_y

        # Calculate WCS RA
        try:
            # NOTE: image function operates on DATA space coords
            image = viewer.get_image()
            if image is None:
                # No image loaded
                return
            ra_txt, dec_txt = image.pixtoradec(fits_x, fits_y,
                                               format='str', coords='fits')
        except Exception as e:
            self.logger.warning("Bad coordinate conversion: %s" % (
                str(e)))
            ra_txt = 'BAD WCS'
            dec_txt = 'BAD WCS'
        if (fits_x > 2048 or fits_x <0) or (fits_y > 2048 or fits_y <0):
            text = "X: Y:  Value:"
            self.readout.setText(text)
        else:
            text = f"X: {int(fits_x)} Y: {int(fits_y)}  Value: {value}"
            self.readout.setText(text)

    def quit(self, *args):
        self.logger.info("Attempting to shut down the application...")
        self.stop_scan()
        time.sleep(2)
        self.threadpool = False
        QtGui.QApplication.instance().quit()

    ##Full frame stuff
    def start_scan(self):
        self.scanning = True
        # hdu = fits.PrimaryHDU()
        # try:
        #     hdu.writeto('procImage.fits')
        # except OSError:
        #     os.remove('procImage.fits')
        #     hdu.writeto('procImage.fits')
        self.cachedFiles = self.walkDirectory()
        print("scan started...")
        scanner = Scanner(self.scan)
        scanner.signals.load.connect(self.load_file)
        self.threadpool.start(scanner)

    def stop_scan(self):
        self.scanning = False
        print('Scanning stopped.')

    def load_file(self, filepath):
        image = load_data(filepath, logger=self.logger)
        self.fitsimage.set_image(image)
        # self.setWindowTitle(filepath)
        try:
            self.fitsimage.get_canvas().get_object_by_tag(self.picktag)
            self.fitsimage.get_canvas().delete_object_by_tag(self.picktag)
        except KeyError:
            pass
        # width, height = image.get_size()
        # data_x, data_y = width / 2.0, height / 2.0
        # # x, y = self.fitsimage.get_canvas_xy(data_x, data_y)
        # radius = float(max(width, height)) / 20
        # self.fitsimage.get_canvas().add(self.compdc(data_x, data_y, radius, color='skyblue',
        #                                fontsize=8))
        name = 'tmp'
        text = f"Image: {name}"
        self.file_info.setText(text)

    def open_file(self):
        res = QtGui.QFileDialog.getOpenFileName(self, "Open FITS file",
                                                str(self.nightpath()))
        print(res)
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.load_file(fileName)

    # def load_sky(self):
    #     res = QtGui.QFileDialog.getOpenFileName(self, "Open Sky file",
    #                                             str(self.nightpath()))
    #     if isinstance(res, tuple):
    #         fileName = res[0]
    #     else:
    #         fileName = str(res)
    #     if len(fileName) != 0:
    #         self.subtract_sky(fileName)

    # def subtract_sky(self, filename):
    #     skyname, skyheader, skyfitsData, skyfilter = self.addWcs(filename)
    #     name, header, fitsData, filter = self.addWcs(self.rawfile)
    #     with_sky = fitsData - skyfitsData
    #     mask = fits.getdata('/kroot/rel/ao/qfix/data/Trick/BadPix_1014Hz.fits', ext=0)
    #     text = f"Sky: {skyname}"
    #     self.sky_info.setText(text)
    #     self.load_file(self.writeFits(header, np.multiply(with_sky, mask)))


    ##Start of image find and processing code

    def scan(self, file_callback):
        while self.scanning:
            hasNewFiles, files, self.cachedFiles = self.updateFileCache(self.cachedFiles)
            if hasNewFiles and ('.fits' in files[0] or '.fits.gz' in files[0]):
                print("New Image Detected!")
                filen = files[0]
                self.waitForFileToBeUnlocked(filen, 1)
                file_callback.emit(filen)
            time.sleep(1)

    def walkDirectory(self):
        directory = self.nightpath()
        return [abspath(join(directory, f)) for f in listdir(directory) if isfile(join(directory, f))]

    def updateFileCache(self, cachedFiles):
        updatedFileList = self.walkDirectory()
        filtered = [i for i in updatedFileList if not i in cachedFiles]
        cachedFiles = updatedFileList
        return len(filtered) > 0, filtered, cachedFiles

    def fileIsCurrentlyLocked(self, filepath):
        locked = None
        hdulist = None
        file_object = None
        if os.path.exists(filepath):
            try:
                print("Trying to open %s." % filepath)
                #time.sleep(15) #place holder if OSError catch doesn't work
                hdulist = fits.open(filepath)

                file_object = np.sum([1 for hdu in hdulist if type(hdu) in
                        	[fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU]
                        	and hdu.data is not None])
                if file_object:
                    locked = False

            except TypeError:
                locked = True

            except OSError:
                locked = True

            finally:
                if file_object:
                    hdulist.close()

        else:
            print("%s not found." % filepath)

        return locked


    #    Checks if the files are ready.
    #    For a file to be ready it must exist and can be opened in append mode.

    def waitForFileToBeUnlocked(self, filename, wait_time):
        # if the file doesn't exist, wait wait_time seconds and try again until it's found
        while not os.path.exists(filename):
            print("%s hasn't arrived. Waiting %s seconds." % (filename, wait_time))
            time.sleep(wait_time)

        # if the file exists but locked, wait wait_time seconds and check
        # again until it's no longer locked by another process
        while self.fileIsCurrentlyLocked(filename):
            # print(self.progress.read())
            time.sleep(wait_time)

    def nightpath(self):
        nightly = Path('/s/sdata1500/nires1/2023sep25')
        # date = datetime.datetime.utcnow()
        # year, month, day = str(date.strftime("%y")), str(date.strftime("%m")), str(date.strftime("%d"))
        # nightly = nightly / year / month / day / 'Trick'
        return nightly

    # def processData(self, filename):
    #     self.rawfile = filename
    #     name, header, fitsData, filter = self.addWcs(filename)
        # mask = fits.getdata('/kroot/rel/ao/qfix/data/Trick/BadPix_1014Hz.fits', ext=0)
        # if filter == 'H':
        #     background = fits.getdata('/kroot/rel/ao/qfix/data/Trick/sky_H.fits')
        #     self.sky_info.setText('Sky: sky_H.fits')
        # else:
        #     background = fits.getdata('/kroot/rel/ao/qfix/data/Trick/sky_Ks.fits')
        #     self.sky_info.setText('sky_Ks.fits')
        # subtracted_data = fitsData-background
        # self.load_file(self.writeFits(header, np.multiply(subtracted_data, mask)))
    #     self.load_file(self.writeFits(header, fitsData))
    #     text = f"Image: {name}"
    #     self.image_info.setText(text)
    #     text = f"Filter: {filter}"
    #     self.filt_info.setText(text)
    #     self.wsky.setEnabled(True)

    # def addWcs(self, filen):
    #     w = wcs.WCS(naxis=2)
    #     fitsData = fits.getdata(filen, ext=0)
    #     header = fits.getheader(filen)
    #     ht, wd = fitsData.shape[:2]
    #     y = ht//2
    #     x = wd//2
    #     name = header['DATAFILE']
    #     ra = float(header['RA'])
    #     dec = float(header['DEC'])
    #     rot = float(header['ROTPOSN'])
    #     filter = header['TRFWNAME']
    #     w.wcs.crpix = [y, x]
    #     w.wcs.cdelt = np.array([-0.05, 0.05])
    #     w.wcs.crota = np.array([0.05, rot])
    #     w.wcs.crval = [ra, dec]
    #     w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    #     pixcrd = np.array([[0, 0], [24, 38], [45, 98]], dtype=np.float64)
    #     world = w.wcs_pix2world(pixcrd, 0)
    #     # Convert the same coordinates back to pixel coordinates.
    #     pixcrd2 = w.wcs_world2pix(world, 0)
    #     # These should be the same as the original pixel coordinates, modulo
    #     # some floating-point error.
    #     assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6
    #     # Now, write out the WCS object as a FITS header
    #     header = w.to_header()
    #     return name, header, fitsData, filter

    # def writeFits(self, headerinfo, image_data):
    #     hdu = fits.PrimaryHDU(header=headerinfo, data=image_data)
    #     filename = 'procImage.fits'
    #     try:
    #         hdu.writeto(filename)
    #     except OSError:
    #         os.remove(filename)
    #         hdu.writeto(filename)
    #     return filename

    ##Find star stuff
    def cutdetail(self, image, shape_obj):
        view, mask = image.get_shape_view(shape_obj)

        data = image._slice(view)

        y1, y2 = view[0].start, view[0].stop
        x1, x2 = view[1].start, view[1].stop

        # mask non-containing members
        mdata = np.ma.array(data, mask=np.logical_not(mask))

        return x1, y1, x2, y2, mdata

    def findstar(self):
        image = self.fitsimage.get_image()
        obj = self.pickbox
        shape = obj
        x1, y1, x2, y2, data = self.cutdetail(image, shape)
        ht, wd = data.shape[:2]
        xc, yc = wd // 2, ht // 2
        radius = min(xc, yc)
        peaks = [(xc, yc)]
        peaks = self.iqcalc.find_bright_peaks(data,
                                              threshold=None,
                                              radius=radius)

        xc, yc = peaks[0]
        xc += 1
        yc += 1
        return int(xc), int(yc), data

    def fitstars(self, y_line):
        x = np.linspace(-30, 30, 60)
        model_gauss = models.Gaussian1D()
        fitter_gauss = fitting.LevMarLSQFitter()
        # gx = fitter_gauss(model_gauss, x, x_line)
        gy = fitter_gauss(model_gauss, x, y_line)
        # amplitude = (gx.amplitude.value+gy.amplitude.value)/2
        amplitude = gy.amplitude.value
        # fwhm = ((gx.stddev.value+gy.stddev.value)/2)*0.118 #2.355*PixelScale
        # fwhm = (gy.stddev.value)*0.118 #2.355*PixelScale
        fwhm = (gy.stddev.value)*2.355 #pixels instead of arcseconds
        return amplitude, fwhm

    def pickstar(self, viewer):
        try:
            self.fitsimage.get_canvas().get_object_by_tag(self.picktag)
            self.fitsimage.get_canvas().delete_object_by_tag(self.picktag)
            self.pickbox = self.recdc(self.xclick-30, self.yclick-30, self.xclick+30, self.yclick+30, color='red')
            self.fitsimage.get_canvas().add(self.pickbox, tag=self.picktag, redraw=True)
        except KeyError:
            self.pickbox = self.recdc(self.xclick-30, self.yclick-30, self.xclick+30, self.yclick+30, color='red')
            self.fitsimage.get_canvas().add(self.pickbox, tag=self.picktag, redraw=True)
        image = self.fitsimage.get_image()
        try:
            xc, yc, data = self.findstar()
            # x_line = data[40-yc, 0:40] doesn't work well for some reason
            y_line = data[0:60, xc]
            # amplitude, fwhm = self.fitstars(x_line, y_line)
            amplitude, fwhm = self.fitstars(y_line)
            text = f"Amplitude: {amplitude:.2f} FWHM: {fwhm:.2f}"
            self.box_readout.setText(text)
        except IndexError:
            text = "Amplitude: N/A FWHM: N/A"
            self.box_readout.setText(text)


    def btndown(self, canvas, event, data_x, data_y):
        # self.fitsimage.set_pan(data_x, data_y)
        self.xclick = data_x
        self.yclick = data_y
        self.pickstar(self.fitsimage)

def main():
    ##Write dummy file so walkDirectory caches it in the beginning

    app = QtGui.QApplication([])

    # ginga needs a logger.
    # If you don't want to log anything you can create a null logger by
    # using null=True in this call instead of log_stderr=True
    logger = log.get_logger("NIRESSlitViewer", log_stderr=True, level=40)

    w = FitsViewer(logger)
    w.show()
    app.setActiveWindow(w)
    w.raise_()
    w.activateWindow()
    # w.start_scan()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
