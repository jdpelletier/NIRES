import os, time, sys, threading, math
# from os import listdir
# from os.path import abspath, isfile, join
from pathlib import Path
# import datetime
# import csv
# import copy

import numpy as np
from astropy.io import fits
# from astropy import wcs
import astropy.units as u
# from astropy.stats import gaussian_sigma_to_fwhm
from astropy.modeling import models, fitting
# import PIL.Image as PILimage

from ginga import Bindings, cmap
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
        # self.slit_filename = ktl.cache('nids', 'FILENAME')
        # self.slit_filename.monitor()
        # self.slit_lastfile = ktl.cache('nids', 'LASTFILE')
        # self.slit_sdiff = ktl.cache('nids', 'LASTFILE')
        # self.go = ktl.cache('nids', 'GO')
        # self.go.monitor()
        # self.test = ktl.cache('nids', 'test')
        # self.test.monitor()
        self.display2 = ktl.cache('nids', 'display2')
        self.display2.monitor()
        self.dispname2 = ktl.cache('nids', 'dispname2')
        self.dispname2.monitor()
        # self.tempsky2 = ktl.cache('nids', 'TEMPSKY2')
        # self.tempsky2.monitor()
        self.dcs = ktl.Service('dcs')

        self.rawfile = ''
        self.mode = ''

        self.img = AstroImage()

        self.threadpool = QtCore.QThreadPool()

        self.iqcalc = iqcalc.IQCalc(self.logger)

        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params('zscale')
        fi.enable_autozoom('off')
        fi.enable_autocenter('off')
        # fi.set_color_map('YlOrBr_r')
        # fi.set_callback('drag-drop', self.drop_file)
        # fi.set_bg(0.2, 0.2, 0.2)
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
        viewer_hbox.setContentsMargins(QtCore.QMargins(4,1,4,1))
        hw = QtGui.QWidget()
        hw.setLayout(viewer_hbox)
        vbox.addWidget(hw)
        click_hbox = QtGui.QHBoxLayout()
        click_hbox.setObjectName("click_hbox")
        self.clickinfo = QtGui.QLabel("Click the image to pan.")
        self.clickinfo.setObjectName("clickinfo")
        click_hbox.addWidget(self.clickinfo)
        click_hbox.setContentsMargins(QtCore.QMargins(4,1,4,1))
        hw = QtGui.QWidget()
        hw.setLayout(click_hbox)
        vbox.addWidget(hw)
        readout_hbox = QtGui.QHBoxLayout()
        readout_hbox.setObjectName("readout_hbox")
        self.readout = QtGui.QLabel("X:                 Y:                    Value:")
        self.readout.setObjectName("readout")
        self.readout.setMinimumSize(QtCore.QSize(350, 0))
        readout_hbox.addWidget(self.readout)
        self.wcmap = QtGui.QComboBox()
        for name in cmap.get_names():
            self.wcmap.addItem(name)
        self.wcmap.currentIndexChanged.connect(self.cmap_change)
        self.wcmap.setCurrentText('gray')
        readout_hbox.addWidget(self.wcmap)
        self.wcut = QtGui.QComboBox()
        for name in fi.get_autocut_methods():
            self.wcut.addItem(name)
        self.wcut.currentIndexChanged.connect(self.cut_change)
        self.wcut.setCurrentText('zscale')
        readout_hbox.addWidget(self.wcut)
        self.wcolor = QtGui.QComboBox()
        for name in fi.get_color_algorithms():
            self.wcolor.addItem(name)
        self.wcolor.currentIndexChanged.connect(self.color_change)
        readout_hbox.addWidget(self.wcolor)
        readout_hbox.setContentsMargins(QtCore.QMargins(4,1,4,1))
        hw = QtGui.QWidget()
        hw.setLayout(readout_hbox)
        vbox.addWidget(hw)
        file_hbox = QtGui.QHBoxLayout()
        file_hbox.setObjectName("file_hbox")
        self.file_info = QtGui.QLabel("File: ")
        self.file_info.setObjectName("file_info")
        file_hbox.addWidget(self.file_info)
        self.box_readout = QtGui.QLabel("Amplitude:                  FWHM: ")
        self.box_readout.setMinimumSize(QtCore.QSize(350, 0))
        self.box_readout.setObjectName("box_readout")
        file_hbox.addWidget(self.box_readout)
        file_hbox.setContentsMargins(QtCore.QMargins(4,1,4,1))
        hw = QtGui.QWidget()
        hw.setLayout(file_hbox)
        vbox.addWidget(hw)
        sky_hbox = QtGui.QHBoxLayout()
        sky_hbox.setObjectName("sky_hbox")
        self.sky_info = QtGui.QLabel("Sky: ")
        self.sky_info.setObjectName("sky_info")
        sky_hbox.addWidget(self.sky_info)
        sky_hbox.setContentsMargins(QtCore.QMargins(8,1,8,1))
        hw = QtGui.QWidget()
        hw.setLayout(sky_hbox)
        vbox.addWidget(hw)
        # sky_hbox = QtGui.QHBoxLayout()
        # file_hbox.setObjectName("sky_hbox")
        # self.sky_info = QtGui.QLabel("Sky: ")
        # self.sky_info.setObjectName("sky_info")
        # sky_hbox.addWidget(self.sky_info)
        # hw = QtGui.QWidget()
        # hw.setLayout(sky_hbox)
        # vbox.addWidget(hw)
        buttons_hbox = QtGui.QHBoxLayout()
        buttons_hbox.setObjectName("buttons_hbox")
        buttons_vbox_left = QtGui.QVBoxLayout()
        buttons_vbox_left.setContentsMargins(QtCore.QMargins(0, 0, 10, 0))
        buttons_vbox_left.setObjectName("buttons_vbox_left")
        self.wrecenter = QtGui.QPushButton("Re-center Image")
        self.wrecenter.setObjectName("wrecenter")
        self.wrecenter.clicked.connect(self.recenter)
        buttons_vbox_left.addWidget(self.wrecenter)
        self.wmovSlitCent = QtGui.QPushButton("Center on Slit")
        self.wmovSlitCent.setObjectName("wmovSlitCent")
        self.wmovSlitCent.clicked.connect(self.movSlitCent)
        buttons_vbox_left.addWidget(self.wmovSlitCent)
        self.wmovObj = QtGui.QPushButton("Move Object")
        self.wmovObj.setObjectName("wmovObj")
        self.wmovObj.clicked.connect(self.movObj)
        buttons_vbox_left.addWidget(self.wmovObj)
        hw = QtGui.QWidget()
        hw.setLayout(buttons_vbox_left)
        buttons_hbox.addWidget(hw)
        # buttons_vbox_cent = QtGui.QVBoxLayout()
        # buttons_vbox_cent.setObjectName("buttons_vbox_cent")
        # self.wsky = QtGui.QPushButton("Load Sky")
        # self.wsky.setObjectName("wsky")
        # self.wsky.clicked.connect(self.load_sky)
        # self.wsky.setEnabled(False)
        # buttons_vbox_cent.addWidget(self.wsky)
        # self.wclearsky = QtGui.QPushButton("Clear Sky")
        # self.wclearsky.setObjectName("wclearsky")
        # self.wclearsky.clicked.connect(self.clearsky)
        # self.wclearsky.setEnabled(False)
        # buttons_vbox_cent.addWidget(self.wclearsky)
        # hw = QtGui.QWidget()
        # hw.setLayout(buttons_vbox_cent)
        # buttons_hbox.addWidget(hw)
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

        self.movSlitCursor = False
        self.autoCenter = False
        self.init_x = 0.0
        self.init_y = 0.0
        self.second_click =False

        fi.set_callback('cursor-changed', self.motion_cb)
        fi.add_callback('cursor-down', self.btndown)

        # self.recdc, self.compdc = self.add_canvas()
        self.recdc = self.add_canvas()
        self.picktag = "pick-box"

        # self.sky = ""
        self.curentfile = ""
        # self.sdiff_done = False

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
        try:
            sky = os.path.realpath(os.readlink("TEMPSKY2.fits"))
        except FileNotFoundError:
            sky = "None"
        self.sky_info.setText(f"Sky: {sky}")

    def add_canvas(self, tag=None):
        # add a canvas to the view
        my_canvas = self.fitsimage.get_canvas()
        RecCanvas = my_canvas.get_draw_class('rectangle')
        # CompCanvas = my_canvas.get_draw_class('compass')
        # return RecCanvas, CompCanvas
        return RecCanvas

    def cmap_change(self):
        self.fitsimage.set_color_map(self.wcmap.currentText())

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
        # self.threadpool = False
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
        print("scan started...")
        scanner = Scanner(self.scan)
        scanner.signals.load.connect(self.load_file)
        self.threadpool.start(scanner)
        # self.wstartscan.setEnabled(False)
        # self.wstopscan.setEnabled(True)

    def stop_scan(self):
        self.scanning = False
        print('Scanning stopped.')
        # self.wstartscan.setEnabled(True)
        # self.wstopscan.setEnabled(False)

    def load_file(self, filepath):
        print(filepath)
        self.currentfile = filepath
        recenter = False
        if self.fitsimage.get_image() == None:
            recenter = True
        image = load_data(filepath, logger=self.logger)
        self.fitsimage.set_image(image)
        # self.setWindowTitle(filepath)
        try:
            self.fitsimage.get_canvas().get_object_by_tag(self.picktag)
            self.fitsimage.get_canvas().delete_object_by_tag(self.picktag)
        except KeyError:
            pass
        if recenter == True:
            self.recenter()
        print(f"Loaded {filepath}")
        self.file_info.setText(f"File: {filepath}")
        # self.wsky.setEnabled(True)
        # self.wsdiff.setEnabled(True)
        # self.sdiff_done = False
        # self.wsdiff.setText("SDiff")
        # if self.sky != "":
        #     self.subtract_sky(self.sky)

    # def load_sky(self):
    #     res = QtGui.QFileDialog.getOpenFileName(self, "Load Sky file",
    #                                             str(self.nightpath()))
    #     if isinstance(res, tuple):
    #         fileName = res[0]
    #     else:
    #         fileName = str(res)
    #     if len(fileName) != 0:
    #         self.sky = fileName
    #         self.wclearsky.setEnabled(True)
    #         self.subtract_sky(self.sky)
    
    # def clearsky(self):
    #     self.sky = ""
    #     image = load_data(self.currentfile, logger=self.logger)
    #     self.fitsimage.set_image(image)
    #     self.wclearsky.setEnabled(False)

    def open_file(self):
        res = QtGui.QFileDialog.getOpenFileName(self, "Open FITS file",
                                                str(self.nightpath()))
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.load_file(fileName)

    # def sdiff(self):
    #     if self.sdiff_done == False:
    #         image = self.fitsimage.get_image()
    #         data = image.get_data()
    #         previous = fits.getdata(str(self.previous_image))
    #         subtracted = data - previous
    #         self.fitsimage.set_data(subtracted)
    #         self.wsdiff.setText("Undo SDiff")
    #         self.sdiff_done = True
    #     else:
    #         image = load_data(self.currentfile, logger=self.logger)
    #         self.fitsimage.set_image(image)
    #         self.wsdiff.setText("SDiff")
    #         self.sdiff_done = False

    # def subtract_sky(self, file):
    #     image = self.fitsimage.get_image()
    #     data = image.get_data()
    #     sky = fits.getdata(file)
    #     try:
    #         subtracted = data - sky
    #         self.fitsimage.set_data(subtracted)
    #     except ValueError:
    #         self.fitsimage.set_data(data)

    ##Start of image find and processing code

    def scan(self, file_callback):
        # self.previous_image = self.slit_lastfile.read() #TODO this is to get first previous image, might remove.
        while self.scanning:
            # if (self.go == 1 or self.test == 1 or self.display2 == 1) and ("v" in self.slit_filename or "TEMP" in self.slit_filename):
            if self.display2 == 1:
                # previm = self.slit_lastfile.read()
                print("Taking image")
                # self.waitForFileToBeUnlocked(0.5)
                file_callback.emit(str(self.dispname2.read()))
                self.waitForZero(0.25)
                # self.previous_image = previm
            time.sleep(0.25)

    def fileIsCurrentlyLocked(self):
        print(f'display2 {self.display2} locked')
        locked = True
        # if int(self.go.read()) == 0 and int(self.test.read()) == 0:
        if int(self.display2.read()) == 0:
            print(f'display2 {self.display2} unlocked')
            locked = False
        return locked
    
    def waitForZero(self, wait_time):
        while self.display2 == 1:
            time.sleep(wait_time)
        print(f'display2 {self.display2}')

    def waitForFileToBeUnlocked(self, wait_time):
        while self.fileIsCurrentlyLocked():
            time.sleep(wait_time)

    def nightpath(self):
        file = self.dispname2
        dir = str(file).split("//")
        # dir = str(file).replace('sdiff.fits', '')
        dir = dir[0] + "/"
        print(dir)
        # path = dir[0]
        nightly = Path(dir)
        return nightly

    def writeFits(self, headerinfo, image_data):
        hdu = fits.PrimaryHDU(header=headerinfo, data=image_data)
        filename = 'procImage.fits'
        try:
            hdu.writeto(filename)
        except OSError:
            os.remove(filename)
            hdu.writeto(filename)
        return filename

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
    
    def recenter(self):
        self.fitsimage.zoom_fit()
    
    def movSlitCent(self):
        if self.movSlitCursor == False:
            self.movSlitCursor = True
            self.autoCenter = True
            self.clickinfo.setText("Click star to center on slit.")
            self.wmovSlitCent.setText("Cancel Slit Center")
        else:
            self.movSlitCursor = False
            self.autoCenter = False
            self.clickinfo.setText("Click image to pan.")
            self.wmovSlitCent.setText("Center on Slit")
    
    def movObj(self):
        if self.movSlitCursor == False:
            self.movSlitCursor = True
            self.clickinfo.setText("Click star you want to move.")
            self.wmovObj.setText("Cancel Object Move")
        else:
            self.movSlitCursor = False
            self.autoCenter = False
            self.clickinfo.setText("Click image to pan.")
            self.wmovObj.setText("Move Object")
        
    def btndown(self, canvas, event, data_x, data_y):
        self.xclick = data_x
        self.yclick = data_y
        if self.movSlitCursor == True:
            if self.autoCenter == True:
                self.movAuto(data_x, data_y)
                self.movSlitCursor = False
                self.autoCenter = False
                self.wmovSlitCent.setEnabled(True)
                self.clickinfo.setText("Click image to pan.")
                self.wmovSlitCent.setText("Center on Slit")
            elif self.second_click == True:
                self.movManual(self.init_x, self.init_y, data_x, data_y)
                self.movSlitCursor = False
                self.autoCenter = False
                self.wmovSlitCent.setEnabled(True)
                self.clickinfo.setText("Click image to pan.")
                self.wmovObj.setText("Move Object")
            else:
                self.init_x = data_x
                self.init_y = data_y
                self.second_click == True
                self.clickinfo.setText("Click where you want to move the object.")
        else:
            self.fitsimage.set_pan(data_x, data_y)
            # self.pickstar(self.fitsimage)
    
    def movAuto(self, data_x, data_y):
        pscale = 0.123
        dx = (data_x - 121.8) * pscale 
        dy = (data_y - 464.3) * pscale
        # self.dcs['instxoff'].write(dx, rel2curr = 't')
        # self.dcs['instyoff'].write(dy, rel2curr = 't')

    def movManual(self, x1, y1, x2, y2):
        pscale = 0.123
        dx = (x1 - x2) * pscale 
        dy = (y1 - y2) * pscale
        # self.dcs['instxoff'].write(dx, rel2curr = 't')
        # self.dcs['instyoff'].write(dy, rel2curr = 't')

def main():

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
    w.start_scan()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
