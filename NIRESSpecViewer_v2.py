import os, time, sys
# from os import listdir
# from os.path import abspath, isfile, join
from pathlib import Path
# import datetime
# import csv
# import copy
from functools import partial
import configparser
# import warnings
# warnings.filterwarnings("error")

import numpy as np
from astropy.io import fits
# from astropy import wcs
import astropy.units as u
# from astropy.stats import gaussian_sigma_to_fwhm
from astropy.modeling import models, fitting
# import PIL.Image as PILimage

from ginga import cmap
from ginga.misc import log
from ginga.qtw.QtHelp import QtGui, QtCore
from ginga.qtw.ImageViewQt import CanvasView
from ginga.util import iqcalc, io_fits, plots
from ginga.util.loader import load_data
from ginga.AstroImage import AstroImage
from ginga.gw import Plot, Widgets

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

class NewFileSignals(QtCore.QObject):
    load = QtCore.Signal()

class NewFile(QtCore.QRunnable):
    '''
    NewFile thread
    '''
    def __init__(self, fn, *args, **kwargs):
        super(NewFile, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = NewFileSignals()
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

##Cuts
class Cuts(Widgets.Box):

    def __init__(self, logger, fitsimage, dispname, coadds, bm):
        super(Cuts, self).__init__(fitsimage)

        bm.reset_mode(fitsimage)

        self.dispname = dispname

        self.logger = logger

        self.layertag = 'cuts-canvas'
        self._new_cut = 'New Cut'
        self.cutstag = self._new_cut
        self.tags = [self._new_cut]
        # get Cuts preferences
        self.fitsimage = fitsimage
        self.coadds = coadds

        self.dc = self.fitsimage.get_canvas().get_draw_classes()
        canvas = self.dc.DrawingCanvas()
        canvas.enable_draw(True)
        canvas.enable_edit(True)
        canvas.set_drawtype('line', color='cyan', linestyle='dash')
        canvas.set_callback('draw-event', self.draw_cb)
        canvas.add_draw_mode('move', down=self.buttondown_cb,
                             move=self.motion_cb, up=self.buttonup_cb,
                             key=self.keydown)
        canvas.set_draw_mode('draw')
        canvas.register_for_cursor_drawing(self.fitsimage)
        canvas.set_surface(self.fitsimage)
        self.canvas = canvas

        self.cuts_image = None

        self.gui_up = False

        vbox = Widgets.VBox()

        self.cuts_plot = plots.CutsPlot(logger=self.logger,
                                        width=700, height=400)
        self.plot = Plot.PlotWidget(self.cuts_plot)
        self.plot.resize(400, 400)
        ax = self.cuts_plot.add_axis()
        ax.grid(True)
        vbox.add_widget(self.plot)
        control_hbox = Widgets.HBox()
        self.freedraw = Widgets.Button("Free Draw")
        self.freedraw.add_callback('activated', self.free_draw_cb)
        control_hbox.add_widget(self.freedraw)
        self.horizontaldraw = Widgets.Button("Horizontal")
        self.horizontaldraw.add_callback('activated', self.horizontal_draw_cb)
        control_hbox.add_widget(self.horizontaldraw)
        self.verticaldraw = Widgets.Button("Vertical")
        self.verticaldraw.add_callback('activated', self.vertical_draw_cb)
        control_hbox.add_widget(self.verticaldraw)
        vbox.add_widget(control_hbox)
        self.closebtn = Widgets.Button("Close")
        self.closebtn.add_callback('activated', self.dismiss)
        vbox.add_widget(self.closebtn)
        self.add_widget(vbox)

        self.cut_mode = "Free"
        self.freedraw.set_enabled(False)

        self.start()
        self.gui_up = True
        self.threadpool = QtCore.QThreadPool()
        self.start_filecheck()

    def free_draw_cb(self, event):
        self.cut_mode = "Free"
        self.freedraw.set_enabled(False)
        self.horizontaldraw.set_enabled(True)
        self.verticaldraw.set_enabled(True)
    
    def horizontal_draw_cb(self, event):
        self.cut_mode = "Horizontal"
        self.freedraw.set_enabled(True)
        self.horizontaldraw.set_enabled(False)
        self.verticaldraw.set_enabled(True)
    
    def vertical_draw_cb(self, event):
        self.cut_mode = "Vertical"
        self.freedraw.set_enabled(True)
        self.horizontaldraw.set_enabled(True)
        self.verticaldraw.set_enabled(False)
    
    def delete_all(self):
        self.canvas.delete_all_objects()
        self.tags = [self._new_cut]
        self.cutstag = self._new_cut
        self.cuts_plot.clear()

    def add_cuts_tag(self, tag):
        if tag not in self.tags:
            self.tags.append(tag)

    def start(self):
        # start line cuts operation
        self.canvas.enable_draw(True)
        self.cuts_plot.set_titles(rtitle="Cuts")

        # insert canvas, if not already
        p_canvas = self.fitsimage.get_canvas()
        try:
            p_canvas.get_object_by_tag(self.layertag)

        except KeyError:
            # Add ruler layer
            p_canvas.add(self.canvas, tag=self.layertag)

        self.resume()

    def pause(self):
        self.canvas.ui_set_active(False)

    def resume(self):
        # turn off any mode user may be in
        # self.modes_off()

        self.canvas.ui_set_active(True, viewer=self.fitsimage)
        self.replot_all()

    def stop(self):
        self.gui_up = False
        # remove the canvas from the image
        p_canvas = self.fitsimage.get_canvas()
        p_canvas.delete_object_by_tag(self.layertag)

    def redo(self):
        """This is called when a new image arrives or the data in the
        existing image changes.
        """

        self.replot_all()

    def _plotpoints(self, obj, color):

        image = self.fitsimage.get_vip()

        # Get points on the line
        points = image.get_pixels_on_line(int(obj.x1), int(obj.y1),
                                                  int(obj.x2), int(obj.y2))

        points = np.divide(np.array(points), float(self.coadds))

        self.cuts_plot.cuts(points, title = f"{self.cut_mode} Cut", xtitle="Line Index", ytitle="ADUs/COADD",
                            color=color)


    def _replot(self, lines):
        for idx in range(len(lines)):
            line= lines[idx]
            self._plotpoints(line, "blue")

        return True

    def replot_all(self):
        self.cuts_plot.clear()
        # self.w.delete_all.set_enabled(False)
        # self.save_cuts.set_enabled(False)

        # idx = 0
        for cutstag in self.tags:
            if cutstag == self._new_cut:
                continue
            obj = self.canvas.get_object_by_tag(cutstag)
            lines = self._getlines(obj)
            self._replot(lines)

        self.cuts_plot.draw()

        self.canvas.redraw(whence=3)

        return True

    def _create_cut(self, x1, y1, x2, y2, color='cyan'):
        self.delete_all()
        text = "cut"
        # if not self.settings.get('label_cuts', False):
        #     text = ''
        line_obj = self.dc.Line(x1, y1, x2, y2, color=color,
                                showcap=False)
        text_obj = self.dc.Text(0, 0, text, color=color, coord='offset',
                                ref_obj=line_obj)
        obj = self.dc.CompoundObject(line_obj, text_obj)
        # this is necessary for drawing cuts with width feature
        obj.initialize(self.canvas, self.fitsimage, self.logger)
        obj.set_data(cuts=True)
        return obj

    def _create_cut_obj(self, cuts_obj, color='cyan'):
        self.delete_all()
        text = "cut"
        # if not self.settings.get('label_cuts', False):
        #     text = ''
        cuts_obj.showcap = False
        cuts_obj.linestyle = 'solid'
        #cuts_obj.color = color
        color = cuts_obj.color
        args = [cuts_obj]
        text_obj = self.dc.Text(0, 0, text, color=color, coord='offset',
                                ref_obj=cuts_obj)
        args.append(text_obj)

        obj = self.dc.CompoundObject(*args)
        obj.set_data(cuts=True)

        return obj

    def _getlines(self, obj):
        return [obj.objects[0]]

    def buttondown_cb(self, canvas, event, data_x, data_y, viewer):
        return self.motion_cb(canvas, event, data_x, data_y, viewer)

    def motion_cb(self, canvas, event, data_x, data_y, viewer):


        if self.cutstag == self._new_cut:
            return True
        obj = self.canvas.get_object_by_tag(self.cutstag)
        # Assume first element of this compound object is the reference obj
        obj = obj.objects[0]
        obj.move_to_pt((data_x, data_y))
        canvas.redraw(whence=3)

        if self.drag_update:
            self.replot_all()
        return True

    def buttonup_cb(self, canvas, event, data_x, data_y, viewer):
        if self.cutstag == self._new_cut:
            return True
        obj = self.canvas.get_object_by_tag(self.cutstag)
        # Assume first element of this compound object is the reference obj
        obj = obj.objects[0]
        obj.move_to_pt((data_x, data_y))

        self.replot_all()
        return True

    def keydown(self, canvas, event, data_x, data_y, viewer):
        return True

    def cut_at(self, cuttype):
        """Perform a cut at the last mouse position in the image.
        `cuttype` determines the type of cut made.
        """
        data_x, data_y = self.fitsimage.get_last_data_xy()
        image = self.fitsimage.get_image()
        wd, ht = image.get_size()

        coords = []
        if cuttype == 'horizontal':
            coords.append((0, data_y, wd - 1, data_y))
        elif cuttype == 'vertical':
            coords.append((data_x, 0, data_x, ht - 1))

        tag = "cut"
        cuts = []
        for (x1, y1, x2, y2) in coords:
            # calculate center of line
            wd = x2 - x1
            dw = wd // 2
            ht = y2 - y1
            dh = ht // 2
            x, y = x1 + dw + 4, y1 + dh + 4

            cut = self._create_cut(x1, y1, x2, y2, color='cyan')
            cuts.append(cut)
            
        cut = cuts[0]

        cut.set_data(count=True)

        self.canvas.delete_object_by_tag(tag)
        self.canvas.add(cut, tag=tag)
        self.add_cuts_tag(tag)

        self.logger.debug("redoing cut plots")
        return self.replot_all()

    def draw_cb(self, canvas, tag):
        if self.cut_mode == "Horizontal":
            return self.cut_at('horizontal')
        elif self.cut_mode == "Vertical":
            return self.cut_at('vertical')
        obj = canvas.get_object_by_tag(tag)
        canvas.delete_object_by_tag(tag)

        tag = "cut"

        cut = self._create_cut_obj(obj, color='cyan')
        cut.set_data(count=True)

        canvas.delete_object_by_tag(tag)
        self.canvas.add(cut, tag=tag)
        self.add_cuts_tag(tag)

        self.logger.debug("redoing cut plots")
        return self.replot_all()
    
    
    def start_filecheck(self):
        self.currentfile = self.fitsimage.get_image().get_header()['DISPNAME']
        self.filechecking = True
        filechecker = NewFile(self.new_file)
        filechecker.signals.load.connect(self.file_compare)
        self.threadpool.start(filechecker)
    
    def new_file(self, file_callback):
        while self.filechecking:
            file_callback.emit()
            time.sleep(0.5)

    def stop_filecheck(self):
        self.filechecking = False
    
    def file_compare(self):
        image = self.fitsimage.get_image()
        name = image.get_header()['DISPNAME']
        if self.currentfile != name:
            self.replot_all()
            self.currentfile = name

    def dismiss(self, event):
        self.stop_filecheck()
        self.stop()
        self.delete()
# END

class MathWindow(Widgets.Box):
    """
    This "window" is a QWidget. If it has no parent, it
    will appear as a free-floating window as we want.
    """
    def __init__(self, logger, fitsimage, loadfile, dispname, sdiff):
        super(MathWindow, self).__init__(fitsimage)

        self.logger = logger
        self.fitsimage = fitsimage
        self.load_file = loadfile
        self.dispname = dispname
        self.sdiff = sdiff

        # self.math_path = "/home/jpelletier/NIRES/"

        vbox = Widgets.VBox()
        math_hbox = Widgets.HBox()
        filebutton_vbox = Widgets.VBox()
        self.wfileone = Widgets.Button("Open File 1")
        self.wfileone.add_callback('activated', self.openfileone)
        filebutton_vbox.add_widget(self.wfileone)
        self.wfiletwo = Widgets.Button("Open File 2")
        self.wfiletwo.add_callback('activated', self.openfiletwo)
        filebutton_vbox.add_widget(self.wfiletwo)
        math_hbox.add_widget(filebutton_vbox)
        filename_vbox = Widgets.VBox()
        self.filenameone = Widgets.Label ("File one")
        filename_vbox.add_widget(self.filenameone)
        self.filenametwo = Widgets.Label ("File two")
        filename_vbox.add_widget(self.filenametwo)
        math_hbox.add_widget(filename_vbox)
        function_vbox = Widgets.VBox()
        self.wsubtract = Widgets.Button("Subtract")
        self.wsubtract.add_callback('activated', self.imageSubtract)
        function_vbox.add_widget(self.wsubtract)
        self.wadd = Widgets.Button("Add")
        self.wadd.add_callback('activated', self.imageAdd)
        function_vbox.add_widget(self.wadd)
        math_hbox.add_widget(function_vbox)
        vbox.add_widget(math_hbox)
        self.wreload = Widgets.Button("Reload Latest Image")
        self.wreload.add_callback('activated', self.reload)
        vbox.add_widget(self.wreload)
        self.closebtn = Widgets.Button("Close")
        self.closebtn.add_callback('activated', self.dismiss)
        vbox.add_widget(self.closebtn)
        self.add_widget(vbox)
        self.resize(500, 0)

    def nightpath(self):
        dir = str(self.dispname)
        if "//" in str(dir):
            dir = str(dir).split("//")
            dir = dir[0] + "/"
            nightly = Path(dir)
        else: 
            dir = Path(dir)
            nightly = dir.parent
        return nightly
    
    def mathFileNames(self, firstfile, secondfile, operation):
        if "//" in firstfile:
            firstfile = firstfile.split("//")
            firstfile = firstfile[-1]
        else: 
            firstfile = firstfile.split("/")
            firstfile = firstfile[-1]
        if "//" in secondfile:
            secondfile = secondfile.split("//")
            secondfile = secondfile[-1]
        else: 
            secondfile = secondfile.split("/")
            secondfile = secondfile[-1]
        return f'{firstfile} {operation} {secondfile}.fits'
    
    def openfileone(self, event):
        res = QtGui.QFileDialog.getOpenFileName(caption="Open FITS file 1",
                                                directory = str(self.nightpath()))
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.filenameone.set_text(fileName)
        return
    
    def openfiletwo(self, event):
        res = QtGui.QFileDialog.getOpenFileName(caption="Open FITS file 1",
                                                directory = str(self.nightpath()))
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.filenametwo.set_text(fileName)
        return

    def imageSubtract(self, event):
        try:
            imageone_data = fits.getdata(self.filenameone.get_text())
        except FileNotFoundError:
            return
        except OSError:
            self.filenameone.set_text("File still exposing")
            return
        try:
            imagetwo_data = fits.getdata(self.filenametwo.get_text())
        except FileNotFoundError:
            return
        except OSError:
            self.filenametwo.set_text("File not loaded")
            return
        image_header = fits.getheader(self.filenameone.get_text())
        subtracted = imageone_data - imagetwo_data
        hdu = fits.PrimaryHDU(header=image_header, data=subtracted)
        filename = self.mathFileNames(self.filenameone.get_text(), self.filenametwo.get_text(), '-')
        # full_path = Path(self.math_path + filename)
        full_path = filename
        try:
            hdu.writeto(full_path)
        except OSError:
            os.remove(full_path)
            hdu.writeto(full_path)
        self.load_file(str(full_path))
        os.remove(full_path)
        self.sdiff_done = False
        return
    
    def imageAdd(self, event):
        try:
            imageone_data = fits.getdata(self.filenameone.get_text())
            imagetwo_data = fits.getdata(self.filenametwo.get_text())
            image_header = fits.getheader(self.filenameone.get_text())
            added = imageone_data + imagetwo_data
            hdu = fits.PrimaryHDU(header=image_header, data=added)
            filename = self.mathFileNames(self.filenameone.get_text(), self.filenametwo.get_text(), '+')
            # full_path = Path(self.math_path + filename)
            full_path = filename
            try:
                hdu.writeto(full_path)
            except OSError:
                os.remove(full_path)
                hdu.writeto(full_path)
            self.load_file(str(full_path))
            os.remove(full_path)
            self.sdiff_done = False
        except FileNotFoundError:
            return
        return

    def reload(self, event):
        self.load_file(str(self.dispname.read()))

    def stop(self):
        self.gui_up = False

    def dismiss(self, event):
        self.stop()
        self.delete()


class FitsViewer(QtGui.QMainWindow):
    # resized = QtCore.Signal()

    def __init__(self, logger):
        super(FitsViewer, self).__init__()

        # self.resized.connect(self.recenter)

        # def resizeEvent(self, event):
        #     self.resized.emit()
        #     return super(FitsViewer, self).resizeEvent(event)
        
        self.logger = logger

        self.cachedFiles = None
        #KTL stuff
        #Cache KTL keywords
        # self.slit_filename = ktl.cache('nids', 'FILENAME')
        # self.slit_filename.monitor()
        self.spec_lastfile = ktl.cache('nsds', 'LASTFILE')
        # self.slit_sdiff = ktl.cache('nids', 'LASTFILE')
        # self.go = ktl.cache('nids', 'GO')
        # self.go.monitor()
        # self.test = ktl.cache('nids', 'test')
        # self.test.monitor()
        self.display = ktl.cache('nsds', 'display')
        self.display.monitor()
        # self.dispname = ktl.cache('nsds', 'dispname')
        self.dispname = ktl.cache('nsds', 'filename')
        self.dispname.monitor()
        self.coadds = ktl.cache('nids', 'coadds')
        self.coadds.monitor()

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

        menubar = self.menuBar()

        filemenu = menubar.addMenu("File")

        item = QtGui.QAction("Open File", menubar)
        item.triggered.connect(self.open_file)
        filemenu.addAction(item)

        sep = QtGui.QAction(menubar)
        sep.setSeparator(True)
        filemenu.addAction(sep)

        item = QtGui.QAction("Quit", menubar)
        item.triggered.connect(self.quit)
        filemenu.addAction(item)

        mathmenu = menubar.addMenu("Math")

        item = QtGui.QAction("SDiff/Undo", menubar)
        item.triggered.connect(self.sdiff)
        mathmenu.addAction(item)

        sep = QtGui.QAction(menubar)
        sep.setSeparator(True)
        mathmenu.addAction(sep)

        item = QtGui.QAction("Math", menubar)
        item.triggered.connect(self.math_popup)
        mathmenu.addAction(item)

        cutmenu = menubar.addMenu("Cuts")

        item = QtGui.QAction("Cut GUI", menubar)
        item.triggered.connect(self.cuts_popup)
        cutmenu.addAction(item)

        colormenu = menubar.addMenu("Colors")
        for cm_name in cmap.get_names():
            item = QtGui.QAction(cm_name, menubar)
            
            item.triggered.connect(partial(self.cmap_change, cm_name))
            colormenu.addAction(item)

        cutmenu = menubar.addMenu("Display Parameters")
        for cut_name in fi.get_autocut_methods():
            item = QtGui.QAction(cut_name, menubar)
            
            item.triggered.connect(partial(self.cut_change, cut_name))
            cutmenu.addAction(item)

        cutmenu = menubar.addMenu("Stretch")
        for stretch_name in fi.get_color_algorithms():
            item = QtGui.QAction(stretch_name, menubar)
            
            item.triggered.connect(partial(self.color_change, stretch_name))
            cutmenu.addAction(item)

        self.bd = fi.get_bindings()
        self.bd.enable_all(False)
        self.bm = fi.get_bindmap()
        self.bm.reset_mode(fi)
        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setObjectName("vbox")
        viewer_hbox = QtGui.QHBoxLayout()
        viewer_hbox.setObjectName("viewer_hbox")
        w = fi.get_widget()
        # w.setMinimumSize(QtCore.QSize(1200, 600))
        viewer_hbox.addWidget(w)
        viewer_hbox.setContentsMargins(QtCore.QMargins(4,1,4,1))
        viewerHB = QtGui.QWidget()
        viewerHB.setLayout(viewer_hbox)
        # vbox.addWidget(hw)
        buttons_vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setObjectName("bvbox")
        click_hbox = QtGui.QHBoxLayout()
        click_hbox.setObjectName("click_hbox")
        self.clickinfo = QtGui.QLabel("Left-drag to manually adjust levels.")
        self.clickinfo.setObjectName("clickinfo")
        # self.clickinfo.setMinimumSize(QtCore.QSize(300, 0))
        click_hbox.addWidget(self.clickinfo)
        self.wzoomin = QtGui.QPushButton("Zoom In")
        self.wzoomin.setObjectName("wzoomin")
        self.wzoomin.clicked.connect(self.zoomIn)
        self.wzoomin.setMaximumSize(QtCore.QSize(75, 40))
        click_hbox.addWidget(self.wzoomin)
        self.wzoomreset = QtGui.QPushButton("Reset Zoom")
        self.wzoomreset.setObjectName("wzoomreset")
        self.wzoomreset.clicked.connect(self.zoomReset)
        self.wzoomreset.setMaximumSize(QtCore.QSize(80, 40))
        click_hbox.addWidget(self.wzoomreset)
        self.wzoomout = QtGui.QPushButton("Zoom Out")
        self.wzoomout.setObjectName("wzoomout")
        self.wzoomout.clicked.connect(self.zoomOut)
        self.wzoomout.setMaximumSize(QtCore.QSize(75, 40))
        click_hbox.addWidget(self.wzoomout)
        self.wsetpan = QtGui.QPushButton("Pan")
        self.wsetpan.setObjectName("wsetpan")
        self.wsetpan.clicked.connect(self.setPan)
        self.wsetpan.setMaximumSize(QtCore.QSize(75, 40))
        click_hbox.addWidget(self.wsetpan)
        self.wrecenter = QtGui.QPushButton("Recenter")
        self.wrecenter.setObjectName("wrecenter")
        self.wrecenter.clicked.connect(self.recenter)
        self.wrecenter.setMaximumSize(QtCore.QSize(75, 40))
        click_hbox.addWidget(self.wrecenter)
        # self.wwindow = QtGui.QPushButton("Window")
        # self.wwindow.setObjectName("wwindow")
        # self.wwindow.clicked.connect(self.window)
        # self.wwindow.setMaximumSize(QtCore.QSize(75, 40))
        # click_hbox.addWidget(self.wwindow)
        click_hbox.setContentsMargins(QtCore.QMargins(4,1,4,1))
        hw = QtGui.QWidget()
        hw.setLayout(click_hbox)
        buttons_vbox.addWidget(hw)
        readout_hbox = QtGui.QHBoxLayout()
        readout_hbox.setObjectName("readout_hbox")
        self.readout = QtGui.QLabel("X:                 Y:                    Value:                   Wavelength:")
        self.readout.setObjectName("readout")
        # self.readout.setMinimumSize(QtCore.QSize(550, 0))
        readout_hbox.addWidget(self.readout)
        readout_hbox.setContentsMargins(QtCore.QMargins(4,1,4,1))
        hw = QtGui.QWidget()
        hw.setLayout(readout_hbox)
        buttons_vbox.addWidget(hw)
        file_hbox = QtGui.QHBoxLayout()
        file_hbox.setObjectName("file_hbox")
        self.file_info = QtGui.QLabel("File: ")
        self.file_info.setObjectName("file_info")
        file_hbox.addWidget(self.file_info)
        file_hbox.setContentsMargins(QtCore.QMargins(4,1,4,1))
        hw = QtGui.QWidget()
        hw.setLayout(file_hbox)
        buttons_vbox.addWidget(hw)
        buttonsVB = QtGui.QWidget()
        buttonsVB.setLayout(buttons_vbox)

        splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        splitter.addWidget(viewerHB)
        splitter.addWidget(buttonsVB)
        splitter.setStretchFactor(0, 1)
        splitter.setSizes([600, 40])
        vbox.addWidget(splitter)

        vw = QtGui.QWidget()
        self.setCentralWidget(vw)
        vw.setLayout(vbox)

        fi.set_callback('cursor-changed', self.motion_cb)
        fi.add_callback('cursor-down', self.btndown)
        fi.add_callback('cursor-move', self.drag_cb)

        self.recdc = self.add_canvas()
        self.picktag = "pick-box"

        # self.sky = ""
        self.curentfile = ""
        self.sdiff_done = False
        self.c = None
        self.m = None
        self.panning = False
        self.base_zoom = 0
        # self.sdiff = False
        current_directory = os.path.dirname(os.path.abspath(__file__))
        config_file = Path(current_directory) / 'config.ini'

        config = configparser.ConfigParser()
        config.read(config_file)
        file_path = config.get('DEFAULT', 'FILE_PATH')
        self.wavelength_data = np.flip((fits.getdata(file_path)), 0)

        self.fitsimage.set_color_map('heat')

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
        nightly = self.nightpath()
        # tempsky = str(nightly) + "/TEMPSKY2.fits"
        # try:
        #     sky = os.path.realpath(os.readlink(tempsky))
        # except FileNotFoundError:
        #     sky = "None"
        # except NotADirectoryError:
        #     sky = "None"
        # self.sky_info.setText(f"Sky: {sky}")

    def add_canvas(self, tag=None):
        # add a canvas to the view
        my_canvas = self.fitsimage.get_canvas()
        RecCanvas = my_canvas.get_draw_class('rectangle')
        # CompCanvas = my_canvas.get_draw_class('compass')
        # return RecCanvas, CompCanvas
        return RecCanvas

    def cmap_change(self, cm_name):
        self.fitsimage.set_color_map(cm_name)

    def cut_change(self, cut_name):
        self.fitsimage.set_autocut_params(cut_name)

    def color_change(self, stretch_name):
        self.fitsimage.set_color_algorithm(stretch_name)

    def motion_cb(self, viewer, button, data_x, data_y):

        # Get the value under the data coordinates
        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            value = viewer.get_data(int(data_x + 0.5), int(data_y + 0.5))

        except Exception:
            value = None

        fits_x, fits_y = data_x, data_y
        try:
            wavelength = int(self.wavelength_data[int(fits_x), int(fits_y)])
        except Exception:
            wavelength = None
            
        if value != None:
            value = round(value, 2)

        if (fits_x > 2048 or fits_x <0) or (fits_y > 2048 or fits_y <0):
            text = "X:    Y:    Value:    Wavelength: "
            self.readout.setText(text)
        else:
            fx = str(int(fits_x))
            fy = str(int(fits_y))
            wavelength = str(wavelength)
            value = str(value)
            text = f"X: {fx : <7} Y: {fy : <7}  Value: {value : <10}  Wavelength: {wavelength}"
            self.readout.setText(text)
        
    def drag_cb(self, viewer, button, data_x, data_y):
        low, high = viewer.get_cut_levels()
        dx = data_x - self.xclick
        dy = data_y - self.yclick
        low = low + dx/10
        high = high + dy/10
        viewer.cut_levels(low, high)
        self.xclick = data_x
        self.yclick = data_y


    def quit(self, *args):
        self.logger.info("Attempting to shut down the application...")
        self.stop_scan()
        time.sleep(1)
        # self.threadpool = False
        QtGui.QApplication.instance().quit()

    def start_scan(self):
        self.scanning = True
        print("scan started...")
        scanner = Scanner(self.scan)
        scanner.signals.load.connect(self.load_file)
        self.threadpool.start(scanner)

    def stop_scan(self):
        self.scanning = False
        print('Scanning stopped.')

    def load_file(self, filepath):
        if '/s/' not in filepath and 'sdata' in filepath:
            filepath = '/s' + filepath
        print(f"Loading {filepath}")
        recenter = False
        if self.fitsimage.get_image() == None:
            recenter = True
        try:
            image = load_data(filepath, logger=self.logger)
            self.fitsimage.set_image(image)
            try:
                self.fitsimage.get_canvas().get_object_by_tag(self.picktag)
                self.fitsimage.get_canvas().delete_object_by_tag(self.picktag)
            except KeyError:
                pass
            if recenter == True:
                self.recenter()
            print(f"Loaded {filepath}")
            self.file_info.setText(f"File: {filepath}")
            self.base_zoom = self.fitsimage.get_zoom()
        except io_fits.FITSError:
            self.file_info.setText(f"File: error loading, empty image")
        except OSError:
            self.file_info.setText(f"File: error loading, possible server crash")

    def open_file(self):
        filters = "Images (s*_*.fits)"
        selected_filter = "Images (s*_*.fits)"
        res = QtGui.QFileDialog.getOpenFileName(self, "Open FITS file",
                                                str(self.nightpath()), filters, selected_filter)
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.load_file(fileName)

    def cuts_popup(self):
        if self.c != None:
            try:
                self.c.dismiss(None)
            except AttributeError:
                pass
        self.c = Cuts(self.logger, self.fitsimage, self.dispname, self.coadds, self.bm)
        self.c.show()

    def sdiff(self):
        if self.sdiff_done == False:
            try:
                ds = self.dispname.read()
                if '/s/' not in ds:
                    ds = '/s' + ds
                image_data = fits.getdata(ds)
                image_header = fits.getheader(ds)
                previous = self.previous_file(ds)
                previous_data = fits.getdata(previous)
            except FileNotFoundError:
                print("No previous file.")
                return
            subtracted = image_data - previous_data
            hdu = fits.PrimaryHDU(header=image_header, data=subtracted)
            filename = self.mathFileNames(ds, previous, '-')
            full_path = filename
            try:
                hdu.writeto(full_path)
            except OSError:
                os.remove(full_path)
                hdu.writeto(full_path)
            self.load_file(str(full_path))
            self.sdiff_done = True
            os.remove(full_path)
        else:
            self.load_file(str(self.dispname.read()))
            self.sdiff_done = False
    
    def previous_file(self, fn):
        prev_dir = str(fn).split("_")[0]
        prev_file = f"{int(str(fn).split('_')[-1].split('.')[0])-1}.fits"
        while len(prev_file)<9:
            prev_file = "0" + prev_file
        return prev_dir + "_" + prev_file
    
    def mathFileNames(self, firstfile, secondfile, operation):
        if "//" in firstfile:
            firstfile = firstfile.split("//")
            firstfile = firstfile[-1]
        else: 
            firstfile = firstfile.split("/")
            firstfile = firstfile[-1]
        if "//" in secondfile:
            secondfile = secondfile.split("//")
            secondfile = secondfile[-1]
        else: 
            secondfile = secondfile.split("/")
            secondfile = secondfile[-1]
        return f'{firstfile} {operation} {secondfile}.fits'
    
    def scan(self, file_callback):
        in_file = str(self.dispname)
        file_callback.emit(in_file)
        while self.scanning:
            cur_file = str(self.dispname)
            if in_file != cur_file:
                print("New Image Detected")
                self.waitForFileToBeUnlocked(cur_file, 1)
                file_callback.emit(cur_file)
                in_file = cur_file
            time.sleep(1)

    # def updateFileCache(self, cachedFiles):
    #     updatedFileList = self.walkDirectory()
    #     if updatedFileList == None:
    #         return False, None, None
    #     filtered = [i for i in updatedFileList if not i in cachedFiles]
    #     cachedFiles = updatedFileList
    #     return len(filtered) > 0, filtered, cachedFiles
        

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
                        	and hdu.data is not None
                            and hdu[0].header.isascii() == True])
                if file_object:
                    locked = False

            except TypeError:
                locked = True

            except OSError:
                locked = True

            # except Warning as w:
            #     print("waiting on header")
            #     locked = True

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
            time.sleep(wait_time)

    # def scan(self, file_callback):
    #     while self.scanning:
    #         if self.display == 1:
    #             print("Taking image")
    #             file_callback.emit(str(self.dispname.read()))
    #             self.waitForZero(0.25)
    #         time.sleep(0.25)

    # def fileIsCurrentlyLocked(self):
    #     print(f'display {self.display} locked')
    #     locked = True
    #     if int(self.display.read()) == 0:
    #         print(f'display {self.display} unlocked')
    #         locked = False
    #     return locked
    
    # def waitForZero(self, wait_time):
    #     while self.display == 1:
    #         time.sleep(wait_time)
    #     print(f'display {self.display}')

    # def waitForFileToBeUnlocked(self, wait_time):
    #     while self.fileIsCurrentlyLocked():
    #         time.sleep(wait_time)

    def nightpath(self):
        dir = str(self.dispname)
        if '/s/' not in dir:
            dir = '/s' + dir
        if "//" in str(dir):
            dir = str(dir).split("//")
            dir = dir[0] + "/"
            nightly = Path(dir)
        else: 
            dir = Path(dir)
            nightly = dir.parent
        return nightly

    def writeFits(self, headerinfo, image_data):
        hdu = fits.PrimaryHDU(header=headerinfo, data=image_data)
        filename = 'procImage.fits'
        full_path = filename
        try:
            hdu.writeto(full_path)
        except OSError:
            os.remove(full_path)
            hdu.writeto(full_path)
        return full_path
    
    def math_popup(self):
        self.m = MathWindow(self.logger, self.fitsimage, self.load_file, self.dispname, self.spec_lastfile)
        self.m.show()

    def window(self):
        self.w = Window(self.logger, self.fitsimage)
        self.w.show()

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

    def zoomIn(self):
        current = self.fitsimage.get_zoom()
        return self.fitsimage.zoom_to(current + 1)
    
    def zoomReset(self):
        return self.fitsimage.zoom_to(self.base_zoom)
    
    def zoomOut(self):
        current = self.fitsimage.get_zoom()
        return self.fitsimage.zoom_to(current - 1)
    
    def setPan(self):
        if self.panning == False:
            self.panning = True
            self.fitsimage.block_callback('cursor-move')
            self.fitsimage.switch_cursor('pan')
            self.clickinfo.setText("Click the image to pan.")
            self.wsetpan.setText("Cancle Pan")
        else:
            self.panning = False
            self.fitsimage.unblock_callback('cursor-move')
            self.fitsimage.switch_cursor('pick')
            self.clickinfo.setText("Left-drag to manually adjust levels.")
            self.wsetpan.setText("Pan")


    def recenter(self):
        self.fitsimage.zoom_fit()
        self.base_zoom = self.fitsimage.get_zoom()
    
        
    def btndown(self, canvas, event, data_x, data_y):
        self.xclick = data_x
        self.yclick = data_y
        if self.panning == True:
            self.fitsimage.set_pan(data_x, data_y)
            self.setPan()
        # self.pickstar(self.fitsimage)


def main():

    app = QtGui.QApplication([])

    # ginga needs a logger.
    # If you don't want to log anything you can create a null logger by
    # using null=True in this call instead of log_stderr=True
    logger = log.get_logger("NIRESSlitViewer", log_stderr=True, level=40)

    w = FitsViewer(logger)
    w.resize(1200,711)
    w.show()
    app.setActiveWindow(w)
    w.raise_()
    w.activateWindow()
    w.start_scan()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
