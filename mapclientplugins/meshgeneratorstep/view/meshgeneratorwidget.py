"""
Created on Aug 29, 2017

@author: Richard Christie
"""


import types
from threeWrapper import BlackfynnGet

from PySide import QtGui, QtCore
from functools import partial

from mapclientplugins.meshgeneratorstep.model.fiducialmarkermodel import FIDUCIAL_MARKER_LABELS
from mapclientplugins.meshgeneratorstep.view.ui_meshgeneratorwidget import Ui_MeshGeneratorWidget
from mapclientplugins.meshgeneratorstep.model.blackfynnECGgraphics import EcgGraphics
from mapclientplugins.meshgeneratorstep.model.blackfynnMesh import Blackfynn_2d_plate

from opencmiss.zinc.node import Node

from opencmiss.utils.maths import vectorops
import time

# imports added for pop up graph
import pyqtgraph as pg
import numpy as np

class MeshGeneratorWidget(QtGui.QWidget):

    def __init__(self, model, parent=None):
        super(MeshGeneratorWidget, self).__init__(parent)
        self._ui = Ui_MeshGeneratorWidget()
        self._ui.setupUi(self)
        self._model = model
        self._model.registerTimeValueUpdateCallback(self._updateTimeValue)
        self._model.registerFrameIndexUpdateCallback(self._updateFrameIndex)
        self._generator_model = model.getGeneratorModel()
        self._plane_model = model.getPlaneModel()

        self._fiducial_marker_model = model.getFiducialMarkerModel()
        self._ui.sceneviewer_widget.setContext(model.getContext())
        self._ui.sceneviewer_widget.setModel(self._plane_model)
        self._model.registerSceneChangeCallback(self._sceneChanged)
        self._doneCallback = None
        self._populateFiducialMarkersComboBox()
        self._marker_mode_active = False
        self._have_images = False
        self.x = 0
        self.y = 0
        # self._populateAnnotationTree()
        meshTypeNames = self._generator_model.getAllMeshTypeNames()
        for meshTypeName in meshTypeNames:
            self._ui.meshType_comboBox.addItem(meshTypeName)
        self._makeConnections()

        self._ui.sceneviewer_widget.foundNode = False
        self._ecg_graphics = model.getEcgGraphics()
        self.blackfynn = BlackfynnGet()
        self.data = {}
        self.blackfynn.loaded = False
        self.y_scaled = 0
        self.pw = None
        self.time = 0

        self._ui.sceneviewer_widget.grid = []

    def _graphicsInitialized(self):
        """
        Callback for when SceneviewerWidget is initialised
        Set custom scene from model
        """
        sceneviewer = self._ui.sceneviewer_widget.getSceneviewer()
        if sceneviewer is not None:
            self._model.loadSettings()
            self._refreshOptions()
            scene = self._model.getScene()
            self._ui.sceneviewer_widget.setScene(scene)
            # self._ui.sceneviewer_widget.setSelectModeAll()
            sceneviewer.setLookatParametersNonSkew([2.0, -2.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
            sceneviewer.setTransparencyMode(sceneviewer.TRANSPARENCY_MODE_SLOW)
            self._autoPerturbLines()
            self._viewAll()

    def _sceneChanged(self):
        sceneviewer = self._ui.sceneviewer_widget.getSceneviewer()
        if sceneviewer is not None:
            if self._have_images:
                self._plane_model.setSceneviewer(sceneviewer)
            scene = self._model.getScene()
            self._ui.sceneviewer_widget.setScene(scene)
            self._autoPerturbLines()

    def _sceneAnimate(self):
        sceneviewer = self._ui.sceneviewer_widget.getSceneviewer()
        if sceneviewer is not None:
            self._model.loadSettings()
            scene = self._model.getScene()
            self._ui.sceneviewer_widget.setScene(scene)
            self._autoPerturbLines()
            self._viewAll()


    def _autoPerturbLines(self):
        """
        Enable scene viewer perturb lines iff solid surfaces are drawn with lines.
        Call whenever lines, surfaces or translucency changes
        """
        sceneviewer = self._ui.sceneviewer_widget.getSceneviewer()
        if sceneviewer is not None:
            sceneviewer.setPerturbLinesFlag(self._generator_model.needPerturbLines())

    def _makeConnections(self):
        self._ui.sceneviewer_widget.graphicsInitialized.connect(self._graphicsInitialized)
        self._ui.done_button.clicked.connect(self._doneButtonClicked)
        self._ui.viewAll_button.clicked.connect(self._viewAll)
        self._ui.meshType_comboBox.currentIndexChanged.connect(self._meshTypeChanged)
        self._ui.deleteElementsRanges_lineEdit.returnPressed.connect(self._deleteElementRangesLineEditChanged)
        self._ui.deleteElementsRanges_lineEdit.editingFinished.connect(self._deleteElementRangesLineEditChanged)
        self._ui.scale_lineEdit.returnPressed.connect(self._scaleLineEditChanged)
        self._ui.scale_lineEdit.editingFinished.connect(self._scaleLineEditChanged)
        self._ui.displayAxes_checkBox.clicked.connect(self._displayAxesClicked)
        self._ui.displayElementNumbers_checkBox.clicked.connect(self._displayElementNumbersClicked)
        self._ui.displayLines_checkBox.clicked.connect(self._displayLinesClicked)
        self._ui.displayNodeDerivatives_checkBox.clicked.connect(self._displayNodeDerivativesClicked)
        self._ui.displayNodeNumbers_checkBox.clicked.connect(self._displayNodeNumbersClicked)
        self._ui.displaySurfaces_checkBox.clicked.connect(self._displaySurfacesClicked)
        self._ui.displaySurfacesExterior_checkBox.clicked.connect(self._displaySurfacesExteriorClicked)
        self._ui.displaySurfacesTranslucent_checkBox.clicked.connect(self._displaySurfacesTranslucentClicked)
        self._ui.displaySurfacesWireframe_checkBox.clicked.connect(self._displaySurfacesWireframeClicked)
        self._ui.displayXiAxes_checkBox.clicked.connect(self._displayXiAxesClicked)
        self._ui.activeModel_comboBox.currentIndexChanged.connect(self._activeModelChanged)
        self._ui.toImage_pushButton.clicked.connect(self._imageButtonClicked)
        self._ui.displayImagePlane_checkBox.clicked.connect(self._displayImagePlaneClicked)
        self._ui.fixImagePlane_checkBox.clicked.connect(self._fixImagePlaneClicked)
        self._ui.timeValue_doubleSpinBox.valueChanged.connect(self._timeValueChanged)
        self._ui.timePlayStop_pushButton.clicked.connect(self._timePlayStopClicked)
        self._ui.frameIndex_spinBox.valueChanged.connect(self._frameIndexValueChanged)
        self._ui.framesPerSecond_spinBox.valueChanged.connect(self._framesPerSecondValueChanged)
        self._ui.timeLoop_checkBox.clicked.connect(self._timeLoopClicked)
        self._ui.displayFiducialMarkers_checkBox.clicked.connect(self._displayFiducialMarkersClicked)
        self._ui.fiducialMarker_comboBox.currentIndexChanged.connect(self._fiducialMarkerChanged)
        self._ui.submitButton.clicked.connect(self._submitClicked)
        self._ui.displayEEGAnimation_checkBox.clicked.connect(self._EEGAnimationClicked)
        self._ui.pushButton.clicked.connect(self._exportWebGLJson)
        # self._ui.treeWidgetAnnotation.itemSelectionChanged.connect(self._annotationSelectionChanged)
        # self._ui.treeWidgetAnnotation.itemChanged.connect(self._annotationItemChanged)

        # currently not able to loop it (will have to do later
        self._ui.LG3.clicked.connect(self._lg3)
        self._ui.LG4.clicked.connect(self._lg4)
        self._ui.LG10.clicked.connect(self._lg10)
        self._ui.LG3.setText('')
        self._ui.LG3.setStyleSheet("background-color: rgba(255, 255, 255, 0);")
        self._ui.LG4.setText('')
        self._ui.LG4.setStyleSheet("background-color: rgba(255, 255, 255, 0);")
        self._ui.LG10.setText('')
        self._ui.LG10.setStyleSheet("background-color: rgba(255, 255, 255, 0);")



    def _fiducialMarkerChanged(self):
        self._fiducial_marker_model.setActiveMarker(self._ui.fiducialMarker_comboBox.currentText())

    def _displayFiducialMarkersClicked(self):
        self._fiducial_marker_model.setDisplayFiducialMarkers(self._ui.displayFiducialMarkers_checkBox.isChecked())

    def _populateFiducialMarkersComboBox(self):
        self._ui.fiducialMarker_comboBox.addItems(FIDUCIAL_MARKER_LABELS)

    def _createFMAItem(self, parent, text, fma_id):
        item = QtGui.QTreeWidgetItem(parent)
        item.setText(0, text)
        item.setData(0, QtCore.Qt.UserRole + 1, fma_id)
        item.setCheckState(0, QtCore.Qt.Unchecked)
        item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsTristate)

        return item

    def _populateAnnotationTree(self):
        tree = self._ui.treeWidgetAnnotation
        tree.clear()
        rsh_item = self._createFMAItem(tree, 'right side of heart', 'FMA_7165')
        self._createFMAItem(rsh_item, 'ventricle', 'FMA_7098')
        self._createFMAItem(rsh_item, 'atrium', 'FMA_7096')
        self._createFMAItem(rsh_item, 'auricle', 'FMA_7218')
        lsh_item = self._createFMAItem(tree, 'left side of heart', 'FMA_7166')
        self._createFMAItem(lsh_item, 'ventricle', 'FMA_7101')
        self._createFMAItem(lsh_item, 'atrium', 'FMA_7097')
        self._createFMAItem(lsh_item, 'auricle', 'FMA_7219')
        apex_item = self._createFMAItem(tree, 'apex of heart', 'FMA_7164')
        vortex_item = self._createFMAItem(tree, 'vortex of heart', 'FMA_84628')

        self._ui.treeWidgetAnnotation.addTopLevelItem(rsh_item)
        self._ui.treeWidgetAnnotation.addTopLevelItem(lsh_item)
        self._ui.treeWidgetAnnotation.addTopLevelItem(apex_item)
        self._ui.treeWidgetAnnotation.addTopLevelItem(vortex_item)

    def getModel(self):
        return self._model

    def registerDoneExecution(self, doneCallback):
        self._doneCallback = doneCallback

    def _updateUi(self):
        if self._have_images:
            frame_count = self._plane_model.getFrameCount()
            self._ui.numFramesValue_label.setText("{0}".format(frame_count))
            self._ui.frameIndex_spinBox.setMaximum(frame_count)
            self._ui.timeValue_doubleSpinBox.setMaximum(frame_count / self._model.getFramesPerSecond())
        else:
            self._generator_model.disableAlignment()
            self._plane_model.disableAlignment()
            self._ui.alignment_groupBox.setVisible(False)
            self._ui.fiducialMarkers_groupBox.setVisible(False)
            self._ui.video_groupBox.setVisible(False)
            self._ui.displayImagePlane_checkBox.setVisible(False)
            self._ui.displayFiducialMarkers_checkBox.setVisible(False)

    def setImageInfo(self, image_info):
        self._plane_model.setImageInfo(image_info)
        self._have_images = image_info is not None
        self._updateUi()

    def _doneButtonClicked(self):
        self._ui.dockWidget.setFloating(False)
        self._model.done()
        self._model = None
        self._doneCallback()

    def _imageButtonClicked(self):
        sceneviewer = self._ui.sceneviewer_widget.getSceneviewer()
        normal, up, offset = self._plane_model.getPlaneInfo()
        _, current_lookat_pos = sceneviewer.getLookatPosition()
        _, current_eye_pos = sceneviewer.getEyePosition()
        view_distance = vectorops.magnitude(vectorops.sub(current_eye_pos, current_lookat_pos))
        eye_pos = vectorops.add(vectorops.mult(normal, view_distance), offset)
        lookat_pos = offset
        sceneviewer.setLookatParametersNonSkew(eye_pos, lookat_pos, up)

    def _updateTimeValue(self, value):
        self._ui.timeValue_doubleSpinBox.blockSignals(True)
        frame_count = self._plane_model.getFrameCount()
        max_time_value = frame_count / self._ui.framesPerSecond_spinBox.value()
        self.time = self._model._current_time

        if value > max_time_value:
            self._ui.timeValue_doubleSpinBox.setValue(max_time_value)
            self._timePlayStopClicked()
        else:
            self._ui.timeValue_doubleSpinBox.setValue(value)
            if self.pw is not None:
                self.line.setValue(round(value, 3)) # adjust time marker
            if self._ui.displayEEGAnimation_checkBox.isChecked() and self.data is not False:
                pass
                # use model to update colours
                #self.updateAllNodes(value)
            #self.updatePlate(value)
        self._ui.timeValue_doubleSpinBox.blockSignals(False)

    def updateAllNodes(self, time):
        colours_at_current_time = []
        for key in self.data['scaled']:
            colours_at_current_time.append(self.data['scaled'][key][self.currentFrame(time)])
        self._ecg_graphics.updateEEGnodeColours(colours_at_current_time)

    def updatePlate(self, time):
        colours_at_current_time = []
        for key in self.data['scaled']:
            colours_at_current_time.append(self.data['scaled'][key][self.currentFrame(time)])
        self.updatePlateColoursTemporary(colours_at_current_time)

    def updatePlateColoursTemporary(self, values):
        reg = self._generator_model._region.findChildByName('ecg_plane')
        fm = reg.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()
        colour = fm.findFieldByName('colour2')
        colour = colour.castFiniteElement()
        nodeset = fm.findNodesetByName('nodes')
        for i in range(10000, 10064):
            node = nodeset.findNodeByIdentifier(i)
            cache.setNode(node)
            colour.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, values[(i % (len(values)-1))])
        fm.endChange()

    def scaleCacheData(self):
        tempDict = {}
        for i, key in enumerate(self.data['cache']):
            tempDict[str(i)] = self.scaleData(key)
        self.data['scaled'] = tempDict

    def scaleData(self, key):
        numFrames = self._plane_model.getFrameCount()
        y = np.array(self.data['cache'][key])
        x = np.linspace(0, 16, len(y))
        xterp = np.linspace(0, 16, numFrames)
        yterp = np.interp(xterp, x, y)
        return yterp

    def initialiseSpectrum(self, data):
        maximum = -1000000
        minimum = 1000000
        for key in data['cache']:
            array_max = max(data['cache'][key])
            array_min = min(data['cache'][key])
            maximum = max(array_max, maximum)
            minimum = min(array_min, minimum)
        scene = self._generator_model._region.findChildByName('ecg_plane').getScene()
        specMod = scene.getSpectrummodule()
        spectrum = specMod.findSpectrumByName('eegColourSpectrum2')
        spectrum_component = spectrum.getFirstSpectrumcomponent()
        spectrum_component.setRangeMaximum(maximum)
        spectrum_component.setRangeMinimum(minimum)

    def _EEGAnimationClicked(self):
        if self.data and self._ecg_graphics.initialised is False:

            self.scaleCacheData()

            self._ecg_graphics.setRegion(self._generator_model._region)

            # create our ecg graphics if we have defined a box
            if len(self._ui.sceneviewer_widget.grid) >= 4:
                self._ecg_graphics.plane_normal = self._ui.sceneviewer_widget.plane_normal
                self._ecg_graphics.createGraphics(new_points=True,
                                                  point1=self._ui.sceneviewer_widget.grid[0],
                                                  point2=self._ui.sceneviewer_widget.grid[1],
                                                  point3=self._ui.sceneviewer_widget.grid[2],
                                                  point4=self._ui.sceneviewer_widget.grid[3])

                self._ui.sceneviewer_widget.grid = []
            else:
                if self._ecg_graphics.settingsLoaded:
                    self._ecg_graphics.createGraphics()
                else:
                    self._ecg_graphics.createGraphics(new_points=True)
            self._ecg_graphics.initialiseSpectrum(self.data)
            self._ecg_graphics.initialised = True

        else:
            self._ecg_graphics.clearAll()
            self._ecg_graphics.__init__()
            self._ecg_graphics.initialised = False




    def currentFrame(self, value):
        frame_count = self._plane_model.getFrameCount()
        frame_vals = np.linspace(0, 16, frame_count)
        currentFrame = (np.abs(frame_vals - value)).argmin()
        return currentFrame

    def find_nearest(array, value):
        # fin_nearets() Find the index of the nearest value in an array
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
            return array[idx - 1]
        else:
            return array[idx]

    def _updateFrameIndex(self, value):
        self._ui.frameIndex_spinBox.blockSignals(True)
        self._ui.frameIndex_spinBox.setValue(value)
        self._ui.frameIndex_spinBox.blockSignals(False)

    def _timeValueChanged(self, value):
        self._model.setTimeValue(value)

    def _timeDurationChanged(self, value):
        self._model.setTimeDuration(value)

    def _timePlayStopClicked(self):
        play_text = 'Play'
        stop_text = 'Stop'
        current_text = self._ui.timePlayStop_pushButton.text()
        if current_text == play_text:
            self._ui.timePlayStop_pushButton.setText(stop_text)
            self._model.play()
        else:
            self._ui.timePlayStop_pushButton.setText(play_text)
            self._model.stop()

    def _timeLoopClicked(self):
        self._model.setTimeLoop(self._ui.timeLoop_checkBox.isChecked())

    def _frameIndexValueChanged(self, value):
        self._model.setFrameIndex(value)

    def _framesPerSecondValueChanged(self, value):
        self._model.setFramesPerSecond(value)
        self._ui.timeValue_doubleSpinBox.setMaximum(self._plane_model.getFrameCount()/value)

    def _fixImagePlaneClicked(self):
        self._plane_model.setImagePlaneFixed(self._ui.fixImagePlane_checkBox.isChecked())

    def _submitClicked(self):
    # submitClicked initialises all the blackfynn functionality and updates login fields.
        if self._ui.api_key.displayText() != 'API Key' and self._ui.api_secret.text() != '***************************':
            self.pw = pg.plot(title='Blackfynn electrode graph',
                              labels={'left': f'EEG value of node', 'bottom': 'time in seconds'})
            self._ui.Login_groupBox.setTitle(QtGui.QApplication.translate("MeshGeneratorWidget", "Login details saved, click on a node to open graphs", None,
                                                                       QtGui.QApplication.UnicodeUTF8))
            self.initialiseBlackfynnData()
            self._ui.api_secret.setText('***************************')
            self.blackfynn.loaded = True


    def initialiseBlackfynnData(self):
        # self.blackfynn.api_key = self._ui.api_key.text()  <- commented so that we do not have to enter key each time
        # self.blackfynn.api_secret = self._ui.api_secret.text()
        self.blackfynn.set_api_key_login()
        self.blackfynn.set_params(channels='LG4', window_from_start=16) # need to add dataset selection
        self.data = self.blackfynn.get()
        self.updatePlot(4)
        self.scaleCacheData()
        self.initialiseSpectrum(self.data)

    def updatePlot(self, key):

        try:
            self.data['cache'][f'LG{key}']
        except KeyError:
            print('ERROR: selected data could not be found')
            self.pw.plot(title='Error in data collection')
            return
        self.pw.clear()

        self.pw.plot(self.data['x'],
                     self.data['cache'][f'LG{key}'],
                     pen='b',
                     title=f'EEG values from {key} LG{key}',
                     )
        self.pw.setTitle(f'EEG values from {key} (LG{key})')
        self.line = self.pw.addLine(x=self.time,
                                    pen='r')  # show current time


    # For linking each EEG node
    def _lg3(self):
        self.updatePlot(3)
    def _lg4(self):
        self.updatePlot(4)
    def _lg10(self):
        self.updatePlot(10)

    def EEGSelectionDisplay(self, key):
        # For selecting EEG (brain) points
        print(f'key {key} clicked!')
        if self.data:
            self.pw.clear()
            self.pw.plot(self.data['x'], self.data['cache'][f'LG{key}'], pen='b',
                         title=f'EEG values from {key} (LG{key})',
                         labels={'left': f'EEG value of node LG{key}', 'bottom': 'time in seconds'})
            self.line = self.pw.addLine(x=self.time, pen='r')  # show current time



    def _displayImagePlaneClicked(self):
        self._plane_model.setImagePlaneVisible(self._ui.displayImagePlane_checkBox.isChecked())

    def _activeModelChanged(self, index):
        if index == 0:
            self._ui.sceneviewer_widget.setModel(self._plane_model)
        else:
            self._ui.sceneviewer_widget.setModel(self._generator_model)

    def _meshTypeChanged(self, index):
        meshTypeName = self._ui.meshType_comboBox.itemText(index)
        self._generator_model.setMeshTypeByName(meshTypeName)
        self._refreshMeshTypeOptions()
        #self._ecg_graphics.createGraphics()



    def _meshTypeOptionCheckBoxClicked(self, checkBox):
        self._generator_model.setMeshTypeOption(checkBox.objectName(), checkBox.isChecked())

    def _meshTypeOptionLineEditChanged(self, lineEdit):
        self._generator_model.setMeshTypeOption(lineEdit.objectName(), lineEdit.text())
        finalValue = self._generator_model.getMeshTypeOption(lineEdit.objectName())
        lineEdit.setText(str(finalValue))

    def _refreshMeshTypeOptions(self):
        layout = self._ui.meshTypeOptions_frame.layout()
        # remove all current mesh type widgets
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
              child.widget().deleteLater()
        optionNames = self._generator_model.getMeshTypeOrderedOptionNames()
        for key in optionNames:
            value = self._generator_model.getMeshTypeOption(key)
            # print('key ', key, ' value ', value)
            if type(value) is bool:
                checkBox = QtGui.QCheckBox(self._ui.meshTypeOptions_frame)
                checkBox.setObjectName(key)
                checkBox.setText(key)
                checkBox.setChecked(value)
                callback = partial(self._meshTypeOptionCheckBoxClicked, checkBox)
                checkBox.clicked.connect(callback)
                layout.addWidget(checkBox)
            else:
                label = QtGui.QLabel(self._ui.meshTypeOptions_frame)
                label.setObjectName(key)
                label.setText(key)
                layout.addWidget(label)
                lineEdit = QtGui.QLineEdit(self._ui.meshTypeOptions_frame)
                lineEdit.setObjectName(key)
                lineEdit.setText(str(value))
                callback = partial(self._meshTypeOptionLineEditChanged, lineEdit)
                lineEdit.returnPressed.connect(callback)
                lineEdit.editingFinished.connect(callback)
                layout.addWidget(lineEdit)

    def _refreshOptions(self):
        self._ui.identifier_label_2.setText('Identifier:  ' + self._model.getIdentifier())
        self._ui.deleteElementsRanges_lineEdit.setText(self._generator_model.getDeleteElementsRangesText())
        self._ui.scale_lineEdit.setText(self._generator_model.getScaleText())
        self._ui.displayAxes_checkBox.setChecked(self._generator_model.isDisplayAxes())
        self._ui.displayElementNumbers_checkBox.setChecked(self._generator_model.isDisplayElementNumbers())
        self._ui.displayLines_checkBox.setChecked(self._generator_model.isDisplayLines())
        self._ui.displayNodeDerivatives_checkBox.setChecked(self._generator_model.isDisplayNodeDerivatives())
        self._ui.displayNodeNumbers_checkBox.setChecked(self._generator_model.isDisplayNodeNumbers())
        self._ui.displaySurfaces_checkBox.setChecked(self._generator_model.isDisplaySurfaces())
        self._ui.displaySurfacesExterior_checkBox.setChecked(self._generator_model.isDisplaySurfacesExterior())
        self._ui.displaySurfacesTranslucent_checkBox.setChecked(self._generator_model.isDisplaySurfacesTranslucent())
        self._ui.displaySurfacesWireframe_checkBox.setChecked(self._generator_model.isDisplaySurfacesWireframe())
        self._ui.displayXiAxes_checkBox.setChecked(self._generator_model.isDisplayXiAxes())
        self._ui.displayImagePlane_checkBox.setChecked(self._plane_model.isDisplayImagePlane())
        self._ui.displayFiducialMarkers_checkBox.setChecked(self._fiducial_marker_model.isDisplayFiducialMarkers())
        self._ui.fixImagePlane_checkBox.setChecked(self._plane_model.isImagePlaneFixed())
        self._ui.framesPerSecond_spinBox.setValue(self._model.getFramesPerSecond())
        self._ui.timeLoop_checkBox.setChecked(self._model.isTimeLoop())
        index = self._ui.meshType_comboBox.findText(self._generator_model.getMeshTypeName())
        self._ui.meshType_comboBox.blockSignals(True)
        self._ui.meshType_comboBox.setCurrentIndex(index)
        self._ui.meshType_comboBox.blockSignals(False)
        index = self._ui.fiducialMarker_comboBox.findText(self._fiducial_marker_model.getActiveMarker())
        self._ui.fiducialMarker_comboBox.blockSignals(True)
        self._ui.fiducialMarker_comboBox.setCurrentIndex(0 if index == -1 else index)
        self._ui.fiducialMarker_comboBox.blockSignals(False)
        self._refreshMeshTypeOptions()

    def _deleteElementRangesLineEditChanged(self):
        self._generator_model.setDeleteElementsRangesText(self._ui.deleteElementsRanges_lineEdit.text())
        self._ui.deleteElementsRanges_lineEdit.setText(self._generator_model.getDeleteElementsRangesText())

    def _scaleLineEditChanged(self):
        self._generator_model.setScaleText(self._ui.scale_lineEdit.text())
        self._ui.scale_lineEdit.setText(self._generator_model.getScaleText())

    def _displayAxesClicked(self):
        self._generator_model.setDisplayAxes(self._ui.displayAxes_checkBox.isChecked())
        # for testing, we delete our ecg nodes and reload the entire mesh
        self._ecg_graphics.deleteAll()
        self._meshTypeChanged(9)
        self._meshTypeChanged(10)

        # prepare data
        self.scaleCacheData()
        self._ecg_graphics.initialiseSpectrum(self.data)
        ECGmatrix = []
        for key in self.data['cache']:
            ECGmatrix.append(self.data['cache'][key][0::10])
        for i in range(len(ECGmatrix)):
            ECGmatrix[i].append(ECGmatrix[i][-1])
        ECGtimes = np.linspace(0, 1, len(ECGmatrix[:][0]))

        # clear all of the current mesh data by going to a mesh with nothing in it
        self._generator_model.deleteAll()
        # self._meshTypeChanged(3)

        # create our new mesh with the Blackfynn_2d_plate class
        pm = Blackfynn_2d_plate(self._generator_model._region, self._ecg_graphics.node_coordinate_list)
        pm.ECGtimes = ECGtimes.tolist()
        pm.ECGcoloursMatrix = ECGmatrix
        pm.generateMesh()
        pm.drawMesh()


    def _exportWebGLJson(self):
        '''
            Export graphics into JSON formats. Returns an array containing the
       string buffers for each export
            '''

        try:
            self.data
            ECGmatrix = []
            for key in self.data['cache']:
                ECGmatrix.append(self.data['cache'][key][0::10])
            for i in range(len(ECGmatrix)):
                ECGmatrix[i].append(ECGmatrix[i][-1])
            ECGtimes = np.linspace(0, 1, len(ECGmatrix[:][0]))

            ecg_region = self._generator_model._region.findChildByName('ecg_plane')
            scene = ecg_region.getScene()
            sceneSR = scene.createStreaminformationScene()
            sceneSR.setIOFormat(sceneSR.IO_FORMAT_THREEJS)
            sceneSR.setInitialTime(ECGtimes[0])
            sceneSR.setFinishTime(ECGtimes[-1])
            sceneSR.setNumberOfTimeSteps(len(ECGtimes))
            sceneSR.setOutputTimeDependentColours(1)

            # Get the total number of graphics in a scene/region that can be exported
            number = sceneSR.getNumberOfResourcesRequired()
            resources = []
            # Write out each graphics into a json file which can be rendered with our
            # WebGL script
            for i in range(number):
                resources.append(sceneSR.createStreamresourceMemory())
            scene.write(sceneSR)
            # Write out each resource into their own file

            buffer = [resources[i].getBuffer()[1] for i in range(number)]

            for i, content in enumerate(buffer):
                f = open(f'webGLExport{i+1}.json', 'w')
                f.write(content)
        except:
            pass

    def _displayElementNumbersClicked(self):
        self._generator_model.setDisplayElementNumbers(self._ui.displayElementNumbers_checkBox.isChecked())


        self.exportWebGLJson([1,2,3,4])
    def _displayLinesClicked(self):
        self._generator_model.setDisplayLines(self._ui.displayLines_checkBox.isChecked())
        self._autoPerturbLines()

    def _displayNodeDerivativesClicked(self):
        self._generator_model.setDisplayNodeDerivatives(self._ui.displayNodeDerivatives_checkBox.isChecked())

    def _displayNodeNumbersClicked(self):
        self._generator_model.setDisplayNodeNumbers(self._ui.displayNodeNumbers_checkBox.isChecked())

    def _displaySurfacesClicked(self):
        self._generator_model.setDisplaySurfaces(self._ui.displaySurfaces_checkBox.isChecked())
        self._autoPerturbLines()

    def _displaySurfacesExteriorClicked(self):
        self._generator_model.setDisplaySurfacesExterior(self._ui.displaySurfacesExterior_checkBox.isChecked())

    def _displaySurfacesTranslucentClicked(self):
        self._generator_model.setDisplaySurfacesTranslucent(self._ui.displaySurfacesTranslucent_checkBox.isChecked())
        self._autoPerturbLines()

    def _displaySurfacesWireframeClicked(self):
        self._generator_model.setDisplaySurfacesWireframe(self._ui.displaySurfacesWireframe_checkBox.isChecked())

    def _displayXiAxesClicked(self):
        self._generator_model.setDisplayXiAxes(self._ui.displayXiAxes_checkBox.isChecked())

    def _annotationItemChanged(self, item):
        print(item.text(0))
        print(item.data(0, QtCore.Qt.UserRole + 1))

    def _viewAll(self):
        """
        Ask sceneviewer to show all of scene.
        """
        if self._ui.sceneviewer_widget.getSceneviewer() is not None:
            self._ui.sceneviewer_widget.viewAll()

    def keyPressEvent(self, event):
        if event.modifiers() & QtCore.Qt.CTRL and QtGui.QApplication.mouseButtons() == QtCore.Qt.NoButton:
            self._marker_mode_active = True
            self._ui.sceneviewer_widget._model = self._fiducial_marker_model
            self._original_mousePressEvent = self._ui.sceneviewer_widget.mousePressEvent
            self._ui.sceneviewer_widget.original_mousePressEvent = self._ui.sceneviewer_widget.mousePressEvent
            self._ui.sceneviewer_widget.plane_model_temp = self._plane_model

            self._ui.sceneviewer_widget._calculatePointOnPlane = types.MethodType(_calculatePointOnPlane, self._ui.sceneviewer_widget)
            self._ui.sceneviewer_widget.mousePressEvent = types.MethodType(mousePressEvent, self._ui.sceneviewer_widget)
            #self._ui.sceneviewer_widget.foundNode = False
            self._model.printLog()

    def keyReleaseEvent(self, event):
        if self._marker_mode_active:
            self._marker_mode_active = False
            self._ui.sceneviewer_widget._model = self._plane_model
            self._ui.sceneviewer_widget._calculatePointOnPlane = None
            self._ui.sceneviewer_widget.mousePressEvent = self._original_mousePressEvent
            if self._ui.sceneviewer_widget.foundNode and len(self._ui.sceneviewer_widget.grid) is 2:
                #self.updatePlot(self._ui.sceneviewer_widget.nodeKey) # updates plot if a node is clicked
                if self._ui.sceneviewer_widget.nodeKey in self._ecg_graphics.node_corner_list:
                    self._ecg_graphics.updateGrid(self._ui.sceneviewer_widget.nodeKey,  self._ui.sceneviewer_widget.grid[1])
                self._ui.sceneviewer_widget.foundNode = False




def mousePressEvent(self, event):
    if self._active_button != QtCore.Qt.NoButton:
        return

    if (event.modifiers() & QtCore.Qt.CTRL) and event.button() == QtCore.Qt.LeftButton:
        point_on_plane = self._calculatePointOnPlane(event.x(), event.y())
        print('Location of click (x,y): (' + str(event.x()) + ', ' + str(event.y()) +')')
        node = self.getNearestNode(event.x(), event.y())
        if node.isValid():
            print(f'node {node.getIdentifier()} was clicked')
            self.foundNode = True
            self.nodeKey = node.getIdentifier()
            self.node = node
            self.grid = []
        # return sceneviewers 'mouspressevent' function to its version for navigation

        self._model = self.plane_model_temp
        self._calculatePointOnPlane = None
        self.mousePressEvent = self.original_mousePressEvent

        if len(self.grid) > 2 and self.foundNode:
            self.grid = []
            self.foundNode = False
        self.grid.append(point_on_plane)

        if len(self.grid) > 4:
            self.grid = []

    return [event.x(), event.y()]




def _calculatePointOnPlane(self, x, y):
    from opencmiss.utils.maths.algorithms import calculateLinePlaneIntersection

    far_plane_point = self.unproject(x, -y, -1.0)
    near_plane_point = self.unproject(x, -y, 1.0)
    plane_point, plane_normal = self._model.getPlaneDescription()
    point_on_plane = calculateLinePlaneIntersection(near_plane_point, far_plane_point, plane_point, plane_normal)
    self.plane_normal = plane_normal
    print(point_on_plane)
    print(f'normal: {plane_normal}')
    return point_on_plane


