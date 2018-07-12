"""
Created on 4th July, 2018 from mapclientplugins.meshgeneratorstep.

"""

import string

from opencmiss.zinc.field import Field
from opencmiss.zinc.glyph import Glyph
from opencmiss.zinc.graphics import Graphics
from opencmiss.zinc.node import Node
from scaffoldmaker.scaffoldmaker import Scaffoldmaker
from scaffoldmaker.utils.zinc_utils import *
import numpy as np

from mapclientplugins.meshgeneratorstep.model.meshalignmentmodel import MeshAlignmentModel

STRING_FLOAT_FORMAT = '{:.8g}'


class EcgGraphics(object):
    """
    Framework for generating meshes of a number of types, with mesh type specific options
    """

    def __init__(self):
        self.initialised = False
        pass

    def setRegion(self, region):
        self._region = region
        self._scene = self._region.getScene()
        self.numberInModel = 0

    def updateEEGcolours(self, value):
        fm = self._region.getFieldmodule()
        scene = self._region.getScene()
        displaySurface = scene.findGraphicsByName('displaySurfaces')
        constant = fm.createFieldConstant(value)
        displaySurface.setDataField(constant)

    def findMeshLocation(self):
        fm = self._region.getFieldmodule()
        cache = fm.createFieldcache()
        nodeset = fm.findNodesetByName('nodes')
        node = nodeset.findNodeByIdentifier(9)
        cache.setNode(node)
        search_value = fm.createFieldFiniteElement(3)
        search_value.assignReal(cache, [1, 1, 1.1])

        mesh = fm.findMeshByName('mesh2d')
        mesh_location = fm.createFieldStoredMeshLocation(mesh)
        found_mesh_location = fm.createFieldFindMeshLocation(search_value, mesh_location, mesh)




    def initialiseSpectrum(self, data):
        maximum = -1000000
        minimum = 1000000
        for key in data['cache']:
            array_max = max(data['cache'][key])
            array_min = min(data['cache'][key])
            maximum = max(array_max, maximum)
            minimum = min(array_min, minimum)
        specMod = self._scene.getSpectrummodule()
        spectrum = specMod.findSpectrumByName('eegColourSpectrum')
        spectrum_component = spectrum.getFirstSpectrumcomponent()
        spectrum_component.setRangeMaximum(maximum)
        spectrum_component.setRangeMinimum(minimum)



    def updateEEGnodeColours(self, values):
        fm = self._region.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()
        colour = fm.findFieldByName('colour')
        colour = colour.castFiniteElement()
        nodeset = fm.findNodesetByName('nodes')
        for i in range(self.eegSize):
            node = nodeset.findNodeByIdentifier(self.numberInModel+1+ i)
            cache.setNode(node)
            colour.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, values[(i % (len(values)-1))])
        fm.endChange()

    def initialiseTimeSequences(self, data):
        fm = self._region.getFieldmodule()
        cache = fm.createFieldcache()
        colour = fm.findFieldByName('colour')


    def createEEGPoints(self, region, eeg_group, eeg_coord, i, cache):
        # createEEGPoints creates subgroups of points that use the 'colour' field to change colour

        fm = region.getFieldmodule()
        coordinates = fm.findFieldByName('coordinates')
        coordinates = coordinates.castFiniteElement()
        colour = fm.findFieldByName('colour')
        colour = colour.castFiniteElement()

        # Create templates
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.defineField(colour)
        nodetemplate.setValueNumberOfVersions(colour, -1, Node.VALUE_LABEL_VALUE, 1)

        # Assign values for the new EEG subset
        eeg_group.removeAllNodes()
        eegNode = nodes.createNode(self.numberInModel*2 + i + 1, nodetemplate)
        cache.setNode(eegNode)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, eeg_coord[i])

        cache.setNode(eegNode)

        # Find the z location of mesh for our new node
        mesh = fm.findMeshByName('mesh3d')
        mesh_location = fm.createFieldStoredMeshLocation(mesh)
        found_mesh_location = fm.createFieldFindMeshLocation(coordinates, coordinates, mesh)
        found_mesh_location.setSearchMode(found_mesh_location.SEARCH_MODE_NEAREST)

        [el, coords] = found_mesh_location.evaluateMeshLocation(cache, 3)
        print(f'node number {self.numberInModel + i + 1},  ({self.numberInModel + i + 1 + 25})')
        print(f'element coordinates: {coords}')
        ef = el.getElementfieldtemplate(coordinates, 1)
        cache.setMeshLocation(el, coords)
        [result, global_coords] = coordinates.evaluateReal(cache, 3)
        print(f'global coordinates: {global_coords}')


        # create new node for projection
        eegNode = nodes.createNode(self.numberInModel + i + 1, nodetemplate)
        cache.setNode(eegNode)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, global_coords)
        eeg_group.addNode(eegNode)



    def generateGridPoints(self, point1, point2, number_on_side):

        elements_on_side = 4

        grid_size_x = abs(point1[0] - point2[0])
        grid_size_y = abs(point1[1] - point2[1])

        #scale sides so they have same number of points
        step_size_x = grid_size_x/number_on_side
        step_size_y = grid_size_y/number_on_side

        eeg_coord = []

        for i in range(elements_on_side):
            for j in range(elements_on_side):
                eeg_coord.append([point1[0] + i * step_size_x, point1[1] + j * step_size_y, -1.15])

        # Add the colour bar node
        eeg_coord.append([1, 1.2, 1.2])
        return eeg_coord

    def createGraphics(self, point1=[0, 0], point2=[1, 1], number_on_side=3):
        # Node numbers are generated here
        fm = self._region.getFieldmodule()
        # make graphics
        scene = self._region.getScene()
        scene.beginChange()
        coordinates = fm.findFieldByName('coordinates')

        # Add EEG nodes
        eeg_coord = self.generateGridPoints(point1,point2,number_on_side)
        self.eegSize = len(eeg_coord)

        # Add Spectrum
        spcmod = scene.getSpectrummodule()
        spec = spcmod.getDefaultSpectrum()
        spec.setName('eegColourSpectrum')

        cache = fm.createFieldcache()

        # Initialise all subgroup parameters
        self.ndsg = []  # (node set group)
        self.pointattrList = []
        self.spectrumList = []
        self.nodeColours = []
        finite_element_field = []
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.numberInModel = nodes.getSize()

        #create all EEG subgroups
        colour = fm.createFieldFiniteElement(1)
        colour.setName('colour')
        colour.setManaged(True)
        for i in range(len(eeg_coord)):
            #create new graphics for our subgroup
            self.nodeColours.append(scene.createGraphicsPoints())
            self.nodeColours[i].setFieldDomainType(Field.DOMAIN_TYPE_NODES)
            self.nodeColours[i].setCoordinateField(coordinates)

            # create new subgroup containing our node
            fng = fm.createFieldNodeGroup(fm.findNodesetByName('nodes'))
            self.ndsg.append(fng.getNodesetGroup())
            self.createEEGPoints(self._region, self.ndsg[i], eeg_coord, i, cache)
            self.nodeColours[i].setSubgroupField(fng)

            #set attributes for our new node
            self.nodeColours[i].setSpectrum(spec)
            self.nodeColours[i].setDataField(colour)
            self.pointattrList.append(self.nodeColours[i].getGraphicspointattributes())
            # self.pointattrList[i].setLabelText(1, f'ECG Node {i}')
            # self.pointattrList[i].setLabelOffset([1.5, 1.5, 0])
            self.pointattrList[i].setGlyphShapeType(Glyph.SHAPE_TYPE_SPHERE)
            self.pointattrList[i].setBaseSize([.05, .05, .05])

        # Add a colour bar for the spectrum
        gm = self._scene.getGlyphmodule()
        colour_bar = gm.createGlyphColourBar(spec)
        colour_bar.setLabelDivisions(6)
        self.pointattrList[-1].setGlyph(colour_bar)
        self.pointattrList[-1].setBaseSize([.1, .8, .1])

        scene.endChange()

        #del self.pointattrList[-1]
