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
    ECG Graphics is used as a home for creating and modifying displays to visualise ECG data on the model
    """

    def __init__(self):
        self.initialised = False
        self.number_of_points_on_grid_side = 8
        self.node_corner_list = [0]*4
        self.node_corner_points = [[0, 0],
                                   [0, 0],
                                   [0, 0],
                                   [0, 0]]
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
            node = nodeset.findNodeByIdentifier(self.numberInModel + 1 + i)
            cache.setNode(node)
            colour.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, values[(i % (len(values)-1))])
        fm.endChange()

    def initialiseTimeSequences(self, data):
        fm = self._region.getFieldmodule()
        cache = fm.createFieldcache()
        colour = fm.findFieldByName('colour')

    def generateGridPoints(self, point1, point2, number_on_side):
        # generateGridPoints generates a rectangular grid on along a plane ( x, y or z )
        # using two points defined on the plan and the number of points per side of the grid

        grid_size_x = abs(point1[0] - point2[0])
        grid_size_y = abs(point1[1] - point2[1])

        #scale sides so they have same number of points
        step_size_x = grid_size_x/number_on_side
        step_size_y = grid_size_y/number_on_side

        eeg_coord = []
        for i in range(number_on_side):
            for j in range(number_on_side):
                eeg_coord.append([point1[0] + i * step_size_x, .65, point1[1] + j * step_size_y])

        # Add the colour bar node
        eeg_coord.append([1, 1.2, 1.2])
        return eeg_coord

    def generateGridPoints4(self, p1, p2, p3, p4, number_on_side):
        # We generate our grid points by having 4 points that we assign weightings to
        # based on how far we are away from them.
        # (1 being on the point 0 being in a region the point does not effect the grid)

        ns = number_on_side
        ns1 = number_on_side - 1

        grid_coord = []
        for i in range(number_on_side):
            for j in range(number_on_side):

                # Create our weightings
                w1 = i*j/(ns1**2)
                w2 = (j/ns1) * (ns1 - i)/ns1
                w3 = (i/ns1) * (ns1 - j)/ns1
                w4 = ((ns1-i)*(ns1-j))/(ns1**2)

                # Use our weightings to find coordinates of our new point
                x = p4[0]*w1 + p3[0]*w2 + p2[0]*w3 + p1[0]*w4
                z = p4[1]*w1 + p3[1]*w2 + p2[1]*w3 + p1[1]*w4

                grid_coord.append([x, .6, z])

        # Add the colour bar node
        grid_coord.append([1, 1.2, 1.2])
        return grid_coord


    def moveNode(self, nodeKey, x, z):
        # createEEGPoints creates subgroups of points that use the 'colour' field to change colour

        # Re-aquire openzinc variables
        fm = self._region.getFieldmodule()
        coordinates = fm.findFieldByName('coordinates')
        coordinates = coordinates.castFiniteElement()
        cache = fm.createFieldcache()

        # Create templates
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)

        # Create our first new node for the search
        y_offset = .1
        old_node = nodes.findNodeByIdentifier(nodeKey)
        cache.setNode(old_node)
        [result, old_coords] = coordinates.evaluateReal(cache, 3)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1,
                                      [x, old_coords[1]+y_offset, z])

        # Create our mesh search
        mesh = fm.findMeshByName('mesh3d')
        mesh_location = fm.createFieldStoredMeshLocation(mesh)
        found_mesh_location = fm.createFieldFindMeshLocation(coordinates, coordinates, mesh)
        found_mesh_location.setSearchMode(found_mesh_location.SEARCH_MODE_NEAREST)

        it = 1
        tol = .0002
        old_coords = [x, 3, z]
        new_coords = [10, 10, 10]
        while abs(abs(new_coords[0]) - abs(x)) > tol \
                or abs(abs(new_coords[2]) - abs(z)) > tol:
            # ^^ test if x and y changes are within tolerence

            # Find nearest mesh location
            [el, coords] = found_mesh_location.evaluateMeshLocation(cache, 3)
            cache.setMeshLocation(el, coords)
            [result, new_coords] = coordinates.evaluateReal(cache, 3)

            # Update our search location
            cache.setNode(old_node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1,
                                          [x,
                                           new_coords[1],
                                           z])

            # Break in case we can not converge
            it += 1
            if it > 20:
                print(f'Could not converge on node {nodeKey}')
                break

    def updateGrid(self, new_point_id, new_point):
        index = self.node_corner_list.index(new_point_id)
        self.node_corner_points[index] = new_point

        eeg_coord = self.generateGridPoints4(self.node_corner_points[0],
                                             self.node_corner_points[1],
                                             self.node_corner_points[2],
                                             self.node_corner_points[3],
                                             number_on_side=8)

        self._scene.beginChange()
        for i in range(len(eeg_coord)-1):
            self.moveNode(i + self.numberInModel + 1, eeg_coord[i][0], eeg_coord[i][2])

        self._scene.endChange()

    def createEEGPoints(self, region, eeg_group, eeg_coord, i, cache):
        # createEEGPoints creates subgroups of points that use the 'colour' field to change colour

        # Re-aquire openzinc variables
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
        eegNode = nodes.createNode(self.numberInModel + i + 1, nodetemplate)
        cache.setNode(eegNode)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, eeg_coord[i])

        # Find the z location of mesh for our new node
        mesh = fm.findMeshByName('mesh3d')
        mesh_location = fm.createFieldStoredMeshLocation(mesh)
        found_mesh_location = fm.createFieldFindMeshLocation(coordinates, coordinates, mesh)
        found_mesh_location.setSearchMode(found_mesh_location.SEARCH_MODE_NEAREST)

        # Use the FindMeshLocation function as a solver to see where the model intersects our gird point in a
        # particular axis ( currently using the y axis as this is the surface of the heart you see in open surgery )
        it = 1
        tol = .002
        global_coords = [3, 3, 3]
        while abs(abs(global_coords[0]) - abs(eeg_coord[i][0])) > tol \
                or abs(abs(global_coords[2]) - abs(eeg_coord[i][2])) > tol:
            # ^^ test if x and y changes are within tolerence

            # if i == 16:
            #     break

            # Find nearest mesh location
            [el, coords] = found_mesh_location.evaluateMeshLocation(cache, 3)
            cache.setMeshLocation(el, coords)
            [result, global_coords] = coordinates.evaluateReal(cache, 3)

            # Update our search location
            # eegNodeSearch = nodes.createNode(1000 + self.numberInModel + i + 1, nodetemplate)
            cache.setNode(eegNode)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1,
                                          [eeg_coord[i][0],
                                           global_coords[1],
                                           eeg_coord[i][2]])
            # cache.setNode(eegNodeSearch)

            # Break in case we can not converge
            it += 1
            if it > 15:
                print(f'Could not converge on node {i}')
                break

        # Create the final node with our search coordinates
        eeg_group.addNode(eegNode)

    def createGraphics(self, point1=[0, 1], point2=[0, 0], point3=[1, 1], point4=[.8, 0], number_on_side=8):
        # createGraphics creates our EEG points and assigns them a spectrum for future animation.

        fm = self._region.getFieldmodule()
        scene = self._region.getScene()
        scene.beginChange()
        coordinates = fm.findFieldByName('coordinates')
        cache = fm.createFieldcache()

        #save points
        self.node_corner_points[0] = point1
        self.node_corner_points[1] = point2
        self.node_corner_points[2] = point3
        self.node_corner_points[3] = point4

        # Add EEG nodes
        eeg_coord = self.generateGridPoints4(point1, point2, point3, point4, number_on_side)
        self.eegSize = len(eeg_coord)

        # Add Spectrum
        spcmod = scene.getSpectrummodule()
        spec = spcmod.getDefaultSpectrum()
        spec.setName('eegColourSpectrum')



        # Initialise all subgroup parameters
        self.ndsg = []  # (node set group)
        self.pointattrList = []
        self.spectrumList = []
        self.nodeColours = []
        finite_element_field = []
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.numberInModel = nodes.getSize()

        # Create all EEG subgroups
        colour = fm.createFieldFiniteElement(1)
        colour.setName('colour')
        colour.setManaged(True)
        for i in range(len(eeg_coord)):
            # Create new graphics for our subgroup
            self.nodeColours.append(scene.createGraphicsPoints())
            self.nodeColours[i].setFieldDomainType(Field.DOMAIN_TYPE_NODES)
            self.nodeColours[i].setCoordinateField(coordinates)

            # create new subgroup containing our node
            fng = fm.createFieldNodeGroup(fm.findNodesetByName('nodes'))
            self.ndsg.append(fng.getNodesetGroup())
            self.ndsg[i].removeAllNodes()
            self.createEEGPoints(self._region, self.ndsg[i], eeg_coord, i, cache)
            self.nodeColours[i].setSubgroupField(fng)

            # Set attributes for our new node
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

        # Create node corner list (used to check which nodes are on corner later
        base_node = self.numberInModel
        self.node_corner_list[0] = base_node + 1
        self.node_corner_list[2] = base_node + self.number_of_points_on_grid_side
        self.node_corner_list[1] = base_node + 1 + self.number_of_points_on_grid_side * (self.number_of_points_on_grid_side - 1)
        self.node_corner_list[3] = base_node + self.number_of_points_on_grid_side ** 2
        #del self.pointattrList[-1]
