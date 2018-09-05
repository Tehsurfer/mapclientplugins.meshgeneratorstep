"""
Created on 4th July, 2018 from mapclientplugins.meshgeneratorstep.

"""

import string

from opencmiss.utils.zinc import createFiniteElementField
from opencmiss.zinc.field import Field
from opencmiss.zinc.glyph import Glyph
import opencmiss.zinc.scenecoordinatesystem as Scenecoordinatesystem
from opencmiss.zinc.graphics import Graphics
from opencmiss.zinc.node import Node
from scaffoldmaker.scaffoldmaker import Scaffoldmaker
from scaffoldmaker.utils.zinc_utils import *
import numpy as np
import time

from mapclientplugins.meshgeneratorstep.model.blackfynnMesh import Blackfynn_2d_plate
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
        self.node_corner_points = [[0, 0, 0],
                                   [0, 0, 0],
                                   [0, 0, 0],
                                   [0, 0, 0]]
        self.settingsLoaded = False
        self.plane_normal = [0, 1, 0]
        self.node_coordinate_list = []
        self._child_region = None
        pass

    def getSettings(self):
        if self.node_corner_list[0] is not 0:
            return self.node_corner_points

    def setSettings(self, settings):
        self.node_corner_points = settings
        self.settingsLoaded = True

    def setRegion(self, region):
        self._region = region
        self._scene = self._region.getScene()
        self._child_region = self._region.createChild('ecg_plane')
        self.numberInModel = 0

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

    def updatePlateColours(self, values):
        self.plateMesh.updatePlateColours(values)

    def initialiseTimeSequences(self, data):
        fm = self._region.getFieldmodule()
        cache = fm.createFieldcache()
        colour = fm.findFieldByName('colour')

    def clearAll(self):
        fm = self._region.getFieldmodule()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        for i in range(self.eegSize):
            node = nodes.findNodeByIdentifier(self.numberInModel+i+1)
            nodes.destroyNode(node)

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
                transformed_x = point1[0] + i * step_size_x

                eeg_coord.append([point1[0] + i * step_size_x, .65, point1[1] + j * step_size_y])

        eeg_coord2 = []
        for i in range(len(eeg_coord)):
            eeg_coord2.append(np.cross(eeg_coord[i], self.plane_normal))

        return eeg_coord2

    def generateGridPoints4(self, number_on_side):
        # We generate our grid points by having 4 points that we assign weightings to
        # based on how far we are away from them.
        # (1 being on the point 0 being in a region the point does not effect the grid)

        p1 = self.node_corner_points[0]
        p2 = self.node_corner_points[1]
        p3 = self.node_corner_points[2]
        p4 = self.node_corner_points[3]

        ns = number_on_side
        ns1 = number_on_side - 1

        plane_normal_offset = .4  # For offsetting the solver to solve from outside the mesh -> on it

        grid_coord = []
        for i in range(number_on_side):
            for j in range(number_on_side):

                # Create our weightings (since we are setting points in a ccwise fashion our diagonal is w3
                w1 = i*j/(ns1**2)
                w2 = (j/ns1) * (ns1 - i)/ns1
                w4 = (i/ns1) * (ns1 - j)/ns1  # The 'bottom left' point, p4
                w3 = ((ns1-i)*(ns1-j))/(ns1**2)  # The diagonal point, p3

                # Use our weightings to find coordinates of our new point
                x = p4[0] * w1 + p3[0] * w2 + p2[0] * w3 + p1[0] * w4
                y = p4[1] * w1 + p3[1] * w2 + p2[1] * w3 + p1[1] * w4
                z = p4[2] * w1 + p3[2] * w2 + p2[2] * w3 + p1[2] * w4

                grid_coord.append([x, y, z])
        plane_norm = np.array(self.plane_normal)
        eeg_coord2 = []
        for i in range(len(grid_coord)):
            #projected_point = np.cross(grid_coord[i], plane_norm)
            shifted_point = grid_coord[i] + plane_norm*plane_normal_offset
            eeg_coord2.append(shifted_point.tolist())

        return eeg_coord2



    def moveNode(self, nodeKey, cache, tol=.01, max_iterations=20):
        # createEEGPoints creates subgroups of points that use the 'colour' field to change colour

        # Re-aquire openzinc variables
        fm = self._region.getFieldmodule()
        coordinates = fm.findFieldByName('coordinates')
        coordinates = coordinates.castFiniteElement()

        # Create templates
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)

        # Create our first new node for the search
        plane_normal_offset = .1
        old_node = nodes.findNodeByIdentifier(nodeKey)
        cache.setNode(old_node)
        [result, old_coords] = coordinates.evaluateReal(cache, 3)

        plane_norm = np.array(self.plane_normal)
        shifted_point = old_coords + plane_norm * plane_normal_offset

        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, shifted_point.tolist())

        # Create our mesh search
        mesh = fm.findMeshByName('mesh3d')
        mesh_location = fm.createFieldStoredMeshLocation(mesh)
        found_mesh_location = fm.createFieldFindMeshLocation(coordinates, coordinates, mesh)
        found_mesh_location.setSearchMode(found_mesh_location.SEARCH_MODE_NEAREST)

        it = 1
        old_coords = shifted_point
        test_coords = [10, 10, 10]
        new_coords = shifted_point
        start3 = time.clock()
        while abs(np.linalg.norm(np.dot(test_coords, plane_norm) - np.dot(new_coords, plane_norm))) > tol:
            # ^^ test if x and y changes are within tolerence
            end3 = time.clock()
            # Find nearest mesh location
            start = time.clock()
            [el, coords] = found_mesh_location.evaluateMeshLocation(cache, 3)
            end = time.clock()
            cache.setMeshLocation(el, coords)
            [result, mesh_coords] = coordinates.evaluateReal(cache, 3)

            # Update our search location
            start2 = time.clock()
            new_coords = old_coords + np.dot(mesh_coords - old_coords, plane_norm)*plane_norm
            cache.setNode(old_node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, new_coords.tolist())
            end2 = time.clock()

            test_coords = old_coords
            old_coords = new_coords

            # Break in case we can not converge
            it += 1
            if it > max_iterations:
                print(f'Could not converge on node {nodeKey}')
                break



            print(f'Find mesh took: {end-start}')
            print(f'Update search took: {end2-start2}')
            print(f'Stop evaluation took: {end3-start3}')
            start3 = time.clock()
        self.node_coordinate_list.append(new_coords)
        print(f'Node {nodeKey} was solved in {it-1} iterations' )

    def nudgeNode(self, nodeKey, eegCoords):
        # createEEGPoints creates subgroups of points that use the 'colour' field to change colour

        tol = .01
        max_iterations = 10

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

        # Create our new node for the search
        plane_normal_offset = .15
        old_node = nodes.findNodeByIdentifier(nodeKey)
        cache.setNode(old_node)

        # Update our nodes location
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, eegCoords)
        [result, old_coords] = coordinates.evaluateReal(cache, 3)

        plane_norm = np.array(self.plane_normal)
        shifted_point = old_coords + plane_norm * plane_normal_offset

        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, shifted_point.tolist())

        # Create our mesh search
        mesh = fm.findMeshByName('mesh3d')
        mesh_location = fm.createFieldStoredMeshLocation(mesh)
        found_mesh_location = fm.createFieldFindMeshLocation(coordinates, coordinates, mesh)
        found_mesh_location.setSearchMode(found_mesh_location.SEARCH_MODE_NEAREST)

        it = 1
        old_coords = shifted_point
        test_coords = [10, 10, 10]
        new_coords = shifted_point
        while abs(np.linalg.norm(np.dot(test_coords, plane_norm) - np.dot(new_coords, plane_norm))) > tol:
            # ^^ test if x and y changes are within tolerence
            # Find nearest mesh location

            [element, local_coords] = found_mesh_location.evaluateMeshLocation(cache, 3)
            cache.setMeshLocation(element, local_coords)
            [result, mesh_coords] = coordinates.evaluateReal(cache, 3)

            # Update our search location
            new_coords = old_coords + np.dot(mesh_coords - old_coords, plane_norm) * plane_norm
            cache.setNode(old_node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, new_coords.tolist())

            test_coords = old_coords
            old_coords = new_coords

            # Break in case we can not converge
            it += 1
            if it > max_iterations:
                print(f'Could not converge on node {nodeKey}')
                break

        self.node_coordinate_list.append(new_coords)
        print(f'Node {nodeKey} was solved in {it-1} iterations')

    def updateGrid(self, new_point_id, new_point):
        index = self.node_corner_list.index(new_point_id)
        self.node_corner_points[index] = new_point

        eeg_coord = self.generateGridPoints4(8)

        self._scene.beginChange()
        self.node_coordinate_list = []
        for i in range(len(eeg_coord)-1):
            self.nudgeNode(i + self.numberInModel + 1, eeg_coord[i])

        self._scene.endChange()

    def createEEGPointsWithNormal(self, region, eeg_group, eeg_coord, i, cache):
        # createEEGPoints creates subgroups of points that use the 'colour' field to change colour

        # Re-aquire openzinc variables
        fm = region.getFieldmodule()
        coordinates = fm.findFieldByName('coordinates')
        coordinates = coordinates.castFiniteElement()
        colour = fm.findFieldByName('colour')
        colour = colour.castFiniteElement()

        eegNode = nodes.createNode(self.numberInModel + i + 1, nodetemplate)
        cache.setNode(eegNode)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, eeg_coord[i])



        #user our solver to find the nodes location
        self.moveNode(self.numberInModel + i + 1, cache)

        # Create the final node with our search coordinates
        eeg_group.addNode(eegNode)

    def createGraphics(self, new_points=False, point1=[0, 1], point2=[0, 0], point3=[1, 1], point4=[.8, 0], number_on_side=8):
        # createGraphics creates our EEG points and assigns them a spectrum for future animation.

        fm = self._region.getFieldmodule()
        scene = self._region.getScene()
        scene.beginChange()
        coordinates = fm.findFieldByName('coordinates')
        coordinates = coordinates.castFiniteElement()
        cache = fm.createFieldcache()

        #save points
        if new_points:
            self.node_corner_points[0] = point1
            self.node_corner_points[1] = point2
            self.node_corner_points[2] = point3
            self.node_corner_points[3] = point4

        # Add EEG nodes
        eeg_coord = self.generateGridPoints4(number_on_side)
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
        self.node_coordinate_list = []
        finite_element_field = []
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self.numberInModel = nodes.getSize()

        # Create all EEG subgroups
        colour = fm.createFieldFiniteElement(1)
        colour.setName('colour')
        colour.setManaged(True)
        # Create new graphics for our subgroup
        nodeColours = self._scene.createGraphicsPoints()
        nodeColours.setFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodeColours.setCoordinateField(coordinates)

        # create new subgroup containing our node
        fng = fm.createFieldNodeGroup(fm.findNodesetByName('nodes'))
        ndsg = fng.getNodesetGroup()
        ndsg.removeAllNodes()

        # Create templates
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.defineField(colour)
        nodetemplate.setValueNumberOfVersions(colour, -1, Node.VALUE_LABEL_VALUE, 1)

        # Assign values for the new EEG subset


        for i in range(len(eeg_coord)):
            eegNode = nodes.createNode(self.numberInModel + i + 1, nodetemplate)
            cache.setNode(eegNode)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, eeg_coord[i])

            # user our solver to find the nodes location
            self.moveNode(self.numberInModel + i + 1, cache)

            # Create the final node with our search coordinates
            ndsg.addNode(eegNode)

        self.ndsg = ndsg

        nodeColours.setSubgroupField(fng)

        # Set attributes for our new node
        nodeColours.setSpectrum(spec)
        nodeColours.setDataField(colour)
        pointattr = nodeColours.getGraphicspointattributes()
        # pointattr.setGlyphShapeType(Glyph.SHAPE_TYPE_SPHERE)
        # pointattr.setBaseSize([.05, .05, .05])

        # Add a colour bar for the spectrum
        check = nodes.findNodeByIdentifier(1000)
        if not check.isValid():
            screen_coords = fm.createFieldFiniteElement(2)
            spectrum_template = nodes.createNodetemplate()
            spectrum_template.defineField(screen_coords)
            spectrum_node = nodes.createNode(1000, spectrum_template)
            cache.setNode(spectrum_node)
            screen_coords.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [-.95, -.78])
            fng = fm.createFieldNodeGroup(nodes)
            spectrum_group = fng.getNodesetGroup()
            spectrum_group.addNode(spectrum_node)


            spectrum_graphics = scene.createGraphicsPoints()
            spectrum_graphics.setScenecoordinatesystem(Scenecoordinatesystem.SCENECOORDINATESYSTEM_NORMALISED_WINDOW_FIT_BOTTOM)
            spectrum_graphics.setFieldDomainType(Field.DOMAIN_TYPE_NODES)
            spectrum_graphics.setCoordinateField(screen_coords)
            spectrum_graphics.setSubgroupField(fng)
            spectrum_graphics.setSpectrum(spec)
            spectrum_point_attr = spectrum_graphics.getGraphicspointattributes()

            gm = self._scene.getGlyphmodule()
            colour_bar = gm.createGlyphColourBar(spec)
            colour_bar.setLabelDivisions(6)

            spectrum_point_attr.setGlyph(colour_bar)
            spectrum_point_attr.setBaseSize([.3, .4,])

        scene.endChange()

        # Create node corner list (used to check which nodes are on corner later
        base_node = self.numberInModel
        self.node_corner_list[0] = base_node + 1
        self.node_corner_list[2] = base_node + self.number_of_points_on_grid_side
        self.node_corner_list[1] = base_node + 1 + self.number_of_points_on_grid_side * (self.number_of_points_on_grid_side - 1)
        self.node_corner_list[3] = base_node + self.number_of_points_on_grid_side ** 2
        #del self.pointattrList[-1]

        #self.generateMesh()

    def deleteAll(self):
        fm = self._region.getFieldmodule()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        for i in range(self.eegSize):
            node = nodes.findNodeByIdentifier(i + self.numberInModel + 1)
            nodes.destroyNode(node)

    def generateMesh(self):


        plateMesh = Blackfynn_2d_plate(self._region, self.node_coordinate_list)
        plateMesh.drawMesh(self._region, self.node_coordinate_list)