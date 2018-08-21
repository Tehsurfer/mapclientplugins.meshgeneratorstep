#!/usr/bin/python
"""
PyZinc examples
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""

import sys
try:
    from PySide import QtGui
except ImportError:
    from PyQt4 import QtGui

from opencmiss.zinc.context import Context
from opencmiss.zinc.graphics import Graphics
from opencmiss.zinc.material import Material
import json



class Export(QtGui.QWidget):
    '''
    This example demonstrates how to read and export a simple mesh
    '''

    def __init__(self, scene, sceneviewer):
        self.scene = scene
        self.sceneviewer = sceneviewer
        self._prefix = 'MeshTest'
        self.Path = 'C:/Users/jkho021/Projects/zinc test/ZincJS-Tutorials/tutorials/examples/basic_mesh/html/'

    def exportViewJson(self):
        '''Export sceneviewer parameters to JSON format'''
        sceneviewer = self.sceneviewer
        farPlane = sceneviewer.getFarClippingPlane()
        nearPlane = sceneviewer.getNearClippingPlane()
        returnValue = sceneviewer.getLookatParameters()
        eyePosition = returnValue[1]
        targetPosition = returnValue[2]
        upVector = returnValue[3]
        viewData = {}
        viewData['farPlane'] = farPlane
        viewData['nearPlane'] = nearPlane
        viewData['eyePosition'] = eyePosition
        viewData['targetPosition'] = targetPosition
        viewData['upVector'] = upVector
        f = open(self.Path + self._prefix + '_view' + '.json', 'w+')
        json.dump(viewData, f)
        f.close()

    def exportWebGLJson(self):
        '''
        Export graphics into JSON format, one json export represents one
        surface graphics.
        '''
        scene = self.scene
        sceneSR = scene.createStreaminformationScene()
        sceneSR.setIOFormat(sceneSR.IO_FORMAT_THREEJS)
        ''' Get the total number of graphics in a scene/region that can be exported'''
        number = sceneSR.getNumberOfResourcesRequired()
        resources = []
        '''Write out each graphics into a json file which can be rendered with our WebGL script'''
        for i in range(number):
            resources.append(sceneSR.createStreamresourceMemory())
        scene.write(sceneSR)
        '''Write out each resource into their own file'''
        for i in range(number):
            f = None
            if i == 0:
                f = open(self.Path + self._prefix + '_' + 'metadata.json', 'w+')
            else:
                f = open(self.Path + self._prefix + '_' + str(i) + '.json', 'w+')
            buffer = resources[i].getBuffer()[1]
            if i == 0:
                for j in range(number - 1):
                    replaceName = '' + self._prefix + '_' + str(j + 1) + '.json'
                    old_name = 'memory_resource' + '_' + str(j + 2)
                    buffer = buffer.replace(old_name, replaceName)
            f.write(buffer)
            f.close()
        ''' the following function exports the camera settings'''
        self.exportViewJson()

