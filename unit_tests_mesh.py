from __future__ import print_function

import salome
import geompy
import GEOM
import smesh
import SMESH

from numpy import array, ndarray, arange, cross
from numpy.linalg import norm
from numpy import float64 as data_type

from MyMesh.Types import *
from MyMesh.Tools import *

class UnitTester(object):

    def testTools(self):
        self.testFindMesh()

    def testFindMesh(self):

        mesh1 = find_mesh('Mesh_1')

    def testTypes(self):
        self.testTria3()

    def testTria3(self):
        mesh = find_mesh('Mesh_1')
        filter_tri = smesh.GetFilter(smesh.FACE, smesh.FT_ElemGeomType, smesh.Geom_TRIANGLE)
        ids_tri = mesh.GetIdsFromFilter(filter_tri)
        tria3 = Tria3(mesh,ids_tri[0])

        print('Tria3 Tests: ',
              tria3.computeNormalOp(),
              tria3.computeArea(),
              tria3.getArea(),
              tria3.computeNormal(),
              tria3.getNormal(),
              )

    def __init__(self):

        self.testTypes()
        self.testTools()


UnitTester()
