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
        self.testQuad4()

    def testTria3(self):
        mesh = find_mesh('Mesh_1')
        filter_tri = smesh.GetFilter(smesh.FACE, smesh.FT_ElemGeomType, smesh.Geom_TRIANGLE)
        ids_tri = mesh.GetIdsFromFilter(filter_tri)
        tria3 = Tria3(mesh,ids_tri[0])

        tria3.computeArea()
        tria3.computeNormal()

        print('Tria3 Tests: ',
              tria3.getNodes(),
              tria3.computeNormalOp(),
              tria3.getArea(),
              tria3.getNormal(),
              )

    def testQuad4(self):
        mesh = find_mesh('Mesh_2')
        filter_quad = smesh.GetFilter(smesh.FACE, smesh.FT_ElemGeomType, smesh.Geom_QUADRANGLE)
        ids_quad = mesh.GetIdsFromFilter(filter_quad)
        quad4 = Quad4(mesh,ids_quad[0])

        quad4.computeArea()
        quad4.computeNormal()

        print('Quad4 Tests: ',
              quad4.getNodes(),
              quad4.computeNormalOp(),
              quad4.getArea(),
              quad4.getNormal(),
              )


    def __init__(self):

        self.testTypes()
        self.testTools()


UnitTester()
