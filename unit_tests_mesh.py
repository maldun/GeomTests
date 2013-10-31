from __future__ import print_function

# import salome
# import geompy
# import GEOM
# import smesh
# import SMESH

import salome
salome.salome_init()
import GEOM
from salome.geom import geomBuilder
geompy = geomBuilder.New(salome.myStudy)

import SMESH, SALOMEDS
from salome.smesh import smeshBuilder
smesh =  smeshBuilder.New(salome.myStudy)

from numpy import array, ndarray, arange, cross
from numpy.linalg import norm
from numpy import float64 as data_type

from MyMesh.Types import *
from MyMesh.Tools import *

from MyGeom.Tools import find_object



class UnitTester(object):

    def testTools(self):
        self.testFindMesh()
        self.testApplyLinearElements()

    def testFindMesh(self):
        pass
        # mesh1 = find_mesh('Mesh_1')

    def testApplyLinearElements(self):
        mesh = find_mesh('Mesh_2')
        node_id = mesh.GetNodesId()[73]
        elems = mesh.GetNodeInverseElements(node_id)

        print('Test apply_linear_elements: ',
              apply_linear_elements(mesh,elems),
              )

### Types

    def testTypes(self):
        self.testTria3()
        self.testQuad4()
        self.testNormalVectorField()

    def testTria3(self):
        mesh = find_object('Mesh_1')
        mesh = smesh.Mesh(mesh)
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

    def testNormalVectorField(self):
        
        mesh = find_mesh('Mesh_1')
        norm_field = NormalVectorField(mesh)

        nodes = mesh.GetNodesId()
        normals = [norm_field.getVectorOnNode(node) for node in nodes]
        truth = array([normal == array((0.0,0.0,1.0)) for normal in normals])

        norm_field.scalarMultiplication(0.5)
        normals2 = [norm_field.getVectorOnNode(node) for node in nodes]
        truth2 = array([normal == array((0.0,0.0,0.5)) for normal in normals2])

        norm_field.setScalar(1.0)
        norm_field2 = 0.5*norm_field
        normals3 = [norm_field2.getVectorOnNode(node) for node in nodes]
        truth3 = array([normal == array((0.0,0.0,0.5)) for normal in normals3])

        mesh6 = smesh.Mesh('Mesh_6')
        #new_ids = [norm_field2.applyVectorOnNode(node,mesh6) for node in nodes]
        faces = mesh.GetElementsByType(FACE)
        #new_ids = [norm_field2.applyVectorFieldOnFace(face,mesh6) for face in faces]
        new_ids = norm_field2.applyVectorFieldOnSurface(mesh6)
        
        
        mesh7 = smesh.CopyMesh( mesh, "Mesh_7")
        boundary = mesh.GetElementsByType(smesh.EDGE)
        boundary = mesh7.MakeGroupByIds("boundary",smesh.EDGE,boundary)
        norm_field3 = 5.0*NormalVectorField(mesh7)
        stuff = norm_field3.computeSurfaceExtrusion()#edge_groups=[boundary])
        new_surf = mesh7.MakeGroupByIds('new_surf',smesh.FACE, stuff[0])
        #new_edge_group = mesh7.MakeGroupByIds('new_boundary',smesh.EDGE, stuff[2][0])
        stuff1 = norm_field3.extrudeSurface(group=new_surf)#,edge_groups=[new_edge_group])

        mesh3 = find_mesh('Mesh_3')
        mesh8 = smesh.CopyMesh( mesh3, "Mesh_8")
        fix1 = find_object('fix1')
        fix2 = find_object('fix2')
        norm_field4 = 2.0*NormalVectorField(mesh8)
        stuff1 = norm_field4.extrudeSurface(edge_groups=[fix1,fix2])

        print('Test normal vector field: ',
              truth.all(),
              truth2.all(),
              truth3.all(),
              new_ids,
              stuff[0],
              stuff[1]
              )

    def __init__(self):

        self.testTypes()
        self.testTools()


UnitTester()
salome.sg.updateObjBrowser(0)
