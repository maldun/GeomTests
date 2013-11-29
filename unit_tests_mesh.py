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
        self.testNormalVectorFieldProbs()
        self.testMeanCurvatureStuff()

    def testTria3(self):
        mesh = find_object('Mesh_1')
        mesh = smesh.Mesh(mesh)
        filter_tri = smesh.GetFilter(smesh.FACE, smesh.FT_ElemGeomType, smesh.Geom_TRIANGLE)
        ids_tri = mesh.GetIdsFromFilter(filter_tri)
        tria3 = Tria3(mesh,ids_tri[0])
        tria_node1 =  tria3.getNodes()[0]
        tria3.computeArea()
        tria3.computeNormal(tria_node1)

        # Test the formula for mean normal curvature
        # center node id = 2; element to test 19
        mesh4 = find_mesh('Mesh_4')
        tria3C = Tria3(mesh4,19)
        
        print('Tria3 Tests: ',
              tria3.getNodes(),
              tria3._computeNormalOp(),
              tria3.getArea(),
              tria3.getNormal(tria_node1),
              tria3.getNormals(),
              norm(tria3C.getCurvatureVector(2)-array([-0.5,0.5,0.0])) < 0.01,
              )

    def testQuad4(self):
        mesh = find_mesh('Mesh_2')
        filter_quad = smesh.GetFilter(smesh.FACE, smesh.FT_ElemGeomType, smesh.Geom_QUADRANGLE)
        ids_quad = mesh.GetIdsFromFilter(filter_quad)
        quad4 = Quad4(mesh,ids_quad[0])
        quad4_node1 =  quad4.getNodes()[0]

        quad4.computeArea()
        quad4.computeNormal(quad4_node1)

        #Test mean curvature vector
        # Center id = 9; element id = 12
        mesh_cq = find_mesh('Mesh_curv_quad')
        quad4C = Quad4(mesh_cq,12)

        print('Quad4 Tests: ',
              quad4.getNodes(),
              quad4._computeNormalOp(),
              quad4.getArea(),
              #quad4.getNormal(quad4_node1),
              quad4.getNormals(),
              norm(quad4C.getCurvatureVector(9)-array([0.0,1.0,0.0])) < 0.01,
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

        mesh9 = smesh.CopyMesh(mesh,"Mesh_9")
        norm_field5 = 1.0*NormalVectorField(mesh9)

        boundary2 = mesh9.GetElementsByType(smesh.EDGE)
        boundary2 = mesh9.MakeGroupByIds("boundary",smesh.EDGE,boundary2)
        
        norm_field5.extrudeSurfaceTimes([0.5,1.0,2.0],edge_groups=[boundary2])
        print('Test normal vector field: ',
              truth.all(),
              truth2.all(),
              truth3.all(),
              new_ids,
              stuff[0],
              stuff[1]
              )

    def testNormalVectorFieldProbs(self):
        # Test problematic case
        meshTA = find_mesh('MA_T')
        edge_groupTA = find_object('G_2233T')

        
        try:
            norm_fieldTA = 0.2*NormalVectorField(meshTA)
        except NotImplementedError:
            print('Test normalvector problematic cases: correct error handling')

    def testMeanCurvatureStuff(self):
        mesh4 = find_mesh('Mesh_4')
        norm_field4 = NormalVectorField(mesh4)
        # get triangles
        filter_tri = smesh.GetFilter(smesh.FACE, smesh.FT_ElemGeomType, smesh.Geom_TRIANGLE)
        ids_tri = mesh4.GetIdsFromFilter(filter_tri)
        tria_elems = [Tria3(mesh4,id_tri) for id_tri in ids_tri]
        print('test mean curvature Formula for tria3: ',
              norm(norm_field4.meanNormalCurvatureFormula(tria_elems,2)) < 1e-10,
            ) 
        
        # quadrangles
        mesh_cq = find_mesh('Mesh_curv_quad')
        norm_field_cq = NormalVectorField(mesh_cq)
        
        filter_q = smesh.GetFilter(smesh.FACE, smesh.FT_ElemGeomType, smesh.Geom_QUADRANGLE)
        ids_q = mesh_cq.GetIdsFromFilter(filter_q)
        quad_elems = [Quad4(mesh_cq,id_q) for id_q in ids_q]
        print('test mean curvature Formula for quad4: ',
              norm(norm_field_cq.meanNormalCurvatureFormula(quad_elems,9)) < 1e-10,
            ) 
        
    def __init__(self):

        self.testTypes()
        self.testTools()


UnitTester()
salome.sg.updateObjBrowser(0)
