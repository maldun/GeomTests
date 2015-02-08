#!/usr/bin/python
# -*- coding: utf-8 -*-

# GeomTests Module - Unit Tests for MyGeom and MyMesh
# unit_tests_mesh.py: Unit tests for MyMesh 
#
# Copyright (C) 2015  Stefan Reiterer - stefan.harald.reiterer@gmail.com
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

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

from numpy import array, ndarray, arange, cross, finfo, float32, zeros
from numpy.linalg import norm
from numpy import float64 as data_type
eps = 10*finfo(float32).eps

from MyMesh.Types import *
from MyMesh.Tools import *

from MyGeom.Tools import find_object, add_to_study

from numpy import pi, sin, cos, array, zeros, cross

import os, sys, inspect

def GetAndChangeCurrentDir():
    # abspath = os.path.abspath(__file__)
    # dname = os.path.dirname(abspath)
    # os.chdir(dname)
    dname = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    return dname

script_dir = GetAndChangeCurrentDir()
print(script_dir)

class UnitTester(object):

### Tools

    def testTools(self):
        self.testFindMesh()
        self.testApplyLinearElements()
        self.testComputeVoroniAreaOfTriangle()
        self.testComputegravityCenter()

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


    def testComputeVoroniAreaOfTriangle(self):

        #case 1: obtuse in x_i
        
        w1 = w2 = pi/6.0
        l1 = array([-1.0,0.0,0.0])
        l2 = array([-cos(4*pi/6),-sin(4*pi/6),0.0])
        case1 = compute_voroni_area_of_triangle(w1,w2,l1,l2)
        area = norm(cross(l1,l2))/2.0

        #case 2: obtuse in different point
        case2 = compute_voroni_area_of_triangle(4*w1,w2,-l2,-l2+l1)

        #case 3: not obtuse
        l3 = array([0.0,-1.0,0.0])
        w3 = w4 = pi/4
        case3 = compute_voroni_area_of_triangle(w3,w4,l1,l3)
        print("test compute Voroni region on Triangle: ",
              abs(case1-area/2.0) < 1e-6,
              abs(case2-area/4.0) < 1e-6,
              abs(case3-1.0/4.0) < 1e-6,
              )
    def testComputegravityCenter(self):

        S = array([10.0,8.793144,132.093783])
        mesh_file = script_dir + '/test_gravity.med'
        test_grav = smesh.CreateMeshesFromMED(mesh_file)[0][0]
        assert norm(S - compute_gravity_center(test_grav))
        mesh_file2 = script_dir + '/test_groups.med'
        test_groups = smesh.CreateMeshesFromMED(mesh_file2)[0][0]
        test_group_group = test_groups.GetGroups()[0]
        assert norm(zeros(3) - compute_gravity_center(test_groups,test_group_group)) < eps
        print("TestComputeGravityCenter: ", True)
### Types

    def testTypes(self):
        self.testTria3()
        self.testQuad4()
        self.testNormalVectorField()
        self.testNormalVectorFieldProbs()
        self.testMovingMethods()
        self.testMeanCurvatureStuff()
        self.testMeanCurvatureNormal()
        self.testMeanCurvatureNormalFiner()
        self.testPlaneProjections()
        self.test_FaceProjectVectorField()

    def testTria3(self):
        mesh = find_object('Mesh_1')
        mesh = smesh.Mesh(mesh)
        filter_tri = smesh.GetFilter(SMESH.FACE, SMESH.FT_ElemGeomType, SMESH.Geom_TRIANGLE)
        ids_tri = mesh.GetIdsFromFilter(filter_tri)
        tria3 = Tria3(mesh,ids_tri[0])
        tria_node1 =  tria3.getNodes()[0]
        tria3.computeArea()
        tria3.computeNormal(tria_node1)

        # Test the formula for mean normal curvature
        # center node id = 2; element to test 19
        mesh4 = find_mesh('Mesh_4')
        tria3C = Tria3(mesh4,19)

        # compute center of gravity
        node_ids = tria3C.getNodes()
        coord_matrix = array([mesh4.GetNodeXYZ(node) for node in node_ids])
        center = apply_along_axis(sum,0,coord_matrix)/3.0

        # test on trialtriangle
        A = array([1.0,2.0,3.0])
        B = array([5.0,3.0,7.0])
        B = array([1.0,-1.0,4.0])
        S = array([7.0/3.0,4.0/3.0,14.0/3.0])
        mesh_file = script_dir + '/test_tria.med'
        test_tria = smesh.CreateMeshesFromMED(mesh_file)[0][0]
        TestTria = Tria3(test_tria,4)
        
        print('Tria3 Tests: ',
              tria3.getNodes(),
              tria3._computeNormalOp(),
              tria3.getArea(),
              tria3.getNormal(tria_node1),
              tria3.getNormals(),
              norm(tria3C.computeCurvatureVector(2)-array([-0.5,0.5,0.0])) < 0.01,
              abs(tria3C.computeCurvatureVector(2,voroni=True)[1]-0.5/8) < 1e-3,
              norm(tria3C.computeGravityCenter() - center) < eps,
              norm(TestTria.computeGravityCenter() - S) < eps,
              )

    def testQuad4(self):
        mesh = find_mesh('Mesh_2')
        filter_quad = smesh.GetFilter(SMESH.FACE, SMESH.FT_ElemGeomType, SMESH.Geom_QUADRANGLE)
        ids_quad = mesh.GetIdsFromFilter(filter_quad)
        quad4 = Quad4(mesh,ids_quad[0])
        quad4_node1 =  quad4.getNodes()[0]

        quad4.computeArea()
        quad4.computeNormal(quad4_node1)

        #Test mean curvature vector
        # Center id = 9; element id = 12
        mesh_cq = find_mesh('Mesh_curv_quad')
        quad4C = Quad4(mesh_cq,12)

        # test center of gravity
        S = array([2.647059,-1.372549,4.745098])
        mesh_file = script_dir + '/test_quad.med'
        test_quad = smesh.CreateMeshesFromMED(mesh_file)[0][0]
        TestQuad = Quad4(test_quad,5)
        

        
        print('Quad4 Tests: ',
              quad4.getNodes(),
              quad4._computeNormalOp(),
              quad4.getArea(),
              #quad4.getNormal(quad4_node1),
              quad4.getNormals(),
              norm(quad4C.computeCurvatureVector(9)-array([0.0,1.0,0.0])) < 0.01,
              (quad4C.computeCurvatureVector(9,voroni=True)[1]-1.0/8.0) < 1e-3,
              norm(TestQuad.computeGravityCenter() - S) < eps,
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
        boundary = mesh.GetElementsByType(SMESH.EDGE)
        boundary = mesh7.MakeGroupByIds("boundary",SMESH.EDGE,boundary)
        norm_field3 = 5.0*NormalVectorField(mesh7)
        stuff = norm_field3.computeSurfaceExtrusion()#edge_groups=[boundary])
        new_surf = mesh7.MakeGroupByIds('new_surf',SMESH.FACE, stuff[0])
        #new_edge_group = mesh7.MakeGroupByIds('new_boundary',SMESH.EDGE, stuff[2][0])
        stuff1 = norm_field3.extrudeSurface(group=new_surf)#,edge_groups=[new_edge_group])

        mesh3 = find_mesh('Mesh_3')
        mesh8 = smesh.CopyMesh( mesh3, "Mesh_8")
        fix1 = find_object('fix1')
        fix2 = find_object('fix2')
        norm_field4 = 2.0*NormalVectorField(mesh8)
        stuff1 = norm_field4.extrudeSurface(edge_groups=[fix1,fix2])

        mesh9 = smesh.CopyMesh(mesh,"Mesh_9")
        norm_field5 = 1.0*NormalVectorField(mesh9)

        boundary2 = mesh9.GetElementsByType(SMESH.EDGE)
        boundary2 = mesh9.MakeGroupByIds("boundary",SMESH.EDGE,boundary2)
        
        norm_field5.extrudeSurfaceTimes([0.5,1.0,2.0],edge_groups=[boundary2],face_groups=[])

        mesh10 = find_mesh('test_faces')
        norm_field6 = 1.0*NormalVectorField(mesh10)
        Floor = find_object('floorF')
        norm_field6.extrudeSurfaceTimes(2,face_groups = [Floor])

        mesh11 = find_mesh('test_groups')
        norm_field7 = 1.0*NormalVectorField(mesh11)
        circle = find_object('circleF')
        norm_field7.applyVectorFieldOnSurface(group = circle)

        print('Test normal vector field: ',
              truth.all(),
              truth2.all(),
              truth3.all(),
              new_ids,
              stuff[0],
              stuff[1]
              )


        
    def testMovingMethods(self):
        # test the methods for moving nodes
        mesh_file = script_dir + '/Mesh_4.med'
        mesh4 = smesh.CreateMeshesFromMED(mesh_file)[0][0]
        norm_field_mov = NormalVectorField(mesh4)
        node_id = 1
        place = array(mesh4.GetNodeXYZ(node_id))
        norm_vec = norm_field_mov.computeVectorOnNode(node_id)
        pot_new_place = norm_field_mov.computeNewPosition(node_id)
        norm_field_mov.moveNode(node_id)
        new_place = array(mesh4.GetNodeXYZ(node_id))
        assert norm((place+norm_vec) - new_place) < eps
        assert norm(pot_new_place - new_place) < eps

        # test moving of groups and meshes
        mesh_file2 = script_dir + '/test_groups.med'
        test_groups = smesh.CreateMeshesFromMED(mesh_file2)[0][0]
        norm_vec2 = 5*NormalVectorField(test_groups)
        group = norm_vec2.mesh.GetGroups()[0]
        norm_vec2.moveSurface(group)
        group_ids = group.GetNodeIDs()
        nodes = norm_vec2.mesh.GetNodesId()
        vectors = [[ids,norm_vec2.mesh.GetNodeXYZ(ids)] for ids in nodes if not (ids in group_ids)] 
        vectors += [[ids,norm_vec2.computeNewPosition(ids)] for ids in group_ids]
        vectors = dict(vectors) 
        norm_vec2.moveSurface(group)
        check = [norm(array(vectors[ids]) - array(norm_vec2.mesh.GetNodeXYZ(ids))) < eps
                 for ids in nodes]
        assert all(check)
        print("Test moving of nodes by vectorfield: ", True)
        
        
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
        filter_tri = smesh.GetFilter(SMESH.FACE, SMESH.FT_ElemGeomType, SMESH.Geom_TRIANGLE)
        ids_tri = mesh4.GetIdsFromFilter(filter_tri)
        tria_elems = [Tria3(mesh4,id_tri) for id_tri in ids_tri]
        print('test mean curvature Formula for tria3: ',
              norm(norm_field4.meanCurvatureNormalFormula(tria_elems,2)) < 1e-10,
            ) 
        
        # quadrangles
        mesh_cq = find_mesh('Mesh_curv_quad')
        norm_field_cq = NormalVectorField(mesh_cq)
        
        filter_q = smesh.GetFilter(SMESH.FACE, SMESH.FT_ElemGeomType, SMESH.Geom_QUADRANGLE)
        ids_q = mesh_cq.GetIdsFromFilter(filter_q)
        quad_elems = [Quad4(mesh_cq,id_q) for id_q in ids_q]
        print('test mean curvature Formula for quad4: ',
              norm(norm_field_cq.meanCurvatureNormalFormula(quad_elems,9)) < 1e-10,
            ) 
        
    def testMeanCurvatureNormal(self):

        mesh = find_mesh('Mesh_Sphere')
        mean = MeanCurvatureNormal(mesh)
        normal = NormalVectorField(mesh)

        nodes = mesh.GetNodesId()

        mean_vecs = [mean.computeVectorOnNode(node) for node in nodes]
        mean_vecs_normed = [vec/norm(vec) for vec in mean_vecs]
        # On sphere: n(x) = x
        normal_vecs = [mesh.GetNodeXYZ(node) for node in nodes]
        
        tester = [norm(mean_vecs[i]*0.5-normal_vecs[i]) for i in range(len(nodes))]
        tester2 = [norm(mesh.GetNodeXYZ(node)-normal.computeVectorOnNode(node)) for node in nodes]
        tester3 = [norm(mean_vecs_normed[i]-normal_vecs[i]) for i in range(len(nodes))]
        
        print('Test mean curvature Normal: ',
              'Mean curvature error: ', max(tester),
              'Standard Normal error: ', max(tester2),
              'Normed mean curvature error: ', max(tester3),
              )

    def testMeanCurvatureNormalFiner(self):

        mesh = find_mesh('Mesh_Sphere2')
        mean = MeanCurvatureNormal(mesh)
        normal = NormalVectorField(mesh)

        nodes = mesh.GetNodesId()

        mean_vecs = [mean.computeVectorOnNode(node) for node in nodes]
        mean_vecs_normed = [vec/norm(vec) for vec in mean_vecs]
        # On sphere: n(x) = x
        normal_vecs = [mesh.GetNodeXYZ(node) for node in nodes]
        
        tester = [norm(mean_vecs[i]*0.5 - normal_vecs[i]) for i in range(len(nodes))]
        tester2 = [norm(mesh.GetNodeXYZ(node)-normal.computeVectorOnNode(node)) for node in nodes]
        tester3 = [norm(mean_vecs_normed[i]-normal_vecs[i]) for i in range(len(nodes))]
        
        print('Test mean curvature Normal on finer grid: ',
              'Mean curvature error: ', max(tester),
              'Standard Normal error: ', max(tester2),
              'Normed mean curvature error: ', max(tester3),              
              )

    def testPlaneProjections(self):
        from numpy import sqrt
        O = array([0.0,0.0,-0.5])
        u = array([1.0,0.0,0.0])
        v = array([0.0,sqrt(0.5),sqrt(0.5)])
        w = array([0.0,-sqrt(0.5),sqrt(0.5)])
        Q = array([u,v,w]).transpose()

        mesh_file = script_dir + '/test_disk2.med'
        mesh = smesh.CreateMeshesFromMED(mesh_file)[0][0]

        # init class
        flachi = PlaneProjectionVectorField(mesh,O,Q,0.5)
        S = compute_gravity_center(mesh)

        # test_correctness of computation
        node2 = array(flachi.mesh.GetNodeXYZ(2))
        traf_vec = flachi.trafo(node2)
        lokal_z = traf_vec.reshape(3)[-1]
        traf_vec[-1] = 0.0
        result = flachi.inv_trafo(traf_vec).reshape(3)
        assert norm(result + lokal_z*w - node2) < eps
        assert norm(result -flachi.computeSingleProjection(2)) < eps
        assert norm(flachi.computeSingleProjection(2) - (flachi.computeVectorOnNode(2)+node2)) < eps
        vectors = flachi.getNodeVectors()
        traf_vecs = flachi.trafo(vectors)
        assert flachi._makeChecks(traf_vecs) is None
        projections = flachi.computeProjections().transpose()
        ids = flachi._internal_ids
        node_ids = flachi.mesh.GetNodesId()
        assert all([norm(projections[ids[node],:] - flachi.computeSingleProjection(node))<eps for 
             node in node_ids])

        flachi._vectors = projections.transpose()
        assert all([norm(flachi.computeVectorOnNode(node) + array(flachi.mesh.GetNodeXYZ(node)) 
                         - flachi.computeSingleProjection(node))<eps for node in node_ids])
        #flachi.moveSurface()
        group = flachi.mesh.GetGroups()[0]
        print("Test PlaneProjections: ", True)
        
    def test_FaceProjectVectorField(self):
        from numpy import sqrt
        O = array([0.0,0.0,-0.5])
        u = array([1.0,0.0,0.0])
        v = array([0.0,sqrt(0.5),sqrt(0.5)])
        w = array([0.0,-sqrt(0.5),sqrt(0.5)])
        Q = array([u,v,w]).transpose()

        mesh_file = script_dir + '/test_disk2.med'
        mesh = smesh.CreateMeshesFromMED(mesh_file)[0][0]


        # init class
        flachi = PlaneProjectionVectorField(mesh,O,Q,0.5)

        # import plane
        plane = find_object("Scale_1")
        mesh2 = smesh.CreateMeshesFromMED(mesh_file)[0][0]
        flachi2 = FaceProjectVectorField(mesh2,plane,0.5)

        nodes = mesh.GetNodesId()
        assert all([norm(flachi.computeVectorOnNode(node) - flachi2.computeVectorOnNode(node)) < eps 
             for node in nodes])
        print("Test FaceProjectVectorField: ", True)
        
    def __init__(self):

        self.testTypes()
        self.testTools()

testi = UnitTester()

