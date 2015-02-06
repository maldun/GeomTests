#!/usr/bin/python
# -*- coding: utf-8 -*-
# MyGeom Module - API for easier Salome geompy usage
# Tools.py: Tool functions for MyGeom module
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

"""
Unit tests for the MyMath module
"""
from __future__ import print_function
import warnings
import numpy as np
from numpy import ones, eye, dot
from numpy import cos, sin, abs
from numpy.linalg import norm
eps = 10*np.finfo(np.float32).eps
warnings.simplefilter('always', UserWarning)

class TypeTests(object):
    """
    Test class for the Type module.
    """

    def checkTests(self,name,passed):
        """
        Checks if all tests of a test with name failed
        """

        if not all(passed):
            raise Exception("Error: " + name + " tests didn't passed!")
        print("All " + name + " tests passed!")
    
    def testMathOperator(self):
        """
        Base Class
        """
        from Types import MathOperator
        passed = [False]
        x = 5
        # init
        
        try:
            with warnings.catch_warnings(record=True) as warn:
                mathi = MathOperator()
                assert self.fallback_warning in str(warn[-1].message)
                try:
                    mathi(x)
                except NotImplementedError:
                    passed[0] = True
        except AssertionError:
            pass


        try:
            mathi._pythonOP(x)
        except NotImplementedError:
            passed += [True]

        try:
            mathi._optimizedOP(x)
        except NotImplementedError:
            passed += [True]

        self.checkTests('MathOperator',passed)

    def testGeometricTransformation(self):
        """
        testing the Geometric Transformation class
        """
        from Types import GeometricTransformation
        passed = [False]
        # Tests for Init
        # test if b is a vector
        try:
            GeometricTransformation(ones((3,3)),ones((3,3)))
        except ValueError:
            passed[0] = True

        # test if Q is orthogonal works
        try:
            GeometricTransformation(ones((3,3)),ones(3))
            passed += [False]
        except ValueError:
            passed += [True]

        # test if fallback works
        with warnings.catch_warnings(record=True) as warn:
            mathi = GeometricTransformation(eye(3),ones(3),method=666)
            assert self.fallback_warning in str(warn[-1].message)
            passed += [True]
        
        # Mathematical correctness:
        # special case
        tester = GeometricTransformation(eye(3),ones(3))
        assert norm(tester(ones(3))-ones(3)*2) < eps
        
        # series of random tests
        for i in range(self.nr_tests):
            QR = np.random.randn(3,3)
            QR = np.linalg.qr(QR)[0]
            bR = np.random.randn(3,1)
            testerR = GeometricTransformation(QR,bR)
            xR = np.random.randn(3,1)
            assert norm(testerR(xR) - dot(QR,xR) - bR) < eps
            testerR_inv = testerR.inv()
            assert norm(testerR_inv(testerR(xR)) - xR) < eps
        
        passed += [True]
        self.checkTests('GeometricTransformation',passed)

    def testPolarCoordinates(self):
        from Types import PolarCoordinates
        pi = np.pi
        from numpy import sqrt, cos, sin
        
        passed = [False]
        with warnings.catch_warnings(record=True) as warn:
            polarDummy = PolarCoordinates(method=666)
            assert issubclass(warn[-1].category, UserWarning)
            assert self.fallback_warning in str(warn[-1].message)
        passed[0] = True
        polar = PolarCoordinates(method=0)
        # test special cases
        v = np.array([0.0,1.0]); w = np.array([1.0,pi/2]) 
        assert norm(polar(v) - w) < eps, "Error: Special case 1 false!"
        v = np.array([0.0,-2.0]); w = np.array([2.0,-pi/2])
        assert norm(polar(v) - w) < eps, "Error: Special case 2 false!"
        passed += [True]
        # first sector 
        v = np.array([sqrt(0.5),sqrt(0.5)]); w = np.array([1.0,pi/4])
        assert norm(polar(v) - w) < eps, "Error: Sector 1 false!"
        v = np.array([3*cos(7*pi/8),3*sin(7*pi/8)]); w = np.array([3.0,7*pi/8])
        assert norm(polar(v) - w) < eps, "Error: Sector 2 false!"
        v = np.array([5*cos(-7*pi/8),5*sin(-7*pi/8)]); w = np.array([5.0,-7*pi/8])
        assert norm(polar(v) - w) < eps, "Error: Sector 3 false!"
        v = np.array([8*cos(-pi/5),8*sin(-pi/5)]); w = np.array([8.0,-pi/5])
        assert norm(polar(v) - w) < eps, "Error: Sector 4 false!"

        nr_tests = self.nr_tests
        R = np.random.rand(nr_tests)*nr_tests
        PHI = 2*pi*(np.random.rand(nr_tests)-0.5)
        for i in range(nr_tests):
            v = np.array([R[i]*cos(PHI[i]),R[i]*sin(PHI[i])]); w = np.array([R[i],PHI[i]])
            assert norm(polar(v) - w) < eps, "Error: Random test" + str(i) + " false!"
        
        self.checkTests("PolarCoordinates",passed)

    def testCartesianCoordinates(self):
        from Types import CartesianCoordinates
        pi = np.pi
        from numpy import sqrt, cos, sin
        
        passed = [False]
        with warnings.catch_warnings(record=True) as warn:
            cartesianDummy = CartesianCoordinates(method=666)
            assert issubclass(warn[-1].category, UserWarning)
            assert self.fallback_warning in str(warn[-1].message)
        passed[0] = True
        cart = CartesianCoordinates(method=0)
        # test wron Input
        try:
            cart(np.array((5,7)))
        except ValueError:
            passed += [True]
        
        # test special cases
        v = np.array([0.0,1.0]); w = np.array([1.0,pi/2]) 
        assert norm(cart(w) - v) < eps, "Error: Special case 1 false!"
        v = np.array([0.0,-2.0]); w = np.array([2.0,-pi/2])
        assert norm(cart(w) - v) < eps, "Error: Special case 2 false!"
        passed += [True]
        # first sector 
        v = np.array([sqrt(0.5),sqrt(0.5)]); w = np.array([1.0,pi/4])
        assert norm(cart(w) - v) < eps, "Error: Sector 1 false!"
        v = np.array([3*cos(7*pi/8),3*sin(7*pi/8)]); w = np.array([3.0,7*pi/8])
        assert norm(cart(w) - v) < eps, "Error: Sector 2 false!"
        v = np.array([5*cos(-7*pi/8),5*sin(-7*pi/8)]); w = np.array([5.0,-7*pi/8])
        assert norm(cart(w) - v) < eps, "Error: Sector 3 false!"
        v = np.array([8*cos(-pi/5),8*sin(-pi/5)]); w = np.array([8.0,-pi/5])
        assert norm(cart(w) - v) < eps, "Error: Sector 4 false!"

        nr_tests = self.nr_tests
        R = np.random.rand(nr_tests)*nr_tests
        PHI = 2*pi*(np.random.rand(nr_tests)-0.5)
        for i in range(nr_tests):
            v = np.array([R[i]*cos(PHI[i]),R[i]*sin(PHI[i])]); w = np.array([R[i],PHI[i]])
            assert norm(cart(w) - v) < eps, "Error: Random test" + str(i) + " false!"
        
        self.checkTests("CartesianCoordinates",passed)

    def testSphericalCoordinates(self):
        from Types import SphericalCoordinates
        pi = np.pi
        from numpy import sqrt, cos, sin
        passed = [False]
        with warnings.catch_warnings(record=True) as warn:
            sphericalDummy = SphericalCoordinates(method=666)
            assert issubclass(warn[-1].category, UserWarning)
            assert self.fallback_warning in str(warn[-1].message)
        passed[0] = True
        
        spher = SphericalCoordinates(method=0)
        nr_tests = self.nr_tests
        R = np.random.rand(nr_tests)*nr_tests
        PHI = 2*pi*(np.random.rand(nr_tests)-0.5)
        PSI = pi*(np.random.rand(nr_tests))
        for i in range(nr_tests):
            v = np.array([R[i]*cos(PHI[i])*sin(PSI[i]),
                          R[i]*sin(PHI[i])*sin(PSI[i]),
                          R[i]*cos(PSI[i])]) 
            w = np.array([R[i],PHI[i],PSI[i]])
            assert norm(spher(v) - w)  < eps, "Error: Random test" + str(i) + " false!"
        passed += [True]
            
        self.checkTests("SphericalCoordinates",passed)


    def testCartesianCoordinates3D(self):
        from Types import CartesianCoordinates3D
        pi = np.pi
        from numpy import sqrt, cos, sin
        
        passed = [False]
        with warnings.catch_warnings(record=True) as warn:
            cartesianDummy = CartesianCoordinates3D(method=666)
            assert issubclass(warn[-1].category, UserWarning)
            assert self.fallback_warning in str(warn[-1].message)
        passed[0] = True
        cart = CartesianCoordinates3D(method=0)
        # test wron Input
        try:
            cart(np.array((5,7,0.1)))
        except ValueError:
            passed += [True]
        try:
            cart(np.array((5,0.1,5)))
        except ValueError:
            passed += [True]
        
        nr_tests = self.nr_tests
        R = np.random.rand(nr_tests)*nr_tests
        PHI = 2*pi*(np.random.rand(nr_tests)-0.5)
        PSI = pi*(np.random.rand(nr_tests))
        for i in range(nr_tests):
            v = np.array([R[i]*cos(PHI[i])*sin(PSI[i]),
                          R[i]*sin(PHI[i])*sin(PSI[i]),
                          R[i]*cos(PSI[i])]) 
            w = np.array([R[i],PHI[i],PSI[i]])
            assert norm(cart(w) - v)  < eps, "Error: Random test" + str(i) + " false!"

        self.checkTests("CartesianCoordinates3D",passed)


        
    def computeGivens(self,i,j,phi,dim):
        c = cos(phi)
        s = sin(phi)
        G = eye(dim)
        if j < i:
            G[i,i] = c
            G[j,j] = c
            G[i,j] = -s
            G[j,i] = s
        else:
            G[i,i] = c
            G[j,j] = c
            G[i,j] = s
            G[j,i] = -s

        return G
        
    def testGivensmatvec(self,i,j,phi,dim):
        from Types import GivensRotator
        G1 = GivensRotator(i,j,phi,dim=dim)
        G2 = self.computeGivens(i,j,phi,dim)
        result = eye(dim)
        result = np.apply_along_axis(G1.matvec,0,result)
        result = norm(G2 - result)
        return result < eps

    def testGivensmatmat(self,i,j,phi,dim):
        from Types import GivensRotator
        G1 = GivensRotator(i,j,phi,dim=dim)
        G2 = self.computeGivens(i,j,phi,dim)
        result = eye(dim)
        result = G1(result)
        result = norm(G2 - result)
        return result < eps

    def testGivensmattrans(self,i,j,phi,dim):
        from Types import GivensRotator
        G1 = GivensRotator(i,j,phi,dim=dim)
        G2 = G1.transpose()
        G3 = G1.inv()
        x = np.random.rand(dim)
        result1 = G1.matvec(x)
        result2 = norm(G2.matvec(result1) - x)
        result3 = norm(G3.matvec(result1) - x)

        return result2 < eps and result3 < eps

    
    
    def testGivensRotator(self):

        from Types import GivensRotator
        pi = np.pi
        from numpy import sqrt, cos, sin, abs
        passed = [False]
        # test creation
        with warnings.catch_warnings(record=True) as warn:
            rotDummy = GivensRotator(0,1,0.0,method=666)
            assert issubclass(warn[-1].category, UserWarning)
            assert self.fallback_warning in str(warn[-1].message)
            rotDummy.matvec
            rotDummy.getPhi()
        passed[0] = True

        # Check creation
        rot1 = GivensRotator(1,0,0.0)
        assert abs(rot1.c-1.0) < eps
        assert abs(rot1.s) < eps
        rot2 = GivensRotator(1,0,pi/2)
        assert abs(rot2.c-0.0) < eps
        assert abs(rot2.s-1.0) < eps
        rot3 = GivensRotator(1,0,sqrt(0.5),sqrt(0.5))
        assert abs(rot3.c-sqrt(0.5)) < eps
        assert abs(rot3.s-sqrt(0.5)) < eps
        passed += [True]
        try:
            GivensRotator(0,1,0.5,0.5)
            passed += [False]
        except ValueError:
            passed += [True]
        try:
            GivensRotator(1,1,0.5)
            passed += [False]
        except ValueError:
            passed += [True]

        try:
            GivensRotator(0,1,0.5,dim=1)
            passed += [False]
        except ValueError:
            passed += [True]
        try:
            GivensRotator(4,1,0.5)
            passed += [False]
        except IndexError:
            passed += [True]
        try:
            GivensRotator(1,4,0.5)
            passed += [False]
        except IndexError:
            passed += [True]


            
        # test if phi is computd correctly 
        assert abs(rot3.getPhi()-pi/4) < eps
        # test set method
        rot4 = GivensRotator(0,1,0.0,copy=False)
        rot4.setCopy(True)
        assert rot4.copy
           
        #test correctness of matrix vector multiplication
        for dim in range(2,self.nr_tests//2):
            phis = (np.random.rand(self.nr_tests)-0.5)*2*pi
            for phi in phis:
                test_vals = [self.testGivensmatvec(i,j,phi,dim) 
                             for i in range(dim)
                             for j in range(dim)
                             if i != j
                             ]
                assert all(test_vals)
                passed += [True]

        # test matrix-matrix multiplication
        for dim in range(2,self.nr_tests//2):
            phis = (np.random.rand(self.nr_tests)-0.5)*2*pi
            for phi in phis:
                test_vals = [self.testGivensmatmat(i,j,phi,dim) 
                             for i in range(dim)
                             for j in range(dim)
                             if i != j
                             ]
                assert all(test_vals)
                passed += [True]

        # test transpose and inverse
        for dim in range(2,self.nr_tests//2):
            phis = (np.random.rand(self.nr_tests)-0.5)*2*pi
            for phi in phis:
                test_vals = [self.testGivensmattrans(i,j,phi,dim) 
                             for i in range(dim)
                             for j in range(dim)
                             if i != j]
                assert all(test_vals)
                passed += [True]

                
        self.checkTests("GivensRotator",passed)

    def testGivensRotationsQR(self,A):

        from Types import GivensRotations
        dim = A.shape[0]
        Q = GivensRotations(dim=dim)
        result = A
        for i in range(dim-1):
            for j in range(i+1,dim):
                result = Q.computeRotation(i,j,result)

        Qmat = Q(eye(dim))
        R = result
        return Q.transpose(), Qmat.transpose(), R 
        
    def testGivensRotations(self):
        from Types import GivensRotator
        from Types import GivensRotations
        pi = np.pi
        from numpy import sqrt, cos, sin, abs
        passed = [False]
        # test creation
        with warnings.catch_warnings(record=True) as warn:
            rotDummy = GivensRotations(dim=5,method=666)
            assert issubclass(warn[-1].category, UserWarning)
            assert self.fallback_warning in str(warn[-1].message)
            rotDummy.matvec
            rotDummy.computeRotation
            rotDummy.computeRotationParameters
            # test if initiated with unit matrix
            test_matrix = np.random.rand(*(rotDummy.shape))
            assert norm(rotDummy(test_matrix)-test_matrix) < eps
        passed[0] = True
        # test list creation
        try:
            liste = [GivensRotator(0,1,2,dim=2),GivensRotator(0,1,2,dim=3)]
            GivensRotations(liste)
            passed += [False]
        except ValueError:
            passed += [True]
        liste = [GivensRotator(0,1,2,copy=False),GivensRotator(0,1,2,copy=False)]
        rot1 = GivensRotations(liste,copy=True)
        liste2 = rot1.getRotations()
        liste2 = [G.copy for G in liste2]
        passed += [all(liste2)]
        # test setCopy method
        rot1.setCopy(False)
        liste2 = rot1.getRotations()
        passed += [not rot1.copy] + [not G.copy for G in liste2]

        #test inversion and transponation (and indirectly matvec mult)
        for dim in range(2,self.nr_tests//2):
            phis = np.random.rand(dim)
            liste = []
            liste = [GivensRotator(i,j,phis[i],dim=dim) for i in range(dim) 
                     for j in range(i,dim) if i != j]
            rot2 = GivensRotations(liste,dim=dim)
            rot2Inv = rot2.inv()
            rot2Tra = rot2.transpose()
            test_matrix = np.random.rand(*(rot2.shape))
            assert norm(rot2Inv(rot2(test_matrix))-test_matrix) < eps
            assert norm(rot2Tra(rot2(test_matrix))-test_matrix) < eps

        # check correctness of compute rotation method
        rot3 = GivensRotations(dim=2)
        Rs = np.random.rand(self.nr_tests)*10 + 1
        Phis = 2*pi*(np.random.rand(self.nr_tests)-0.5)
        for i in range(self.nr_tests):
            a = Rs[i]*cos(Phis[i])
            b = Rs[i]*sin(Phis[i])
            c,s,r = rot3.computeRotationParameters(a,b)
            assert abs(r - Rs[i]) < eps
            assert abs(c - cos(Phis[i]))<eps
            assert abs(s - sin(Phis[i]))<eps
            assert norm(dot(np.array([[c,s],[-s,c]]),np.array([a,b]))-np.array([r,0.0]))<eps
        passed += [True]

        for dim in range(2,self.nr_tests//2):
            for tester in range(self.nr_tests):
                A = np.random.rand(dim,dim)

                Q,Qmat,R = self.testGivensRotationsQR(A)
                
                assert norm(Q(R)-A) < eps
        self.checkTests("GivensRotations",passed)

    def testGivensQR(self):
        from Types import GivensQR
        givens_message =  "Warning: Dimension < 2! Q is scalar not GivensRotation!"
        passed = [False]
        with warnings.catch_warnings(record=True) as warn:
            dummy_givens_qr = GivensQR(method=666)
            assert issubclass(warn[-1].category, UserWarning)
            assert self.fallback_warning in str(warn[-1].message)
            passed[0] = True   

        givens_qr = GivensQR()
        # test special case for dim = 1
        with warnings.catch_warnings(record=True) as warn:
            A1 = np.array([[1.0]])
            Q1,R1 = givens_qr(A1)
            assert issubclass(warn[-1].category, UserWarning)
            assert givens_message in str(warn[-1].message)
            assert np.all(Q1 == np.array([[1.0]]))
            assert np.all(R1 == A1)
        
        with warnings.catch_warnings(record=True) as warn:
            A1 = np.array([[1.0,1.0,1.0]])
            Q1,R1 = givens_qr(A1)
            assert issubclass(warn[-1].category, UserWarning)
            assert givens_message in str(warn[-1].message)
            assert np.all(Q1 == np.array([[1.0]]))
            assert np.all(R1 == A1)

         
        for i in range(2,self.nr_tests):
            for j in range(2,self.nr_tests):
                A = np.random.rand(i,j)
                Q,R = givens_qr(A)
                assert norm(Q(R)-A) < eps
                
        self.checkTests("GivensQR",passed)

        
    def __init__(self,nr_tests):
        """
        Method for executing tests.
        """
        self.fallback_warning = "Warning: Make fallback since no other version is implemented!"
        self.nr_tests = nr_tests

        print("Start Type Tests...")
        self.testMathOperator()
        self.testGeometricTransformation()
        self.testPolarCoordinates()
        self.testCartesianCoordinates()
        self.testSphericalCoordinates()
        self.testCartesianCoordinates3D()
        self.testGivensRotator()
        self.testGivensRotations()
        self.testGivensQR()


# class ToolTests(object):

#     def checkTests(self,name,passed):
#         if not all(passed):
#             raise Exception("Error: " + name + " tests didn't passed!")
#         print("All " + name + " tests passed!")

    

#     def __init__(self,nr_tests):

#         self.nr_tests = nr_tests
        
#         print("Start Tool Tests ...")


nr_tests = 10
TypeTests(nr_tests)
# ToolTests(nr_tests)

