
## deals with formatting a triangulation for snappea.
import sys
import copy

"""
In SnapPy a positively oriented tetrahedra is one where the edges (30, 31, 32)
form a right handed coordinate system. The tetrahedra shape parameter in SnapPy
corresponds to the edge 23 (and 01).
23 <-> z
02 <-> 1/(1-z)
12 <-> (z-1)/z
"""

class Triangulation:
    # tetrahedra = list of Tetrahedron
    # cusps = number of cusps
    
    def __init__(self, n, tri_title=None):
        # tri_title is the title of the triangulation as written in the triangulation file
        self.tetrahedra = []
        self.cusps = 0
        for i in range(n):
            self.tetrahedra.append(Tetrahedron(i))
        self.tri_title = tri_title
            
    def addTetrahedron(self, num):
        tet = Tetrahedron(num)
        self.tetrahedra.append(tet)
        return tet
        
    def removeTetrahedron(self, num):
        # Removes Tetrahedron with getTetNo() == num
        for i in range(self.numTetrahedra()):
            if self.tetrahedra[i].getTetNo() == num:
                del self.tetrahedra[i]
                break
                
    def renumberTriangulation(self):
        for i in range(self.numTetrahedra()):
            self.tetrahedra[i].n = i
                
    def renumberTetrahedron(self, old, new):
        # Since tet.getNeighbour(i) returns a reference to a tetrahedron,
        # rather than a tetrahedron number, this is simple.
        self.getTetrahedronByNum(old).n = new
        
        
    def threeToTwo(self, tet_nos, t1_faces):
        # Assumes that the three tetrahedra do not externally glue to each other.
        # Sets all the new cusp nos to 0. This may not work well if there are multiple cusps.
        # t1_faces 
        nTri = self.copy()
        
        tet_nos[1], face = nTri.gluedFace(tet_nos[0], [t1_faces[0][0], t1_faces[1][0], t1_faces[1][1]])
        opVert = nTri.oppositeVertex(face)
        t2_faces = [[face[0], face[2], opVert], [face[1], opVert, face[2]]]
        print 't2_faces', t2_faces
        
        tet_nos[2], face = nTri.gluedFace(tet_nos[1], [t2_faces[0][0], t2_faces[1][0], t2_faces[1][1]])
        opVert = nTri.oppositeVertex(face)
        t3_faces = [[face[0], face[2], opVert], [face[1], opVert, face[2]]]
        print 't3_faces', t3_faces
        
        tet1 = nTri.getTetrahedronByNum(tet_nos[0])
        tet2 = nTri.getTetrahedronByNum(tet_nos[1])
        tet3 = nTri.getTetrahedronByNum(tet_nos[2])
        
        n1 = nTri.nextTetIndex()
        ntet1 = nTri.addTetrahedron(n1)
        
        n2 = nTri.nextTetIndex()
        ntet2 = nTri.addTetrahedron(n2)  
        
        nTri.switchFace(face)
        ntet1.joinTo(0, ntet2, [0,1,3,2])
        ntet2.joinTo(0, ntet1, [0,1,3,2])
        # def switchFace(self, tri):
        # def gluedFace(self, tet_no, tri):
        
        ntet1.setCuspNos([0]*4)
        ntet2.setCuspNos([0]*4)
        
        print 'z_new_tet1(01) = ', 'z_' + str(tet_nos[1]), str(t2_faces[0][0]) + str(t2_faces[0][2]), '* z_' + str(tet_nos[2]), str(t3_faces[0][0]) + str(t3_faces[0][1])                
        print 'z_new_tet2(01) = ', 'z_' + str(tet_nos[1]), str(t2_faces[1][0]) + str(t2_faces[1][1]), '* z_' + str(tet_nos[2]), str(t3_faces[1][0]) + str(t3_faces[1][2])                  
        opVert = nTri.oppositeVertex(t1_faces[0])
        g = tet1.getGluing(opVert)
        print 'neighbour:', tet1.getNeighbour(opVert).n
        ntet1.joinTo(1, tet1.getNeighbour(opVert),
                           [g[t1_faces[0][0]],
                           g[opVert],
                           g[t1_faces[0][1]],
                           g[t1_faces[0][2]]])
                           
        opVert = nTri.oppositeVertex(t2_faces[0])
        g = tet2.getGluing(opVert)
        print 'neighbour:', tet2.getNeighbour(opVert).n
        ntet1.joinTo(2, tet2.getNeighbour(opVert),
                           [g[t2_faces[0][0]],
                           g[t2_faces[0][2]],
                           g[opVert],
                           g[t2_faces[0][1]]])        

        opVert = nTri.oppositeVertex(t3_faces[0])
        g = tet3.getGluing(opVert)
        print 'neighbour:', tet3.getNeighbour(opVert).n
        ntet1.joinTo(3, tet3.getNeighbour(opVert),
                           [g[t3_faces[0][0]],
                           g[t3_faces[0][1]],
                           g[t3_faces[0][2]],
                           g[opVert]])
                           
        # tet 2                   
        opVert = nTri.oppositeVertex(t1_faces[1])
        g = tet1.getGluing(opVert)
        print 'neighbour:', tet1.getNeighbour(opVert).n
        ntet2.joinTo(1, tet1.getNeighbour(opVert),
                           [g[t1_faces[1][0]],
                           g[opVert],
                           g[t1_faces[1][1]],
                           g[t1_faces[1][2]]])
                                   
        opVert = nTri.oppositeVertex(t2_faces[1])
        g = tet2.getGluing(opVert)
        print 'neighbour:', tet2.getNeighbour(opVert).n
        ntet2.joinTo(3, tet2.getNeighbour(opVert),
                           [g[t2_faces[1][0]],
                           g[t2_faces[1][1]],
                           g[t2_faces[1][2]],
                           g[opVert]])        

        opVert = nTri.oppositeVertex(t3_faces[1])
        g = tet3.getGluing(opVert)
        print 'neighbour:', tet3.getNeighbour(opVert).n
        ntet2.joinTo(2, tet3.getNeighbour(opVert),
                           [g[t3_faces[1][0]],
                           g[t3_faces[1][2]],
                           g[opVert],
                           g[t3_faces[1][1]]])
               
        print tet_nos
        for i in range(3):
            nTri.removeTetrahedron(tet_nos[i])

        nTri.renumberTriangulation()
        return nTri
        
            
    def twoToThree(self, tet1_no, tet2_no, vert1, vert2, face1, face2):
        # Assumes that the pair of tetrahedra only glue along one face, i.e. no external gluings.
        # Returns a new triangulation object
        nTri = self.copy()
        
        tet1 = nTri.getTetrahedronByNum(tet1_no)
        tet2 = nTri.getTetrahedronByNum(tet2_no)
        
        n1 = nTri.nextTetIndex()
        ntet1 = nTri.addTetrahedron(n1)
        
        n2 = nTri.nextTetIndex()
        ntet2 = nTri.addTetrahedron(n2)     

        n3 = nTri.nextTetIndex()
        ntet3 = nTri.addTetrahedron(n3)
        
        # def joinTo(self, myFace, you, gluing):
        
        # set cusp numbers
        ntet1.setCuspNos([tet2.getCuspNo(vert2), tet1.getCuspNo(vert1), tet1.getCuspNo(face1[1]), tet1.getCuspNo(face1[0])])
        ntet2.setCuspNos([tet2.getCuspNo(vert2), tet1.getCuspNo(vert1), tet1.getCuspNo(face1[2]), tet1.getCuspNo(face1[1])])
        ntet3.setCuspNos([tet2.getCuspNo(vert2), tet1.getCuspNo(vert1), tet1.getCuspNo(face1[0]), tet1.getCuspNo(face1[2])])
        
        
        print 'z_new_tet1(23) = ', 'z_' + str(tet1_no), str(face1[0]) + str(face1[1]), '* z_' + str(tet2_no), str(face2[0]) + str(face2[1])                
        print 'z_new_tet2(23) = ', 'z_' + str(tet1_no), str(face1[1]) + str(face1[2]), '* z_' + str(tet2_no), str(face2[1]) + str(face2[2])     
        print 'z_new_tet3(23) = ', 'z_' + str(tet1_no), str(face1[2]) + str(face1[0]), '* z_' + str(tet2_no), str(face2[2]) + str(face2[0])
        
        # face pairings of first new tetrahedron
        g = tet1.getGluing(face1[2])
        ntet1.joinTo(0, tet1.getNeighbour(face1[2]), 
                    [g[face1[2]], 
                    g[vert1], 
                    g[face1[1]], 
                    g[face1[0]]]
                    )
                    
        g = tet2.getGluing(face2[2]) 
        ntet1.joinTo(1, tet2.getNeighbour(face2[2]), 
                    [g[vert2], 
                    g[face2[2]], 
                    g[face2[1]], 
                    g[face2[0]]]
                    )
        ntet1.joinTo(2, ntet3, [0,1,3,2])
        ntet1.joinTo(3, ntet2, [0,1,3,2])
        
        
        # face pairings of second new tetrahedron
        g = tet1.getGluing(face1[0])

        ntet2.joinTo(0, tet1.getNeighbour(face1[0]), 
                    [g[face1[0]], 
                    g[vert1], 
                    g[face1[2]], 
                    g[face1[1]]]
                    )
        g = tet2.getGluing(face2[0])     

        ntet2.joinTo(1, tet2.getNeighbour(face2[0]), 
                    [g[vert2], 
                    g[face2[0]], 
                    g[face2[2]], 
                    g[face2[1]]]
                    )
        ntet2.joinTo(2, ntet1, [0,1,3,2])
        ntet2.joinTo(3, ntet3, [0,1,3,2])
        
        # face pairings of third new tetrahedron
        g = tet1.getGluing(face1[1])

        ntet3.joinTo(0, tet1.getNeighbour(face1[1]), 
                    [g[face1[1]], 
                    g[vert1], 
                    g[face1[0]], 
                    g[face1[2]]]
                    )
        g = tet2.getGluing(face2[1])     

        ntet3.joinTo(1, tet2.getNeighbour(face2[1]), 
                    [g[vert2], 
                    g[face2[1]], 
                    g[face2[0]], 
                    g[face2[2]]]
                    )
        ntet3.joinTo(2, ntet2, [0,1,3,2])
        ntet3.joinTo(3, ntet1, [0,1,3,2])

        nTri.removeTetrahedron(tet1_no)
        nTri.removeTetrahedron(tet2_no)
        nTri.renumberTriangulation()
        return nTri
   
    def nextTetIndex(self):
        return max(map(lambda tet: tet.getTetNo(), self.tetrahedra))+1
     
    def copy(self):
        # Requires the tetrahedra in self to be numbered 0..(n-1)
        nTri = Triangulation(self.numTetrahedra())
        nTri.setNoCusps(self.getNoCusps())
        
        for i in range(self.numTetrahedra()):
            nTet = nTri.getTetrahedron(i)
            oldTet = self.getTetrahedron(i)
            nTet.setGluings(oldTet.getGluings())
            nTet.setCuspNos(oldTet.getCuspNos())
                
        for i in range(self.numTetrahedra()):
            nTet = nTri.getTetrahedron(i)
            oldTet = self.getTetrahedron(i)
            for j in range(4):
                neighbour_no = oldTet.getNeighbour(j).getTetNo()
                nTet.setNeighbour(j, nTri.getTetrahedronByNum(neighbour_no))
        return nTri
            
    def getTetrahedronByNum(self, n):
        for tet in self.getTetrahedraList():
            if tet.getTetNo() == n:
                return tet
                
    def getTetrahedraList(self):
        return self.tetrahedra
        
    # return i-th Tetrahedron in list of tetrahedra
    # this may be different to Tetrahedron.getTetNo() == i
    def getTetrahedron(self, i):
        return self.tetrahedra[i]
        
    def numTetrahedra(self):
        return len(self.tetrahedra)
        
    def setNoCusps(self, n):
        self.cusps = n
        
    def getNoCusps(self):
        return self.cusps
        
    def oppositeVertex(self, tri):
        for i in range(4):
            if i not in tri:
                return i
        
    def switchFace(self, tri):
        """ if tri = [a,b,c], then want [a,b,d] """
        return tri[:2] + [self.oppositeVertex(tri)]
                
    def gluedFace(self, tet_no, tri):
        tet = self.tetrahedra[tet_no]
        face = self.oppositeVertex(tri)
        g = tet.getGluing(face)
        return (tet.getNeighbour(face).n, [g[tri[0]], g[tri[1]], g[tri[2]]])
    
    def rotateFace(self, tet_no, tri):
        """ tri = [a,b,c], then a-b is the special edge """
        return self.gluedFace(tet_no, self.switchFace(tri))
        
    def writeSnapPea(self, foutname):
        f = open(foutname, 'w')
        f.write('% Triangulation\n')
        if self.tri_title is not None:
            f.write(self.tri_title + '\n')
        else:
            f.write('Traintrack_Triangulation\n')
        f.write('not_attempted 0.0\n')
        f.write('unknown_orientability\n')
        f.write('CS_unknown\n')
        f.write('{0} 0\n'.format(self.cusps))
        for i in range(self.cusps):
            f.write('\ttorus\t0.0\t0.0\n')
        f.write(str(self.numTetrahedra()) + '\n')
        
        for i in range(self.numTetrahedra()):
            tet = self.getTetrahedron(i)
            f.write('\t{0}\t{1}\t{2}\t{3}\n'.format(tet.getNeighbour(0).getTetNo(),tet.getNeighbour(1).getTetNo(),
                                                    tet.getNeighbour(2).getTetNo(),tet.getNeighbour(3).getTetNo()))
            
            for j in range(4):
                g = tet.getGluing(j)
                f.write('\t{0}{1}{2}{3}'.format(g[0],g[1],g[2],g[3]))
            f.write('\n')
            f.write('\t{0}\t{1}\t{2}\t{3}\n'.format(tet.getCuspNo(0),tet.getCuspNo(1),tet.getCuspNo(2),tet.getCuspNo(3)))
            for k in range(4):
                f.write('\t0'*16 + '\n')
            f.write('0.0\t0.0\n')
            


class Tetrahedron:
    # Tetrahedron neighbour[0-3]
    # Gluings gluing[0-3][4] where gluing[i][j] pairs face i and vertex j to tet neighbour[i] vertex gluing[i][j]
    # Cusp numbers corresponding to vertices cusp_no[0-3]
    
    def __init__(self, n = 0):
        self.neighbour = [0]*4
        self.gluing = [[]]*4
        self.cusp_no = [-1]*4
        self.n = n
        
    def joinTo(self, myFace, you, gluing):
        self.setNeighbour(myFace, you)
        self.setGluing(myFace, gluing)
        
        you.setNeighbour(gluing[myFace], self)
        
        you_gluing = [0]*4
        for i in range(4):
            you_gluing[gluing[i]] = i
            
        you.setGluing(gluing[myFace], you_gluing)
        
    def setNeighbour(self, myFace, you):
        self.neighbour[myFace] = you
        
    def setGluings(self, gluings):
        self.gluing = copy.deepcopy(gluings)
        
    def setGluing(self, myFace, gluing):
        self.gluing[myFace] = gluing
        
    def setCuspNo(self, vert, n):
        self.cusp_no[vert] = n
        
    def setCuspNos(self, cusps):
        self.cusp_no = copy.deepcopy(cusps)
        
    def getNeighbour(self, myFace):
        if type(self.neighbour[myFace]) == int:
            print 'returning int, tet:', self.n
        return self.neighbour[myFace]
        
    def getGluing(self, myFace):
        return self.gluing[myFace]
        
    def getGluings(self):
        return copy.deepcopy(self.gluing)
        
    def getCuspNos(self):
        return self.cusp_no[:]
        
    def getCuspNo(self, vert):
        return self.cusp_no[vert]
        
    def getTetNo(self):
        return self.n

def SnappyToEquation(vecs):
    """ takes a snappy gluing equation (form="rect") and writes it out as an actual equation."""
    eqns = []
    for vec in vecs:
        A, B, c = vec
        eqn = ""
        for i in range(len(A)):
            if A[i] != 0:
                eqn += "* z{0}^{1}".format(i,A[i],B[i])
            if B[i] != 0:
                eqn += "* (1 - z{0})^{2}".format(i,A[i],B[i])
        eqn = eqn[2:] # remove leading '* '
        eqn += " = " + str(c)
        eqns.append(eqn)
    return eqns
    
    