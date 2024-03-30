
# triangulation_builder
import sys
import math
import os
import string

sys.path.append('/mount/autofs/home_stude/aiissa/Research/bigdecimal version')
sys.path.append('C:/pr2011/Research/Agol Construction/bigdecimal version')
from snappea_triangulation import *
from train_track_parser import *

# Requires xlwt package to handle excel files


def tri_around_switch(tt, s_no):
    switch = tt.switches[s_no]
    if len(switch.branches[0]) == 1:
        return [abs(x) for x in switch.branches[0] + switch.branches[1]]
    else:
        return [abs(x) for x in switch.branches[1] + switch.branches[0]]
        
def oriented_tri_around_switch(tt, s_no):
    switch = tt.switches[s_no]
    if len(switch.branches[0]) == 1:
        return [x for x in switch.branches[0] + switch.branches[1]]
    else:
        return [x for x in switch.branches[1] + switch.branches[0]]

def dual_triangles(tt):
    tris = []
    for switch in tt.switches.values():
        tris.append(tri_around_switch(tt, switch.n))
    return tris
    
def signed_dual_triangles(tt):
    tris = []
    for switch in tt.switches.values():
        tris.append(oriented_tri_around_switch(tt, switch.n))
    return tris
    
def orient_triangle(signed_tris, tri):
    for t in signed_tris:
        if equal_tri([abs(x) for x in t], tri):
            return t
    
def diag_exch_tris(tt, b_no):
    init_no, end_no = tt.branches[b_no].b_ends # branch terminal switch numbers
    return [oriented_tri_around_switch(tt, init_no), oriented_tri_around_switch(tt, end_no)]
    
def equal_tri(tri_1, tri_2):
    ctri_1 = tri_1[:]
    for i in range(3):
        if ctri_1 == tri_2:
            return True
        ctri_1.append(ctri_1.pop(0))
    return False
    
    
class Layer:
    def __init__(self):
        # each element of tris in form: layer_tri, (tet_no, tet)
        self.tris = []
        
    def add_tri(self, tri):
        self.tris.append([tri])
        
    # returns None if can't find
    def get_tri(self, base_tri):
        for tri in self.tris:
            if equal_tri(base_tri, tri[0]):
                return tri
        return None
        
    def remove_tri(self, tri):
        for i in range(len(self.tris)):
            if equal_tri(tri, self.tris[i][0]):
                del self.tris[i]
                return None
        
    # gets tetrahedra associated with tri
    # returns None if no tetrahedra associated with tri
    # returns None if more than 1 tetrahedra
    def get_tet(self, tri):
        actual_tri = self.get_tri(tri)
        if actual_tri == None or len(actual_tri) != 2:
            return None
        layer_tri, (tet_no, tet) = actual_tri
        c_tet = tet[:]
        c_layer_tri = layer_tri[:]
        for i in range(3):
            if tri == c_layer_tri:
                return (tet_no, c_tet)
            else:
                c_tet.append(c_tet.pop(0))
                c_layer_tri.append(c_layer_tri.pop(0))
        
    def associate(self, base_tri, tet_no, tet_tri):
        equiv_base_tri = self.get_tri(base_tri)[0]
        c_base_tri = base_tri[:] # create copies
        c_tet_tri = tet_tri[:] # copy
        for i in range(3):
            if c_base_tri != equiv_base_tri:
                c_base_tri.append(c_base_tri.pop(0))
                c_tet_tri.append(c_tet_tri.pop(0))
            else:
                for t in self.tris:
                    if c_base_tri == t[0]:
                        t.append((tet_no, c_tet_tri))
                        return None
        print 'error in associating faces'
                
    def __str__(self):
        return 'layer: ' + str(self.tris)
        
# returns a list of the other two triangles after the diagonal exchange
def flipped_tris(tris, tt, b_no):
    t1, t2 = tris

    init_no, end_no = tt.branches[b_no].b_ends # branch terminal switch numbers
    b1 = tt.switches[init_no].small_branches()
    b2 = tt.switches[end_no].small_branches()
    m = max(tt.branches[abs(b1[0])].weight, tt.branches[abs(b1[1])].weight, 
            tt.branches[abs(b2[0])].weight, tt.branches[abs(b2[1])].weight)
            
    # the rest of this is to decide whether to put t1[0] or t2[0] in the first slot of the tri
    # these are now SIGNED, so it matters
    if tt.branches[abs(b1[0])].weight == m or tt.branches[abs(b2[0])].weight == m:
        return  [[t1[0], t2[2], t1[1]], [t2[0],t1[2],t2[1]]]
    else:
        return  [[t2[0], t2[2], t1[1]], [t1[0],t1[2],t2[1]]]
        
# inserts a new tetrahedron and associates faces with layer1 and layer2 triangles
# orientation convention of the tetrahedra is determined here.
# We're sticking with the old convention.
def update_associations(layer1, layer2, tt, b_no, tet_no):
    tris = diag_exch_tris(tt, b_no)
    
    # old convention: (02, 01, 03) form right-handed coordinate system
    layer1.associate(tris[0], tet_no, [1,2,3])
    layer1.associate(tris[1], tet_no, [3,0,1])
    
    # new convention: (01, 02, 03) form right-handed coordinate system
    # agrees with snappy
    #layer1.associate(tris[0], tet_no, [1,3,2])
    #layer1.associate(tris[1], tet_no, [2,0,1])
    
    f_tris = flipped_tris(tris, tt, b_no)
    
    #old convention
    layer2.associate(f_tris[0], tet_no, [0,1,2])
    layer2.associate(f_tris[1], tet_no, [2,3,0])
    
    # new convention
    #layer2.associate(f_tris[0], tet_no, [0,1,3])
    #layer2.associate(f_tris[1], tet_no, [3,2,0])
        
# base and final are layers
def monod_associations(layers, punc_bij):
    for t in layers[0].tris:
        cor_tri = map(lambda x: punc_bij[x], t[0])
        for i in range(len(layers)-1):
            tet = layers[len(layers) - i - 1].get_tet(cor_tri)
            if tet != None:
                t.append(tet)
                break
                
def opposite_vertex(tri):
    for i in range(4):
        if i not in tri:
            return i
            
def perm_bw_tris(tri1, tri2):
    perm = [0,0,0,0]
    perm[opposite_vertex(tri1)] = opposite_vertex(tri2)
    for i in range(3):
        perm[tri1[i]] = tri2[i]
    return perm
    
def inverse_perm(p):
    p_inv = [0]*len(p)
    for i in range(len(p)):
        p_inv[abs(p[i])] = i if p[i] >= 0 else -i
    return p_inv
                

def build_triangulation(layers, emap, split_seq, sing_to_cusp, foutname=None, tri_title=None):
    associated_layer = [] # contains [tri-in-layer-zero, above-tet, below-tet]
    no_tet = len(layers)-1
    triang = Triangulation(no_tet, tri_title=tri_title)
    triang.setNoCusps(1+max(sing_to_cusp)) # set the number of cusps, used by snappea
    joins = 0
        
    for i in range(len(layers)-1):
        layer = layers[i]

        for t in layer.tris:
            ind = i
            found = False
            tris_layer0 = []
                
            if len(t) == 3: # associates two faces together
                (tet1_no, tri1), (tet2_no, tri2) = t[1:]
                if i == 0: # connects through layer 0
                    associated_layer.append(t)
                found = True
            elif len(t) == 2: # it's got 1 tetrahedron
                tet1_no, tri1 = t[1]
                # loop through the rest of the layers to find an associated triangle
                found = False
                
                # two situations, we either search forwards (upwards) or backwards (downwards)
                # if tet1_no == i then we search backwards
                
                if tet1_no == i:
                    t0 = t[0]
                    ind = i
                    while True: # going down
                        if ind == 0:
                            tris_layer0.append(t0)
                            
                        ind -= 1
                        
                        if ind == -1: 
                            # factor by monodromy
                            t1 = t0[:]
                            for j in range(len(t0)):
                                if t0[j] >= 0:
                                    t1[j] = emap[t0[j]]
                                else:
                                    t1[j] = -emap[-t0[j]]
                            t0 = t1
                            ind += len(layers)

                        tet2 = layers[ind].get_tet(t0)
                        if tet2 != None: # found the matching pair
                            tet2_no, tri2 = tet2
                            t.append((tet2_no, tri2))
                            layers[ind].remove_tri(t0)
                            found = True
                            for tri0 in tris_layer0:
                                associated_layer.append([tri0, (tet1_no, tri1), (tet2_no, tri2)])
                            break
                else:
                    t0 = t[0]
                    ind = i
                    if ind == 0:
                        tris_layer0.append(t0)
                    while True:
                        # going up
                        
                        ind += 1
                        if ind == len(layers):
                            inv_bij = inverse_perm(emap)
                            t1 = t0[:]
                            for j in range(len(t0)):
                                if t0[j] >= 0:
                                    t1[j] = inv_bij[t0[j]]
                                else:
                                    t1[j] = -inv_bij[abs(t0[j])]
                            t0 = t1
                            ind = 0
                        if ind == 0:
                            tris_layer0.append(t0)
                            
                        tet2 = layers[ind].get_tet(t0)
                        if tet2 != None: # found the matching pair
                            tet2_no, tri2 = tet2
                            t.insert(1, (tet2_no, tri2))
                            layers[ind].remove_tri(t0)
                            found = True
                            for tri0 in tris_layer0:
                                associated_layer.append([tri0, (tet2_no, tri2), (tet1_no, tri1)])
                            break
            else:
                continue
            # glue (tet1_no, tri1) to (tet2_no, tri2)
            if found:
                p = perm_bw_tris(tri1, tri2)
                # 'associating:', (tet1_no, tri1), 'to', (tet2_no, tri2)
                triang.getTetrahedron(tet1_no).joinTo(opposite_vertex(tri1),
                                triang.getTetrahedron(tet2_no), [p[0],p[1],p[2],p[3]])
                branches = t[0] # e.g. [15,10,4]
                tt = split_seq[i][0] # get the train track

                oriented_branches = branches
                tet1 = triang.getTetrahedron(tet1_no)
                tet2 = triang.getTetrahedron(tet2_no)
                for m in range(3):
                    cusp_no = sing_to_cusp[tt.left_puncture(oriented_branches[m])]
                    tet1.setCuspNo(tri1[m], cusp_no)
                    tet2.setCuspNo(tri2[m], cusp_no)
                joins += 1
                                
    if foutname is not None:
        triang.writeSnapPea(foutname)
        print 'written to file:', foutname
    return (triang, associated_layer)
              
def trains_to_snappea_tri(fname=None, foutname='triang.tri', lines=False, debug=False, tri_title=None):
    # lines is a list of lines of Trains output, otherwise fname is the filename of trains output
    # tri_title is the title of the triangulation as saved in the snappea file
    # it's convenient to have this include the surface type and monodromy
    if lines:
        tt, growth_rate, punc_bij = parse_train_track(lines=lines, debug=debug)
    else:
        tt, growth_rate, punc_bij = parse_train_track(fname, debug=debug)
    initial_tt = tt
    
    if debug:
        tt.print_traintrack()
        print 'number of strands originally:', tt.strands
        
    sing_to_cusp = tt.cusp_info(punc_bij)
    
    if debug:
        print 'cusp info:', sing_to_cusp
        
    split_seq = []
    res = splitting_sequence(tt, growth_rate, punc_bij)
    emap = res[1][0] # take the first equivalence
    
    if debug:
        print 'edge map:', emap
        
    for ss in res[0]:
        split_seq.extend(ss)

    layers = []
    for i in range(len(split_seq)):
        layers.append(Layer())
        for tri in signed_dual_triangles(split_seq[i][0]):
            layers[i].add_tri(tri)
        
    for i in range(len(split_seq)-1):
        tt, b_no = split_seq[i]
        update_associations(layers[i], layers[i+1], tt, b_no, i)

    build_triangulation(layers, emap, split_seq, sing_to_cusp, foutname, tri_title=tri_title)
    
    return (initial_tt, growth_rate, punc_bij, res)
    
def all_trains_to_snappea_tris(fname, foutname, lines=False, debug=False):
    # foutname shouldn't have the 'tri' extension
    # if lines is true, then fname is a list of lines in memory which represents trains output
    if lines:
        tt, growth_rate, punc_bij = parse_train_track(lines=fname, debug=debug)
    else:
        tt, growth_rate, punc_bij = parse_train_track(fname, debug=debug)
    initial_tt = tt
    
    if debug:
        tt.print_traintrack()
        print 'number of strands originally:', tt.strands
        
    sing_to_cusp = tt.cusp_info(punc_bij)
    
    if debug:
        print 'cusp info:', sing_to_cusp
        
    split_seq = []
    res = splitting_sequence(tt, growth_rate, punc_bij, debug=debug)

    count = 1
    for emap in res[1]:
        split_seq = []
        if debug:
            print 'edge map:', emap
        for ss in res[0]:
            split_seq.extend(ss)

        layers = []
        for i in range(len(split_seq)):
            layers.append(Layer())
            for tri in signed_dual_triangles(split_seq[i][0]):
                layers[i].add_tri(tri)
        
        for i in range(len(split_seq)-1):
            tt, b_no = split_seq[i]
            update_associations(layers[i], layers[i+1], tt, b_no, i)

        build_triangulation(layers, emap, split_seq, sing_to_cusp, foutname + '_' + str(count) + '.tri')
        count += 1
    return (initial_tt, growth_rate, punc_bij, res)
    
def rotate(tri):
    tri.append(tri.pop(0))
    
def surTriangles(tt, n):
    """ returns a list of signed triangles surrounding puncture n in traintrack tt in anti-clockwise order,
        each of the triangles is oriented so that if [a,b,c] is the oriented triangle, then
        the puncture is to the left of edge c.
    """
    sub_tris = []
    b_nos = tt.sur_branches(n) # clockwise around puncture
    
    b_nos.reverse()
    b_nos = [-x for x in b_nos]
    # b_nos now anti-clockwise around puncture
    for b_no in b_nos:
        switch = tt.get_terminal_switch(b_no)
        tri = switch.all_branches()
        while tri[0] != -b_no:
            rotate(tri)
        sub_tris.append(tri)
    return sub_tris
    

def tetAboveTris(tt, growth_rate, punc_bij, punc, ss_res=None, emap_index=0):
    """ returns a list of pairs (tet_no, [e1,e2,e3]) 
        ei = 0,1,2,3
        where e1 corresponds to the puncture we're focused on
        e1-e2-e3 is a face of the tetrahedron in the fibre, ccw order looking above
        we'll rotate around the edge e which isn't in [e1,e2,e3] so that
        [e1,e2,e] gets rotated to [e1,e3,e]
        
        Caution: punc is the puncture number IN THE SURFACE of the train track
        NOT the puncture i.e. cusp of the manifold.
    """
    sing_to_cusp = tt.cusp_info(punc_bij)
    split_seq = []
    if ss_res is not None:
        res = ss_res
    else:
        res = splitting_sequence(tt, growth_rate, punc_bij)
    emap = res[1][emap_index]

    for ss in res[0]:
        split_seq.extend(ss)

    layers = []
    for i in range(len(split_seq)):
        layers.append(Layer())
        for tri in signed_dual_triangles(split_seq[i][0]):
            layers[i].add_tri(tri)
             
    for i in range(len(split_seq)-1):
        tt, b_no = split_seq[i]
        update_associations(layers[i], layers[i+1], tt, b_no, i)
        
    sur_tris = surTriangles(split_seq[0][0], punc)
    # print 'sur_tris:', sur_tris # [7,4,-6] would mean puncture interested in to the left of -6
    triang, assoc_layer = build_triangulation(layers, emap, split_seq, sing_to_cusp)
    sur_tet = []
    
    for sur_tri in sur_tris:
        for assoc in assoc_layer:
            tri, (tet1_no, tri1), (tet2_no, tri2) = assoc
            c_tri = tri[:]
            c_tri1 = tri1[:]
            c_tri2 = tri2[:]
            if equal_tri(sur_tri, tri):
                while sur_tri[0] != c_tri[0]:
                    c_tri.append(c_tri.pop(0))
                    c_tri1.append(c_tri1.pop(0))
                    c_tri2.append(c_tri2.pop(0))
                c_tri1 = c_tri1[-1:] + c_tri1[:2] # we want the main edge we're going to rotate around to be the first two elements
                sur_tet.append((tet1_no, c_tri1))
                
    return (triang, sur_tet)
    
def meridianCurve(triang, sur_tet):
    e_inv = [] # sequence of edge invariants
    for i in range(len(sur_tet)):
        tet_no, layer_face = sur_tet[i]
        e_inv.append((tet_no, [layer_face[0], triang.oppositeVertex(layer_face)], -1))
        face = [layer_face[0], triang.oppositeVertex(layer_face), layer_face[1]]
        ntet_no, nlayer_face = sur_tet[(i+1) % len(sur_tet)]
        
        rot_tet_no, rot_face = triang.rotateFace(tet_no, face)
        rot_face = [rot_face[0], rot_face[2], rot_face[1]]

        while True:
            if (ntet_no == rot_tet_no and triang.oppositeVertex(rot_face) == nlayer_face[2]
                and rot_face[:2] == nlayer_face[:2]):
                break
            else:
                e_inv.append((rot_tet_no, rot_face[:2], 1))
                rot_tet_no, rot_face = triang.rotateFace(rot_tet_no, rot_face)

    return e_inv
    
def toSnappyVector(num_tet, e_invs):
    vec = [0] * (3*num_tet)
    A = [[2,3],[3,2],[1,0],[0,1]]
    B = [[0,2],[2,0],[1,3],[3,1]]
    C = [[1,2],[2,1],[0,3],[3,0]]
    for e_inv in e_invs:
        # e looks like (5, [3,0], -1)
        tet_no, edge, sign = e_inv
        if edge in A:
            ind = 0
        elif edge in B:
            ind = 1
        elif edge in C:
            ind = 2
        if sign == 1:
            vec[tet_no*3+ind] += 1
        else:
            vec[tet_no*3+ind] -= 1
    return vec
    
def computeMeridian(tt, growth_rate, punc_bij, punc, ss):
    triang, sur_tet = tetAboveTris(tt, growth_rate, punc_bij, punc, ss, 0)
    return toSnappyVector(triang.numTetrahedra(), meridianCurve(triang, sur_tet))
    
    
#trains_to_snappea_tri('C:/pr2011/veering_code/two_punc.txt', 
#                    'C:/pr2011/veering_code/two_punc.tri', debug=False)

# def trains_to_snappea_tri(fname, foutname, lines=False, debug=False):