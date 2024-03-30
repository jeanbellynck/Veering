
import string
import subprocess # needed to run a Trains background process
import os, sys
from train_track import TrainTrack

class GeneralGraphMap:
    # n_periph, number of peripheral loops
    # n_edges, number of edges
    # n_verts, number of vertices
    # ribbon = [[edges around vert 1], [edges around vert 2], ...]
    # edge_map = [[image of edge 1], [image of edge 2], ...]
    # vert_map = [image of vert 1, image of vert 2, ...]
    # periph = [None if edge i is not periph otherwise any integer value, ...]
    
    # edge_puncs = [[left punc/region, right punc/right of edge 1], ...]
    #               Negative values for the punc/region indicates that it's a region with no puncture
    # n_puncs = number of punctures
    # punc_perm = [image of punc 1, image of punc 2, ...]
    # edge_term = [[initial vert of edge 1, final vert of edge 1], ...] (NOT SURE IF WE NEED THIS YET)
    
    # we want routines which takes a graph map then "fills in" edge_puncs and punc_perm for us
    
    def shorten_isotopy(self):
        for i, r in enumerate(self.ribbon):
            if len(r) <= 1:
                continue
            shortened = True
            while shortened:
                shortened = False
                found_empty = False
                for e_no in r:
                    if len(self.getEdgeImage(e_no)) == 0:
                        found_empty = True
                        break
                if found_empty:
                    break
                    
                starts = [self.getEdgeImage(e_no)[0] for e_no in r]
                if len(set(starts)) == 1: # if init edge image all the same
                    self.vert_map[i] = self.getTermVert(starts[0])
                    for e_no in r:
                        im = self.getEdgeImage(e_no)
                        del im[0]
                        self.setEdgeImage(e_no, im)
                    shortened = True
        
    
    def inferVertexMap(self):
        """ Use the images of edges to deduce the vertex map. """
        for i in range(self.n_verts):
            e = self.ribbon[i][0]
            e_image = self.getEdgeImage(e)[0]
            self.vert_map[i] = self.getInitVert(e_image)
            #print 'vert',i,'has edge',e, 'maps to',e_image,'with init vert', self.vert_map[i] 
    
    def inferPeripheral(self):
        """ Uses the permutations of punctures to determine appropriate peripheral
            loops for the BH algorithm. All but one orbit of punctures must consist of
            peripheral loops.
        """
        self.periph = [None for i in range(self.n_edges)]
        orbit = 0
        
        # we're going to label elements of self.periph with orbit numbers
        for i in range(len(self.punc_perm)):
            if self.periph[i] is not None:
                continue
            j = i
            while True:
                self.periph[j] = orbit
                j = self.punc_perm[j]
                if j == i:
                    orbit += 1
                    break

        for i in range(len(self.periph)):
            if self.periph[i] is not None and self.periph[i] == 0:
                self.periph[i] = None
                
        self.n_periph = sum([1 for x in self.periph if x is not None])
        
    def surLoops(self):
        """ returns a dictionary dic indexed by puncture numbers, so that
            dic[punc_no] gives a set of edges forming a loop going counter clockwise around
            the puncture.
        """
        dic = {}
        for e_no in range(1, self.n_edges+1):
            for i in range(2):
                punc = self.edge_puncs[e_no-1][i]
                if punc not in dic:
                    dic[punc] = self.followEdge(e_no*(-1)**i)
        return dic
                
        
    def inferEdgePuncs(self, respect_original_puncs=False):
        """Assigns puncture numbers to each complementary region and updates edge_puncs.
           If respect_original_puncs is true then we only set complementary regions as
           regions and not punctures, and we start from -1.
        """
        cur_punc_no = 0
        if not respect_original_puncs:
            self.edge_puncs = [[None,None] for i in range(self.n_edges)] # reset them
            increment = 1
        else:
            cur_punc_no = -1
            increment = -1

        for e_no in range(1,self.n_edges+1):
            if self.edge_puncs[e_no-1][0] is None: # if it hasn't been set, then set it
                sur = self.followEdge(e_no)

                for e in sur:
                    self.setLeftEdgePunc(e, cur_punc_no)

                cur_punc_no += increment
                
            if self.edge_puncs[e_no-1][1] is None:
                sur = self.followEdge(-e_no)
                for e in sur:
                    self.setRightEdgePunc(-e, cur_punc_no)

                cur_punc_no += increment
                
        if not respect_original_puncs:
            self.n_puncs = cur_punc_no
            self.punc_perm = [i for i in range(self.n_puncs)]
        if self.punc_perm == []:
            self.n_puncs = 1+max([max(p[0], p[1]) for p in self.edge_puncs])
            self.punc_perm = range(self.n_puncs)
        
    def getPuncturePerm(self):
        return self.punc_perm[:]
        
    def inferPuncPerm(self):
        """Uses the graph map i.e. edge maps to determine the permutation of punctures
            contained in complementary regions."""
            
        # for each puncture enumerate the surrounding edges
        regions = [] # we assume these are punctures labelled 0...(n_puncs-1)
        all_sur_edges = [] # list of (lists of surrounding edges), all_sur_edges[i] surrounds regions[i]
        
        for e_no in range(1, self.n_edges+1):
            left,right = self.edge_puncs[e_no-1]
            if left not in regions and left >= 0:
                all_sur_edges.append(self.followEdge(e_no))
                regions.append(left)
            if right not in regions and right >= 0:
                all_sur_edges.append(self.followEdge(-e_no))
                regions.append(right)
        #print 'regions:',regions
        image_regions = []
        all_sur_edges = [TrainTrack.cyclic_pulltight(x) for x in all_sur_edges]
        #print 'all_sur_edges:', all_sur_edges
        for i in range(len(all_sur_edges)):
            image = []
            for e_no in all_sur_edges[i]:
                image.extend(self.getEdgeImage(e_no))
            image = TrainTrack.cyclic_pulltight(image)
            found = False
            for j in range(len(all_sur_edges)):
                #print 'image',image,'all_sur_edges[j]',all_sur_edges[j]
                if TrainTrack.cyclic_equality(image, all_sur_edges[j]):
                    self.punc_perm[regions[i]] = regions[j]
                    found = True
                    break
            if not found:
                raise Exception('WARNING: Could not find image puncture of region '+str(i))
    
    def followEdge(self, e):
        """returns a list of the edges surrounding a region in ccw order
            starting at edge e"""
        lst = []
        e_puncs = self.getEdgePuncs()
                
        def nextEdge(cur):
            # cur is a signed edge number
            # returns the edge number at the (signed) terminal
            # edge of cur, bounding the region to the left of cur

            for r in self.ribbon:
                for j in range(len(r)):
                    if -cur == r[j]:
                        return r[j-1]
                        
        lst.append(e)
        while nextEdge(lst[-1]) != lst[0]:
            lst.append(nextEdge(lst[-1]))
        return lst

    def dehnTwist(self, edges, right=False):
        """ edges is the sequence of edges which form a loop in the graph map,
            in which a left twist is done. Here by left twist we mean the left side
            of the oriented loop formed by the edges is twisted in the direction
            of the oriented loop.
            
            If right = True then performs a twist to the right.
        """
        def between(e1, e2):
            # returns edges between e2 and e1 in cw order
            v = self.getInitVert(e2)
            r = self.getRibbon(v)
            for i in range(len(r)):
                if r[0] == e2:
                    for j in range(len(r)):
                        if r[j] == e1:
                            return r[1:j]
                else:
                    r.append(r.pop(0))
        
        for i in range(-1,len(edges)-1):
            edges_inv = between(-edges[i], edges[i+1])

            loop = edges[i+1:]+edges[:i+1]
            if not right:
                loop.reverse()
                loop = [-x for x in loop]
            loop_im = []
            for e in loop:
                loop_im.extend(self.getEdgeImage(e))
            for e in edges_inv:
                im = self.getEdgeImage(e)
                self.setEdgeImage(e, loop_im + im)
            
    
    def surEdges(self, region):
        """returns a list of the edges surrounding region in ccw order"""
        lst = []
        
        e_puncs = self.getEdgePuncs()
        # find a first edge, put its number in e
        for i in range(self.n_edges):
            if e_puncs[i][0] == region:
                e = i+1
                break
            elif e_puncs[i][1] == region:
                e = -(i+1)
                break
                
        return self.followEdge(e)
    
    def removeEdge(self, n):
        ind = n-1
        if self.periph[ind] is not None:
            self.n_periph -= 1
        del self.periph[ind]
            
        self.n_edges -= 1

        # delete edge from ribbon
        for r in self.ribbon:
            for i in reversed(range(len(r))):
                if abs(r[i]) == n:
                    del r[i]
                elif r[i] > n: # shift indices to reflect deletion of edge
                    r[i] -= 1
                elif r[i] < -n:
                    r[i] += 1
                   
        # shift indices in image of edges
        del self.edge_map[ind]
        for e_image in self.edge_map:
            for i in range(len(e_image)):
                if e_image[i] > n:
                    e_image[i] -= 1
                elif e_image[i] < -n:
                    e_image[i] += 1
        
        del self.edge_puncs[ind]
        
    def removeEmptyVerts(self):
        for i in reversed(range(len(self.ribbon))):
            if len(self.ribbon[i]) == 0:
                del self.ribbon[i]
                del self.vert_map[i]
                self.n_verts -= 1
    
    def __init__(self, n_periph=0, n_edges=0, n_verts=0, ribbon=None, edge_map=None, 
                        vert_map=None, periph=None, edge_puncs=None, n_puncs=0, punc_perm=None):
        self.n_periph = n_periph
        self.n_edges = n_edges
        self.n_verts = n_verts
        
        if edge_puncs is not None:
            self.edge_puncs = edge_puncs
        else:
            self.edge_puncs = []
            
        self.n_puncs = n_puncs
        if punc_perm is not None:
            self.punc_perm = punc_perm
        else:
            self.punc_perm = []
            for i in range(self.n_puncs):
                self.punc_perm.append(i)
        
        if ribbon is not None:
            self.ribbon = ribbon
        else:
            self.ribbon = []
            for i in range(n_verts):
                self.ribbon.append([])
                
        if edge_map is not None:
            self.edge_map = edge_map
        else:
            self.edge_map = []
            for i in range(n_edges):
                self.edge_map.append([i+1])
                
        if vert_map is not None:
            self.vert_map = vert_map
        else:
            self.vert_map = []
            for i in range(n_verts):
                self.vert_map.append(i+1)
                
        if periph is not None:
            self.periph = periph
        else:
            self.periph = [None]*n_edges
            
    def setEdgeImage(self, n, edge_image):
        if n > 0:
            self.edge_map[n-1] = edge_image
        else:
            edge_image_c = [-x for x in edge_image]
            edge_image_c.reverse()
            self.edge_map[(-n)-1] = edge_image_c
        
    def getEdgePuncs(self):
        return self.edge_puncs
        
    def setLeftEdgePunc(self, n, punc_no):
        # note a negative punc_no is treated as a region
        # allows n negative
        if n > 0:
            self.edge_puncs[n-1][0] = punc_no
        elif n < 0:
            self.edge_puncs[(-n)-1][1] = punc_no
            
    def setRightEdgePunc(self, n, punc_no):
        self.setLeftEdgePunc(-n, punc_no)
        
    def getEdgeImage(self, n):
        if n > 0:
            return self.edge_map[n-1][:]
        elif n < 0:
            im = self.edge_map[abs(n)-1][:]
            im.reverse()
            im = [-x for x in im]
            return im
        
    def setRibbon(self, n, cyc):
        self.ribbon[n] = cyc
        
    def getRibbon(self, vert_no=None):
        if vert_no is None:
            return self.ribbon
        else:
            return self.ribbon[vert_no-1][:]
        
    def getInitVert(self, edge):
        for i in range(len(self.ribbon)):
            if edge in self.ribbon[i]:
                return (i+1)
        return None # failed
        
    def getTermVert(self, edge):
        return self.getInitVert(-edge)
        
    def retract(self):
        """ Uses left/right region information of each edge to turn this graph map
            into one in which every complementary region is homeomorphic to a
            punctured disk. This is done by successively retracting edges which bound
            different regions to its left and right (at least one of which is not a punctured disk).
        """
        for i in reversed(range(self.n_edges)):
            if self.periph[i] is not None: # we don't want to retract peripheral loops
                continue
            e_no = i+1
            left = self.edge_puncs[i][0]
            right = self.edge_puncs[i][1]
            #print 'left, right', left, right
            if left != right and (left < 0 or right < 0): # different regions to its left and right
                # which side should we retract onto? the one without a puncture
                if left < 0: # negative value for region mean it's not a puncture
                    # retract left
                    sur_edges_ccw = self.surEdges(left)
                    for j in range(len(sur_edges_ccw)):
                        if sur_edges_ccw[0] == e_no:
                            sur_edges_ccw.pop(0)
                            break
                        else:
                            sur_edges_ccw.append(sur_edges_ccw.pop(0))
                    sur_edges_cw = sur_edges_ccw[:]
                    sur_edges_cw.reverse()
                    sur_edges_cw = [-x for x in sur_edges_cw]
                    
                    # replace all occurences of e_no with sur_edges_cw, and -e_no by sur_edges_ccw
                    # in the image of edges
                    for j in range(1,1+self.n_edges):
                        if j == e_no:
                            continue
                        image = self.getEdgeImage(j)
                        for k in reversed(range(len(image))):
                            if image[k] == e_no:
                                image = image[:k] + sur_edges_cw + image[k+1:]
                            elif image[k] == -e_no:
                                image = image[:k] + sur_edges_ccw + image[k+1:]
                        self.setEdgeImage(j, image)
                        
                    #print 'sur_edges_ccw', sur_edges_ccw
                    
                    # update regions to the left/right of edges
                    for e in sur_edges_ccw:
                        self.setLeftEdgePunc(e, right)
                    
                    #print 'about to remove edge', e_no
                    #print '\n', self
                    # delete the edge now
                    self.removeEdge(e_no)
                    #print '\nafter removal\n', self
                    
                    self.retract()
                    return
                elif right < 0:
                    # retract right
                    sur_edges_ccw = self.surEdges(right)
                    
                    for j in range(len(sur_edges_ccw)):
                        if sur_edges_ccw[0] == -e_no:
                            sur_edges_ccw.pop(0)
                            break
                        else:
                            sur_edges_ccw.append(sur_edges_ccw.pop(0))
                            
                    sur_edges_cw = sur_edges_ccw[:]
                    sur_edges_cw.reverse()
                    sur_edges_cw = [-x for x in sur_edges_cw]
                    
                    # replace all occurences of e_no with sur_edges, and -e_no by sur_edges_rev
                    # in the image of edges
                    for j in range(1,1+self.n_edges):
                        if j == e_no:
                            continue
                        image = self.getEdgeImage(j)
                        for k in reversed(range(len(image))):
                            if image[k] == e_no:
                                image = image[:k] + sur_edges_ccw + image[k+1:]
                            elif image[k] == -e_no:
                                image = image[:k] + sur_edges_cw + image[k+1:]
                        self.setEdgeImage(j, image)
                    
                    # update regions to the left/right of edges
                    for e in sur_edges_cw:
                        self.setRightEdgePunc(e, left)
                        
                    #print 'about to remove edge', e_no
                    #print '\n', self
                    # delete the edge now
                    self.removeEdge(e_no)
                    #print '\nafter removal\n', self
                    
                    self.retract()
                    return
                    
                    
    def __str__(self):
        return ("verts: " + str(self.ribbon) + "\nvert map: " + str(self.vert_map)
                        + "\nedge puncs: " + str(self.edge_puncs) +
                        "\nedge_paths: " + str(self.edge_map) + 
                        "\npunc_perm: " + str(self.punc_perm) + "\nperiph: " + str(self.periph))
        
    def pullTight(self, path=None):
        """ returns a tightened copy of a path. """
        path = path[:]
        i = 0
        while i < len(path)-1 and i >= 0:
            if path[i] == -path[i+1]:
                del path[i:i+2]
                i -= 1
                i = max(0, i)
            else:
                i += 1
        return path
        
    def tighten(self):
        new_edge_map = []
        for e_path in self.edge_map:
            new_edge_map.append(self.pullTight(e_path))
        self.edge_map = new_edge_map    
    
        
    def saveAs(self, foutname):
        """
            The Trains graph map file format is undocumented. I worked out (guessed?) 
            the format based on examples and email correspondence with Toby Hall.
            
            At the moment this has only been tested to work for once punctured surfaces,
            i.e. no peripheral loops.
        """
        f = open(foutname, 'w')
        f.write('V 4.3\n')
        # i think the 20 20 below allocates memory for something, may need to increase
        f.write('{0} {1} {2} 0 {3} {4}\n'.
                    format(self.n_periph, self.n_edges, self.n_verts, self.n_verts+1, self.n_edges+1))
        for i in range(self.n_verts):
            f.write('{0} {1}\n'.format(i+1, self.vert_map[i]))
            for j in self.ribbon[i]:
                f.write('{0}\n'.format(j))
            f.write('0\n')
        for i in range(self.n_edges):
            f.write('{0} {1} {2}\n'.format(self.getInitVert(i+1), self.getTermVert(i+1), i+1))
            if self.periph[i] == 0:
                f.write('0\n')
            else:
                f.write('1 {0}\n'.format(self.periph[i]))
            e_image = self.edge_map[i]
            for j in range(len(e_image)):
                f.write('{0} '.format(e_image[j]))
            f.write('0\n')
        f.write('5\n0\n0\n') # this might need to be 0 0 0 for once punctured g=2 surface.
                             # this might need to be 5 0 0 for once punctured g=3 surface.
                             
    """
        Runs a background trains process and sends the graphmap data to trains. Obtains
        invariant train track information as string, which is returned.
        
        Editing may need to be done to handle homeomorphisms with more than one puncture, i.e.
        peripheral curve information needs to be sent to trains but this is not currently
        handled in sendToTrains().
    """
    def sendToTrains(self, debug=False):
        cur_path = os.path.dirname(os.path.abspath(__file__))
        platform = sys.platform.lower()
        filename = "trains.exe"
        if platform in ['win32', 'cygwin']:
            filename = "trains.exe"
        elif platform in ['linux2']:
            filename = "trains.out"
        elif platform in ['darwin']:
            filename = "trains"
            
        if debug:
            stdin = sys.stdout
        else:
            trains_proc = subprocess.Popen(cur_path + "/" + filename, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=None)
            stdin = trains_proc.stdin
        
        stdin.write("input\n")
        stdin.write("{0} {1} {2}\n".format(self.n_periph, self.n_edges, self.n_verts))
        
        for i in range(self.n_verts):
            stdin.write("{0}\n".format(self.vert_map[i]))
            stdin.write(str(self.ribbon[i])[1:-1].replace(",", "") + " 0\n") # cyclic order of edges around vertex i
            
        punc = 1
        for i in range(self.n_edges):
            if self.n_periph > 0:
                if self.periph[i] is not None:
                    stdin.write("1\n") # yes peripheral edge
                    stdin.write(str(punc) + "\n") # peripheral about puncture i+1
                    punc += 1
                else:
                    stdin.write("0\n") # no not peripheral edge
            stdin.write(str(self.edge_map[i])[1:-1].replace(",", "") + " 0\n")
            
        # add loops to keep track of punctures
        dic = self.surLoops()
        for i in sorted(dic.keys()):
            stdin.write('addloop\n')
            stdin.write(str(dic[i])[1:-1].replace(",", "") + " 0\n")
            stdin.write("Punc " + str(i-1) + "\n") # Trains name of loop, first punc is 0
            
        stdin.write("train\n") # ask Trains to compute invariant train track
        stdin.write("print\n") # print train track
        stdin.write("loops\n")
        stdin.write("exit\n")
        
        train_info = []
        start_saving = False
        for line in trains_proc.stdout:
            #if "Now have an efficient" in line:
            if "Isotopy class is Pseudo-Anosov" in line:
                start_saving = True
            if start_saving:
                train_info.append(line.replace(">","").strip())
                #print line
            #sys.stdout.flush()
        train_info.insert(0, "Non-braid")
        return (train_info, start_saving) # start_saving is True if map is pA

    def saveTriangulation(self, filename, lines=None, tri_title=None):
        # if lines is not None then it's expected to be a list of lines of trains output
        if not lines:
            lines, is_pA = self.sendToTrains()
        else:
            is_pA = True
        output = ()
        if is_pA:
            from triangulation_builder import trains_to_snappea_tri
            output = trains_to_snappea_tri(fname=None, foutname=filename, lines=lines, debug=False, tri_title=tri_title)
        return is_pA, output
        
    def splitting_info(self):
        lines, is_pA = self.sendToTrains()
        if is_pA:
            from train_track_parser import parse_train_track, splitting_sequence
            tt, growth_rate, punc_bij = parse_train_track(lines=lines)
            ss = splitting_sequence(tt, growth_rate, punc_bij)
            return (ss, growth_rate, punc_bij)
        else:
            return False

def torusTwoPuncs(twists):
    """Convention: If twists = [1,2,3] then twist in 3 then 2 then 1."""
    gm = GeneralGraphMap(n_periph=2, n_edges=6, n_verts=2, ribbon=[[3,2,-3,-1,5,-5], [-2,6,-6,4,1,-4]],
                         edge_map=None, vert_map=None, periph=[0, 0, 0, 0, 1, 2], 
                         edge_puncs=[[-1,-1], [-2,-2], [-2,-1], [-1,-2], [1,-1], [2,-2]],
                         n_puncs=2, punc_perm=None)
                         
    def dehnTwist(gm, curve):
        a_1 = gm.getEdgeImage(1)
        a_2 = gm.getEdgeImage(2)
        b_1 = gm.getEdgeImage(3)
        b_2 = gm.getEdgeImage(4)
        p_1 = gm.getEdgeImage(5)
        p_2 = gm.getEdgeImage(6)
        
        A_1 = gm.getEdgeImage(-1)
        A_2 = gm.getEdgeImage(-2)
        B_1 = gm.getEdgeImage(-3)
        B_2 = gm.getEdgeImage(-4)
        P_1 = gm.getEdgeImage(-5)
        P_2 = gm.getEdgeImage(-6)
        
        if curve == 1: # left twist in loop (edge_1 edge_2)
            gm.setEdgeImage(3, b_1 + a_2 + a_1)
            gm.setEdgeImage(4, b_2 + a_1 + a_2)
        elif curve == 2: # left twist in loop edge_3
            gm.setEdgeImage(2, B_1 + a_2)
        elif curve == 3: # left twist in loop edge_4
            gm.setEdgeImage(1, B_2 + a_1)
        elif curve == -1:
            gm.setEdgeImage(3, b_1 + A_1 + A_2)
            gm.setEdgeImage(4, b_2 + A_2 + A_1)
        elif curve == -2:
            gm.setEdgeImage(2, b_1 + a_2)
        elif curve == -3:
            gm.setEdgeImage(1, b_2 + a_1)
            
    for edge in twists:
        dehnTwist(gm, edge)
        gm.tighten()
    return gm
    
def surfaceGM(g, p, twists):
    """ g is the genus of the surface.
        p is the number of punctures
        twists is a string of twists
        
        p1, p2, ... pp are right half-twists in punctures (1,2), (2,3), etc
        if twists = [1,2,3] then twist in 1 followed by twist in 2 followed by twist in 3 etc.
    """
    if g == 0:
        n_edges = 3*p
    elif g == 1:
        n_edges = 3*p+2
    else:
        n_edges=3*p+9+(g-2)*8
        
    n_verts = max(2*g+p, 1+p) # if g == 0 then there's 1 vert
    ribbon = []
    # there are 2*g vertices [with (g-1) of them behind the surface]
    for i in range(n_verts):
        ribbon.append([])
    
    # first vertex with surrounding punctures
    for i in range(p):
        ribbon[0].extend([n_edges-p+i+1,i+1+p])
        
    if g > 0:
        ribbon[0].append(2*p+1)
    for i in range(p): # the other side of rightmost puncture
        ribbon[0].append(-(2*p-i))
    if g > 0:
        ribbon[0].append(-(2*p+2))
    
    if g > 0:
        ribbon[1].extend([2*p+2,-(2*p+1)])
    if g>1:
        ribbon[1].insert(1,2*p+3)
        ribbon[1].extend([2*p+5, 2*p+7, -(2*p+4), -(2*p+8), -(2*p+6)])
        
    # the middle front vertices
    for i in range(g-2):
        ribbon[i+2] = [2*p+11+8*i, 2*p+9+8*i, -(2*p+7+8*i), -(2*p+9+8*i),
                       2*p+13+8*i, 2*p+15+8*i, -(2*p+12+8*i), -(2*p+16+8*i), 
                       -(2*p+14+8*i), -(2*p+10+8*i), 2*p+8+8*i, 2*p+10+8*i]
                       
    # leftmost vertex  
    if g>1:
        i = 2*p+7+(g-2)*8
        ribbon[g] = [-i, -(i+2), i+1, i+2]
        
    # back vertices
    for i in range(g-1):
        j = 2*p+3+8*i
        ribbon[g+1+i] = [-j, j+3, j+1, -(j+2)]
        
    # peripheral vertices
    for i in range(p):
        ribbon[n_verts-p+i] = [(i+1), -(i+1), -(n_edges-p+i+1)]
    
    # set peripheral edges, these are the edges which loop around a puncture.
    
    periph = [None for i in range(n_edges)]
    for i in range(p):
        periph[i] = i+1
        
    edge_puncs = [[None,None] for i in range(n_edges)]
    for i in range(p):
        edge_puncs[i] = [i+1, None]
    
    gm = GeneralGraphMap(n_periph=0, n_edges=n_edges, n_verts=n_verts,
                        ribbon=ribbon, edge_map=None, vert_map=None,
                        periph=periph, edge_puncs=edge_puncs, n_puncs=p, punc_perm=None)
                        
    gm.inferEdgePuncs(True)
    
    def dehnTwist(w):
        if len(w) != 2:
            raise Exception('Invalid input: ' + str(w))
        i = int(w[1:])
        twist_loop = []
        if w[0].lower() == "a":
            twist_loop = [p+i]
        elif w[0].lower() == "b":
            if i == 1:
                twist_loop = [2*p+1, 2*p+2]
            else:
                twist_loop = [2*p-1+8*(i-1), 2*p+8*(i-1)]
        elif w[0].lower() == "c":
            if g == 1: # in this case c1 is the same as a ap twist
                twist_loop = [2*p]
            elif i == g:
                twist_loop = [2*p+3+8*(g-2)+6]
            else:
                twist_loop = [2*p+3+8*(i-1), 2*p+4+8*(i-1)]
        elif w[0].lower() == "d":
            if i % 2 == 1:
                twist_loop = [2*p+9+4*(i-1)]
            else:
                twist_loop = [2*p+2+4*i]
        elif w[0].lower() == "e":
            twist_loop = [2*p+5+8*(i-1), 2*p+6+8*(i-1)]
        elif w[0].lower() == "p":
            a = gm.getEdgeImage(i+p)
            b = gm.getEdgeImage(i)
            c = gm.getEdgeImage(i+1)
            d = gm.getEdgeImage(n_edges-p+i)
            e = gm.getEdgeImage(n_edges-p+i+1)
            
            A = gm.getEdgeImage(-(i+p))
            B = gm.getEdgeImage(-i)
            C = gm.getEdgeImage(-(i+1))
            D = gm.getEdgeImage(-(n_edges-p+i))
            E = gm.getEdgeImage(-(n_edges-p+i+1))

            if w[0] == "p": # clockwise half-twist
                gm.setEdgeImage(i+p, e+C+E+d+b+D+a) # a
                gm.setEdgeImage(i, c) # b
                gm.setEdgeImage(i+1, b) # c
                gm.setEdgeImage(n_edges-p+i, e) # d
                gm.setEdgeImage(n_edges-p+i+1, e+C+E+d) # e
            else: # ccw half-twist
                gm.setEdgeImage(i+p, d+b+D+e+C+E+a) # a
                gm.setEdgeImage(i, c) # b
                gm.setEdgeImage(i+1, b) # c
                gm.setEdgeImage(n_edges-p+i, d+b+D+e) # d
                gm.setEdgeImage(n_edges-p+i+1, d) # e
                
            gm.punc_perm[i-1], gm.punc_perm[i] = gm.punc_perm[i], gm.punc_perm[i-1]
            
            # Update vertex images.
            # Permute the image of these vertices.
            gm.vert_map[n_verts-p-1 + i], gm.vert_map[n_verts-p + i] = \
              gm.vert_map[n_verts-p + i], gm.vert_map[n_verts-p-1 + i]
        else:
            raise Exception("Invalid input: " + str(w))
        if w[0].lower() != "p":
            if w[0] != w[0].lower(): # if not lower case
                gm.dehnTwist(twist_loop, right=True)
            else:
                gm.dehnTwist(twist_loop, right=False)
        
    for w in reversed(twists.split(" ")):
        sys.stdout.flush()
        dehnTwist(w)
        gm.tighten()
    #print 'finished twisting'
    sys.stdout.flush()
        
    #print gm
    gm.inferPeripheral()
    #print '\nbefore retract:\n', gm
    #print 'size:',sum([len(x) for x in gm.edge_map])
    
    
    gm.retract()
    
    
    
    #print '\nafter retract:\n', gm
    #print gm
    #gm.inferEdgePuncs(False)
    #print '\nafter inferEdgePuncs:', gm
    #gm.inferPuncPerm()
    #gm.removeEmptyVerts()
    
    gm.tighten()
    
    
    #print 'size:',sum([len(x) for x in gm.edge_map])
    #print '\nafter inferPuncPerm:', gm
    #gm.removeEmptyVerts()
    #gm.inferVertexMap()
    #print '\n',gm
    return gm

def generators(g, p):
    # returns a list of all Dehn twists and half-twists on a
    # genus g surface with p punctures.
    gens = []
    
    for i in range(1, g+1):
        gens.append("b" + str(i))
        gens.append("c" + str(i))
        
    for i in range(1, g):
        gens.append("e" + str(i))
        
    for i in range(1, g-1):
        gens.append("d" + str(2*i-1))
        gens.append("d" + str(2*i))
    
    for i in range(1, p+1):
        gens.append("a" + str(i))
        
    for i in range(1, p):
        gens.append("p" + str(i))
        
    for w in gens[:]:
        gens.append(w.upper())
        
    return gens
    
def toTwister(twists):
    # this is for S(2,1) and S(3,1) surfaces
    # 1) reverse order of twists since we read from right to left, twister reads from left to right
    # Conventions: C1 -> a_0, C2 -> b_1, C3 -> b_2, C4 -> b_3, C5 -> c
    # Twister convention: Figure 13 of "Lebruere and Paris -- Presentations for the punctured
    # mapping class groups."
    t = twists[:]
    mapping = {-1:'a_0', -2:'b_1', -3:'b_2', -4:'b_3', -5:'c',-6:'b_4',-7:'b_5',
                1:'A_0',2:'B_1',3:'B_2',4:'B_3',5:'C',6:'B_4',7:'B_5'}

    t.reverse()
    return map(lambda x: mapping[x], t)
    
