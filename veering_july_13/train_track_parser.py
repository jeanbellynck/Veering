

from train_track import *

import re, math, copy

from collections import deque, defaultdict

import perron_frobenius

import os

from general_graph_map import GeneralGraphMap



# need to set the precision here

#getcontext().prec = 100

    

def parse_train_track(fname=None, lines=None, debug=False, weight_dec=40):

    if lines is None:

        f = open(fname, 'r')

        lines = f.readlines()



    pv1 = re.compile('Vertex number (\d+) with image vertex (\d+):')

    pv2 = re.compile('Edges at vertex are: (.*)')

    

    pe1 = re.compile('Edge number (\d+) from vertex (\d+) to vertex (\d+):')

    pe2 = re.compile('Type: Peripheral about puncture number (\d+)')

    pe3 = re.compile('Image is: (.*)')

    pe4 = re.compile('Path \(\d+ -> \d+\):(.*)')

    

    pg1 = re.compile('Vertex (\d+):')

    

    pp1 = re.compile('Punc (\d+): (.*)')

    

    punc_pat = re.compile('Graph on surface with (\d+) peripheral loops:')

    

    

    # edge format:

    # [edge number, [vert1, vert2], [left puncture, right puncture], [edge images], [isotopy path]]

    # this changes later

    

    edges = []

    edge = []

    

    # we've decided to keep the image paths

    image_paths = {}

    

    # we store the edge numbers of infinitesimal edges

    # this is needed to compute the transition matrix

    inf_edge_nos = []

    

    # keeps track of which vertices the extra punctures correspond to

    # we need this to work out where a puncture is sent

    punc_to_vert = {}

    

    # are we dealing with a braid?

    is_braid = True

    

    """

    when we add infinitesimal edges we relabel the vertices

    we need to keep track of the old vertices to obtain the new image paths

    """

    new_to_old_verts = {} # this might not be needed, take a closer look

                

    # vert format:

    # [vertex number, surround edge in cyclic (counter-clockwise) order, region vertex is in]

    vert_bij = {}

    verts = []

    vert = []

    

    n_punc = 0

    punc_bij = {}

    punc_bij[0] = 0 # might need to remove this in the future, ASSUMES 1 puncture in g>0 surface.

    braid_pres = []

    

    growth_rate = 0

    perron_eigenvector = []

    

    cur_vert = 0 # keep track of vertex number for infinitesimal edges

    gates = {}

    inf_edges = {}



    for line in lines:

        line = line.strip()

        if line.startswith('Vertex number'):

            vert_ends = [int(s) for s in pv1.match(line).groups()]

            vert_bij[vert_ends[0]] = vert_ends[1]

            vert = [vert_ends[0]]

        elif line.startswith('Non-braid') or line.startswith('surface homeomorphism'):

            is_braid = False

        elif line.startswith('Edges at vertex'):

            vert.append([int(s) for s in pv2.match(line).group(1).split()])

            if not is_braid:

                vert.append(0) # fake 'Region 0'

                verts.append(vert)

        elif line.startswith('Region'):

            vert.append(int(line[7:]))

            verts.append(vert)

            

        # edges information

        elif line.startswith('Edge number'):

            grps = pe1.match(line).groups()

            edge = []

            edge.append(int(grps[0])) # edge number

            edge.append([int(grps[1]), int(grps[2])]) # terminal vertices

        elif line.startswith('Type: '):

            m = pe2.match(line)

            if m != None: # got a match, so edge is peripheral

                puncture = int(m.group(1))

                edge.append([puncture, 0]) # puncture on left

            else: # non-peripheral, surrounding puncture is 0

                edge.append([0, 0])

        elif line.startswith('Image is'):

            edge.append([int(s) for s in pe3.match(line).group(1).split()])

            # fake 'Path (1 -> 1): 1 -1...'

            if not is_braid:

                edge.append([]) # empty isotopy info, not keeping track

                edges.append(edge)

        elif line.startswith('Path'):

            # add isotopy information

            edge.append([int(s) for s in pe4.match(line).group(1).split()])

            edges.append(edge)

            

            

            

     ### If we're working on other surfaces we'll need to set braid_pres in some other way

        elif line.startswith('Braid: '):

            braid_pres = [int(s) for s in line[7:].strip().split()]

            

        # add permutation of punctures information as though it were a braid

        elif line.startswith('Puncture permutation: '):

            is_braid = False

            perm = [int(s) for s in line[22:].strip().split()]

            for j in range(len(perm)):

                punc_bij[j] = perm[j]

        elif line.startswith('Graph on surface with'):

            n_punc = int(punc_pat.match(line).group(1)) # this is true for braids, not in general

    

        # infinitesimal edge information

        elif line.startswith('Vertex'): # note, must be done after 'Vertex number..', for infinitesimal edge info

            cur_vert = int(line.strip()[7:-1])

        elif line.startswith('Gates are:'):

            partition = line[11:].split('}, {')

            partition[0] = partition[0][1:]

            partition[-1] = partition[-1][:-1]

            # partition will look like  ['{-2}', '{9, -8}', '{2}']

            # now to remove the braces

            partition = [map(int, s.split(', ')) for s in partition]

            gates[cur_vert] = partition

        elif line.startswith('Infinitesimal edges join '):

            # eg Infinitesimal edges join 1 to 4, -1 to 4

            inf_edges[cur_vert] = [map(int, w.split(' to ')) for w in line[25:].split(', ')]

        elif line.startswith("Punc"):

            punc, loop = pp1.match(line).groups()

            punc = int(punc)

            loop = [int(x) for x in loop.split()]

            for i in loop:

                if i > 0:

                    edges[abs(i)-1][2][0] = punc

                else:

                    edges[abs(i)-1][2][1] = punc

            

        

    # want to reorder elements of each gate so that we get the cyclic order

    def order_gate():

        for i in range(0, len(verts)):

            vert_cyc = verts[i][1]

            def cyc_key(g):

                return min([vert_cyc.index(w) for w in g])

            gates[verts[i][0]].sort(key=cyc_key)

                

    order_gate()

    # WARNING: assuming that within each gate the edges are oriented correctly

    # which I assume is a Trains convention.

    

    # determine missing infinitesimal edge

    # returns None if only 2 gates (therefore 1 infinitesimal edge)

    # rotates gate so that the last element doesn't link to the first

    def missing_edge(gate, edge_pairs):

        if len(gate) == 2:

            return None

        

        def connects(g1, g2):

            for p in edge_pairs:

                if g1 in p and g2 in p:

                    return True

            return False

            

        for i in range(len(gate)):

            if not connects(gate[0][0], gate[1][0]):

                gate.append(gate.pop(0))

                break

            gate.append(gate.pop(0))

                

    for i in range(1, 1+len(gates)):

        missing_edge(gates[i], inf_edges[i])

    #now missing edge os 0:-1

    

    # compute permutation of punctures

    if is_braid:

        for i in range(0,n_punc+1):

            punc_bij[i] = i

        

        #print 'braid presentation:', braid_pres

        puncs = range(1,n_punc+1)

        for twist in braid_pres:

            pos_twist = abs(twist)

            puncs[pos_twist-1], puncs[pos_twist] = puncs[pos_twist], puncs[pos_twist-1]

        for i in range(1, n_punc+1):

            punc_bij[puncs[i-1]] = i

    # if not braid then punc_bij needs to be set separately.



    # now to create transition matrix and image graph, should have weights

    # create transition matrix

    mat = [[0]*len(edges) for i in range(len(edges))]

    for edge in edges:

        image_edges = edge[3]

        for i in image_edges:

            abs_i = 0

            if (i >= 0):

                abs_i = i

            else:

                abs_i = -i

            mat[abs_i-1][edge[0]-1] += 1

        

    # edge format:

    # [edge number, [vert1, vert2], [left puncture, right puncture], [edge images], [isotopy path]]

    # vert format:

    # [vertex number, surround edge in cyclic (counter-clockwise) order, region vertex is in]

    

    verts_2 = copy.deepcopy(verts)

    for i in range(len(verts)):

        for key in vert_bij.keys():

            if vert_bij[key] == verts_2[i][0]:

                verts_2[i][0] = key

                break

    

    edges_2 = []

    for edge in edges:

        new_edge = []

        new_edge.append(edge[0]) # edge number

        # use bijection to determine vertex at ends of edge

        

        new_edge.append(edge[1][:])

        # use permutation of punctures to determine left/right punctures

        if is_braid:

            new_edge.append( [ punc_bij[edge[2][0]], punc_bij[edge[2][1]] ] )

        else:

            new_edge.append(edge[2])

        isotopy = []

        edge_images = edge[3]

        for e in edge_images:

            abs_e = abs(e)

            for ed in edges:

                if ed[0] == abs_e:

                    ed_c = ed[4][:]

                    if e < 0:

                        ed_c.reverse()

                        ed_c = [-x for x in ed_c]

                    isotopy.extend(ed_c)

                    break

        new_edge.append(isotopy)

        # add weights

        new_edge.append(Decimal(0))

        edge.append(Decimal(0))

        

        edges_2.append(new_edge)

    

        

    for edge in edges:

        # this is clumsy but we've decided to keep this information

        # we need it to compute the transition matrix with infinitesimal edges

        # which is needed to work out the weights of infinitesimal edges

        # when the infinitesimal edges form a polygon at a vertex

        image_paths[edge[0]] = edge[3]

        del(edge[3]) # delete image information

        

    # NEW edge format:

    # [image edges] has been removed, no longer needed

    # [edge number, [vert1, vert2], [left puncture, right puncture], [isotopy path], weight]

    

    # vert format:

    # [vertex number, surround edge in cyclic (counter-clockwise) order, region vertex is in]



    if not is_braid:

        verts.sort(key=lambda vert: vert[0])

        edges.sort(key=lambda e: e[0])

        gm = GeneralGraphMap(n_edges=len(edges), n_verts=len(verts), ribbon=[v[1] for v in verts],

                        edge_map=[image_paths[edges[i][0]] for i in range(len(edges))], 

                        vert_map=vert_bij.values(), edge_puncs=[e[2] for e in edges], n_puncs=0)

                        

        gm.inferEdgePuncs(respect_original_puncs=True)

        gm.inferPuncPerm()

        n_punc = gm.n_puncs

        for i in range(len(edges)):

            edges[i][2] = gm.edge_puncs[i]

        #print 'puncture perm from gm:', gm.getPuncturePerm()

        for i in range(n_punc):

            punc_bij[i] = gm.getPuncturePerm()[i]

        #print gm

    

    """

    At this point we're now ready to create the train track. The train track will have

    switches of valence greater than 3 in general. We'll need to fix that later.

    """

    # gates

    # inf_edges

    

    def get_edge(edges, n): # returns edge numbered...

        for edge in edges:

            if edge[0] == n:

                return edge

        raise Exception('error in get_edge: '+str(n))

        return None # should error here

    

    def weight(n):

        if n < 0:

            return weight(-n)

        for edge in edges:

            if edge[0] == n:

                return edge[4] # weight

        return 0

        

    def left_puncture(edges, gate):

        left_punc = 0

        if gate[0] >= 0:

            left_punc = get_edge(edges,gate[0])[2][0] # left of edge

        else:

            left_punc = get_edge(edges,-gate[0])[2][1] # right of edge

        return left_punc

            

    def right_puncture(edges, gate):

        right_punc = 0

        if gate[-1] >= 0:

            right_punc = get_edge(edges,gate[-1])[2][1] # right of edge

        else:

            right_punc = get_edge(edges,-gate[-1])[2][0] # left of edge

        return right_punc

        

    # for the edge numbered e_no (with sign) set its terminal vertex to be v_no

    def set_head_vertex(edges, e_no, v_no):

        e = get_edge(edges, abs(e_no))

        if e_no < 0:

            e[1][0] = v_no

        else:

            e[1][1] = v_no

    

    tt = TrainTrack()

    

    """Code originally for braids, in which case n_punc is the number of punctures in the disk.

        This is one less than the actual number of punctures in the surface. For general surfaces

        we'll decrease n_punc, execute the code below, then increment it back to get the actual

        number of punctures"""

    if not is_braid:

        n_punc -= 1

        

    tt.no_of_puncs = n_punc # this will actually need to be changed

    tt.strands = n_punc # this remains fixed

    cur_vert = 1 # keep track of how many vertices we've added to our train track

    cur_edge = len(edges)+1

    

    # add current

    """

    step 1

    work out if any of the two ending vertices on each infinitesimal circle

    are not really new vertices

    """

    for i in range(len(verts)):

        v = verts[i]

        gate = gates[v[0]]

        inf_edge = inf_edges[v[0]] # no longer need this, all information contained in gates



        left_punc = left_puncture(edges, gate[0])

        

        if len(gate) == 2: # only two gates

            edges.append([cur_edge, [cur_vert, cur_vert+1], 

                            [left_puncture(edges,gate[0]),right_puncture(edges,gate[0])], [], sum(map(weight, gate[0]))])

            inf_edge_nos.append(cur_edge)

                            

            tt.switches[cur_vert] = Switch(n=cur_vert, region=v[2], branches=[gate[0],[cur_edge]])

            new_to_old_verts[cur_vert] = v[0]

            # each of the edges connected to this gate needs to have its terminal vertex updated

            for e_no in gate[0]:

                set_head_vertex(edges, -e_no, cur_vert)

            

            cur_vert += 1

            tt.switches[cur_vert] =  Switch(n=cur_vert, region=v[2], branches=[gate[1],[-cur_edge]])

            new_to_old_verts[cur_vert] = v[0]

                

            # each of the edges connected to this gate needs to have its terminal vertex updated

            for e_no in gate[1]:

                set_head_vertex(edges, -e_no, cur_vert)



            cur_vert += 1

            cur_edge += 1

        elif len(gate) == len(inf_edge): # infinitesimal edges form polygon and bound singularity

            # need to update punc_bij and tt.no_of_puncs

            n_punc += 1

            no_of_inf_edges = len(inf_edge)



            """

            the new punctures correspond to vertices, to see where they get mapped to:

            new puncture --> vertex ---> image vertex ---> corresponding puncture

            """

            punc_to_vert[n_punc] = v[0]

            # need to do something special to compute infinitesimal weights in this case.

            

            left_punc = n_punc # we know left punc of each inf. edge is the new puncture

            left_edge = (cur_edge + no_of_inf_edges - 1)



            for k in range(len(gate)):

                g = gate[k]

                if k == len(gate) - 1: # last gate

                    # the terminal vertex is the same as the initial vertex of the first gate

                    edges.append([cur_edge, [cur_vert, cur_vert+1-no_of_inf_edges], 

                                [left_punc,right_puncture(edges,g)], [], 0])

                else:

                    edges.append([cur_edge, [cur_vert, cur_vert+1], 

                            [left_punc,right_puncture(edges,g)], [], 0])

                          

                inf_edge_nos.append(cur_edge)



                tt.switches[cur_vert] = Switch(n=cur_vert, region=v[2], branches=[g, [cur_edge, -left_edge]])

                new_to_old_verts[cur_vert] = v[0]

                

                for e_no in g:

                    set_head_vertex(edges, -e_no, cur_vert)

                    

                cur_vert += 1

                left_edge = cur_edge

                cur_edge += 1

            

        else: # this is the only other possibility, polygon missing an edge

            # have to check first and last vertex on inf. circle

            left_edge = -gate[0][0]

            edges.append([cur_edge, [cur_vert, cur_vert+1], 

                        [left_puncture(edges,[left_edge]),right_puncture(edges,[left_edge])], [], sum(map(weight, gate[0]))])

            inf_edge_nos.append(cur_edge)

                

            tt.switches[cur_vert] = Switch(n=cur_vert, region=v[2], branches=[gate[0],[cur_edge]])

            new_to_old_verts[cur_vert] = v[0]

                

            for e_no in gate[0]:

                set_head_vertex(edges, -e_no, cur_vert)

                    

            cur_vert += 1

            left_edge = cur_edge

            cur_edge += 1

                

            for g in gate[1:-2]:

                edges.append([cur_edge, [cur_vert, cur_vert+1], 

                        [left_puncture(edges,[left_edge]),right_puncture(edges,g)], [], sum(map(weight, g)) - get_edge(edges, left_edge)[4]])

                inf_edge_nos.append(cur_edge)

                

                tt.switches[cur_vert] = Switch(n=cur_vert, region=v[2], branches=[g, [cur_edge, -left_edge]])

                new_to_old_verts[cur_vert] = v[0]

                

                # each of the edges connected to this gate needs to have its terminal vertex updated

                for e_no in g:

                    set_head_vertex(edges, -e_no, cur_vert)

                    

                cur_vert += 1

                left_edge = cur_edge

                cur_edge += 1



            edges.append([cur_edge, [cur_vert, cur_vert+1], 

                    [left_puncture(edges,[left_edge]),right_puncture(edges,gate[-2])], [], sum(map(weight, gate[-1]))])

            inf_edge_nos.append(cur_edge)

                        

            tt.switches[cur_vert] = Switch(n=cur_vert, region=v[2], branches=[gate[-2], [cur_edge, -left_edge]])

            new_to_old_verts[cur_vert] = v[0]

            for e_no in gate[-2]:

                set_head_vertex(edges, -e_no, cur_vert)

                    

            cur_vert += 1

            tt.switches[cur_vert] = Switch(n=cur_vert, region=v[2], branches=[gate[-1], [-cur_edge]])

            new_to_old_verts[cur_vert] = v[0]

            for e_no in gate[-1]:

                set_head_vertex(edges, -e_no, cur_vert)

            cur_vert += 1

            left_edge = cur_edge

            cur_edge += 1

                

    # need to add branches from 'edges' to tt

    

    if not is_braid:

        n_punc += 1

        

    # due to extra singularities we may have increased number of punctures

    tt.no_of_puncs = n_punc

    

    def vert_to_punc(v_no):

        for key in punc_to_vert.keys():

            if punc_to_vert[key] == v_no:

                return key

    

    if is_braid:

        for i in range(1, n_punc+1):

            if i not in punc_bij.keys():

                punc_bij[i] = vert_to_punc(vert_bij[punc_to_vert[i]])

    else:

        for i in range(n_punc):

            if i not in punc_bij.keys():

                punc_bij[i] = vert_to_punc(vert_bij[punc_to_vert[i]])

 

    for e in edges:

        edge_no = e[0]

        term_verts = e[1]

        LR_puncs = e[2]

        #isotopy = e[3]

        #isotopy = defaultdict(int)

        #isotopy[edge_no] = 1

        #isotopy = [edge_no]

        isotopy = []

        weight = e[4]

        tt.branches[edge_no] = Branch(edge_no, LR_puncs, isotopy, term_verts, weight)

    

    if debug:

        print 'punc bij:', punc_bij



    # compute the transition matrix with infinitesimal edges here.



    # need a helper function. takes 2 vertices and returns inbetween infinitesimal edge

    # returns unsigned edge number

    def bw_inf_edge(v1_no, v2_no):

        for e_no in inf_edge_nos:

            b = tt.branches[abs(e_no)].b_ends

            if v1_no == b[0] and v2_no == b[1]:

                return e_no

            elif v1_no == b[1] and v2_no == b[0]:

                return -e_no

    

    # returns the vertex at the end of edge e_no (with sign)

    def term_vert(e_no):

        b = tt.branches[abs(e_no)].b_ends

        if e_no > 0:

            return b[1]

        else:

            return b[0]

            

    # returns the image vertex of v_no

    def image_vert(v_no):

        for b in tt.switches[v_no].branches:

            if abs(b[0]) not in inf_edge_nos:

                if b[0] >= 0:

                    return term_vert(-image_paths[b[0]][0])

                else:

                    return term_vert(image_paths[-b[0]][-1])

                            

    # determines the edge number of the infinitesimal edge which

    # the edge numbered e_no maps to 

    def inf_to_inf_edge(e_no):        

        v1_no, v2_no = tt.branches[abs(e_no)].b_ends

        im_v1_no = image_vert(v1_no)

        im_v2_no = image_vert(v2_no)

        return bw_inf_edge(im_v1_no, im_v2_no)

            

    

    # the lists in image_paths

    trans_mat = [[0]*len(edges) for i in range(len(edges))]

    for edge in edges:

        e_no = edge[0]

        if e_no not in inf_edge_nos:

            for k in range(len(image_paths[abs(e_no)])):

                i = image_paths[abs(e_no)][k]

                trans_mat[abs(i)-1][abs(e_no)-1] += 1

                tt.branches[abs(e_no)].isotopy.append(i)

                if k < len(image_paths[abs(e_no)])-1: # check for inbetween infinitesimal edges

                    v1_no = term_vert(i)

                    v2_no = term_vert(-image_paths[abs(e_no)][k+1])

                    inf_no = bw_inf_edge(v1_no, v2_no)

                    if inf_no != None:

                        trans_mat[abs(inf_no)-1][abs(e_no)-1] += 1

                        tt.branches[abs(e_no)].isotopy.append(inf_no)

                    elif v1_no != v2_no:

                        raise Exception('An error occurred in computing extended transition matrix.')

        else:

            trans_mat[abs(inf_to_inf_edge(e_no))-1][abs(e_no)-1] += 1

            tt.branches[abs(e_no)].isotopy.append(inf_to_inf_edge(e_no))



    perron_eigenvector = perron_frobenius.perron_eigenvector(trans_mat)

    growth_rate = perron_frobenius.vector_norm(perron_eigenvector)

    #perron_eigenvector = perron_frobenius.normalize(perron_eigenvector)

    if debug:

        print 'growth rate:', growth_rate

    

    for i in range(len(perron_eigenvector)):

        tt.branches[i+1].weight = Decimal(str(perron_eigenvector[i])).quantize(Decimal('10e-' + str(weight_dec)))

    getcontext().prec = weight_dec

    

    # turn our train track into a trivalent one

    tt.make_trivalent()

    

    """

    for x in tt.fundamental_cycles():

        image = reduce(lambda x,y: x+y, map(lambda i: tt.getIsotopy(i), x))

        print image

        if len(image) == 0:

            continue

        if tt.init_vert(image[0]) != tt.term_vert(image[-1]):

            print "ERROR IN ISOTOPY PATHS"

        for i in range(len(image)-1):

            if tt.term_vert(image[i]) != tt.init_vert(image[i+1]):

                print 'ERROR IN ISOTOPY PATHS'

    print "PASSED"

    """

    #tt.print_traintrack()

    #print 'growth_rate', growth_rate

    return [tt, growth_rate, punc_bij]

    



def splitting_sequence(tt, growth_rate, punc_bij, debug=True):

    # returns a list of (list of [train track, branch], emaps)

    tts = [tt]

    split_seq = []

    found = False

    for i in xrange(10000):

        for j in range(i):

            sys.stdout.flush()

            #if i<220:

            #    continue

            emap = tts[j].equivalent(tts[i], punc_bij, growth_rate)

            if emap != False:

                no_tets = sum([len(x) for x in split_seq[j:i]])

                #if debug:

                #    print 'periodic sequence found: traintracks', j, 'and', i

                #    print 'number of tetrahedra:', no_tets

                #    print 'track:', j



                emaps = tts[j].all_equivalent(tts[i], punc_bij, growth_rate)

                tt_new_copy = tts[j].copy()

                for b in tt_new_copy.branches.values():

                    b.isotopy = [b.n]

                for k in range(i-j):

                    tt_new_copy.max_split(image_iso=False)

                #tts[j].print_traintrack()

                #tt_new_copy.print_traintrack()

                

                emaps = tts[j].all_equivalent(tt_new_copy, punc_bij, growth_rate)

                #if debug:

                #    print 'found edge maps:', emaps

                #    print 'fundamental cycles:', tts[j].fundamental_cycles()

                correct_emaps = tts[j].filter_by_isotopy(tt_new_copy, emaps)

                if debug:

                    print 'periodic sequence found: traintracks', j, 'and', i

                    print 'Train track bijection map of edges (i -> i-th element of list):', correct_emaps

                    iterator = i
                    while(iterator >= j):
                        print '\n\nTrain track', iterator, ':'

                        print tts[iterator].print_traintrack()

                        iterator = iterator - 1

                #print 'correct_emaps:', correct_emaps

                if len(correct_emaps) > 0:

                    return (split_seq[j:i] + [[[tts[i]]]], correct_emaps)

        ttc = tts[i].copy()

        #if debug:

        #	print '\ntrain track: ', i, ' isotopy info size:', ttc.isotopy_size()

        #	print '\ntrain track: ', i

        #	ttc.print_traintrack()

        #if debug:

        #    ttc.check_paths()

        split_seq.append(ttc.max_split())

        tts.append(ttc)

    raise Exception('Could not detect periodic splitting sequence within 10000 splits.')



"""

tt, growth_rate, punc_bij = parse_train_track(fname='C:/pr2011/veering_code/two_punc.txt', debug=True)

#tt, growth_rate, punc_bij = parse_train_track(

#        fname='C:/pr2011/Research/Agol Construction/bigdecimal version/rand3_braid37.txt', debug=True)

print punc_bij

splitting_sequence(tt, growth_rate, punc_bij, debug=True)

"""