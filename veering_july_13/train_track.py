
# train_track

# NEW edge format:
# [image edges] has been removed, no longer needed
# [edge number, [vert1, vert2], [left puncture, right puncture], [isotopy path], weight]

# vert format:
# [vertex number, surround edge in cyclic (counter-clockwise) order, region vertex is in]

import math
import sys
import copy
from collections import defaultdict
try:
    from cdecimal import *
except ImportError:
    from decimal import *
    
def reverse(w):
    w = w[:]
    w.reverse()
    return [-x for x in w]

class Branch:
    #n = 0
    #b_ends = [0, 0]
    #puncs = [0, 0]
    #isotopy = []
    #weight = 0.0
    
    def __init__(self, n=0, puncs=[0,0], isotopy=[], b_ends=[0,0], weight=0.0):
        self.n=n
        self.puncs=puncs
        self.isotopy=isotopy
        self.b_ends=b_ends
        self.weight=weight
        
    def copy(self):
        bb = Branch(n=self.n, puncs=self.puncs[:], 
                   isotopy=self.isotopy[:], b_ends=self.b_ends[:], weight=self.weight)
        return bb
        
    def __str__(self):
        return ('edge_no=' + str(self.n) + ', puncs=' + str(self.puncs) +
               ', terminal_switches=' +
               str(self.b_ends) + ', w=' + str(self.weight) + ', isotopy=' + str(self.isotopy))

    def pull_tight(self):
        length = len(self.isotopy)
        for i in range(length-1):
            if self.isotopy[length-i-1] == self.isotopy[length-i-2]:
                del self.isotopy[length-i-1]
        return None

class Switch:
    #n = 0 # vertex number
    #region = 0
    #branches = []
    
    def __init__(self, n=0, region=0, branches=[], strands=0):
        self.branches = branches
        self.n = n
        self.region = region
        
    def is_trivalent(self):
        if len(self.branches[0]) + len(self.branches[1]) == 3:
            return True
        else:
            return False
            
    def small_branches(self):
        if len(self.branches[0]) == 2:
            return self.branches[0]
        else:
            return self.branches[1]
            
    def all_branches(self):
        return self.branches[0] + self.branches[1]
            
    def copy(self):
        ss = Switch(n=self.n, region=self.region, branches=copy.deepcopy(self.branches))
        return ss
        
    def __str__(self):
        return 'vert_no='+str(self.n)+', region='+str(self.region)+', branches='+str(self.branches)
    
class TrainTrack:
    #switches = {}
    #branches = {}
    #no_of_puncs = 0
    
    def isotopy_size(self):
        return sum([len(x.isotopy) for x in self.branches.values()])
    
    def tighten(self):
        for branch in self.branches.values():
            branch.isotopy = self.pulltight(branch.isotopy)
            
    def getIsotopy(self, b_no):
        b = self.branches[abs(b_no)]
        iso = b.isotopy[:]
        if b_no < 0:
            iso.reverse()
            iso = [-x for x in iso]
        return iso
    
    def shorten_isotopy(self):
        def getIsotopy(b_no):
            b = self.branches[abs(b_no)]
            iso = b.isotopy[:]
            if b_no < 0:
                iso.reverse()
                iso = [-x for x in iso]
            return iso
                
        def setIsotopy(b_no, iso):
            iso = iso[:]
            if b_no < 0:
                iso.reverse()
                iso = [-x for x in iso]
            self.branches[abs(b_no)].isotopy = iso
        
        for switch in self.switches.values():
            branches = switch.all_branches()
            isos = map(getIsotopy, branches)
            shortened = True
            while shortened:
                shortened = False
                if (len(isos[0]) > 0 and len(isos[1]) > 0 and len(isos[2]) > 0
                    and isos[0][0] == isos[1][0] and isos[1][0] == isos[2][0]):
                    shortened = True
                    del isos[0][0], isos[1][0], isos[2][0]
                    
                    for i in range(3):
                        if branches[i] == -branches[(i+1)%3]:
                            del isos[i][-1], isos[(i+1)%3][-1]
                    
            for i in range(3):
                setIsotopy(branches[i], isos[i])    
        self.tighten()
    
    
    def print_traintrack(self):
        print 'Switch information:'
        for switch in self.switches.values():
            print switch
        print 'Branch information:'
        for branch in self.branches.values():
            print branch
            


    def path_isotopy(self, path):
        iso = []
        for b_no in path:
            sub_iso = self.branches[abs(b_no)].isotopy[:]
            if b_no < 0:
                sub_iso.reverse()
                sub_iso = [-x for x in sub_iso]
            iso.extend(sub_iso)
        return iso
            
    # performs a splitting on branch number branch_no
    def __init__(self, switches=None, branches=None, no_of_puncs=0, strands=0):
        """ no_of_puncs one less than actual number for g>0 surface """
        if branches is not None:
            self.branches = branches
        else:
            self.branches = {}
        
        if switches is not None:
            self.switches = switches
        else:
            self.switches = {}
            
        self.no_of_puncs = no_of_puncs
        self.strands = strands # number of strands before extra singularities
        
        self.memo_sur_branches = {}
        self.uptodate = False # it's false if memo_sur_branches needs to be updated
        
    def get_terminal_switch(self, b_no):
        b = self.branches[abs(b_no)]
        if b_no >= 0:
            s_no = b.b_ends[1]
        else:
            s_no = b.b_ends[0]
        return self.switches[s_no]
    
    def copy(self):
        tt = TrainTrack()
        tt.no_of_puncs = self.no_of_puncs
        tt.strands = self.strands
        for key in self.switches.keys():
            tt.switches[key] = self.switches[key].copy()
            
        #tt.branches = {}
        for key in self.branches.keys():
            tt.branches[key] = self.branches[key].copy()

        return tt
        
    # return a list of (train track, branch)
    # where branch belongs to train track and splits
    def max_split(self, tol=Decimal('10e-17'), image_iso=True, debug=False):
        # see split() for image_iso info
        max_w = max([branch.weight for branch in self.branches.values()])
        
        splitting_branches = []
        for branch in self.branches.values():
            if abs(branch.weight - max_w) < tol:
                splitting_branches.append([self.copy(), branch.n])
                self.split(branch.n, image_iso=image_iso, debug=debug)
                #print 'SPLIT',branch.n
                #self.print_traintrack()
        return splitting_branches
        
        
    """  all edges oriented outwards for e1_no, e2_no, e3_no, e4_no
       \                          /
     e1  \                      /  e3
           \         E        /
            v1 ------>------- v2
           /                  \
     e2  /                      \  e4
       /                          \
    """
    def split(self, branch_no, tolerance = Decimal('10e-20'), image_iso=True, debug=False):
        # image_iso is a flag, if true then .isotopy of edges
        # are computed in terms of the SPLIT track, rather than initial track.
        self.uptodate = False
        E = self.branches[branch_no] # branch we're going to split
        v1, v2 = self.switches[E.b_ends[0]], self.switches[E.b_ends[1]]
        e1_no, e2_no = self.small_branch(v1)
        e4_no, e3_no = self.small_branch(v2)
        e1, e2, e3, e4 = [self.branches[abs(b_no)] for b_no in [e1_no, e2_no, e3_no, e4_no]]
        max_w = max(e1.weight, e2.weight, e3.weight, e4.weight)

        E.weight = abs(e1.weight - e3.weight)
        if E.weight < tolerance: # FIX THIS TOLERANCE LATER
            raise Exception('Small weight encountered. Consider increasing precision.')
            
        case = 0 
        # case will be 1 if new edge goes from top-left to bottom-right
        # case will be 2 if bottom-left to top-right
        if e1.weight == max_w or e4.weight == max_w:
            case = 1
            if debug:
                print 'splitting:', branch_no, ' from ', e1_no, 'to', e4_no
            E.puncs = [self.left_puncture(e4_no), self.left_puncture(e1_no)]
            v1.branches = [[e1_no], [E.n, e3_no]]
            v2.branches = [[e4_no], [-E.n, e2_no]]

            self.extend(-e2_no, E.n) # extend -e2 by E
            self.extend(-e3_no, -E.n)
            
            #keep_cond = [e1_no,-e1_no,e4_no,-e4_no] # for isotopy updating below
        else:
            #print 'e2 or e3'
            case = 2
            if debug:
                print 'splitting:', branch_no, ' from ', e2_no, 'to', e3_no
            E.puncs = [self.right_puncture(e2_no), self.right_puncture(e3_no)]
            v1.branches = [[e2_no], [e4_no, E.n]]
            v2.branches = [[e3_no], [e1_no, -E.n]]
            
            self.extend(-e1_no, E.n)
            self.extend(-e4_no, -E.n)
            
            #keep_cond = [e2_no,-e2_no,e3_no,-e3_no] # for isotopy updating below
            
        if not image_iso:
            return None
            
        self.shorten_isotopy()
        for b in self.branches.values():
            iso = b.isotopy[:]
            n_changes = 0
            new_iso = iso[:]
            for i in range(len(iso)):
                if case == 1:
                    if iso[i] == e3_no:
                        new_iso.insert(i+n_changes, -E.n)
                        n_changes += 1
                    elif iso[i] == e2_no:
                        new_iso.insert(i+n_changes, E.n)
                        n_changes += 1
                    if iso[i] == -e2_no:
                        new_iso.insert(i+n_changes+1, -E.n)
                        n_changes += 1
                    elif iso[i] == -e3_no:
                        new_iso.insert(i+n_changes+1, E.n)
                        n_changes += 1
                elif case == 2:
                    if iso[i] == e1_no:
                        new_iso.insert(i+n_changes, E.n)
                        n_changes += 1
                    elif iso[i] == e4_no:
                        new_iso.insert(i+n_changes, -E.n)
                        n_changes += 1
                    if iso[i] == -e1_no:
                        new_iso.insert(i+n_changes+1, -E.n)
                        n_changes += 1
                    elif iso[i] == -e4_no:
                        new_iso.insert(i+n_changes+1, E.n)
                        n_changes += 1
            if n_changes > 0:
                b.isotopy = self.pulltight(new_iso)
            
    def check(self):
        for switch in self.switches.values():
            sub = abs(sum([self.branches[abs(b_no)].weight for b_no in switch.branches[0]])
                      - sum([self.branches[abs(b_no)].weight for b_no in switch.branches[1]]))
            if sub > Decimal('10e-30'):
                print 'error at switch:', switch.n
        
    # returns a small branch end of switch
    def small_branch(self, switch):
        if len(switch.branches[0]) == 2:
            return switch.branches[0]
        else:
            return switch.branches[1]
        
    
    """        b1_no
    -------->--------
                      \
                        \   b2_no
                         >\
                            \
    """
    # if b1_no < 0 then extends extends by b2_no on b1's tail
    def extend(self, b1_no, b2_no, calc_isotopy=True):
        b1 = self.branches[abs(b1_no)]
        b2 = self.branches[abs(b2_no)]
        
        if b1_no > 0:
            init_switch = b1.b_ends[0]
            if b2_no > 0:
                end_switch = b2.b_ends[1]
            else:
                end_switch = b2.b_ends[0]
        else:
            end_switch = b1.b_ends[1]
            if b2_no > 0:
                init_switch = b2.b_ends[1]
            else:
                init_switch = b2.b_ends[0]

        b1.b_ends = [init_switch, end_switch]
        # fix b1's new isotopy class
        
        if not calc_isotopy: # we don't want to calculate isotopy
            return
            
        #print 'extending:', b1_no, 'by', b2_no
        if b1_no > 0:
            if b2_no > 0:
                b1.isotopy.extend(b2.isotopy)
            else:
                rev_isotopy = b2.isotopy[:] # make a copy
                rev_isotopy.reverse() # reverse order
                rev_isotopy = [-x for x in rev_isotopy]
                b1.isotopy = b1.isotopy + rev_isotopy
        else:
            if b2_no > 0:
                rev_isotopy = b2.isotopy[:] # make a copy
                rev_isotopy.reverse() # reverse order
                rev_isotopy = [-x for x in rev_isotopy]
                b1.isotopy = rev_isotopy + b1.isotopy
            else:
                b1.isotopy = b2.isotopy + b1.isotopy
        
    
    # make_trivalent deletes from branches
    # this function relabels the train track so that
    # the branches are labelled from 1 upto n
    def simplify_branch_labels(self):    
        b_nos = sorted([b.n for b in self.branches.values()])
        branch_bij = {}
        for i in range(1, len(self.branches)+1):
            branch_bij[b_nos[i-1]] = i
            
        # two things need to be updated: 
        # i) number of each branch
        # ii) and labels incident branches to each switch


        # update every branch's isotopy so that branch.n is replaced by branch_bij[branch.n]
        for b in self.branches.values():
            b.isotopy = [branch_bij[x] if x > 0 else -branch_bij[-x] for x in b.isotopy]
            b.n = branch_bij[b.n]
            
        new_branches = {}
        for branch in self.branches.values():
            new_branches[branch.n] = branch
            
        self.branches = new_branches
        
        for switch in self.switches.values():
            for j in range(2):
                for i in range(len(switch.branches[j])):
                    if switch.branches[j][i] > 0:
                        switch.branches[j][i] = branch_bij[switch.branches[j][i]]
                    else:
                        switch.branches[j][i] = -branch_bij[abs(switch.branches[j][i])]
        
    
    """
        Used in make_trivalent.
            init - initial oriented edge
            mid - middle oriented edge
            fin - final oriented edge
    """
    def repair_isotopy(self, init, mid, fin):
        for branch in self.branches.values():
            iso = branch.isotopy
            iso_c = iso[:]
            n_changes = 0
            for i in range(len(iso)-1):
                if iso[i] == init and iso[i+1] == fin:
                    iso_c.insert(i+1+n_changes, mid)
                    n_changes += 1
            branch.isotopy = iso_c
                        
    """
     a train track coming from Trains may have switches which are not trivalent
     in which case one side of the switch has more than 2 branches
     pick any two adjacent branches, attach them to a new branch
     which is connected .
     Now we're dealing with bivalent switches by removing the switch
     and turning the two branches into a single branch by deleting one of the branches
     and relabelling the terminal vertex of the other branch
    """
    def make_trivalent(self):
        switch = self.get_non_trivalent()
        while switch != None:
            if len(switch.branches[0]) + len(switch.branches[1]) == 2:
                # bivalent switch case, turn it into a single branch and remove the vertex
                # this is similar to a valence 2 isotopy in Bestvina-Handel's paper
                b0_no, b1_no = switch.branches[0][0], switch.branches[1][0]
                branch_0 = self.branches[abs(b0_no)]
                branch_1 = self.branches[abs(b1_no)]
                
                # we're going to delete branch_1 so let's reassign branch_0's terminal vertex
                if b1_no > 0:
                    term_vert = branch_1.b_ends[1]
                else:
                    term_vert = branch_1.b_ends[0]

                if b0_no > 0:        
                    branch_0.b_ends[0] = term_vert
                else:
                    branch_0.b_ends[1] = term_vert
                    
                # update isotopy of branch_0, 4 cases
                if b0_no > 0 and b1_no > 0:
                    branch_1.isotopy.reverse()
                    branch_1.isotopy = [-e for e in branch_1.isotopy]
                    branch_0.isotopy = branch_1.isotopy + branch_0.isotopy
                elif b0_no > 0 and b1_no < 0:
                    branch_0.isotopy = branch_1.isotopy + branch_0.isotopy
                elif b0_no < 0 and b1_no > 0:
                    branch_0.isotopy = branch_0.isotopy + branch_1.isotopy
                elif b0_no < 0 and b1_no < 0:
                    branch_1.isotopy.reverse()
                    branch_1.isotopy = [-e for e in branch_1.isotopy]
                    branch_0.isotopy = branch_0.isotopy + branch_1.isotopy
                    
                # also, reassign that terminal vertex's incident branches
                for j in range(2):
                    for i in range(len(self.switches[term_vert].branches[j])):
                        if self.switches[term_vert].branches[j][i] == -b1_no:
                            self.switches[term_vert].branches[j][i] = b0_no
                
                # delete the switch and branch_1
                del self.switches[switch.n]
                del self.branches[abs(b1_no)]
                
                for edge in self.branches.values():
                    edge.isotopy = filter(lambda e: e not in [b1_no,-b1_no], edge.isotopy)
                    
                #self.print_traintrack()
            else:
                # we're in the n-valent case, where n > 3
                #print '\nn-valent case:'
                #self.print_traintrack()
                new_branch_no = max([b.n for b in self.branches.values()])+1
                new_switch_no = max([s.n for s in self.switches.values()])+1
                #print '\nnew_branch_no', new_branch_no
                if len(switch.branches[0]) > 2:
                    branch = switch.branches[0]
                    # branch are the edges which will emanate from the NEW switch
                    op_branch = switch.branches[1]
                else:
                    branch = switch.branches[1]
                    op_branch = switch.branches[0]
                    
                for b3 in self.branches.values():
                    n_changes = 0
                    iso = b3.isotopy[:]
                    for i in range(len(b3.isotopy)):
                        x = b3.isotopy[i]
                        if x in branch[:2]:
                            iso.insert(i+n_changes, new_branch_no)
                            n_changes += 1
                        if x in [-y for y in branch[:2]]:
                            iso.insert(i+n_changes+1, -new_branch_no)
                            n_changes += 1
                    b3.isotopy = iso
                
                new_ends = [switch.n, new_switch_no] # from old switch to new switch
                new_weight = self.branches[abs(branch[0])].weight + self.branches[abs(branch[1])].weight
                new_puncs = [self.left_puncture(branch[1]), self.right_puncture(branch[0])]
                
                iso_info = []
                
                new_branch = Branch(n=new_branch_no, puncs=new_puncs, isotopy=iso_info, b_ends=new_ends, weight=new_weight)
                
                new_switch = Switch(n=new_switch_no, region=switch.region, 
                                    branches=[[-new_branch_no], branch[:2]])            
                
                self.branches[new_branch_no] = new_branch
                self.switches[new_switch_no] = new_switch
                
                if len(switch.branches[0]) > 2:
                    switch.branches[0] = [new_branch_no] + switch.branches[0][2:] # update branch 
                else:
                    switch.branches[1] = [new_branch_no] + switch.branches[1][2:] # update branch 
                
                # the ends of the 2 branches that we pulled get changed
                for i in [0,1]:
                    if branch[i] > 0:
                        self.branches[abs(branch[i])].b_ends[0] = new_switch_no
                    else:
                        self.branches[abs(branch[i])].b_ends[1] = new_switch_no
            #self.print_traintrack()
            self.simplify_branch_labels()
            switch = self.get_non_trivalent()
        return None
        
    def get_non_trivalent(self):
        for switch in self.switches.values():
            if not switch.is_trivalent():
                return switch
        return None
        
        
    # the branch_no may be negative, which means to the left of branch |branch_no| with opposite
    # orientation. Errors if can't find a branch with branch_no
    def left_puncture(self, branch_no):
        if branch_no > 0:
            return self.branches[abs(branch_no)].puncs[0]
        else:
            return self.branches[abs(branch_no)].puncs[1]
        
    def right_puncture(self, branch_no):
        return self.left_puncture(-branch_no)
        
    # returns a list of labels of oriented branches, going clockwise around
    # puncture number p. a branch is positively oriented if it
    # points clockwise around the puncture
    def sur_branches(self, p):
        #if self.memo_sur_branches.get(p, None) is not None and self.uptodate:
        #    return self.memo_sur_branches.get(p)
        e_no = 0
        for branch in self.branches.values():
            if branch.puncs[1] == p:
                e_no = branch.n
                break
            elif branch.puncs[0] == p:
                e_no = -branch.n
                break
                
        if e_no == 0:
            raise Exception('Did not find edge surrounding puncture.')
            return None
        
        def next_branch(b_no):
            if b_no > 0:
                switch = self.switches[self.branches[b_no].b_ends[1]]
            else:
                switch = self.switches[self.branches[-b_no].b_ends[0]]

            branches = switch.branches[0][:] + switch.branches[1][:]

            for i in range(len(branches)): # should loop 3 times, trivalent tt
                if branches[i] == -b_no:
                    return branches[(i+1)%len(branches)]        
            print 'failed in next_branch'
            # something failed if we reach here
            
        sur_b = []
        next_b = e_no
        while True:
            sur_b.append(next_b)
            next_b = next_branch(next_b)

            if next_b == e_no: # we've looped around the puncture
                break
                
        self.memo_sur_branches[p] = sur_b
        self.uptodate = True
        return sur_b
    
    def equivalent(self, t2, perm, growth, p=0, emap=[], tol=Decimal('10e-5')):
        # checks if growth*self = t2, up to combinatorial equivalence
        
        # has to respect permutation of punctures
        t1 = self
        if emap == []:
            emap = [0] * (len(t1.branches) + 1)
            
        if p >= t1.no_of_puncs: # for braids this should be > not >=
            return emap
            
        sur1 = t1.sur_branches(p)
        sur2 = t2.sur_branches(perm[p])
        
        if len(sur1) != len(sur2):
            return False
            
        found = False
        recip_growth = Decimal(1)/growth
        for i in range(len(sur1)):
            found = True
            emapc = emap[:]
            for j in range(len(sur1)):
                if abs(t1.branches[abs(sur1[j])].weight*recip_growth - t2.branches[abs(sur2[j])].weight) > tol: # FIX TOL
                    found = False
                    break
                elif emapc[abs(sur1[j])] != 0 and ((sur1[j] >= 0 and emapc[sur1[j]] != sur2[j]) or 
                        (sur1[j] < 0 and emapc[-sur1[j]] != -sur2[j])): # should be part of above if
                    found = False
                    break
                elif emapc[abs(sur1[j])] == 0:
                    if sur1[j]*sur2[j] >= 0:
                        if t1.branch_type(abs(sur1[j])) != t2.branch_type(abs(sur2[j])):
                            found = False
                            break
                        emapc[abs(sur1[j])] = abs(sur2[j])
                    else:
                        type1 = t1.branch_type(abs(sur1[j]))
                        type1.reverse()
                        if type1 != t2.branch_type(abs(sur2[j])):
                            found = False
                            break
                        emapc[abs(sur1[j])] = -abs(sur2[j])
                
            if found:
                edge_map = t1.equivalent(t2, perm, growth, p+1, emapc, tol)
                if edge_map != False:
                    return edge_map
                        
            # rotate sur1
            sur1.append(sur1.pop(0))
        return False
        
    def branch_type(self, b_no):
        # returns [a, b] where a is True if init vert end is large, else False. similar for b but with term vert
        b = self.branches[b_no]
        type = []
        for i in [0,1]:
            vert = b.b_ends[i]
            small = self.switches[vert].small_branches()
            if b_no in small or -b_no in small:
                type.append(True)
            else:
                type.append(False)
        return type
    
    def all_equivalent_unweighted(self, t2, perm, p=0, emap=[]):
        t1 = self
        all_emaps = []
        if emap == []:
            emap = [0] * (len(t1.branches) + 1)
            
        if p >= t1.no_of_puncs:
            return [emap]
            
        sur1 = t1.sur_branches(p)
        sur2 = t2.sur_branches(perm[p])

        if len(sur1) != len(sur2):
            return False
            
        found = False
        for i in range(len(sur1)):
            found = True
            emapc = emap[:]

            for j in range(len(sur1)):
                if emapc[abs(sur1[j])] != 0 and ((sur1[j] >= 0 and emapc[sur1[j]] != sur2[j]) or 
                        (sur1[j] < 0 and emapc[-sur1[j]] != -sur2[j])):
                    found = False
                    break
                elif emapc[abs(sur1[j])] == 0:
                    if sur1[j]*sur2[j] >= 0:
                        if t1.branch_type(abs(sur1[j])) != t2.branch_type(abs(sur2[j])):
                            found = False
                            break
                        emapc[abs(sur1[j])] = abs(sur2[j])
                    else:
                        type1 = t1.branch_type(abs(sur1[j]))
                        type1.reverse()
                        if type1 != t2.branch_type(abs(sur2[j])):
                            found = False
                            break
                        emapc[abs(sur1[j])] = -abs(sur2[j])
                
            if found:
                edge_maps = t1.all_equivalent_unweighted(t2, perm, p+1, emapc)
                if edge_maps != False:
                    all_emaps.extend(edge_maps)
                        
            # rotate sur1
            sur1.append(sur1.pop(0))
            
        if len(all_emaps) == 0:
            return False
        else:
            return all_emaps
            
    def all_equivalent(self, t2, perm, growth, p=0, emap=[], tol=Decimal('10e-5')):
        t1 = self
        all_emaps = []
        if emap == []:
            emap = [0] * (len(t1.branches) + 1)
            
        if p >= t1.no_of_puncs:
            return [emap]
            
        sur1 = t1.sur_branches(p)
        sur2 = t2.sur_branches(perm[p])

        if len(sur1) != len(sur2):
            return False
            
        found = False
        recip_growth = Decimal(1)/growth
        for i in range(len(sur1)):
            found = True
            emapc = emap[:]

            for j in range(len(sur1)):
                if abs(t1.branches[abs(sur1[j])].weight*recip_growth - t2.branches[abs(sur2[j])].weight) > tol: # FIX TOL
                    found = False
                    break
                elif emapc[abs(sur1[j])] != 0 and ((sur1[j] >= 0 and emapc[sur1[j]] != sur2[j]) or 
                        (sur1[j] < 0 and emapc[-sur1[j]] != -sur2[j])): # should be part of above if
                    found = False
                    break
                elif emapc[abs(sur1[j])] == 0:
                    if sur1[j]*sur2[j] >= 0:
                        if t1.branch_type(abs(sur1[j])) != t2.branch_type(abs(sur2[j])):
                            found = False
                            break
                        emapc[abs(sur1[j])] = abs(sur2[j])
                    else:
                        type1 = t1.branch_type(abs(sur1[j]))
                        type1.reverse()
                        if type1 != t2.branch_type(abs(sur2[j])):
                            found = False
                            break
                        emapc[abs(sur1[j])] = -abs(sur2[j])

            if found:
                edge_maps = t1.all_equivalent(t2, perm, growth, p+1, emapc, tol)
                if edge_maps != False:
                    all_emaps.extend(edge_maps)
                        
            # rotate sur1
            sur1.append(sur1.pop(0))
            
        if len(all_emaps) == 0:
            return False
        else:
            return all_emaps
           
       
    def sing_info(self):
        """ returns a list of length the number of punctures in the surface
            the i-th entry of the list is the number of cusps surrounding
            singularity number i"""
        info = []
        for p in range(self.no_of_puncs):
            branches = self.sur_branches(p)
            n_cusps = 0
            n_branches = len(branches)
            for i in range(n_branches):
                b_no = branches[i]
                if b_no < 0:
                    b_next = branches[(i+1)%n_branches]
                    s_no = self.branches[-b_no].b_ends[0]
                    s = self.switches[s_no].small_branches()
                    if -b_no in s and b_next in s:
                        n_cusps += 1
                else:
                    b_next = branches[(i+1)%n_branches]
                    s_no = self.branches[b_no].b_ends[1]
                    s = self.switches[s_no].small_branches()
                    if -b_no in s and b_next in s:
                        n_cusps += 1           
            info.append(n_cusps)
        return info
        
    def cusp_info(self, punc_bij):
        # returns a pair: a bijection from singularity number to cusp number
        # number of cusps 
            
        sing_to_cusp = [-1]*len(punc_bij)
        cusp = 0
        for i in range(len(punc_bij)):
            if sing_to_cusp[i] == -1:
                orbit = [i]
                while True:
                    j = punc_bij[orbit[-1]]
                    if j not in orbit:
                        orbit.append(j)
                    else:
                        break
                for j in orbit:
                    sing_to_cusp[j] = cusp
                cusp += 1
        return sing_to_cusp
        
    def magnitude(self):
        """returns the L^1 norm of the vector of weights of the branches of this train track"""
        return sum([b.weight for b in self.branches.values()])
        
    def init_vert(self, e):
        if e > 0:
            return self.branches[e].b_ends[0]
        else:
            return self.branches[-e].b_ends[1]
            
    def term_vert(self, e):
        return self.init_vert(-e)
        
        
    def fundamental_cycles(self, vert=None):
        """ returns a list of lists of edges, where each list of edges
            cyclically defines a loop. The set of such loops forms a generating set
            for the fundamental group of the train track (treated as a CW-complex).
            The loops are based at a single point.
            
            vert is the vertex we grow a spanning tree from.
        """
        if vert is None:
            vert = self.switches.values()[0].n # choose any vertex

        path_from_vert = defaultdict(list)
        v_visited = [vert] # list of vertices visited
        tree_edges = [] # list of edges in the spanning tree
        curr_vert = vert
        v_unvisited = self.switches.keys()
        v_unvisited = filter(lambda v: v != vert, v_unvisited)
        
        while len(v_unvisited) > 0:
            for curr_vert in v_visited:
                for b in self.switches[curr_vert].all_branches():
                    if abs(b) in tree_edges:
                        continue
                    s = self.get_terminal_switch(b).n
                    if s not in v_visited:
                        path_from_vert[s] = path_from_vert[curr_vert] + [b]
                        tree_edges.append(abs(b))
                        v_visited.append(s)
                        v_unvisited = filter(lambda v: v != s, v_unvisited)
                        
        all_edges = self.branches.values()
        fundamental_cycles = []
        for edge in all_edges:
            if edge.n in tree_edges:
                continue
                
            rev = path_from_vert[edge.b_ends[1]][:]
            rev.reverse()
            rev = [-x for x in rev]
            cycle = path_from_vert[edge.b_ends[0]] + [edge.n] + rev
            #cycle = TrainTrack.cyclic_pulltight(cycle)
            cycle = TrainTrack.pulltight(cycle)
            fundamental_cycles.append(cycle)
            
        return fundamental_cycles
       
    @staticmethod
    def pulltight(path):
        """ returns a tightened copy of a path. """
        path = path[:] # make a copy
        i = 0
        while i < len(path)-1 and i >= 0:
            if path[i] == -path[i+1]:
                del path[i:i+2]
                i -= 1
                i = max(0, i)
            else:
                i += 1
        return path 

        
    @staticmethod
    def cyclic_pulltight(path):
        """ returns a tightened copy of a cyclic path. """
        path = path[:] # make a copy
        i = 0
        while i < len(path) and i >= 0:
            if i == len(path)-1:
                if len(path) > 1 and path[i] == -path[0]:
                    del path[i]
                    del path[0]
                    i -= 2
                else:
                    i += 1
            elif path[i] == -path[i+1]:
                del path[i+1]
                del path[i]
                i -= 1
                i = max(0, i)
            else:
                i += 1
        return path
        

    @staticmethod
    def cyclic_equality(path1, path2):
        for i in range(len(path1)):
            if path1 == path2:
                return True
            path1.append(path1.pop(0)) # rotate
        return False
        
    def filter_by_isotopyz(self, tt2, all_emaps):
        cycles = self.fundamental_cycles()
        correct_emaps = []
        for emap in all_emaps:
            isotopic = True
            for cycle in cycles:
                p1 = self.cyclic_pulltight(self.path_isotopy(cycle))
                p2 = self.cyclic_pulltight(tt2.path_isotopy(map(lambda a: emap[a] if a > 0 else -emap[-a], cycle)))
                #print 'p1', p1
                #print 'p2', p2
                if not TrainTrack.cyclic_equality(p1, p2):
                    isotopic = False
                    #print 'failed on:', cycle, ' or (image) ', map(lambda a: emap[a] if a > 0 else -emap[-a], cycle)
                    #print p1
                    #print p2
                    break
            if isotopic:
                correct_emaps.append(emap[:])
        return correct_emaps
        
        
    def filter_by_isotopy(self, tt2, all_emaps, debug=False):
        def boundary_tighten(w):
            # returns (y,x) where w = xyx^-1 and y is cyclically reduced
            w = w[:]
            if len(w) <= 1:
                return (w, [])
            x = []
            while len(w) > 1 and w[0] == -w[-1]:
                x.append(w[0])
                del w[0]
                del w[-1]
            return (w, x)
    
        cycles = self.fundamental_cycles()
        correct_emaps = []
        for emap in all_emaps:
            if debug:
                print 'emap:',emap
            isotopic = True
            conj_factor = None
            conj_factors = []
            p2s = []
            common_conj = None # the common conjugating element
            for j, cycle in enumerate(cycles):
                p1 = self.pulltight(self.path_isotopy(cycle))
                p2 = self.pulltight(tt2.path_isotopy(map(lambda a: emap[a] if a > 0 else -emap[-a], cycle)))
                if common_conj != None:
                    if p2 != TrainTrack.pulltight(common_conj + p1 + reverse(common_conj)):
                        isotopic = False
                        break
                    else:
                        continue
                t_p1, x = boundary_tighten(p1)
                t_p2, y = boundary_tighten(p2)
                #print 't_p1',t_p1
                #print 't_p2',t_p2
                z = []
                found = False
                for i in range(len(t_p1)+1):
                    if t_p1 == t_p2:
                        found = True
                        break
                    else:
                        z.append(t_p1[0])
                        t_p1.append(t_p1.pop(0))
                conj_f = TrainTrack.pulltight(y + reverse(z) + reverse(x))
                # p2 = conj_f + p1 + (conj_f)^-1 (reduced)
                conj_factors.append(conj_f)
                p2s.append(p2)
                
                if conj_factor is None:
                    conj_factor = conj_f
                #if conj_f != conj_factor or not found:
                if not found:
                    if debug:
                        print 'did not find conjugating element'
                    isotopic = False
                    break
                
                # A common conjugating element x would be of the form
                # x = p2s[k]^m conj_factors[k] for all k>=0 and some m in Z (depending on k).
                # Thus we have:
                # p2s[0]^n conj_factors[0] = p2s[1]^m' conj_factors[1]
                # p2s[1]^m p2s[0]^n = conj_factors[1] (conj_factors[0])^-1, where m = -m'
                
                # 1) There is a unique pair (m,n) solving the above equation. This is because
                # {p2s[i] : all i} is a basis for the fundamental group i.e. the free group.
                
                # 2) If the RHS reduces to the empty list, i.e. the conj_factors are the same,
                # then m = n = 0.
                if j > 0:
                    if conj_factors[0] == conj_factors[1]:
                        if debug:
                            print 'conj factors were the same to start with'
                        common_conj = conj_factors[0][:]
                        continue
                
                # 3) From (2) the only problem is when conj_factors are all different.
                # The right hand side is of finite length, which we can compute.
                    RHS = TrainTrack.pulltight(conj_factors[1] + reverse(conj_factors[0]))

                    def solvePowers(w1, w2, RHS):
                    # Solves w1^m w2^n (reduced) = RHS returning (m,n) or None if no solution.
                        if not RHS:
                            return (0, 0)
                        t_p1, x = boundary_tighten(w1)
                        t_p2, y = boundary_tighten(w2)

                        def wpow(w, x, n):
                            if n == 0:
                                return []
                            if n < 0:
                                return x + reverse(w)*(-n) + reverse(x)
                            if n > 0:
                                return x + w*n + reverse(x)
               
                        m = 0
                        while True:
                            n = 0
                            while True:
                                all_greater = True
                                for sgn1 in [-1,1]:
                                    for sgn2 in [-1,1]:
                                        z1 = wpow(t_p1, x, sgn1*m)
                                        z2 = wpow(t_p2, y, sgn2*n)
                                        LHS = TrainTrack.pulltight(z1+z2)
                                        #print 'z1',z1,'z2',z2,'LHS', LHS, 'n=',n
                                        if LHS == RHS:
                                            return (sgn1*m,sgn2*n)
                                        if len(LHS) <= len(RHS) and len(LHS) > 0:
                                            all_greater = False
                                if all_greater and n > 0:
                                    break
                                else:
                                    n += 1
                            # We have to compute how much cancellation can occur
                            # in w1^m * w2^n for m fixed (possibly negated), and n in {-1,0,1}
                            c = min(len(TrainTrack.pulltight(w1*m)),
                                    len(TrainTrack.pulltight(w1*m + w2)),
                                    len(TrainTrack.pulltight(w1*m + reverse(w2))),
                                    len(TrainTrack.pulltight(reverse(w1)*m + reverse(w2))),
                                    len(TrainTrack.pulltight(reverse(w1)*m + w2)))
                            if c > len(RHS):
                                break
                            else:
                                m += 1
                        return None
                        
                    sol = solvePowers(p2s[1], p2s[0], RHS)
                    if not sol:
                        if debug:
                            print "Didn't find any solution to inner automorphism conjugation equation."
                        isotopic = False
                        break
                    else:
                        m, n = sol
                        if debug:
                            print 'm,n=',m,n
                        common_conj = TrainTrack.pulltight(p2s[0]*n + conj_factors[0])
                    
            if isotopic:
                if debug:
                    print 'conjugating element:', conj_factor
                correct_emaps.append(emap[:])
        return correct_emaps
        
    def check_paths(self):
        tt = self
        for x in tt.fundamental_cycles():
            image = reduce(lambda x,y: x+y, map(lambda i: tt.getIsotopy(i), x))
            #print image
            if len(image) == 0:
                continue
            if tt.init_vert(image[0]) != tt.term_vert(image[-1]):
                print "ERROR IN ISOTOPY PATHS"
            for i in range(len(image)-1):
                if tt.term_vert(image[i]) != tt.init_vert(image[i+1]):
                    print 'ERROR IN ISOTOPY PATHS'
                    print 'Edge path:', image[i],'to', image[i+1]
                    print 'fundamental cycle:', x
                    raise Exception('Error in isotopy paths')
                
        
        