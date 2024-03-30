
from general_graph_map import *  
from triangulation_builder import *
import change_periph_curves as cpc

def enum_bij(lst):
    bijs = []
    
    if len(lst) == 0:
        return [[]]
        
    for i in range(len(lst)):
        for bij in enum_bij(lst[:i] + lst[(i+1):]):
            bijs.append([lst[i]] + bij)
    return bijs

def conjugate(tt1_c, tt2_c, period, growth_rate, punc_bij1, punc_bij2, emap1, emap2, strings1, strings2):
    tt1 = tt1_c.copy()
    tt2 = tt2_c.copy()
        
    if len(punc_bij1) != len(punc_bij2):
        return False
    if len(emap1) != len(emap2):
        return False
    if strings1 != strings2:
        return False
        
    """
        Are are initial punctures on our surface, and additional punctures.
        For two homeomorphisms to be conjugate, the conjugating map must preserve
        the initial punctures on our surface.
        The puncture numbered 0 must be preserved for braids and once-punctured surfaces.
    """
    init_punc_bij = enum_bij(range(0,strings1))
    addit_punc_bij = enum_bij(range(strings1,len(punc_bij1)))
    bijs = [a+b for a in init_punc_bij for b in addit_punc_bij]

    def image(edge_map, signed_edge):
        if signed_edge < 0:
            return -edge_map[-signed_edge]
        else:
            return edge_map[signed_edge]
            

    for i in range(period):
        for bij in bijs:
            # want high tolerance, might make program do a bit more work, but it'll ensure correct result
            #equivs = tt1.all_equivalent(tt2, bij, tt1.magnitude()/tt2.magnitude(), tol=Decimal('10e50'))
            equivs = tt1.all_equivalent_unweighted(tt2, bij)
            if equivs != False:
                for equiv in equivs:
                    #print 'found equiv:', equiv
                    equality_holds = True
                    for j in range(1,len(emap2)):
                        # check psi(j) = f o phi o f^-1 (j), where f is equiv, phi is emap1, psi emap2
                        if image(equiv, image(emap1,j)) != image(emap2, image(equiv, j)):
                            equality_holds = False
                            break
                    if equality_holds:
                        #print 'Found to be conjugate.'
                        return bij
        tt1.max_split()
    return False

class ConjClass:
    def __init__(self):
        self.conj_classes = defaultdict(list)
        self.n = 0
        # [tt, emap, punc_bij, growth]
        
    def add(self, tt2, emap2, punc_bij2, growth2, period2, word2, strings2):
        i = int(growth2*(10**3))
        for (tt1, emap1, punc_bij1, growth1, period1, word1, c_no, strings) in self.conj_classes[i]:
            if conjugate(tt1, tt2, period1, growth1, punc_bij1, punc_bij2, emap1, emap2, strings, strings2):
                return (c_no, False)
                    
        # not conjugate to anything in our list
        self.n += 1
        self.conj_classes[i].append([tt2.copy(), emap2, punc_bij2, growth2, period2, word2, self.n, strings2])
        return (self.n, True)

class IsomClass:
    def __init__(self):
        self.n = 0 # number of manifolds added
        self.mfolds = {} # dictionary indexed by volume upto 5 decimal places?
        # each value is a pair (manifold, isom class)
            
    def add(self, M, vol):
        """ returns a pair (n, True/False) where True indicates
            that M represents a new isometry class, n is class no. """
        i = int(vol*(10**5))
        if i not in self.mfolds.keys():
            self.n += 1
            self.mfolds[i] = [(M.copy(), self.n)]
            return (self.n, True)
        else:
            for N, m in self.mfolds[i]:
                if N.is_isometric_to(M):
                    return (m, False)
            self.n += 1
            self.mfolds[i].append((M.copy(), self.n))
        return (self.n, True)

def main():
    print("Veering. July 2013.\n\n"
          "Constructs Agol's veering triangulations for punctured orientable\n"
          "surfaces. Relies on the software Trains by Toby Hall.\n\n"
          "Type 'help' for a list of commands")
    
    sys.stdout.flush()
    TOLERANCE = 0.01
    while True:
        cmd = sys.stdin.readline()
        if not cmd:
            break
        cmd = cmd.strip()
        if cmd.startswith("h"): # help
            print("Type 'build' to construct a veering triangulation, and save it to a file.\n" +
                  "Type 'conjugacy' to test whether two pA mapping classes are conjugate in the mapping class group.\n" +
                  "Type 'exit' to quit.")
        elif cmd.startswith("b"): # build
            genus = int(raw_input("Surface genus:\n"))
            puncs = int(raw_input("Number of punctures:\n"))
            if genus < 0 or puncs <= 0:
                raise Exception("Invalid surface type. Genus must be greater " +
                                "than or equal to 0. Number of punctures must" +
                                " be strictly greater than 0.")
            word = raw_input("Enter pA map as product of Dehn twists and half twists:\n")
            gm = surfaceGM(genus, puncs, word)
            lines, is_pA = gm.sendToTrains()
            if not is_pA:
                print "Mapping class is not pseudo-Anosov."
            else:
                fileaddress = raw_input("Enter filename of triangulation file to be saved:\n")
                gm.saveTriangulation(fileaddress, lines=lines, tri_title='g='+str(genus)+', p='+str(puncs)+', monodromy='+word)
        elif cmd.startswith("c"): # conjugacy
            genus = int(raw_input("Surface genus:\n"))
            puncs = int(raw_input("Number of punctures:\n"))
            if genus < 0 or puncs <= 0:
                raise Exception("Invalid surface type. Genus must be greater " +
                                "than or equal to 0. Number of punctures must" +
                                " be strictly greater than 0.")
            word1 = raw_input("Enter first pA map as product of Dehn twists:\n")
            word2 = raw_input("Enter second pA map as product of Dehn twists:\n")
            gm1 = surfaceGM(genus, puncs, word1)
            si1 = gm1.splitting_info()
            gm2 = surfaceGM(genus, puncs, word2)
            si2 = gm2.splitting_info()
            if not si1:
                print 'First map is not pA.',
                if not si2:
                    print 'Second map is not pA. Unable to determine conjugacy.'
                else:
                    print 'Second map is pA. Not conjugate'
                continue
            if not si2:
                print 'Second map is not pA.',
                if not si1:
                    print 'First map is not pA. Unable to determine conjugacy.'
                else:
                    print 'First map is pA. Not conjugate.'
                continue
            ss1, growth1, punc_bij1 = si1
            ss2, growth2, punc_bij2 = si2
            # length of periodic splitting sequence is len(ss[0])-1
            if len(ss1[0]) != len(ss2[0]):
                print 'Different periodic splitting sequence length. Not conjugate.'
                continue
            elif abs(growth1 - growth2) > Decimal(str(TOLERANCE)):
                print 'Different growth rates. Not conjugate.'
                print growth1
                print growth2
                continue
            tt1 = ss1[0][0][0][0]
            tt2 = ss2[0][0][0][0]

            is_conj = conjugate(tt1, tt2, len(ss1[0])-1, growth1, punc_bij1, punc_bij2, 
                                ss1[1][0], ss2[1][0], puncs, puncs)
            #conjugate(tt1_c, tt2_c, period, growth_rate, punc_bij1, punc_bij2, emap1, emap2, strings1, strings2)
            if is_conj:
                print 'Found to be conjugate.'
            else:
                print 'Mapping classes are not conjugate.'
        elif cmd == "exit" or cmd == "quit":
            sys.exit(0)

def toNewConvention(lst):
    word = ''
    for i, w in enumerate(lst):
        if i > 0:
            word += ' '
        if w < 0:
            word += 'P' + str(abs(w))
        else:
            word += 'p' + str(abs(w))
    return word
                 
conjclasses = ConjClass()

counter = 0
def check_all(g, p, cur, poss, depth):
    if depth != 0:
        for w in poss:
            global counter
            counter += 1
            if counter % 200 == 0:
                print 'counter:', counter, ' '.join(cur+[w])
            gm = surfaceGM(g, p, " ".join(cur+[w]))
            si = gm.splitting_info()
            if si: # if pA
                ss, growth, punc_bij= si
                tt = ss[0][0][0][0]
                period = len(ss[0])
                emap = ss[1][0]
                c_no, new_conj = conjclasses.add(tt, emap, punc_bij, growth, period, cur+[w], p)
                if new_conj:
                    gm.saveTriangulation("tri_S_2_1_len5/S"+str(g)+"_"+str(p)+"_"+''.join(cur+[w])+".tri",
                                        tri_title='g='+str(g)+', p='+str(p)+', monodromy='+' '.join(cur+[w]))            
            check_all(g, p, cur + [w], poss, depth-1)
            
def generate(words, commute, word, i, n):
    if i==len(word):
        words.append(word)
        return
        
    lastGen = n
    if (i>0):
        lastGen = word[i-1]
        
    for j in range(n,0,-1):
        if ((j>math.fabs(lastGen) and (j,math.fabs(lastGen)) in commute) or 
                (i>0 and j>word[0])):
            continue
            
        if (j != -lastGen):
            word_copy1 = word[:]
            word_copy1[i] = j
            generate(words, commute, word_copy1, i+1, n)
            
        if (j != lastGen):
            word_copy2 = word[:]
            word_copy2[i] = -j
            generate(words, commute, word_copy2, i+1, n)        

            
def genDTorGMs(n=7):
    lst = []
    commute = [(1,3),(1,4),(1,5),(2,4),(2,5),(3,1),(3,5),(4,1),(4,2),(5,1),(5,2),(5,3)]

    for i in range(3,n+1):
        generate(lst, commute, [0]*i, 0, 5)

    #print len(lst)

    lst = filter(lambda x: ((-2 in x or 2 in x) and
                  (-3 in x or 3 in x) and
                  (-4 in x or 4 in x) and
                  x[0] != -x[-1]),
            lst)
    print len(lst)
    
    global counter
    conv = {}
    conv[1] = 'c2'
    conv[2] = 'b2'
    conv[3] = 'c1'
    conv[4] = 'b1'
    conv[5] = 'a1'
    conv[-1] = 'C2'
    conv[-2] = 'B2'
    conv[-3] = 'C1'
    conv[-4] = 'B1'
    conv[-5] = 'A1'
    g=2
    p=1
    for w in lst:
        w_old = w
        w = convert_to_new(w)
        counter += 1
        if counter % 200 == 0:
            print 'counter:', counter, ' '.join(w)
        gm = surfaceGM(2, 1, " ".join(w))
        si = gm.splitting_info()
        if si: # if pA
            ss, growth, punc_bij= si
            tt = ss[0][0][0][0]
            period = len(ss[0])
            emap = ss[1][0]
            c_no, is_new = conjclasses.add(tt, emap, punc_bij, growth, period, w, p)
            if is_new:
                filename = "S" + str(g) + "_" + str(p) + " " + str(w_old).replace(' ', '') + ".tri"
                gm.saveTriangulation("C:/pr2011/veering_code/tri_S_2_1_len7/" + filename,
                                tri_title='g='+str(g)+', p='+str(p)+', monodromy='+' '.join(w)) 

def convert_to_old(w):
    conv = {}
    conv[1] = 'c2'
    conv[2] = 'b2'
    conv[3] = 'c1'
    conv[4] = 'b1'
    conv[5] = 'a1'
    conv[-1] = 'C2'
    conv[-2] = 'B2'
    conv[-3] = 'C1'
    conv[-4] = 'B1'
    conv[-5] = 'A1'
    res = []
    for x in w:
        for key in conv:
            if conv[key] == x:
                res.append(key)
    return res
    
def convert_to_new(w):
    conv = {}
    conv[1] = 'c2'
    conv[2] = 'b2'
    conv[3] = 'c1'
    conv[4] = 'b1'
    conv[5] = 'a1'
    conv[-1] = 'C2'
    conv[-2] = 'B2'
    conv[-3] = 'C1'
    conv[-4] = 'B1'
    conv[-5] = 'A1'
    return [conv[x] for x in w]
               


def write_info_to_spreadsheet(words, excel_file, g=2, p=1, out_folder = "tri_S_2_1_len7/"):
    import xlwt, snappy
    wbk = xlwt.Workbook()
    sheet = wbk.add_sheet('S({0}, {1})'.format(g, p))
    isomclasses = IsomClass()
    conjclasses = ConjClass()

    columns = ['twists', 'twist_old', 'word_length', 'growth', 'entropy', 'vol_unf', 'vol_addit_fil',
                'n_sing', 'num_cusps', 'num_tet', 'simp_num_tet', 'solution_type', 'sing_info', 'basis',
                'punc_bij', 'sing_to_cusp', 'isom_twister', 'conj_class', 'isom_class']
    col_index = {}
    for i in range(len(columns)):
        col_index[columns[i]] = i
                
    sheet.write(0, col_index['twists'], 'Dehn twist word')
    sheet.write(0, col_index['twist_old'], 'Dehn twist word (old format)')
    sheet.write(0, col_index['word_length'], 'Word length')
    sheet.write(0, col_index['vol_unf'], 'Volume (unfilled)')
    sheet.write(0, col_index['growth'], 'Growth rate')
    sheet.write(0, col_index['entropy'], 'Entropy')
    sheet.write(0, col_index['basis'], 'Change of basis for peripheral curves')
    sheet.write(0, col_index['n_sing'], 'No. singularities in foliation')
    sheet.write(0, col_index['punc_bij'], 'Permutation of singularities')
    sheet.write(0, col_index['sing_to_cusp'], 'Cusp no. of singularities')
    sheet.write(0, col_index['num_tet'], 'No. of tetrahedra in Agol triangulation')
    sheet.write(0, col_index['num_cusps'], 'No. of cusps')
    sheet.write(0, col_index['vol_addit_fil'], 'Volume (with additional cusps filled)')
    sheet.write(0, col_index['solution_type'], 'Solution type')
    sheet.write(0, col_index['sing_info'], 'Prongs at nth singular point')
    sheet.write(0, col_index['isom_twister'], 'Verified upto isometric (using Twister):')
    sheet.write(0, col_index['conj_class'], 'Conjugacy class:')
    sheet.write(0, col_index['simp_num_tet'], 'Simplified no. tet:')
    sheet.write(0, col_index['isom_class'], 'Isom class:')

    conjclasses = ConjClass()
    k = -1
    conj_class = 0
    for w in words:
        k += 1
        gm = surfaceGM(g, p, " ".join(w))
        si = gm.splitting_info()
        if si: # if pA
            ss, growth, punc_bij= si
            tt = ss[0][0][0][0]
            period = len(ss[0])-1
            emap = ss[1][0]
            sing_info = tt.sing_info()
            strings = tt.strands
            sing_to_cusp = tt.cusp_info(punc_bij)
            n_sing = len(sing_to_cusp)
            conj_class, is_new = conjclasses.add(tt, emap, punc_bij, growth, period, w, p)
            if is_new: # take care of inversions
                old_format = convert_to_old(w)
                old_format.reverse()
                old_format = [-x for x in old_format]
                w2 = convert_to_new(old_format)
                gm2 = surfaceGM(g, p, " ".join(w2))
                si2 = gm2.splitting_info()
                ss2, growth2, punc_bij2 = si2
                tt2 = ss2[0][0][0][0]
                period2 = len(ss2[0])-1
                emap2 = ss2[1][0]
                conjclasses.add(tt2, emap2, punc_bij2, growth2, period2, w2, p)
            else:
                print 'word: ', w, 'is not new'
                k -= 1
                continue
            twist_info = ' '.join(w)
            filename = "S"+str(g)+"_"+str(p)+"_"+str(convert_to_old(w)).replace(' ','')+".tri"
            gm.saveTriangulation(out_folder + filename, tri_title='g='+str(g)+', p='+str(p)+', monodromy='+' '.join(w))
            change_of_basis_matrices = periph_curves(out_folder + filename, tt, growth, punc_bij, ss)
            M = snappy.Manifold(out_folder + filename)
            
            sheet.write(k+1, col_index['conj_class'], k+1)
            sheet.write(k+1, col_index['twists'], twist_info)
            sheet.write(k+1, col_index['sing_info'], str(sing_info))
            sheet.write(k+1, col_index['twist_old'], str(convert_to_old(w)).replace(' ',''))
            N = M.copy()
            N.simplify()
            sheet.write(k+1, col_index['vol_unf'], N.volume())
            sheet.write(k+1, col_index['isom_class'], isomclasses.add(N, N.volume())[0])
            sheet.write(k+1, col_index['word_length'], len(w))

            filling_coeffs = [(0,0)]*M.num_cusps()
            for j in range(p, n_sing):
                filling_coeffs[sing_to_cusp[j]] = (1,0)
            N = M.copy()
            N.set_peripheral_curves(change_of_basis_matrices)
            N.dehn_fill(filling_coeffs)
            N = N.filled_triangulation()
            N.simplify()
            twister_word = '*'.join(toTwister(convert_to_old(w)))
            R = snappy.twister.twister(surface=(g,p), monodromy=twister_word)
            for j in range(100):
                if j == 99:
                    raise Exception('isometry check failed')
                try:
                    twister_verified = N.is_isometric_to(R)
                    break
                except:
                    N.randomize()
                    R.randomize()
                    N.simplify()
                    R.simplify()

            # If this is false then we weren't able to verify any were isometric to the twister manifold
            sheet.write(k+1, col_index['isom_twister'], str(twister_verified))
            
            sheet.write(k+1, col_index['vol_addit_fil'], N.volume())
            
            sheet.write(k+1, col_index['growth'], float(growth))
            sheet.write(k+1, col_index['entropy'], math.log(float(growth)))
            sing_perm = []
            for i in sorted(punc_bij.keys()):
                sing_perm.append(punc_bij[i])
            sheet.write(k+1, col_index['punc_bij'], str(sing_perm))
            
            sheet.write(k+1, col_index['num_cusps'], M.num_cusps())
            sheet.write(k+1, col_index['num_tet'], M.num_tetrahedra())
            sheet.write(k+1, col_index['solution_type'], M.solution_type())
            K = M.copy()
            K.simplify()
            sheet.write(k+1, col_index['simp_num_tet'], K.num_tetrahedra())
            sheet.write(k+1, col_index['basis'], str(change_of_basis_matrices))
            sheet.write(k+1, col_index['n_sing'], n_sing)
            sheet.write(k+1, col_index['sing_to_cusp'], str(sing_to_cusp))
    wbk.save(out_folder + excel_file)

def periph_curves(filename, tt, growth_rate, punc_bij, ss):
    try:
        import snappy
    except:
        print 'Could not load snappy module. Make sure snappy python is installed.'
        return None
    sing_to_cusp = tt.cusp_info(punc_bij)
    M = snappy.Manifold(filename)
    N = M.copy()
    g_eqns = M.gluing_equations()
    g_eqns = [[g_eqns[x,y] for y in range(g_eqns.shape[1])] for x in range(g_eqns.shape[0])]
    comp_eqns = g_eqns[-2*M.num_cusps():] # completeness equations in vector form
    change_of_basis_matrices = []
    for j in range(0, M.num_cusps()):
        punc = -1
        for cp in range(len(sing_to_cusp)):
            if sing_to_cusp[cp] == j:
                punc = cp
                break
                    
        new_curve = computeMeridian(tt, growth_rate, punc_bij, punc, ss)
        a,c,b,d = cpc.change_of_basis_matrix(new_curve, comp_eqns[2*j], comp_eqns[2*j+1])
        change_of_basis_matrices.append([(a,b), (c,d)])
    return change_of_basis_matrices

#out_folder = "C:/pr2011/veering_code/"
#write_info_to_spreadsheet([['a1','a1','b1','c1','c2','B2','B2'], ['a1','B1','B1','c1','b2','c1','c2']], 'homeo.xls', g=2, p=1, out_folder=out_folder)

if __name__=='__main__':
    sys.exit(main())