
import os
import timeit
import sys, re
import math
from graph_map import *
from optparse import OptionParser

"""
    This file contains methods to 
        (i) create lists of words (Dehn twist words and braid words),
        (ii) then produce:
            (a) graph map files for each of those words (in the case of non-braids)
                which can be loaded into Trains.
            (b) a Trains batch file to process all homeomorphisms corresponding
                to the list of words. Trains will then output a single text file
                containing train tracks for every pA in the list.
                This file can be used as input to produce triangulations.
"""
    

# stores our list of braids in here
braid_lst = []

# return true if b is in braid_lst
# assumes that braid_lst is sorted
# return index of b in braid_lst, -1 otherwise
def binary_search(b, braid_lst):
    lo, hi = 0, len(braid_lst)
    while lo < hi:
        mid = (lo+hi)/2
        midval = braid_lst[mid]
        if midval < b:
            lo = mid+1
        elif midval > b:
            hi = mid
        else:
            return mid
    return -1


# check if the given braid is cyclically reduced
# assumes braid length is greater than 1
def is_reduced(b):
    for i in range(len(b)-1):
        if b[i] == -b[i+1]:
            return False
    if b[0] == -b[-1]:
        return False
    return True
    
# gives a list of 2*len(b) braids which are
# conjugate to, or inverses of b
def equiv_braids(b):
    e = []
    for i in range(len(b)):
        rot_b = b[i:] + b[:i]
        inv_rot_b = rot_b[:]
        inv_rot_b.reverse()
        inv_rot_b = [-x for x in inv_rot_b]
        e.append(rot_b)
        e.append(inv_rot_b)
    return e

    
# n, number of strands
# depth is length of braid word
# stores in braid_lst
def braids(n, depth, braid_lst=[], b=[]):
    if depth == 0:
        if is_reduced(b):
            for equiv_b in equiv_braids(b):
                if binary_search(equiv_b, braid_lst) != -1:
                    return None
            braid_lst.append(b)
        return None
        
    for i in range(-n+1,n):
        if i != 0:
            braids(n, depth-1, braid_lst, b[:]+[i])
    return None
    

def homeoString(allowed, depth, braid_lst=[], b=[]):
    if depth == 0:
        if is_reduced(b):
            for equiv_b in equiv_braids(b):
                if binary_search(equiv_b, braid_lst) != -1:
                    return None
            braid_lst.append(b)
        return None
        
    for i in allowed:
        if i != 0:
            homeoString(allowed, depth-1, braid_lst, b[:]+[i])
    return None    
    
    
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

            
def genDTorGMs(bfilename='S(2,1)_batch', tfilename='tracks.txt', path='', n=4):
    lst = []
    commute = [(1,3),(1,4),(1,5),(2,4),(2,5),(3,1),(3,5),(4,1),(4,2),(5,1),(5,2),(5,3)]

    for i in range(3,n+1):
        generate(lst, commute, [0]*i, 0, 5)

    #print len(lst)

    # The filter below excludes cases which are trivially not pA
    lst = filter(lambda x: ((-2 in x or 2 in x) and
                  (-3 in x or 3 in x) and
                  (-4 in x or 4 in x) and
                  x[0] != -x[-1] and
                  x[0] != x[1]),
            lst)
    print 'Number of words: ', len(lst)
    dTorSave(lst, bfilename, tfilename, path)
    
def dTorSave(lst, bfilename='S(2,1)_batch', tfilename='tracks.txt', path=''):
    f = open(path + '/' + bfilename, "w")
    f.write('TO ' + tfilename + '\n')
    f.write('IFPA\n')
    f.write('OUT g/d/s\n')
    f.write('PRINT Genus 2 surface with 1 puncture\n')
    for dehn_twist in lst:
        dTorusGraphMap(dehn_twist).saveAs(path + '/dtor'
                                            + str(dehn_twist).replace(' ','') + '.grm')
        f.write('PRINT Surface homeomorphism: ' + str(dehn_twist).replace(' ','') + '\n')
        f.write('LOAD dtor' + str(dehn_twist).replace(' ','') + '.grm' + '\n')
    
    f.flush()
    f.close()
    
def genTTorGMs(bfilename='S(3,1)_batch', tfilename='tracks.txt', path='', n=5):
    lst = []
    commute = [(1,3),(1,4),(1,5),(1,6),(1,7),
                (2,4),(2,5),(2,6),(2,7),
                (3,1),(3,5),(3,6),(3,7),
                (4,1),(4,2),(4,7),
                (5,1),(5,2),(5,3),(5,6),(5,7),
                (6,1),(6,2),(6,3),(6,5),
                (7,1),(7,2),(7,3),(7,4),(7,5)]

    for i in range(3,n+1):
        generate(lst, commute, [0]*i, 0, 7)

    lst = filter(lambda x: ((-2 in x or 2 in x) and
                        (-3 in x or 3 in x) and
                        (-4 in x or 4 in x) and
                        (-6 in x or 6 in x) and
                        (-7 in x or 7 in x) and
                        ((-1 in x or 1 in x) or (-5 in x or 5 in x)) and
                  x[0] != -x[-1] and
                  x[0] != x[1]),
            lst)
    print len(lst)
    tTorSave(lst, bfilename, tfilename, path)

def tTorSave(lst, bfilename='S(3,1)_batch', tfilename='tracks.txt', path=''):
    f = open(path + '/' + bfilename, "w")
    f.write('TO ' + tfilename + '\n')
    f.write('IFPA\n')
    f.write('OUT g/d/s\n')
    f.write('PRINT Genus 3 surface with 1 puncture\n')
    for dehn_twist in lst:
        g3GraphMap(dehn_twist).saveAs(path +
                                          '/ttor' + str(dehn_twist).replace(' ','') + '.grm')
        f.write('PRINT Surface homeomorphism: ' + str(dehn_twist).replace(' ','') + '\n')
        f.write('LOAD ttor' + str(dehn_twist).replace(' ','') + '.grm' + '\n')
    f.flush()
    f.close()



def str_braids(braids, individual_files):
    s = ''
    for b in braids:
        if individual_files:
            s += 'TO braids/braids[' + ','.join(map(str, b)) + '].txt\n'
        s += 'BR '+' '.join(map(str, b)) + ' 0\n'
    return s

# number of strings
def batch(strings, depth):
    lst = []
    braids(strings, depth, lst)
    print 'strings:',strings,'depth:',depth,'number:',len(lst)
    sys.stdout.flush()
    f = open('C:/pr2011/Research/Braids/braids.btt', 'a')
    f.write('OUT b/g\n')
    f.write('IFPA\n')
    f.write('PREC 14\n')
    f.write('STR ' + str(strings) + '\n')
    f.write('TO braids(' + str(strings) + ',' + str(depth) + ').txt\n')
    f.write(str_braids(lst, False))
    f.flush()
    f.close()
    

# run trains here

# this creates a btt file with the pseudo-anosovs, so that trains
# will output individual files for each pseudo-anosov
def batch_pa(strings, folder, filename):
    fpa = open(folder + '/' + filename + '.txt', 'r')
    f2 = open(folder + '/pa_braids.btt', 'a')

    f2.write('OUT b/g/d/s\n')
    f2.write('SS\n')
    f2.write('IFPA\n')
    f2.write('PREC 14\n')
    f2.write('STR '  + str(strings) + '\n')

    n_pa = 0
    for s in fpa:
        n_pa += 1
        f2.write('TO C:/pr2011/Research/Braids/braids/braids_' + str(strings) + '[' + s[7:].replace(' ', ',').strip() + '].txt\n')
        f2.write('BR '+ s[7:].strip() + ' 0\n')
    print "Number of pseudo-anosovs:", n_pa
    f2.flush()
    f2.close()
    
# like the above function but outputs to a single file
def batch_pa_one_file(strings, folder, filename):
    fpa = open(folder + '/' + filename + '.txt', 'r')
    f2 = open(folder + '/pa_braids.btt', 'a')

    f2.write('OUT b/g/d/s\n')
    f2.write('SS\n')
    f2.write('IFPA\n')
    f2.write('PREC 14\n')
    f2.write('STR '  + str(strings) + '\n')
    f2.write('TO C:/pr2011/Research/Braids/pa_train_info.txt\n')
    n_pa = 0
    for s in fpa:
        n_pa += 1
        f2.write('BR '+ s[7:].strip() + ' 0\n')
    print "Number of pseudo-anosovs:", n_pa
    f2.flush()
    f2.close()
    
def braid_batch_gen(path, bfilename='pa_batch.btt', tfilename='tracks.txt', strings=5, wlen=5):
    braid_list = []
    for d in range(2,wlen+1):
        braids(strings, d, braid_list)
    print 'Number of words:', len(braid_list)
    braid_gen_save(braid_list, path, bfilename, tfilename, strings)
    
def braid_gen_save(braid_list, path, bfilename='pa_batch.btt', tfilename='tracks.txt', strings=5):
    f = open(path + '/' + bfilename, 'w')
    f.write('OUT b/g/d/s\n')
    f.write('IFPA\n')
    f.write('PREC 14\n')
    f.write('STR '  + str(strings) + '\n')
    f.write('TO '+tfilename+'\n')    
    
    for s in braid_list:
        f.write('BR ' + str(s)[1:-1].replace(',','') + ' 0\n')
        
    f.flush()
    f.close()
    
    
def main():
    parser = OptionParser()
    parser.add_option('--surface',
                      type='choice',
                      dest='surface',
                      choices=['S(2,1)', 'S(3,1)', 'braid'],
                      help="Type of surface. (Options: 'S(2,1)', 'S(3,1)', 'braid')")
                      
    parser.add_option('--strings',
                      dest='strings',
                      default='6',
                      help='For braid: disk with this many punctures.')
                      
    parser.add_option('--wlen',
                      dest='wlen',
                      default='6',
                      help='Generates words in generators up to this length.')
                      
    parser.add_option('--bfilename',
                      dest='bfilename',
                      default='batch',
                      help="Trains batch file name e.g. 'batch.btt'.")
                      
    parser.add_option('--tfilename',
                      dest='tfilename',
                      default='tracks.txt',
                      help='File name of text file of train tracks to be produced by Trains.')
                      
    parser.add_option('--path',
                      dest='path',
                      help="Folder where produced graph map files and batch file are saved, e.g. 'C:/files'")
                      
    parser.add_option('--word',
                      dest='word',
                      default='',
                      help="Specify a single MCG word e.g. '[3,-4,-2,-1]' (no spaces).")
                      
    (options, args) = parser.parse_args()
    options = vars(options)
    word = options['word']
    bfilename = options['bfilename']
    tfilename = options['tfilename']
    path = options['path']
    wlen = int(options['wlen'])
    surface = options['surface']
    strings = int(options['strings'])
    if word != '':
        pword = map(int, word[1:-1].split(','))
    
    if surface == 'S(2,1)':
        if word != '':
            dTorSave([pword], bfilename, tfilename, path)
        else:
            genDTorGMs(bfilename, tfilename, path, n=wlen)
    elif surface == 'S(3,1)':
        if word != '':
            tTorSave([pword], bfilename, tfilename, path)
        else:
            genTTorGMs(bfilename, tfilename, path, n=wlen)
    elif surface == 'braid':
        if word != '':
            braid_gen_save([pword], path, bfilename, tfilename, strings)
        else:
            braid_batch_gen(bfilename, tfilename, path, wlen, strings=strings)
    
if __name__=='__main__':
    sys.exit(main())
    


