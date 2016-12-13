

# STEPS=10000

from searchHelpers import *
import time

#NDATA = 6  # ?
TOP = 20
STEPS = 5000
BEST_HYP = 10
update_prior = False
analogy = False



def makeLsts1(start, nTimes, thrsh=0.1, totNum=10):
    lst = [start for i in xrange(nTimes)]
    for i in xrange(max(0, totNum - nTimes)):
        tmp = ""
        for j in xrange(len(start)):
            if random.random() < thrsh:
                tmp += str((start[j] == "0") * 1)
            else:
                tmp += str(start[j])


        lst.append(tmp)
    return lst


def transfExp():
    stims = ["010101010101", "000100010001", "000001000001", "101111111110",
             "110110110110", "011111101111", "001111110011", "000100010001",
             "000000011111", "111000000000", "111111111111", "111110000000",
             "001100011100", "100111000011", "111110000111", "011110011100",
             "111110101010", "010101011111", "100100111111", "111111110101",
             "100101000111", "110110110010", "011010111000", "001001001101"]

    endings = ["111111111110", "010010001000", "110110110010", "011000111100"]
    random.shuffle(stims)
    return endings + stims


def someStims():
    stims = ["001100011100", "100111000011", "111110000111", "011110011100","010010001000"]

    #endings = ["111111111110", "010010001000", "110110110010", "011000111100"]
    random.shuffle(stims)
    return stims

def someStims2():
    stims = ["11111110"]

    #endings = ["111111111110", "010010001000", "110110110010", "011000111100"]
    random.shuffle(stims)
    return stims

def run(lst):
    if LOTlib.SIG_INTERRUPTED:
        return ""
    data = [FunctionData(input=(), output={lst: len(lst)})]
    h0 = MyHypothesis()
    tn = TopN(N=TOP)
    # run the sampler
    counter = Counter()
    for h in MHSampler(h0, data, steps=STEPS, acceptance_temperature=1.0, likelihood_temperature=1.0):#, likelihood_temperature=10.0):
        # counter[h] += 1
        tn.add(h)

    z = logsumexp([h.posterior_score for h in tn])
    sort_post_probs = [(h, exp(h.posterior_score - z)) for h in tn.get_all(sorted=True)][::-1]
    return sort_post_probs

def proportionalize(all_vals, totGen):
    z = 0
    props = [0.0 for _ in xrange(totGen)]
    length = len(all_vals.keys())
    for key in all_vals.keys():
        z += 1
        val = all_vals[key]
        for i in xrange(len(key)):
            v = int(key[i]) * val
            props[i] += v

    return props



def main():

    all_predictions = []
    strings = [""]

    #to_run = allStrs(strings=[""], length=0, maxLength=8)
    #to_run = to_run[:len(to_run)/2]
    to_run = someStims2()

    orig_terminals = store_orig_terminals()
    update_ints_grammar(maxInt = len(to_run[0]), integers=True, positives=False, weight=1.0)

    count = 0
    totGen = 8
    for r in to_run:
        t1 = time.time()
        res = run(r)
        all_vals = getReturnVals(r, res, lenth = len(r), numGen=totGen, numHypProd=5)
        #turn outputs, with likelihoods, into likelihood
        #in each position over entire output sequence
        props = proportionalize(all_vals, totGen)
        all_predictions.append((r, copy.copy(props)))

        if update_prior:
            update_grammar_rules(grammar, res, exclude = ["INT", "BOOL", "START"])
        if analogy:
            analogize(grammar, res, mass=0.9, mult=2)
            get_nTERM()
            print_terminals(orig_terminals, thresh=0.05)
        print time.time() - t1
        print
    output_multiple(all_predictions, totGen)



main()

#x.fullprint
