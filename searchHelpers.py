
FLIP_STR = '1*flip_' # the name of "flip" in our grammar

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define the function we will map
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import LOTlib
from model_analogy import *
from LOTlib.DataAndObjects import FunctionData
from LOTlib.Inference.Samplers.MetropolisHastings import MHSampler
from LOTlib.Hypotheses.Hypothesis import Hypothesis
from LOTlib.TopN import TopN
from collections import Counter
from LOTlib.Miscellaneous import logsumexp, qq
from math import exp, log
import copy
from LOTlib.GrammarInference.Precompute import create_counts
import numpy as np

def get_nTERM():
    rules = get_rule_sigs()
    count = len(rules["TERM"])
    return float(count)

def get_rule_sigs():
    rule_sigs = {}
    for i in grammar.get_all_rules():
        sigI = i.get_rule_signature()
        if sigI[0] not in rule_sigs:
            rule_sigs[sigI[0]] = []
        rule_sigs[sigI[0]].append((sigI[1]))
    return rule_sigs


def update_ints_grammar(maxInt, integers=True, positives=False, weight=1.0):
    #maxInt = NDATA
    rule_sigs = get_rule_sigs()
    # grammar.get_all_rules()
    if integers:
        for i in xrange(0,maxInt):
            poss = str(i)
            if ("INT" not in rule_sigs) or (poss not in rule_sigs["INT"]):
                if ("INT" not in rule_sigs):
                    probUse = weight
                else:
                    probUse = weight/float(len(rule_sigs["INT"]))
                grammar.add_rule('INT', str(i), None, probUse)

    grammar.renormalize()

def invert(s):
    inverse = ""
    for j in s:
        if j == "0":
            l = "1"
        else:
            l = "0"
        inverse += l
    return inverse


def toAB(s):
    ret = ""
    for j in s:
        if j == "0":
            l = "a"
        else:
            l = "b"
        ret += l
    return ret

def allStrs(strings=[""],length=0, maxLength=8):
     if length == maxLength:
         return strings
     else:
         copyStrings = copy.copy(strings)
         return (allStrs([string + "0" for string in copyStrings] + [string + "1" for string in copyStrings], length + 1, maxLength))


def update_grammar_rules(grammar, res, delete=False, threshold=0.01, exclude = []):
    hyps = [res[h][0] for h in xrange(len(res))]
    probs = [res[h][1] for h in xrange(len(res))]

    cc = create_counts(grammar, hyps)
    # print cc[1]
    sigs = [i.get_rule_signature() for i in grammar]
    dif_sigs = []
    for s in sigs:
        if s[0] not in dif_sigs:
            dif_sigs.append(s[0])


    for g in grammar:
        sig = g.get_rule_signature()
        if sig[0] not in exclude:
            indx = cc[1][sig]
            g.p = 0
            for h in xrange(len(hyps)):
                rule_count = cc[0][sig[0]][h][indx]
                hyp = hyps[h]
                prob = probs[h]
                cplx = np.sum([np.sum(cc[0][s][h]) for s in dif_sigs])
                new_prior = (1.0 + rule_count) / (grammar.nrules() + cplx)
                g.p += new_prior * prob
            if delete and g.p < threshold:
                deleteRule(hyp)

    grammar.renormalize()

def analogize(grammar, res, mass=0.2, mult=1):
    totMass = 0.0
    hyps = [res[h][0] for h in xrange(len(res))]
    probs = [res[h][1] for h in xrange(len(res))]
    h = 0
    while totMass < mass:
        prob = res[h][1]
        totMass += prob
        rule_sigs = get_rule_sigs()
        hyp = res[h][0]
        s_hyp = str(hyp)
        s_hyp = s_hyp[len("lambda lst, index: "):]

        if (s_hyp not in rule_sigs["TERM"]):
            num_terms = len(rule_sigs["TERM"])

            #print s_hyp[len("lambda lst, index: "):]
            grammar.add_rule('TERM', s_hyp, None, mult * prob/float(num_terms))

        h += 1

    grammar.renormalize()


def output(all_predictions, length):
    strs = ""
    for i in all_predictions:
        str1 = ""
        for j in i[0]:
            if j == "0":
                str1 += "a"
            else:
                str1 += "b"
        strs += str1 + "," + str(i[1]) + "," + str(i[2])  + "\n"


    out = open("lotOut-posterior-%d.csv" % length, "w+")
    out.write(strs)
    out.close()


def output_multiple(all_predictions, length, file="lotOut-posterior-all.csv"):
    strs = ""
    for pred in all_predictions:
        inp = pred[0]
        probs = pred[1]
        strs += "%s, " % toAB(inp)
        strs +=  ("%f," * length) % tuple(pred[1])
        strs += "\n"



    out = open(file, "w+")
    out.write(strs)
    out.close()


def store_orig_terminals():
    orig_terminals = []
    for i in grammar:
        rule_sig = tuple(i.get_rule_signature())
        orig_terminals.append(rule_sig)
    return orig_terminals


def print_terminals(orig_terminals, typeT="TERM", thresh=0.0):

    for i in grammar:
        rule_sig = tuple(i.get_rule_signature())
        if typeT in rule_sig[0] and (rule_sig not in orig_terminals) and (i.p > thresh):
            print rule_sig[1], i.p


def getBinProb(p, n, allVals, prb):
    for out_dem in allStrs([""], 0, n):
        k = out_dem.count("1")
        if out_dem not in allVals:
            allVals[out_dem] = 0.0
        pNum = (p**k) * (1 - p)**(n - k)
        allVals[out_dem] += pNum * prb
        #allVals[out_dem] = pNum
    return allVals



def getFlipOut(res, hyp, prob, allVals, numGen, seq1, lenth):
    s_hyp = str(hyp)
    flip = s_hyp.find("flip_")
    prob_flip = 0.0
    if len(s_hyp) < len("lambda lst, index: 1*flip_(0.1)") + 2:
        prob_flip = s_hyp[flip + len("flip_"):]
        end = s_hyp[flip + len("flip_"):].find(")")
        prob_flip = prob_flip[1:end]
        out_dem = float(prob_flip)
        allVals = getBinProb(out_dem, numGen, allVals, prob)


    else:
        out_dem = 0.0
        sample = 50
        for _ in xrange(sample):

            out_dem = hyp.__call__(noise=0.0, start=lenth, max_length=(len(seq1) + numGen), origSeq=copy.copy(seq1))
            out_dem = out_dem[len(out_dem) - numGen:]
            if out_dem not in allVals:
                allVals[out_dem] = 0.0
            allVals[out_dem] += prob/float(sample)
    return allVals


def getReturnVals(r, res, lenth, numGen=1, numHypProd=5):
    all_pred = 0.0
    all_prob = 0.0
    allVals = dict()
    seq1 = [0 for _ in xrange(lenth)]
    for i in r:
        seq1.append(int(i))
    for i in xrange(len(res)):
        totP = 0.0
        hyp = res[i][0]
        prob = res[i][1]
        s_hyp = str(hyp)
        flipFind = s_hyp.find("flip_")
        if flipFind != -1:
            flip = 1.0
            allVals = getFlipOut(res, hyp, prob, allVals, numGen, seq1, lenth)

        else:
            out_dem = hyp.__call__(noise=0.0, start = lenth,
                                   max_length = len(seq1) + numGen, origSeq=copy.copy(seq1))
            out_dem = out_dem[len(out_dem) - numGen:]

            if out_dem not in allVals:
                allVals[out_dem] = 0.0
            allVals[out_dem] += prob
            flip = 0.0
        #totP += float(out_dem)
        if i < numHypProd:
            z2 = hyp.__call__(noise=0.0, start = lenth,
                                   max_length = len(seq1) + numGen, origSeq=copy.copy(seq1))
            print hyp, z2[len(z2) - numGen:], prob

    vs = sorted([(i, allVals[i]) for i in allVals.keys()], key= lambda tup: 1 - tup[1])
    sumVals = sum([allVals[i] for i in allVals.keys()])
    minVal = min([allVals[i] for i in allVals.keys()])
    print r
    for i in xrange(min(len(vs), 10)):
        print vs[i][0], vs[i][1]/float(sumVals)

    return allVals
