from LOTlib.Primitives import primitive
import math
import random
import copy
@primitive
def mod2_(x,y):
    if math.isnan(x) or math.isnan(y):
        return float("nan")
    elif y < 1:
        return 0
    else:
        return x % y

@primitive
def plusmod2_(x,y,z):
    if math.isnan(x) or math.isnan(y) or math.isnan(z):
        return float("nan")
    elif y < 1:
        return 0
    else:
        return (x + z) % y

@primitive
def timesmod2_(x,y,z):
    if math.isnan(x) or math.isnan(y) or math.isnan(z):
        return float("nan")
    elif y < 1:
        return 0
    else:
        return (x * z) % y

@primitive
def div2_(x,y):
    if math.isnan(x) or math.isnan(y):
        print x, y
        return float("nan")
    elif y < 1:
        return 0
    else:
        return x / y


@primitive
def not2_(x):
    if x != True and x != False:
        return float("nan")
    else:
        if x == False:
            return 1
        else:
            return 0

@primitive
def and2_(x, y):
    if x == 0 or y == 0:
        return 0
    else:
        return 1


@primitive
def or2_(x, y):
    if x == 1 or y == 1:
        return 1
    else:
        return 0


@primitive
def bool2term_(b):
    if b == True:
        return 1
    else:
        return 0

@primitive
def istrue_(x):
    if math.isnan(x):
        return float("nan")
    if x > 0:
        return 1
    else:
        return 0

@primitive
def flip2term_(b):
    return b



"""
    A simple example to explore sequence learning and the tradeoff between deterministic and stochastic hypotheses.
    For strings that have simple descriptions as the output of computations, we will find them.
    If not, we may throw up our hands and go with a stochastic definition.
"""

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define the grammar
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from LOTlib.Grammar import Grammar

grammar = Grammar()

grammar.add_rule('START', '', ['TERM'], 1.0)


grammar.add_rule('TERM', '1', None, 1.0)
grammar.add_rule('TERM', '0', None, 1.0)


grammar.add_rule('TERM', '1*flip_', ['PROB'], 3.0)
grammar.add_rule('TERM', 'if_',  ['BOOL', 'TERM', 'TERM'], 1.0)

grammar.add_rule('TERM', 'lst[-%s]', ['INT'], 1.0) # refer to our list argument
grammar.add_rule('TERM', 'not2_', ['BOOL'], 1.0)
grammar.add_rule('TERM', 'bool2term_', ['BOOL'], 1.0)
grammar.add_rule('TERM', 'istrue_', ['INT'], 1.0)

grammar.add_rule('BOOL', 'if_',  ['BOOL', 'BOOL', 'BOOL'], 0.25)
#grammar.add_rule('BOOL', '(%s == %s)',  ['TERM', 'INT'], 0.5)
#grammar.add_rule('BOOL', '(%s == %s)',  ['INT', 'INT'], 0.5)
grammar.add_rule('BOOL', '(%s == %s)',  ['index', 'INT'], 1.0)
grammar.add_rule('BOOL', '(%s == %s)',  ['INT', 'INT'], 1.0)
grammar.add_rule('BOOL', '%s > %s', ['index', 'INT'], 1.0)
grammar.add_rule('BOOL', '%s > %s', ['INT', 'index'], 1.0)
grammar.add_rule('BOOL', '%s > %s', ['INT', 'INT'], 1.0)

#grammar.add_rule('BOOL', 'flip_', ['PROB'], 1.0)
#grammar.add_rule('BOOL', '(%s == %s)',  ['index', 'INT'], 0.5)

grammar.add_rule('BOOL', 'or_',  ['BOOL', 'BOOL'], 0.5)
#grammar.add_rule('BOOL', 'xor_',  ['BOOL', 'BOOL'], 0.5)
grammar.add_rule('BOOL', 'and_', ['BOOL', 'BOOL'], 0.5)
grammar.add_rule('BOOL', 'not_', ['BOOL'], 1.0)


#grammar.add_rule('TERM', 'not2_', ['TERM'], 0.5)
grammar.add_rule('INT', 'lst[-%s]', ['INT'], 1.0) # refer to our list argument
grammar.add_rule('INT', 'mod2_', ['index', 'INT'], 1.0)
grammar.add_rule('INT', 'plusmod2_', ['index', 'INT', 'INT'], 1.0) # we can't use plus if it goes over the max number
#grammar.add_rule('INT', 'timesmod2_', ['index', 'INT', 'INT'], 0.5)# we can't use times if it goes over the max number
grammar.add_rule('INT', 'div2_', ['index', 'INT'], 1.0)



maxInt = 8
for p in xrange(1,maxInt):
    grammar.add_rule('PROB', str(p/float(maxInt)), None, 0.5)


grammar.renormalize()
#update_grammar = copy.deepcopy(grammar)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define the hypothesis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from LOTlib.Hypotheses.LOTHypothesis import LOTHypothesis
from LOTlib.Hypotheses.Likelihoods.NoisyStochasticLikelihood import NoisyStochasticLikelihood
#from LOTlib.Miscellaneous import attrmem, Infinity


#from LOTlib.Hypotheses.Priors.PCFGPrior import *

class MyHypothesis(NoisyStochasticLikelihood, LOTHypothesis):
    def __init__(self, **kwargs):
        LOTHypothesis.__init__(self, grammar=grammar,  maxnodes=25,
                               display="lambda lst, index: %s", prior_temperature=1.0,
                               likelihood_temperature=1.0, noise=0.1, **kwargs)

    def __call__(self, max_length, start, noise=0.1, origSeq = None):
        # calling a hypothesis will use my fvalue to "unwrap" a sequence
        # here the fvalue is something that will take a sequence and give back
        # the next elements
        #noise = 0.03
        seq = []
        if origSeq == None:
            seq = [0 for _ in xrange(start)] #randomize start
        else:
            seq = copy.copy(origSeq)
        retseq = copy.copy(seq)
        while len(seq) < max_length:
            #if max_length == start:
            ret = self.fvalue(seq, len(seq) - start)
            #else:
                #ret = self.fvalue(seq, len(seq) - (start + max_length)/2)
            seq.append(ret)
            if "flip" in str(self.value) and random.random() < noise:
                ret = (ret == 0) * 1
            retseq.append(ret)

        #print grammar.nrules()
        #rand = random.random()
        #print retseq
        #return ''.join(map(lambda x: str((rand >= noise) * str(x)) + str((rand < noise) * str((x == 0) * 1)), seq[start:]))
        return ''.join(map(str, retseq[start:]))





if __name__ == "__main__":

    from LOTlib.Miscellaneous import qq, display_option_summary

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process Options
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--steps", dest="STEPS", type="int", default=1000000, help="Number of samples to run")
    parser.add_option("--data", dest="DATA", type="str", default='1100110011', help="What data to run on?")
    parser.add_option("--ndata", dest="NDATA", type="int", default=10, help="How much do we have?")
    (options, args) = parser.parse_args()

    display_option_summary(options)

    ## FOR NOW -- CHECK SIZE
    assert len(options.DATA) == 10, "*** Everything is set up for fixed size=10!"

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Set up the data
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    from LOTlib.DataAndObjects import FunctionData

    data = [
        FunctionData(input=(), output={options.DATA: options.NDATA } )
    ]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Run
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    from LOTlib import break_ctrlc

    from LOTlib.Inference.Samplers.MetropolisHastings import MHSampler
    from LOTlib.SampleStream import *

    h0 = MyHypothesis()
    for h in SampleStream(break_ctrlc(MHSampler(h0, data, steps=options.STEPS))): # >> Unique():
        print h.posterior_score, h.prior, h.likelihood, h.value.contains_function("1*flip_"), qq(h)
        print h()

    # for _ in xrange(100):
    #     h = MyHypothesis()
    #     print h
    #     print h()
    #     print h.compute_likelihood(data)







