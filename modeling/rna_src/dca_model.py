"""DCA model/Sampling/Statistics
"""

from sklearn.linear_model import Lasso
from numpy.random import choice, uniform
from numpy import exp, array, matmul, zeros, log
from numpy import sum as npsum
from copy import deepcopy
from itertools import product

from .utils import encode_conf_single, encode_conf_full, read_data,\
    read_data_mat, print_params, encode_conf_triplet

import argparse

ENCODE = encode_conf_full
TYPES = [0, 1, 2]


def init_entropy(h_terms, exp_prob, size):
    "The model diagonal parameters can be initialized using entropy positions."
    for i in range(size*3):
        if exp_prob[i] > 0:
            h_terms[i] = -log(exp_prob[i])
        else:
            h_terms[i] = 0.
    return h_terms


def mutate(conf):
    "Single segment is mutated"
    new_conf = [el for el in conf]
    pos = choice(range(len(conf)))
    new_el = choice(TYPES)
    while conf[pos] == new_el:
        new_el = choice(TYPES)
    new_conf[pos] = new_el
    return new_conf


def update_H_terms(H_terms, exp_prob, cconf, inc):
    "Update the model parameters h = (mc - obs) * step"
    H_terms += (cconf - exp_prob) * inc
    return H_terms


def sample_dca(conf, H_terms, nb_steps, KT=1.0, conf_c=True):
    "Simple MC sampling for segments."
    traj = []
    cconf = array(ENCODE(conf))

    for _ in range(nb_steps):
        nconf = mutate(conf)
        cnconf = array(ENCODE(nconf))

        delta_nrj = sum((cnconf - cconf) * H_terms)

        if delta_nrj < 0 or exp(-delta_nrj/KT) >= uniform(0, 1):
            conf, cconf = nconf, cnconf
        if conf_c:
            traj += [cconf]
        else:
            traj += [conf]

    return cconf, traj


def bl_dca(exp_prob, all_samp, nb_steps=10, nb_mc=1000):
    "Boltzmann learning procedure"
    size = len(all_samp[0])
    h_terms = zeros(shape=int(size * 3 + (size * (size - 1)/2) * 9))
    h_terms = init_entropy(h_terms, exp_prob, size)

    inc = 0.1
    limit = 10**-2
    beta_l = float(limit/inc)**(1.0/float(nb_steps))
    traj_parms = []
    KT = 1.0

    for _ in range(nb_steps):
        conf = all_samp[choice(range(len(all_samp)))]
        _, encoded_traj = sample_dca(conf, h_terms, nb_mc, KT)

        pred_prob = npsum(encoded_traj, axis=0)/nb_mc
        h_terms += (pred_prob - exp_prob) * inc
        traj_parms += [abs(deepcopy(pred_prob) - exp_prob)]
        inc *= beta_l

    return h_terms, array(traj_parms)


def fit_parms(all_samp, prob, nb_steps=500, nb_mc=500, temp=1.0):
    "returns parms, true_nrj, fit_nrj, exp_prob, fit_prob, couplings"

    encoded_samp = array([ENCODE(comb) for comb in all_samp])

    # for each configuration, compute the prob with the random model
    prob /= sum(prob)
    # compute fi fij
    exp_prob = npsum(array([p * el for p, el in zip(prob, encoded_samp)]), axis=0)

    # learn the parameters
    h_terms, traj_parms = bl_dca(exp_prob, all_samp, nb_steps, nb_mc)
    return h_terms


def fit_parms_lasso(all_samp, prob, reg=0.02, temp=1.0):
    "returns parms, true_nrj, fit_nrj, exp_prob, fit_prob, couplings"
    encoded_samp = array([ENCODE(comb) for comb in all_samp])

    # for each configuration, compute the prob with the random model
    prob /= sum(prob)
    lreg = Lasso(alpha=reg, fit_intercept=False, max_iter=1000)
    reg = lreg.fit(encoded_samp, -temp*log(prob/prob[0]))

    return reg.coef_


def parse_arguments():
    """Parsing command line
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input')
    parser.add_argument('-ns', '--nb_steps', type=int, default=500)
    parser.add_argument('-nm', '--nb_mc', type=int, default=1000)
    parser.add_argument('-pm', '--prod_mc', type=int, default=10000)
    parser.add_argument('-r', '--reg', type=float, default=0.0)
    parser.add_argument('-t', '--temp', type=float, default=1.0)
    parser.add_argument('-d', '--drop', type=int, default=0)
    parser.add_argument('-md', '--min_drop', type=int, default=0)
    parser.add_argument('-o', '--out_file')
    parser.add_argument('--sample_in', action="store_true")
    parser.add_argument('--lasso', action="store_true")
    parser.add_argument('--single', action="store_true")
    parser.add_argument('--samp_no_drops', action="store_true")
    parser.add_argument('--count', action="store_true")
    parser.add_argument('--theo', action="store_true")
    parser.add_argument('--triplets', action="store_true")
    parser.add_argument('--back')
    parser.add_argument('--no_h', action="store_true", help="don't use h terms from the energy")
    return parser.parse_args()


def main_dca():
    """Learn the parameters for the model
    """
    args = parse_arguments()
    global ENCODE
    # test using only the single segment model
    if args.single:
        ENCODE = encode_conf_single
    else:
        ENCODE = encode_conf_full

    # read the input configuration
    all_samp, prob = read_data(args.input, nb_drops=args.drop, min_nb_drops=args.min_drop)

    # can use to train the model using a regression approach
    if args.lasso:
        h_parms = fit_parms_lasso(all_samp, prob, args.reg, args.temp)
    else:
        h_parms = fit_parms(all_samp, prob, args.nb_steps, args.nb_mc, args.temp)

    # print the trained parameters in the standard output
    print_params(h_parms)


def main_sample():
    # sample configuration usin
    args = parse_arguments()
    h_terms_val = read_data_mat(args.input, args.no_h)

    size = 8
    h_terms = zeros(shape=int(size * 3 + (size * (size - 1)/2) * 9))
    if args.samp_no_drops:
        global TYPES
        TYPES = [0, 1]

    positions = list(range(size))
    pid = 0
    for i in positions:
        for ii in [0, 1, 2]:
            h_terms[pid] = h_terms_val[(i, i)][(ii, ii)]
            pid += 1

    for i in positions:
        for j in positions[i+1:]:
            for ii in [0, 1, 2]:
                for jj in [0, 1, 2]:
                    h_terms[pid] = h_terms_val[(i, j)][(ii, jj)]
                    pid += 1

    KT = args.temp
    conf = [choice([0, 1, 2]) for i in range(8)]
    if args.theo:
        all_samp = [list(comp) for comp in product(*([TYPES]*8))]
        encoded_samp = [ENCODE(comp) for comp in all_samp]
        fit_nrj = matmul(encoded_samp, h_terms)
        pred_prob = exp(-fit_nrj)/sum(exp(-fit_nrj))

        all_el = list(zip(all_samp, pred_prob))
        all_el.sort(key=lambda el: -el[1])

        for conf, pp in all_el:
            if args.back is None:
                print("{} {}".format(" ".join([str(el) for el in conf]), pp))
            elif args.back == "H1N1" and conf.count(0) > 4:
                print("{} {}".format(" ".join([str(el) for el in conf]), pp))
            elif args.back == "H3N2" and conf.count(1) > 4:
                print("{} {}".format(" ".join([str(el) for el in conf]), pp))
    else:
        _, encoded_traj = sample_dca(conf, h_terms, args.prod_mc, KT, False)
        for conf in encoded_traj:
            print(" ".join([str(el) for el in conf]))


def main_stat():
    "compute pairwise/triplet statistics"
    args = parse_arguments()
    global ENCODE
    if args.triplets:
        ENCODE = encode_conf_triplet

    if args.sample_in:
        traj = []
        for line in open(args.input):
            val = line.strip().split()
            if val.count("2") <= args.drop:
                traj += [ENCODE([int(el) for el in val])]
        if args.count:
            pred_prob = npsum(array(traj), axis=0)
        else:
            pred_prob = npsum(array(traj), axis=0)/len(traj)
    elif args.theo:
        traj = []
        prob = []
        for line in open(args.input):
            el = line.strip().split()
            val, p = el[:8], float(el[8])

            if val.count("2") <= args.drop:
                traj += [ENCODE([int(el) for el in val])]
                prob += [p]
        prob = array(prob)
        prob /= sum(prob)
        pred_prob = npsum(array([array(el)*p for el, p in zip(traj, prob)]), axis=0)
    else:
        confs, prob, counts = read_data(args.input, nb_drops=args.drop, count=True)
        encoded_samp = []
        for el in confs:
            encoded_samp += [ENCODE(el)]
        if args.count:
            pred_prob = npsum(array([array(p) * el for p, el in zip(counts, encoded_samp)]), axis=0)
        else:
            pred_prob = npsum(array([array(p) * el for p, el in zip(prob, encoded_samp)]), axis=0)

    print_params(pred_prob, args.triplets)
