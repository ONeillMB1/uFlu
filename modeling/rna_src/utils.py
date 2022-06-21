from numpy import array
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import Normalize
from scipy.stats import chisquare, chi, fisher_exact, barnard_exact, binom_test
from scipy.stats import pearsonr


def read_seq_data(infile):
    traj = []
    for line in open(infile):
        val = line.strip().split()
        if val.count("2") <= 0:
            traj += [[int(el) for el in val]+[2]*20]
    return traj


def print_params(params, triplets=False):
    size = 8
    positions = range(size)
    pid = 0
    for i in positions:
        print(i, i, 0, 0, params[i*3])
        print(i, i, 1, 1, params[i*3+1])
        print(i, i, 2, 2, params[i*3+2])
        pid += 3

    for i in positions:
        for j in positions[i+1:]:
            for ii in [0, 1, 2]:
                for jj in [0, 1, 2]:
                    print(i, j, ii, jj, params[pid])
                    pid += 1

    if triplets:
        for i in positions:
            for j in positions[i+1:]:
                for k in positions[j+1:]:
                    for ii in [0, 1, 2]:
                        for jj in [0, 1, 2]:
                            for kk in [0, 1, 2]:
                                print(i, j, k, ii, jj, kk, params[pid])
                                pid += 1


def read_data_mat(infile, no_h=False):
    h_term = {}
    for l in open(infile):
        if len(l.strip().split()) == 5:
            pi, pj, ai, aj, v = l.strip().split()
            val = float(v)
            pair = (int(pi), int(pj))
            pair_a = (int(ai), int(aj))
        elif len(l.strip().split()) == 7:
            pi, pj, pk, ai, aj, ak, v = l.strip().split()
            val = float(v)
            pair = (int(pi), int(pj), int(pk))
            pair_a = (int(ai), int(aj), int(ak))

        if no_h and pi == pj:
            # pass if h terms
            pass

        if pair not in h_term:
            h_term[pair] = {pair_a: val}
        else:
            h_term[pair][pair_a] = val
    return h_term


def encode_conf_single(recombination):
    "Encode a configuration in [h0, h1..., hN, J00, J01,...JN-1N]"
    nb_el = len(recombination)
    conf = [0 for i in range(int(nb_el * 3 + (nb_el * (nb_el - 1)/2.0) * 9))]

    for i, pi in enumerate(recombination):
        conf[i*3+pi] = 1
    return conf


def encode_conf_triplet(recombination):
    "Encode a configuration in [h0, h1..., hN, J00, J01,...JN-1N]"
    nb_el = len(recombination)
    conf = [0 for i in range(int(nb_el * 3 + (nb_el * (nb_el - 1)/2.0) * 9 +
                                 (nb_el * (nb_el - 1)/2.0 * (nb_el - 2)/3.0) * 27))]

    pid, pit = 0, 0
    pairs = {}
    triplets = {}
    for ii in [0, 1, 2]:
        for jj in [0, 1, 2]:
            pairs[(ii, jj)] = pid
            pid += 1
            for kk in [0, 1, 2]:
                pairs[(ii, jj, kk)] = pid
                pit += 1

    pid = 0
    for i, pi in enumerate(recombination):
        conf[i*3+pi] = 1
        pid += 3

    for i, pi in enumerate(recombination):
        for j, pj in enumerate(recombination[i+1:]):
            # id_pair = (nb_el*3) + (i*(nb_el-1) + j - int(0.5*i*(i-1))) * 9
            conf[pid + pairs[(pi, pj)]] = 1
            pid += 9

    for i, pi in enumerate(recombination):
        for j, pj in enumerate(recombination[i+1:], start=i+1):
            for k, pk in enumerate(recombination[j+1:], start=j+1):
                conf[pid + pairs[(pi, pj, pk)]] = 1
                pid += 27
    return conf


def encode_conf_full(recombination):
    "Encode a configuration in [h0, h1..., hN, J00, J01,...JN-1N]"
    nb_el = len(recombination)
    conf = [0 for i in range(int(nb_el * 3 + (nb_el * (nb_el - 1)/2.0) * 9))]

    pid = 0
    pairs = {}
    for ii in [0, 1, 2]:
        for jj in [0, 1, 2]:
            pairs[(ii, jj)] = pid
            pid += 1

    pid = 0
    for i, pi in enumerate(recombination):
        conf[i*3+pi] = 1
        pid += 3

    for i, pi in enumerate(recombination):
        for j, pj in enumerate(recombination[i+1:]):
            conf[pid + pairs[(pi, pj)]] = 1
            pid += 9
    return conf


def read_data(infile, count=False, nb_drops=0, min_nb_drops=0, header=False):
    results, prob = [], []
    convert = {"H1N1": 0, "H3N2": 1, "NEITHER": 2}
    with open(infile) as input_dat:
        if header:
            input_dat.readline()    # read header
        for l in input_dat:
            val = l.strip().split()
            segments = val[1:9]
            ndrop = val.count("NEITHER")
            if ndrop <= nb_drops:
                prob += [int(val[9])]
                results += [[convert[el] for el in segments]]

    prob = array(prob)
    if count:
        return results, prob/sum(prob), prob
    else:
        return results, prob/sum(prob)


def read_theo(infile, nb_drops=0, min_nb_drops=0, header=False):
    traj = []
    prob = []
    with open(infile) as sample:
        if header:
            sample.readline()
        for line in sample:
            line_s = line.strip().split()
            val, p = [int(el) for el in line_s[:8]], float(line_s[8])
            if val.count(2) <= nb_drops and len(set(val)) > 1:
                traj += [val]
                prob += [p]
    prob = array(prob)
    return traj, prob/sum(prob)


def read_sample(infile, nb_drops=7):
    traj = []
    results_dic = {}
    for line in open(infile):
        val = line.strip().split()
        if val.count("2") <= nb_drops and len(set(val)) > 1:
            conf = tuple([int(el) for el in val])
            if conf in results_dic:
                results_dic[conf] += 1
            else:
                results_dic[conf] = 1
            traj += [[int(el) for el in val]]
    confs = [list(el) for el in results_dic]
    confs.sort(key=lambda el: -results_dic[tuple(el)])
    tot = sum(results_dic.values())
    prob = [float(results_dic[tuple(el)])/tot for el in confs]
    return confs, prob


def read_new_data(infile, count=False, header=False):
    results, prob, background = [], [], []
    convert = {"H1N1": 0, "H3N2": 1, "NEITHER": 2}
    with open(infile) as input_dat:
        if header:
            input_dat.readline()    # read header
        for l in input_dat:
            val = l.strip().split()
            segments = val[7:15]
            ndrop = val.count("NEITHER")
            prob += [int(val[15])]
            results += [[convert[el] for el in segments]]
            background += [val[4]]
    prob = array(prob)
    if count:
        return results, prob/sum(prob), background, prob
    else:
        return results, prob/sum(prob), background
