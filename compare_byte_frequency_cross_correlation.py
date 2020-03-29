import json
import os
import sys
import ast
import numpy as np
from scipy.stats import chisquare

from matplotlib import pyplot as plt


def compute_bytes_distribution(filename):
    frequencies = 256 * [0]
    with open(filename, "rb") as file:
        data = file.read()
    for byte in data:
        frequencies[int(byte)] += 1
    div = max(frequencies)
    if div == 0:
        return 256 * [0]
    for i, x in zip(range(256), frequencies):
        if frequencies[i] != 0:
            frequencies[i] = frequencies[i] / div
    return frequencies


def import_signatures(imported, filename):
    with open(filename) as file:
        frequency_signature = []
        frequency_strength = []
        char = file.read(1)
        frequency_signature.append(char)
        while char != "]":
            char = file.read(1)
            frequency_signature.append(char)
        frequency_signature = ast.literal_eval("".join(frequency_signature))
        file.readline()
        char = file.read(1)
        frequency_strength.append(char)
        while char != "]":
            char = file.read(1)
            frequency_strength.append(char)
        frequency_strength = ast.literal_eval("".join(frequency_strength))
        file.readline()
        cross_correlation = ast.literal_eval(file.read())
        imported.update({
            "frequencies": frequency_signature,
            "strengths": frequency_strength,
            "cross_correlations": cross_correlation
        })
        # count = 0
        # total = 0
        # for x in range(256):
        #     for y in range(x+1, 256):
        #         count += 1
        #         total += abs(cross_correlation[x][y])
        # print(total/count)
        # exit()


def import_classifications(classifications, filename, mode):
    with open(filename) as file:
        for line in file:
            line = line.split(" ")
            line[1] = line[1].replace("\n", "")
            if mode == "strict":
                if line[1] == "False":
                    classifications.update({line[0]: False})
                else:
                    classifications.update({line[0]: True})
            else:
                if line[1] == "False" or line[1] == "Maybe":
                    classifications.update({line[0]: False})
                else:
                    classifications.update({line[0]: True})


def compute_frequency_score(fingerprint, signature, strengths):
    tmp_f = []
    tmp_s = []
    score = 0
    for i in range(256):
        if strengths[i] > 0.0000001:
            tmp_f.append(fingerprint[i] * 10)
            tmp_s.append(signature[i] * 10)
    score = chisquare(np.array(tmp_f), f_exp=np.array(tmp_s))
    return score.pvalue


def compute_cross_correlation_score(fingerprint, signature):
    score = 0
    for i in range(256):
        for u in range(i+1, 256):
            score += abs(abs(signature[i][u]) - abs(fingerprint[i] - fingerprint[u]))
    return score


if __name__ == "__main__":
    imported = {}
    classifications = {}
    in_file = sys.argv[1]
    signatures = sys.argv[2]
    import_signatures(imported, signatures)
    if len(sys.argv) < 4:
        frequency = compute_bytes_distribution(in_file)
        # score_frequency = compute_frequency_score(frequency, imported["frequencies"], imported["strengths"])
        score_cross_correlation = compute_cross_correlation_score(frequency, imported["cross_correlations"])
    else:
        hits = 0
        misses = 0
        f_pos = 0
        f_neg = 0
        files = os.listdir(in_file)
        import_classifications(classifications, sys.argv[3], sys.argv[5])
        threshold = float(sys.argv[4])
        for file in files:
            frequency = compute_bytes_distribution(os.path.join(in_file, file))
            # score_frequency = compute_frequency_score(frequency, imported["frequencies"], imported["strengths"])
            score_cross_correlation = compute_cross_correlation_score(frequency, imported["cross_correlations"])
            #print(score_cross_correlation)
            if score_cross_correlation < threshold:
                result = False
            else:
                result = True
            if result == classifications[file]:
                hits += 1
            else:
                misses += 1
                if result:
                    f_pos += 1
                else:
                    f_neg += 1
        print("Hits: {}, Misses: {}. False positives: {}, False Negatives: {}".format(str(100 / len(files) * hits),
                                                                                      str(100 / len(files) * misses),
                                                                                      str(100 / len(files) * f_pos),
                                                                                      str(100 / len(files) * f_neg)))
