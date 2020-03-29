import os
import sys
import multiprocessing as mp

from matplotlib import pyplot as plt


def compute_bytes_distribution(filename):
    frequencies = 256 * [0]
    with open(filename, "rb") as file:
        data = file.read()
    for byte in data:
        frequencies[int(byte)] += 1
    div = max(frequencies)
    if div == 0:
        return None
    for i, x in zip(range(256), frequencies):
        if frequencies[i] != 0:
            frequencies[i] = frequencies[i] / div
    return frequencies


def combine(old, new, cnt, size=1):
    return [w / (cnt + size) for w in [sum(y) for y in zip([cnt * x for x in old], [size * z for z in new])]]


def compute_strengths(sums, new, cnt):
    e = 2.71828
    sigma = 0.0375
    tmp = []
    for s, n in zip(sums, new):
        x = n - s
        tmp.append(e ** ((-x ** 2) / (2 * (sigma ** 2))))
    return tmp


def thread_function(mappings, index, queue):
    # Create table for cross correlations
    cross_correlations = [256 * [0] for _ in range(256)]

    # Compute frequencies over first file. Needed for frequency score function and correlation
    old_frequencies = compute_bytes_distribution(mappings[str(index)]["files"][0])
    for x in range(256):
        for y in range(1 + x, 256):
            cross_correlations[x][y] = abs(old_frequencies[x] - old_frequencies[y])
    sums = old_frequencies
    cnt = 1

    # Compute frequency over second file. Needed for correlation
    new_frequencies = compute_bytes_distribution(mappings[str(index)]["files"][1])
    old_correlations = compute_strengths(sums, new_frequencies, cnt)
    for i, y in zip(range(256), new_frequencies):
        sums[i] += y
    for x in range(256):
        for y in range(1 + x, 256):
            cross_correlations[x][y] = abs(new_frequencies[x] - new_frequencies[y])
        if x != 255:
            ind = x + 1
            cross_correlations[x][ind:] = combine(cross_correlations[x][ind:], [r[x] for r in cross_correlations[ind:]],
                                                  cnt)

    # Combine frequencies as needed
    old_frequencies = combine(old_frequencies, new_frequencies, cnt)
    cnt += 1

    # Iterate over files
    for file in mappings[str(index)]["files"][2:]:
        new_frequencies = compute_bytes_distribution(file)
        if new_frequencies:
            new_correlations = compute_strengths(sums, new_frequencies, cnt)
            old_correlations = combine(old_correlations, new_correlations, cnt)
            for i in range(len(new_frequencies)):
                sums[i] = (sums[i] + new_frequencies[i])
            old_frequencies = combine(old_frequencies, new_frequencies, cnt)
            for x in range(256):
                for y in range(1 + x, 256):
                    cross_correlations[x][y] = abs(new_frequencies[x] - new_frequencies[y])
                    if x != 255:
                        ind = x + 1
                        cross_correlations[x][ind:] = combine(cross_correlations[x][ind:],
                                                              [r[x] for r in cross_correlations[ind:]],
                                                              cnt)
            cnt += 1

    # Write output to multiprocessing queue
    queue.put((index, old_frequencies, old_correlations, cross_correlations))
    return


def compute_filetype_signature(directory, threads):
    ts = []
    queue = mp.Queue()
    files = [os.path.join(directory, file) for file in os.listdir(directory)]
    mappings = {}
    for x in range(threads):
        mappings[str(x)] = {"result": None, "strength": None, "files": []}
    size = len(files) // threads
    for thread in range(threads):
        mappings[str(thread)]["files"] = files[:size]
        files = files[:size]
    mappings[str(threads - 1)]["files"] = mappings[str(threads - 1)]["files"] + files
    for i in range(threads):
        t = mp.Process(target=thread_function, args=(mappings, i, queue))
        ts.append(t)
        t.start()
    for t in ts:
        t.join()
    if threads == 1:
        _, freq, strn, corr = queue.get()
    else:
        frequencies = []
        correlations = []
        cross_correlations = []
        for _ in range(threads):
            index, f, c, cc = queue.get()
            frequencies.append((index, f))
            correlations.append((index, c))
            cross_correlations.append((index, cc))
        frequencies.sort()
        correlations.sort()
        cross_correlations.sort()
        frequencies = [x for _, x in frequencies]
        correlations = [x for _, x in correlations]
        cross_correlations = [x for _, x in cross_correlations]
        freq = frequencies[0]
        strn = correlations[0]
        corr = cross_correlations[0]
        cont = len(mappings["0"]["files"])
        for x in range(1, threads):
            freq = combine(freq, frequencies[x], cont, len(mappings[str(x)]["files"]))
            strn = combine(strn, correlations[x], cont, len(mappings[str(x)]["files"]))
            for r in range(255):
                ind = r + 1
                corr[r][ind:] = combine(corr[r][ind:], cross_correlations[x][r][ind:], cont,
                                        len(mappings[str(x)]["files"]))
            cont += len(mappings[str(x)]["files"])
    return freq, strn, corr


if __name__ == "__main__":
    if sys.argv[1] == "single_file":
        filename = sys.argv[2]
        frequencies = compute_bytes_distribution(filename)
    elif sys.argv[1] == "directory":
        directory = sys.argv[2]
        threads = int(sys.argv[3])
        frequencies, strengths, cross_correlations = compute_filetype_signature(directory, threads)
        with open(sys.argv[4], "w") as out_file:
            out_file.write(str(frequencies))
            out_file.write("\n")
            out_file.write(str(strengths))
            out_file.write("\n")
            out_file.write(str(cross_correlations))
