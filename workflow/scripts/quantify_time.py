import numpy as np

def quantify_time(alpha, mu, rho, seq_speed, unblock_path):
    fh = open(unblock_path, 'rt')
    # loop over all reads in the fastq file
    for desc, name, seq, qual in readfq(fh):
        total += 1
        base_total +=len(seq)
        if len(seq) == mu:
            unb += 1
    fh.close()

    time_control += total_bases
    time_control += (batchsize * alpha)
    # BR: all accepted bases and ((mu + rho) * number of unmapped/rejected reads) + (alpha * batch_size)
    bases_br = np.sum([r[0].qlen for r in paf_dict.values()])
    time_boss += bases_br
    time_boss += (n_unmapped * (mu + rho))
    time_boss += (n_reject * (mu + rho))
    time_boss += (batchsize * alpha)


if __name__ == "__main__":
    data = quantify_time(fq=sys.argv[1], mu=int(sys.argv[2]))
    print(data)