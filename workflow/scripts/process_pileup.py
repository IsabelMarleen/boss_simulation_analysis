import numpy as np
import scipy.integrate as integrate
import subprocess
import sys



def parse_pup(pup_path: str) -> np.array:
    '''
    Parse the pileup file to extract the coverage values as an array
    :param pup_path: path to the pileup file
    :return: array of coverage
    '''
    # extract 4th column only
    running2 = subprocess.Popen(
        args=f"cut -f2 {pup_path}",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf-8',
        shell=True
    )
    stdout2, stderr2 = running2.communicate()
    running4 = subprocess.Popen(
        args=f"cut -f4 {pup_path}",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf-8',
        shell=True
    )
    stdout4, stderr4 = running4.communicate()
    # parse to array of coverage
    pos = np.array([i for i in stdout2.split("\n") if i], dtype="int")
    coverage = np.array([i for i in stdout4.split("\n") if i], dtype="int")

    cov = np.array([pos, coverage], dtype="int")
    
    return cov


def calc_lowcov(cov: np.array, otu_size: int, lowcov_threshold: int = 5) -> float:
    '''
    Calculate the proportion of low coverage sites in the OTU.
    Attention the coverage array does not contain sites without any coverage

    :param cov: array of coverage counts
    :param otu_size: genome size
    :return: proportion of sites with low coverage
    '''
    lowcov_n = np.where(cov < lowcov_threshold)[0].shape[0]
    # to account for sites without coverage
    # calc difference in genome length and cov array
    nocov_n = otu_size - cov.shape[0]
    if nocov_n > 0:
        lowcov_n += nocov_n
    lowcov_p = lowcov_n / otu_size
    return lowcov_p

def calc_evenness(cov: np.array, otu_size: int) -> float:
    '''
    Calculate the evenness of coverage sites in the OTU
    as described in Mokry et al, 2010, Nucleic Acids Research.

    :param cov: array of coverage counts
    :param otu_size: genome size
    :return: evenness of coverage
    '''
    def _calc_evn(i, cov, otu_size) -> float:
        avg_cov = np.average(cov)
        cov_norm = cov/avg_cov
        p_i = np.count_nonzero(cov_norm >= (i/avg_cov))

        evn = p_i/otu_size
        return evn

    evenness = integrate.quad(_calc_evn, 0, 1, args=(cov, otu_size))
    
    return evenness[0]



def process_pup(pup: str, otu_size: int, pup_full: str) -> str:
    '''
    Process the pileup file to extract the mean coverage and proportion of low coverage sites
    :param pup: path of pileup file
    :param otu_size: genome size
    :return: csv data
    '''
    # get metadata
    meta = pup.split("/")[-1].split(".")[0].split('_')
    cond, time, otu = meta[0], meta[1], meta[2]
    # get an array of coverage from the pileup file
    cov = parse_pup(pup_path=pup)
    np.savetxt(pup_full, cov)
    # calculate the mean coverage of the OTU
    mean_cov = np.sum(cov[1]) / otu_size
    # calculate the proportion of low coverage sites
    lowcov_p = calc_lowcov(cov=cov[1], otu_size=otu_size)
    # calculate the evenness of coverage of the OTU
    evenness = calc_evenness(cov=cov[1], otu_size=otu_size)
    output = f'{cond},{time},{otu},{mean_cov},{lowcov_p},{evenness}'
    return output



if __name__ == "__main__":
    # input to this script is the output of samtools pileup
    # also pass in the genome size of the OTU to calculate the mean coverage
    data = process_pup(pup=sys.argv[1], otu_size=int(sys.argv[2]), pup_full = sys.argv[3])
    print(data)







