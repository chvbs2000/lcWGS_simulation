import subprocess
import argparse
import multiprocessing as mp
import numpy as np
import random
import gzip
import allel
from util import *
from decimal import Decimal
import pandas as pd
import os
from sklearn import mixture
from datetime import date
import sklearn

def print_versions():
    print(f"subprocess: {subprocess.__version__ if hasattr(subprocess, '__version__') else 'N/A'}")
    print(f"argparse: {argparse.__version__ if hasattr(argparse, '__version__') else 'N/A'}")
    print(f"multiprocessing: {mp.__version__ if hasattr(mp, '__version__') else 'N/A'}")
    print(f"numpy: {np.__version__}")
    print(f"random: {random.__version__ if hasattr(random, '__version__') else 'N/A'}")
    print(f"gzip: {gzip.__version__ if hasattr(gzip, '__version__') else 'N/A'}")
    print(f"allel: {allel.__version__}")
    print(f"decimal: {Decimal.__version__ if hasattr(Decimal, '__version__') else 'N/A'}")
    print(f"pandas: {pd.__version__}")
    print(f"os: {os.__version__ if hasattr(os, '__version__') else 'N/A'}")
    print(f"sklearn: {sklearn.__version__}")
    print(f"datetime: {date.__version__ if hasattr(date, '__version__') else 'N/A'}")

def flatten(mydata):
    # subjects, SNP, REF/ALT counts
    if len(mydata.shape) == 3:
        mydata = np.reshape(mydata, (mydata.shape[0], -1))
    else:  # do one hot encoding, depth=3 because missing (-1) is encoded to all zeroes
        mydata = tf.one_hot(indices=mydata, depth=3)
        mydata = tf.layers.flatten(mydata)  # flattening to simplify calculations later (matmul, add, etc)
    return mydata

def convert_GL2PL(GL_arr):
    def get_PL(n):
        if n < 0.000000001:
            n = 0.000000001
        dn = Decimal(abs(n))
        log_n = (-10)*(dn.log10())
        #log_n = np.log10(n)
        return log_n
    
    raw_PL_arr = list(map(get_PL, GL_arr))
    nor_PL_arr = [x - min(raw_PL_arr) for x in raw_PL_arr]
    nor_PL_arr = list(map(int, nor_PL_arr))
    return nor_PL_arr

def simulate_lowpass_genotype(y_pred_sample, j, is_debug_on):
    # Homo Ref = Ref * (1-ALT)
    # Het = Ref * Alt
    # Homo Alt = (1-Ref) * Alt
    Pa = y_pred_sample[j] * (1 - y_pred_sample[j + 1])
    Pab = y_pred_sample[j] * y_pred_sample[j + 1]
    Pb = (1 - y_pred_sample[j]) * y_pred_sample[j + 1]
    P = [Pa, Pab, Pb]  
    
    choices = ['missing','non-missing']
    weights = [0.36, 0.64]
    outcome = random.choices(choices, weights=weights, k=1)[0]
    
    if outcome == "missing":
        gen = "./."
        return gen
    
    elif outcome == "non-missing":
        # 2 versions of dosage
        # Df0=(1-Pa)+Pb
        PL = convert_GL2PL(P)
        Df1 = Pab + (2 * Pb)
        if Df1 != 0:
            Psum = Pa + Pab + Pb
            Df1 = Df1 / Psum
        D = np.argmax(P)
        Df1 = np.round(Df1, 4)

        if D == 2:
            gen = "1/1:"
        elif D == 0:
            gen = "0/0:"
        elif D == 1:
            gen = "0/1:"
        else:
            print("ERROR dosage = %d  %f%" % (D, P))
            return None
        gen += str(Df1) + ":" + str(PL[0]) + "," + str(PL[1]) + "," + str(PL[2])
        if is_debug_on is True:
            gen += ":" + str(y_pred_sample[j]) + ":" + str(y_pred_sample[j + 1])
        return gen
    
    else:
        print ("ERROR in generating missingness: %s" % (outcome))
        return None

def simulate_gen_low(y_pred_sample, is_debug_on):
    j = 0
    ret_out = []
    while j < len(y_pred_sample):
        gen = None
        gen = simulate_lowpass_genotype(y_pred_sample, j, is_debug_on)
        j += 2
        if gen is None:
            print("ERROR, genotype NULL. y_pred[i][j] %f" % y_pred_sample[j])
        ret_out.append(gen)
    return ret_out

def get_header_line(input_file):
    if(str(input_file).endswith('.gz')):
        with gzip.open(input_file,'rt') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    return line
    else:
        with open(input_file) as f:
            for line in f:
                if line.startswith("#CHROM"):
                    return line
    return None

def export_vcf(posfile, pred_out, infile, output_file):
    my_header = get_header_line(infile)
    if not my_header:
        raise Exception("Could not extract header line from genotype_array file")

    refpos = pd.read_csv(posfile, sep='\t', comment='#', header=None)

    pred_out = np.transpose(pred_out)
    refpos = np.asarray(refpos.values)
        #"##FORMAT=<ID=Paa,Number=1,Type=Float,Description=\"Imputation probability for homozigous reference : Pa=y_pred[i][j]*(1-y_pred[i][j+1])\">\n" + \
        #"##FORMAT=<ID=Pab,Number=1,Type=Float,Description=\"Imputation probability for heterozigous : Pab=y_pred[i][j]*y_pred[i][j+1]\">\n" + \
        #"##FORMAT=<ID=Pbb,Number=1,Type=Float,Description=\"Imputation probability for homozigous alternative : Pb=(1-y_pred[i][j])*y_pred[i][j+1]\">\n" + \
        #"##FORMAT=<ID=AP,Number=1,Type=Float,Description=\"Predicted presence of reference allele (autoencoder raw output)\">\n" + \
        #"##FORMAT=<ID=BP,Number=1,Type=Float,Description=\"Predicted presence of alternative allele (autoencoder raw output)\">\n" + \
    comments = \
        "##fileformat=VCFv4.1\n##filedate=" + str(date.today()) + \
        "\n##source=Imputation_autoencoder\n##contig=<ID=" + str(refpos[0][0])+">\n" + \
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" + \
        "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">\n" + \
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n" + \
        my_header

    columns = ['.', '.', '.', 'GT:DS:PL']*len(refpos)
    columns = np.reshape(columns, (len(refpos), 4))

    refpos = np.concatenate((refpos, columns), axis=1)
    vcf = np.concatenate((refpos, pred_out), axis=1)

    with open(output_file, 'w') as f:
        f.write(comments)  # python will convert \n to os.linesep
    with open(output_file, 'ab') as bf:  # open as binary
        np.savetxt(bf, vcf, delimiter="\t", fmt='%s')
    
class simulator:
    def __init__(self):
        """
        reload pre-trained 1x gmm model for each genotype
        """
        #model directory
        self.model_dir = os.path.join(os.getcwd(), "trained_gmm")
        
        # homozygous reference model
        self.means_homoref = np.load(self.model_dir + "/" + "gmm_1x_homozygous_reference_means.npy")
        self.covar_homoref = np.load(self.model_dir + "/" + "gmm_1x_homozygous_reference_covariances.npy")
        self.gmm_homoref = mixture.GaussianMixture(n_components = len(self.means_homoref), covariance_type='full')
        self.gmm_homoref.precisions_cholesky_ = np.linalg.cholesky(np.linalg.inv(self.covar_homoref))
        self.gmm_homoref.weights_ = np.load(self.model_dir + "/" + "gmm_1x_homozygous_reference_weights.npy")
        self.gmm_homoref.means_ = self.means_homoref
        self.gmm_homoref.covariances_ = self.covar_homoref

        # heterozygous model
        self.means_het = np.load(self.model_dir + "/" + "gmm_1x_heterozygous_means.npy")
        self.covar_het = np.load(self.model_dir + "/" + "gmm_1x_heterozygous_covariances.npy")
        self.gmm_het = mixture.GaussianMixture(n_components = len(self.means_het), covariance_type='full')
        self.gmm_het.precisions_cholesky_ = np.linalg.cholesky(np.linalg.inv(self.covar_het))
        self.gmm_het.weights_ = np.load(self.model_dir + "/" + "gmm_1x_heterozygous_weights.npy")
        self.gmm_het.means_ = self.means_het
        self.gmm_het.covariances_ = self.covar_het
        
        # homozygous alternative model
        self.means_homoalt = np.load(self.model_dir + "/" + "gmm_1x_homozygous_alternative_means.npy")
        self.covar_homoalt = np.load(self.model_dir + "/" + "gmm_1x_homozygous_alternative_covariances.npy")
        self.gmm_homoalt = mixture.GaussianMixture(n_components = len(self.means_homoalt), covariance_type='full')
        self.gmm_homoalt.precisions_cholesky_ = np.linalg.cholesky(np.linalg.inv(self.covar_homoalt))
        self.gmm_homoalt.weights_ = np.load(self.model_dir + "/" + "gmm_1x_homozygous_alternative_weights.npy")
        self.gmm_homoalt.means_ = self.means_homoalt
        self.gmm_homoalt.covariances_ = self.covar_homoalt
    
    def chunks(self, n_data, n_threads):
        for i in range(0, len(n_data), n_threads):
            yield n_data[i:i + n_threads]
    
    def multithread_sampling(self, gt_arr):
        gt_flat = gt_arr.flatten()
        nrow, ncol = gt_arr.shape
        lowpass_gl = np.empty([len(gt_flat),2])
        homoref_idx = np.where(gt_flat == 0)[0]
        het_idx = np.where(gt_flat == 1)[0]
        homoalt_idx = np.where(gt_flat == 2)[0]
        sample_homoref = self.gmm_homoref.sample(len(homoref_idx))
        sample_het = self.gmm_het.sample(len(het_idx))
        sample_homoalt = self.gmm_homoalt.sample(len(homoalt_idx))
        lowpass_gl[homoref_idx] = sample_homoref[0]
        lowpass_gl[het_idx] = sample_het[0]
        lowpass_gl[homoalt_idx] = sample_homoalt[0]
        lowpass_gl = lowpass_gl.reshape((nrow,ncol,2))
        return lowpass_gl

    def simulate_allele(self, highpass_arr): 
        nproc = mp.cpu_count()
        pool = mp.pool.ThreadPool(nproc)
        results = pool.map(self.multithread_sampling, self.chunks(highpass_arr, nproc))
        pool.close()
        pool.join()
        results = [val for sublist in results for val in sublist]
        lowpass_gl = np.array(results)
        return lowpass_gl
    
if __name__ == "__main__":
    print_versions()

    parser = argparse.ArgumentParser(description='Simulate a VCF file.')
    parser.add_argument('--input_vcf', help='Input VCF file path')
    parser.add_argument('--output_dir', help='Output VCF file directory')
    args = parser.parse_args()
    
    if not args.input_vcf:
        parser.error('no input file found!')

    input_vcf = args.input_vcf
    pos_name = os.path.basename(input_vcf)
    file_name = os.path.basename(input_vcf).split(".")[:-2]
    
    tbi_dir = "{}.tbi".format(input_vcf)
    posfile = "{}/{}.position".format(args.output_dir,pos_name)
    print("posfile:", posfile)
    output_file = os.path.join(args.output_dir, ".".join(file_name) + ".simulated.vcf")

    extract_pos_cmd = "zcat {} | grep -v '#' | cut -f1-5 > {}/{}.position".format(input_vcf, args.output_dir, pos_name)
    result_extract = subprocess.run(extract_pos_cmd, shell=True, capture_output=True, text=True)

    if result_extract.returncode != 0:
        print(f"Error: {result_extract.stderr}")
    else:
        print(f"Positions extracted from {input_vcf} and written to {posfile}")

    low_pass_data = allel.read_vcf(input_vcf, tabix=tbi_dir, fields=['calldata/GT'])
    low_pass_gt = low_pass_data['calldata/GT']
    low_pass_gt = allel.GenotypeArray(low_pass_gt,dtype='i1')
    low_pass_gt = low_pass_gt.to_n_alt(fill=-1)
    data=low_pass_gt.T
    lowpass_simulator = simulator()
    simulated_data = lowpass_simulator.simulate_allele(data)
    xs = flatten_data(simulated_data.copy())
    pool_low = mp.Pool(mp.cpu_count())
    out_low = pool_low.starmap(
        simulate_gen_low,
        [(xs[i], False) for i in range(len(xs))]
    )
    pool_low.close()

    export_vcf(posfile, out_low, input_vcf, output_file)

    bgzip_command = "bgzip -c {} > {}.gz; tabix -p vcf {}.gz".format(output_file,output_file,output_file)
    result_bgzip = subprocess.run(bgzip_command, shell=True, capture_output=True, text=True)

    print("RESULT: {}.gz".format(output_file))

    delete_command = "rm {}/{}.position; rm {}".format(args.output_dir, pos_name, output_file)
    result_delete = subprocess.run(delete_command, check=True, shell=True, capture_output=True, text=True)
