import numpy
import pandas
import subprocess
import os

def get_and_make_finemap_dir(working_dir):
    finemap_dir = working_dir + '/finemap_files'
    os.makedirs(finemap_dir, exist_ok=True)
    return finemap_dir

def write_zfiles(sumstat_dfs, finemap_dir):
    def get_zfile_path(clump):
        return finemap_dir + '/clump' + str(clump) + '.z'
    for i in range(len(sumstat_dfs)):
        sumstat_dfs[i].to_csv(get_zfile_path(i), sep=' ', index=False)

def write_ldfiles(ld_matrices, finemap_dir):
    def get_ldfile_path(clump):
        return finemap_dir + '/clump' + str(clump) + '.ld'
    for i in range(len(ld_matrices)):
        numpy.savetxt(get_ldfile_path(i), ld_matrices[i])

def write_masterfile(finemap_dir, num_clumps, num_samples):
    def get_masterfile_path():
        return finemap_dir + '/master'
    def get_masterfile_line(clump):
        extensions = ['.z', '.ld', '.snp', '.config', '.cred', '.log']
        line = ""
        for ext in extensions:
            line += finemap_dir + '/clump' + str(clump) + ext + ';'
        line += str(num_samples)
        return line
    header = "z;ld;snp;config;cred;log;n_samples"
    with open(get_masterfile_path(), 'w') as masterfile:
        masterfile.write(header + '\n')
        for i in range(num_clumps):
            masterfile.write(get_masterfile_line(i) + '\n')

def get_output_file_exts():
    return ['.snp', '.config']

def read_output_files(num_clumps, finemap_dir, exts):
    def get_output_file_path(clump, ext):
        return finemap_dir + '/clump' + str(clump) + ext
    return [[pandas.read_csv(get_output_file_path(i, ext), sep=' ') for i in range(num_clumps)] for ext in exts]

def read_log_files(num_clumps, finemap_dir):
    def get_log_file_path(clump):
        return finemap_dir + '/clump' + str(clump) + '.log_sss'
    logs = []
    for i in range(num_clumps):
        with open(get_log_file_path(i), 'r') as log_file:
            logs.append(log_file.read())
    return logs

def get_post_num_causal_snps(logs, max_causal_snps_per_clump):
    posts = []
    for log in logs:
        lines = log.split('\n')
        clump_posts = []
        for i in range(max_causal_snps_per_clump+1):
            clump_posts.append(float(lines[50+max_causal_snps_per_clump+i][8:].replace(')', '')))
        posts.append(clump_posts)
    return posts

def get_finemap_args(max_causal_snps_per_clump, prior_k0, num_threads, finemap_dir):
    return ['./finemap', '--sss', '--log', '--n-causal-snps', str(max_causal_snps_per_clump), '--prior-k0', str(prior_k0), '--n-threads', str(num_threads), '--in-files', finemap_dir + '/master']

def get_num_threads():
    return os.cpu_count()

def run_finemap(sumstat_dfs, ld_matrices, num_samples, max_causal_snps_per_clump, prior_k0, working_dir):
    finemap_dir = get_and_make_finemap_dir(working_dir)

    num_clumps = len(sumstat_dfs)
    write_masterfile(finemap_dir, num_clumps, num_samples)
    write_zfiles(sumstat_dfs, finemap_dir)
    write_ldfiles(ld_matrices, finemap_dir)

    subprocess.run(get_finemap_args(max_causal_snps_per_clump, prior_k0, get_num_threads(), finemap_dir))

    snp_dfs, config_dfs = read_output_files(num_clumps, finemap_dir, get_output_file_exts())
    logs = read_log_files(num_clumps, finemap_dir)
    posts = get_post_num_causal_snps(logs, max_causal_snps_per_clump)
    return snp_dfs, config_dfs, posts
