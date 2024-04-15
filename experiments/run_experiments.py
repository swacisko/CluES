import threading
from multiprocessing import Pool
import itertools as iter
import subprocess
import os
import time
from concurrent.futures import ThreadPoolExecutor


def runClues( filt ):
    outpath = os.path.join('results', 'jsons', 'res_' + filt[0] + '_' + str(filt[1]) + ext)
    cmd = f'./CluES_exp {filt[0]} {filt[1]} > {outpath}'
    print(f'running command: {cmd} in thread with id: {threading.get_ident()}')
    subprocess.call(cmd, shell=True)
    time.sleep(1) # this is just for testing if calling work indeed independently in threads



if __name__ == '__main__':

    delim = '__'
    ext = '.json'

    os.environ["BENCHMARK_FORMAT"] = "json" #setting environment variable to make json as the standard output format

    indep_benchmars = []
    for (psid,A,B) in iter.product(
                    [ 5,4,3,2,1 ],
                    ['recursion', 'fixed_time', 'algorithm', 'ls_iteration', 'perturbation', 'neg_map', 'granularity'],
                    ['A', 'B', 'C', 'D'],
                    ):
        if psid < 9 or (psid == 9 and A == 'granularity'):
            indep_benchmars.append( (A+delim+B,psid) )

    executor = ThreadPoolExecutor(max_workers=4)
    executor.map(runClues, indep_benchmars)

