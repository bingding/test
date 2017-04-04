#!/usr/bin/python
"""
Compare result tree and tell which files have an unmatched hash or size.
Also apply hooks to allow comparison of files which match except some criteria

Intented for Python 2.7 (dev environment)

"""
__author__ = 'Thibault North <north@newoncology.de>'
__license__ = 'Public Domain'
__version__ = '1.0'

import os
import sys
import hashlib

try:
    from colorama import Fore, Style
except ImportError:
    class Fore:
        RED=''
        GREEN=''
    class Style:
        RESET_ALL=''


from collections import defaultdict, OrderedDict

def filehash(filepath):
        """return the sha256 sum of filepath"""
        blocksize = 64*1024
        sha = hashlib.sha256()
        with open(filepath, 'rb') as fp:
            while True:
                data = fp.read(blocksize)
                if not data:
                    break
                sha.update(data)
        return sha.hexdigest() 

def find_files(root_dir, flist, results_dir, ignore_list=[]):
    print "Comparing files in folders:"
    print root_dir
    print results_dir
    seen_files = []
    for root, dirnames, filenames in os.walk(root_dir):
            for filename in filenames:
                for suffix, folder in flist.items():
                    if filename.endswith(suffix):
                        seen_files.append(suffix)
                        if suffix in ignore_list:
                            print "Ignoring {}".format(filename)
                            continue
                        # hash and look in the results
                        res_path = os.path.sep.join([results_dir, folder, filename])
                        if os.path.exists(res_path):
                            res_sha = filehash(res_path)
                        else:
                            print Fore.RED,
                            print "Couldn't find file {}".format(res_path)
                            print Style.RESET_ALL,
                            continue
                        full_path = os.path.sep.join([root, filename])
                        orig_sha = filehash(full_path)
                        if (orig_sha != res_sha):
                            print Fore.RED,
                            print "Mismatch hash: {} ({})\n{} ({})".format(full_path, orig_sha, res_path, res_sha)
                            print Style.RESET_ALL,
                        else:
                            print Fore.GREEN,
                            print "OK: {} ({})".format(filename, res_sha)
                            print Style.RESET_ALL,
                            
    if set(flist.keys()) != set(seen_files):
        diff = set(flist.keys()).symmetric_difference(set(seen_files))
        diff_whitelist = diff.difference(set(ignore_list))
        if (diff_whitelist):
            print "Not seen: {}".format(", ".join(diff_whitelist))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("d", help="The path to the pipeline workspace2")

    
    
    args = parser.parse_args()

    root_path = args.d
    
    for root, dirnames, filenames in os.walk(root_path):
        for dirname in dirnames:
            if dirname.startswith("."):
                continue
            comp = dirname 
            break # Keep first directory
            
    
    dv_name = comp.split("-")[1].lower()
    
    results_folder = "/mnt/pipeline-results/"
    filemaker_folder = "/mnt/pipeline-results-filemaker"
    bam_folder = "/mnt/bam-{}".format(dv_name)
    
    
    # Check files:_
    fdr = {"MUTS":"MUTS", "ENSPAN":"ENSPAN", "QC":"QC", "COV":"COV",
               "CN":"CN", "CCN":"CCN", "PDF":"PDF", "MSI":"MSI", "PARAMS": "PARAMS"}
    
    files_in_analysis = {"_muts.txt": fdr["MUTS"],
                     "_WFEmuts.txt" : fdr["MUTS"],
                     "_enspan.txt": fdr["ENSPAN"],
                     "_qc.txt" : fdr["QC"],
                     "_cov.txt" : fdr["COV"],
                     "_cn.txt" : fdr["CN"],
                     "_ccn.txt" : fdr["CCN"],
                     "_cn.pdf": fdr["PDF"], # Only DS
                     "_gcn.pdf" : fdr["PDF"],
                     "_mfilter.pdf" : fdr["PDF"],
                     "_CNfocality.pdf" : fdr["PDF"], # Only DS
                     "_msi_del.txt": fdr["MSI"],
                     "_params.txt" : fdr["PARAMS"],
                     "_families_count_ontarget.txt" : fdr["QC"], # Only DS
                     ".log": fdr["QC"],
                     "_cov.pdf" : fdr["PDF"],
                     }
    files_in_bam = {"_raw.seg" : "", 
                   "_cn.seg":"", 
                   "_enspan.bam": "", 
                   "_enspan.bam.bai": ""}
    
    print "Checking folder {}".format(results_folder)
    find_files(root_path, files_in_analysis, results_folder)
    print "Checking BAM folder {}".format(bam_folder)
    find_files(root_path, files_in_bam, bam_folder)
    
    
    ignore_list = ["_families_count_ontarget.txt", ".log", "_WFEmuts.txt"]
    for fmf in ["{}{}".format(filemaker_folder, i) for i in ('', '2')]:
        print "Checking folder {}".format(fmf)
        find_files(root_path, files_in_analysis, fmf, ignore_list)
    
    
