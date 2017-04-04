#!/usr/bin/python
"""
Compare result tree and tell which files have an unmatched hash or size.
Also apply hooks to allow comparison of files which match except some criteria

Intended for Python 2.7 (dev environment)

"""
__author__ = 'Thibault North <north@newoncology.de>'
__license__ = 'Public Domain'
__version__ = '1.0'

import os
import hashlib
import fnmatch
import re
import tempfile


from collections import defaultdict, OrderedDict
from hook_pdf import compare_pdfs

class FileCompare:
    """
    Host class for comparing files by pairs
    """
    rep_name = "RUN_NAME"
    hooks = None
    whitelist = None
    accepted = None
    will_write_temp = True
    """
    Constructor setting paths to the files as well as hooks and
    whitelist
    """
    def __init__(self, path1, path2, will_write_temp=True):
        self.paths = {0:{"path":path1.strip(), "run_name":None},
                      1:{"path":path2.strip(), "run_name":None}}
        self.f_lists = {0:None, 1:None}
        
        # Define hooks and file whitelists
        self.hooks = {"_premut_data.txt":self.hook_premut_data, "_msi_del.txt":self.hook_msi_del,
                 "_qc.txt":self.hook_qc, "_params.txt": self.hook_params, 
                 "_snp_list.txt":self.hook_snp_list, "_muts.txt":self.hook_muts,
                 "_allmuts.txt":self.hook_allmuts, "_allmuts.vcf":self.hook_allmuts,
                 "_mutcall.vcf":self.hook_mutcall, "_mutcall_annotated.vcf":self.hook_mutcall,
                 "_mutcall_filtered.vcf":self.hook_mutcall, 
                 "_mutcall_emfiltered.vcf":self.hook_mutcall,
                 "_pileup_T.txt":self.hook_pileup, "_pileup_N.txt":self.hook_pileup,
                 "_cnum.txt": self.hook_cnum, "_WFEmuts.txt": self.hook_WFTmuts,
                 "_cn.seg":self.hook_cn_seg, "_raw.seg":self.hook_raw_seg,
                 ".pdf":compare_pdfs
                 }
        
        self.whitelist = [] #[".pdf"]
        self.will_write_temp = will_write_temp
        
    def print_whitelist(self):
        """Just print out the contents of the whitelist"""
        for kk in self.whitelist:
            print("Files ending with {} are whitelisted".format(kk))
            
    def print_accepted(self, sha=True, hook=True):
        print Fore.GREEN,
        if sha:
            print("Accepted because their checksum match:")
            for ff in sorted(self.accepted["sha"]):
                print ff
        if hook:
            print("Accepted because they match within a hook:")
            for ff in sorted(self.accepted["hook"]):
                print ff
        print Style.RESET_ALL,
        
    def compare(self):
        """Retrieve the file lists and abstract the file paths"""
        for i in range(2):
            self.f_lists[i] = self.list_recursive(self.paths[i]["path"], 
                                                  path_num=i)
            print "Reading path ", self.paths[i]["path"]
    def check_f_lists(self):
        """Make sure that each set contains the same files."""
        d1_set = set(self.f_lists[0].keys())
        d2_set = set(self.f_lists[1].keys())
        if d2_set != d1_set:
            print d1_set.symmetric_difference(d2_set)
            return False
        return True
    
    def size_sha_comp(self):
        """Compute file size and sha for each file pair"""
        dd = defaultdict(list)
        for d in self.f_lists.values(): 
            for key, value in d.iteritems():
                dd[key].append(value)
        
        # Check number of values
        for key, val in dd.iteritems():
            if len(val) != 2:
                print "Warning: incorrect number of items ({}) in ".format(len(val)), key
                
        # Now hash
        hash_list = []
        size_list = []
        accepted = {"sha":[], "hook":[]}
        for key, val in dd.iteritems():
            size = [os.path.getsize(fpath) for fpath in val]
            sha = [self.filehash(fpath) for fpath in val]
            if len(size) != 2:
                # Only one item, comparison doesn't make sense
                continue
            if len(set(sha)) != 1:
                # Check hook
                if self.apply_hook(*val):
                    accepted["hook"].append(key)
                    continue
                else:
                    print Fore.RED,
                    print "Despite hooks, comparison failed!",
                    print Style.RESET_ALL,
                print Fore.RED,
                print "Checksum differ for files: "
                for ff in val:
                    print ff,
                print "\nFile sizes: ",
                for fs in size:
                    print fs, " ",
                print
                print Style.RESET_ALL,
            else:
                #accepted["sha"].append('\n'.join([': '.join((ss, vv)) for vv,ss in zip(val, sha)]) + '\n')
                accepted["sha"].append(' '.join(val) + " : {}".format(sha[0]))
        
        self.accepted = accepted

    def filehash(self, filepath):
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

    def list_recursive(self, root_dir, path_num=None):
        """
        List file under the root_dir, recursively.
        
        A dictonary is returned. Its values are the complete file paths,
        and its keys contain the file path where the run name has been replaced
        by rep_name
        """
        rep_name = self.rep_name
        run_name = [kk for kk in root_dir.split(os.path.sep) if kk!=''][-1]
        if "_ANALYSIS" in run_name or "_BAM" in run_name:
            spl = run_name.split('_')
            run_name = '_'.join(spl[:-1])
            
        run_date_time = "-".join(run_name.split("-")[-1:])
        run_name_excl = "-".join(run_name.split("-")[:-1])
        
        run_name_alt = "_".join([run_name_excl, run_date_time])
        # Run with _L<something>_
        run_name_alt2 = "_".join([run_name_excl, "L1234", run_date_time])
        run_name_alt3 = '_'.join(run_name.split("-")[-1:])


        
        #if "_ANALYSIS" in run_name_alt2 or "_BAM" in run_name_alt2:
            #run_name_alt2 = '_'.join('_'.split(run_name_alt2)[:-1])
            
        #if "_ANALYSIS" in run_name_alt3 or "_BAM" in run_name_alt3:
            #run_name_alt3 = '_'.join('_'.split(run_name_alt3)[:-1])
            
       
        
        matches = {}
        # remove constant part : keep only the path before the last folder
        root_dir_container = os.path.sep.join(root_dir.split(os.path.sep)[:-2])
        self.root_dir_container = root_dir_container
        
        # Also include the path prefix that may be different for two runs
        run_name_alts = [run_name, run_name_alt, run_name_alt2, run_name_alt3,
                         root_dir_container]
        
                
        for root, dirnames, filenames in os.walk(root_dir):
            for filename in filenames: #fnmatch.filter(filenames, '*.c'):
                curr_filepath = os.path.join(root, filename)
                # Explode to replace file name
                # Turns out the date_time separator prefix is _ or -
                rep_path = curr_filepath[len(root_dir_container):]

                for rname in run_name_alts:
                    rep_path = os.path.sep.join([kk.replace(rname, rep_name) for kk in rep_path.split(os.path.sep)])

                if matches.has_key(rep_path):
                    print("Warning: file already seen: ", curr_filepath, 
                        " processed as ", rep_path)
                matches[rep_path] = curr_filepath
        self.paths[path_num]["run_name"] = run_name_alts
        return matches

    def read_file_pair(self, file1, file2):
        """Open each file and put the content as string in an array"""
        f_content = []
        with open(file1) as f1:
            f_content.append(f1.read())
        with open(file2) as f2:
            f_content.append(f2.read())
        return f_content
    
    def write_temp(self, contents, n1=".txt", n2=".txt"):
        if not self.will_write_temp: return
        dirpath = tempfile.mkdtemp()
        dests = ["f1_" + n1, "f2_" + n2]
        for dest, ct in zip(dests, contents):
            with open(os.path.join(dirpath, dest), 'w') as of:
                of.write(ct)
        print "Hooked files written in {} {}".format(os.path.join(dirpath, dests[0]),
                                                     os.path.join(dirpath, dests[1]))
        
    
    def hook_premut(self, file1, file2):
        """
        Just the second line (t_file path) and the 4th line (out_name path) are different.
        Replace them with self.rep_name
        """
        f_content = self.read_file_pair(file1, file2)
        for ii, _ in enumerate(f_content):
            for kk,_ in enumerate(self.paths[ii]["run_name"]):
                f_content[ii] = f_content[ii].replace(self.paths[ii]["run_name"][kk], self.rep_name)
        
        if (f_content[0] == f_content[1]):
            return True
        #print f_content[0]
        #print f_content[1]
        self.write_temp(f_content)
        return False
    
    def hook_premut_data(self, file1, file2):
        """only the header changes"""
        f_content = self.read_file_pair(file1, file2)
        for ii, _ in enumerate(f_content):
            f_content[ii] = re.sub(r"t_file\t(.*?)\n", "t_file\tT_FILE\n", 
                                   f_content[ii])
            f_content[ii] = re.sub(r"n_file\t(.*?)\n", "b_file\tN_FILE\n", 
                                   f_content[ii])
            f_content[ii] = re.sub(r"out_name\t(.*?)\n", "out_name\tOUT_NAME\n", 
                                   f_content[ii])
        if (f_content[0] == f_content[1]):
            return True
        self.write_temp(f_content)
        return False
    
    def hook_params(self, file1, file2):
        """
        Replace t_file, in_name, out_name, inp_file, out_file as for hook_premut
        """
        f_content = self.read_file_pair(file1, file2)
        names = ('fname', 't_file', 'n_file', 'in_name', 'out_name', 
                 'input_file', 'f_name', 'in_file', 'out_file', 'inp_file', 
                 'bam_file', 'mut_data',)
        pairs = [(kk, kk.upper()) for kk in names]
        
        for ii, _ in enumerate(f_content):
            for name, rep in pairs:
                f_content[ii] = re.sub("{}\t(.*?)\n".format(name), 
                                       "{}\t{}\n".format(name, rep), 
                                       f_content[ii])
        if (f_content[0] == f_content[1]):
            return True
        self.write_temp(f_content)
        return False
    
    def hook_msi_del(self, file1, file2):
        """
        replace each sample name with self.rep_name
        sample name has format: 
        A4180-DV16-7-P19-BC89-DS_S1_L1234_160912_130536
        """
        f_content = self.read_file_pair(file1, file2)
        for ii, _ in enumerate(f_content):
            for kk,_ in enumerate(self.paths[ii]["run_name"]):
                f_content[ii] = f_content[ii].replace(self.paths[ii]["run_name"][kk], self.rep_name)
            # remove root_dir_container
            #f_content[ii] = f_content[ii].replace(self.root_dir_container, '')
        
        if (f_content[0] == f_content[1]):
            return True
        self.write_temp(f_content)
        #print f_content[0] 
        #print f_content[1] 
        return False
    
    def hook_mutcall(self, file1, file2):
        """only the filedate changes"""
        f_content = self.read_file_pair(file1, file2)
        for ii, _ in enumerate(f_content):
            f_content[ii] = re.sub("##fileDate=(\d+)", "##fileDate=DATE", 
                                   f_content[ii])
        if (f_content[0] == f_content[1]):
            return True
        self.write_temp(f_content)
        return False
    
    def hook_qc(self, file1, file2, ext="_qc.txt"):
        """ Same as for premut"""
        f_content = self.read_file_pair(file1, file2)
        for ii, (_, ff) in enumerate(zip(f_content, (file1, file2))):
            path = re.sub("production\/|staging\/", "", ff)
            f_content[ii] = f_content[ii].replace(path.split(ext)[0], self.rep_name)
            # remove root_dir_container
            #f_content[ii] = f_content[ii].replace(self.root_dir_container, '')
        
        if (f_content[0] == f_content[1]):
            return True
        #print f_content[0] 
        #print f_content[1] 
        self.write_temp(f_content)
        return False
            
    def hook_snp_list(self, file1, file2):
        """as for msi hook"""
        return self.hook_pileup(file1, file2)
    
    def hook_cnum(self, file1, file2):
        """as for msi hook"""
        return self.hook_pileup(file1, file2)
    
    def hook_WFTmuts(self, file1, file2):
        """as for msi hook"""
        return self.hook_msi_del(file1, file2)
    
    def hook_muts(self, file1, file2):
        """as for msi hook"""
        return self.hook_msi_del(file1, file2)
    
    def hook_allmuts(self, file1, file2):
        """as for msi hook"""
        return self.hook_msi_del(file1, file2)
    
    def hook_pileup(self, file1, file2):
        """as for msi hook"""
        f_content = self.read_file_pair(file1, file2)
        for ii, _ in enumerate(f_content):
            f_content[ii] = re.sub("BAM_TUMOR=(.*?);", "BAM_TUMOR=FILE;", 
                                   f_content[ii])
        if (f_content[0] == f_content[1]):
            return True
        return self.hook_premut(file1, file2)

        return False
    
    def hook_cn_seg(self, file1, file2):
        """as for premut"""
        return self.hook_qc(file1, file2, ext="_cn.seg")
    
    def hook_raw_seg(self, file1, file2):
        """as for premut"""
        return self.hook_qc(file1, file2, ext="_raw.seg")
    

    def apply_hook(self, file1, file2):
        """
        When the comparison of two files has failed (different hash or size)
        hooks can be applied to remove items such as file names, which differ from
        run to run. 
        Then, the file are compared again, and the test passes if they are similar
        after this modification
        """
        
        hooks = self.hooks
        whitelist = self.whitelist
        for hook in hooks.keys():
            if file1.endswith(hook) and file2.endswith(hook):
                if hooks.has_key(hook):
                    return hooks[hook](file1, file2)
                
        for ff in whitelist:
            if file1.endswith(ff) and file2.endswith(ff):
                return True
            
        print "No hook for files ending with .{}".format(file1.split(".")[-1])
        return False
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("d1", help="The path to the pipeline results at rev0")
    parser.add_argument("d2", help="The path to the pipeline results at revX")
    parser.add_argument("--no-color", help="Disable colored output", action="store_true")
    parser.add_argument("--no-tempfile", help="Disable writing of temporaries when hooks fail", action="store_true")
    
    args = parser.parse_args()
    
    # Fake coloring class
    class Fore:
        RED=''
        GREEN=''
    class Style:
        RESET_ALL=''
        
    if not args.no_color:
        try:
            from colorama import Fore, Style
        except ImportError:
            print "No color available; install python-colorama"

    

    if args.d1 == args.d2:
        print("Error: comparing the same locations")
        sys.exit(1)
    
    comp = FileCompare(args.d1, args.d2, will_write_temp=not args.no_tempfile)
    comp.compare()
    
    if not comp.check_f_lists():
        print("Sets difference found!")
    
    comp.print_whitelist()
    comp.size_sha_comp()
    comp.print_accepted()
