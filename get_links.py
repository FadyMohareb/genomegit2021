#!/usr/bin/env python

"""
INSTRUCTIONS: The script will align all the fasta files with each other  Arguments:
update_dependent_datasets.py -new_file- -threads- -file_size- -template_length-
"""

# Make the imports
import os
import sys
import datetime
import shutil
from subprocess import Popen, check_output
from update_functions import detect_updates
from reconstruct_functions import reconstruct_dataset
from ObtainAlignment_functions import obtain_alignment
from StoreAlignment import store_variables, load_variables, obtain_alignment_pickle
# is this too heavy? also creates commits without fastas
threads = sys.argv[1]
commit_list = int(sys.argv[2]) #cannot do hash thing for now
#print(commit_list)
if(os.path.isdir("./.git/info/temporary_directory")):
    shutil.rmtree("./.git/info/temporary_directory")
os.mkdir("./.git/info/temporary_directory")
if (commit_list == 0):
    # Get list of commits
    print("get commits")
    commit_list = check_output(
        "git log --format=%H", shell=True).decode("utf-8").split()
    print(commit_list)
    # check if alignments are needed
    if (len(commit_list) < 3):
        print("\n\t***All the Assembly files are aligned!***")
        sys.exit()
    #
    head = 0
    print(len(commit_list))
    for x in range(head, (len(commit_list)-1)):
        if (x==(len(commit_list)-2)):
            break
        else:
            #print("X: ",x)
            #print(commit_list[x])
            print("\n\t***Finding Alignment Pickle {}***".format(str(datetime.datetime.now())))
            print("\n\t Warning: If you have added any data into the repository but did not commit, these changes will be lost.")
            ShellCommand = Popen(
                    "git checkout " + commit_list[x] + " 2> /dev/null", shell=True).wait()
            reconstruct_dataset(size=60, directory="./Genome",
                    output_file="./.git/info/temporary_directory/assembly_1.fa", mode="Genome")
            next_commit = x+2
            for y in range(next_commit, len(commit_list)):
                #print("Y: ",y)
                #print(commit_list[y])
                # Reconstruct the second version of the assembly
                ShellCommand = Popen(
                    "git checkout " + commit_list[y] + " 2> /dev/null", shell=True).wait()
                reconstruct_dataset(size=60, directory="./Genome",
                            output_file="./.git/info/temporary_directory/assembly_2.fa", mode="Genome")
                # Obtain alignment pickle of both assemblies
                alignment_pickle1 = obtain_alignment_pickle(
                    "./.git/info/temporary_directory/assembly_1.fa", "./.git/info/temporary_directory/assembly_2.fa")
                alignment_pickle2 = obtain_alignment_pickle(
                    "./.git/info/temporary_directory/assembly_2.fa", "./.git/info/temporary_directory/assembly_1.fa")
                ShellCommand = Popen("git checkout master 2> /dev/null", shell=True).wait()

                # Check if the alignemnt containing the information between these two assemblies already exists
                if(os.path.isdir(alignment_pickle1) | os.path.isdir(alignment_pickle2)):
                    #print("already", alignment_pickle2)
                    continue
                else:
                    #print(alignment_pickle1)
                    print("\n\t***Making Alignment {}***".format(str(datetime.datetime.now())))
                    ToUpdate = detect_updates("./RepoMap.txt")
                    # Obtain the variables
                    variables = obtain_alignment(old_assembly="./.git/info/temporary_directory/assembly_1.fa",
                                                new_assembly="./.git/info/temporary_directory/assembly_2.fa",
                                                directory="./.git/info/temporary_directory",
                                                threads=threads, ToUpdate=ToUpdate,
                                                alignment_pickle=alignment_pickle1, aligner_switch=2,
                                                percent_identity=95, kmer=15, segLength=5000, c_flag=2000, b_flag=1, ms_flag=1)

                    # Store the pickle [tabix_queries,OldNewID_Dict,alignment_pickle,summary_Dict]
                    store_variables(variables=variables, alignment_pickle=alignment_pickle1)
ShellCommand = Popen("git checkout master 2> /dev/null", shell=True).wait()
print("\n\t***All the Assembly files are aligned!***")
        

            
