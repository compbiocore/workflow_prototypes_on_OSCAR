# OUTLINE
# 1. Assuming we are working in the work directory, go through each of the session identifiers (SI) folders
# 2. For each SI, loop through each unique session ID and check last modified
# 3. Return last modified, and sort by modified date
# 4. Give option to clean files greater than n days

import os
import shutil
import datetime
from collections import defaultdict
import sys

n = 30
current_files = defaultdict()

# Our directory hierachy looks like this /work/a9/eaioaec109231/
# want to remove /eaioaec109231
try:
    os.chdir("./work")
    #print("Should be in the work directory", os.getcwd())
except FileNotFoundError:
    sys.exit("./work directory not found!")
    
for file in os.listdir("./"):
    #print("Current file looked at: ", file)
    #print("current WD", os.getcwd())

    if os.path.isdir(file):
        #print(file)
        os.chdir(file)
        #print("Should be in the next directory!", os.getcwd())

        for item in os.listdir("./"):
            #print(item)
            if os.path.isdir(item):
                #print("changed into: ", item)
                #print(os.getcwd())
                t = os.path.getmtime(item)
            #   print(datetime.datetime.fromtimestamp(t))
            #    print(datetime.datetime.now())
                last_modified = datetime.datetime.now() - datetime.datetime.fromtimestamp(t)

                current_files[os.getcwd() + "/" + item] = last_modified
                
                #print("last mod: ", last_modified.days)
                #print("current item", item)
        # print("After adding items to mod folderes", os.getcwd())
        os.chdir("../")
    continue
    

#print(current_files)
sorted_folders = dict(reversed(sorted(current_files.items(), key=lambda x: x[1])))
#print(sorted_folders)
for folder in sorted_folders:
    if sorted_folders[folder].days >= n:
        print(folder, "Last modified: " + str(sorted_folders[folder].days) + " days ago")
    
    
confirmation = input("Do you want to remove these files? (yes/no) ").lower()
while True:
    if confirmation == "yes":
        for folder in sorted_folders:
            if sorted_folders[folder].days >= n:
                shutil.rmtree(folder)

        print("Files have been removed.")
        break
    elif confirmation == "no":
        print("no")
        break
    else:
        confirmation = input("Input not recognized. Do you want to remove these files? (yes/no) ").lower()
        
