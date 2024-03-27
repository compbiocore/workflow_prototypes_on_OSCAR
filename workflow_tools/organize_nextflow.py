from collections import defaultdict
import argparse
import time

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, required=True)
    args = parser.parse_args()
    return args

def main(args):
    file_path = args.file
    processes = []

    # preprocessing the nextflow file to find process names
    with open(file_path, "r") as f:
        for line in f:
            if "process" in line:
                processes.append(line.split(" ")[1])

    """
        ln -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_star.sh
        ln -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_star_slurm.stderr
        ln -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_star_slurm.stdout    
    """

    # storing the commands for each of the processes- simlink before running and copy after running
    # we are also traversing backwards- so we insert each set of cmds first
    process_cmds_sim = []
    for process in processes:
        cmds = []
        cmds.append("\t\tln -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_" + process + ".sh\n")
        cmds.append("\t\tln -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stderr\n")
        cmds.append("\t\tln -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stdout\n")
        process_cmds_sim.insert(0, cmds)
    
    process_cmds_cp = []
    for process in processes:
        cmds = []
        cmds.append("\t\tcp -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_" + process + ".sh\n")
        cmds.append("\t\tcp -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stderr\n")
        cmds.append("\t\tcp -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stdout\n")
        process_cmds_cp.insert(0, cmds)
    

    with open(file_path, "r+") as f:
        contents = f.readlines()
        script_idx = defaultdict(list)
        process_idx = 0

        
        # traversing backwards through the file, so added lines don't mess up the line numbering (repeated numbers)
        for file_idx, line in reversed(list(enumerate(contents))):
            
            if "\"\"\"" in line:
                script_idx[0].append(file_idx)
                # since we know that the script """ comes in pairs, we can add all to a list and check if it's the first in the pair

                if len(script_idx[0]) % 2 == 0:
                    # iterate through each of the unique processes
                    for idx, cmd in enumerate(process_cmds_sim[process_idx]):
                        # begin inserting after the """
                        contents.insert(script_idx[0][-1] + idx + 1, cmd)
                    process_idx += 1

                elif len(script_idx[0]) % 2 == 1:
                    for idx, cmd in enumerate(process_cmds_cp[process_idx]):
                        # begin inserting before the """
                        contents.insert(script_idx[0][-1] + idx, cmd)

            elif "$/" in line or "/$" in line:
                print(file_idx)
                script_idx[1].append(file_idx)
                
                if len(script_idx[1]) % 2 == 0:
                    # iterate through each of the unique processes
                    # begin inserting after the $/
                    for idx, cmd in enumerate(process_cmds_sim[process_idx]):
                        contents.insert(script_idx[1][-1] + idx + 1, cmd)
                    process_idx += 1

                elif len(script_idx[1]) % 2 == 1:
                    # begin inserting before the /$
                    for idx, cmd in enumerate(process_cmds_cp[process_idx]):
                        contents.insert(script_idx[1][-1] + idx, cmd)

                process_idx += 1

    with open(file_path, "w") as f:
        f.writelines(contents)

def main_v2(args):
    file_path = args.file
    processes = []
    process_cmds_sim = defaultdict()
    process_cmds_cp = defaultdict()

    # preprocessing the nextflow file to find process names
    with open(file_path, "r") as f:
        for line in f:
            if "process" in line:
                processes.append(line.split(" ")[1])
                process_cmds_sim[line.split(" ")[1]] = 0
                process_cmds_cp[line.split(" ")[1]] = 0

    """
        ln -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_star.sh
        ln -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_star_slurm.stderr
        ln -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_star_slurm.stdout    
    """

    # storing the commands for each of the processes- simlink before running and copy after running
    # we are also traversing backwards- so we insert each set of cmds first
    for process in process_cmds_sim:
        cmds = []
        cmds.append("\t\tln -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_" + process + ".sh\n")
        cmds.append("\t\tln -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stderr\n")
        cmds.append("\t\tln -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stdout\n")
        process_cmds_sim[process] = cmds
    
    #print(process_cmds_sim)
    for process in process_cmds_cp:
        cmds = []
        cmds.append("\t\tcp -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_" + process + ".sh\n")
        cmds.append("\t\tcp -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stderr\n")
        cmds.append("\t\tcp -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stdout\n")
        process_cmds_cp[process] = cmds
    
    with open(file_path, "r+") as f:
        contents = f.readlines()
        script_idx = defaultdict(list)
        current_process = ""

        # traversing backwards through the file, so added lines don't mess up the line numbering (repeated numbers)
        for file_idx, line in enumerate(contents):
            if "process" in line:
                current_process = line.split(" ")[1]
                print("current process", current_process)
                
            if "\"\"\"" in line:
                script_idx[0].append(file_idx)
                # since we know that the script """ comes in pairs, we can add all to a list and check if it's the first in the pair

                if len(script_idx[0]) % 2 == 1:
                    # iterate through each of the unique processes
                    for idx, cmd in enumerate(process_cmds_sim[current_process]):
                        # begin inserting after the """
                        contents.insert(script_idx[0][-1] + idx + 1, cmd)
                    

                elif len(script_idx[0]) % 2 == 0:
                    for idx, cmd in enumerate(process_cmds_cp[current_process]):
                        # begin inserting before the """
                        contents.insert(script_idx[0][-1] + idx, cmd)
                

            elif "$/" in line:
                print(file_idx)
                script_idx[1].append(file_idx)
                
                if len(script_idx[1]) % 2 == 0:
                    # iterate through each of the unique processes
                    # begin inserting after the $/
                    for idx, cmd in enumerate(process_cmds_sim[current_process]):
                        contents.insert(script_idx[1][-1] + idx + 1, cmd)
                    

                elif len(script_idx[1]) % 2 == 1:
                    # begin inserting before the /$
                    for idx, cmd in enumerate(process_cmds_cp[current_process]):
                        contents.insert(script_idx[1][-1] + idx, cmd)


    with open(file_path, "w") as f:
        f.writelines(contents)

def main_v3(args):
    file_path = args.file
    processes = ["start"]
    process_cmds_sim = defaultdict()
    process_cmds_cp = defaultdict()

    # preprocessing the nextflow file to find process names
    with open(file_path, "r") as f:
        for line in f:
            if "process" in line:
                processes.append(line.split(" ")[1])
                process_cmds_sim[line.split(" ")[1]] = 0
                process_cmds_cp[line.split(" ")[1]] = 0

    """
        ln -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_star.sh
        ln -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_star_slurm.stderr
        ln -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_star_slurm.stdout    
    """

    # storing the commands for each of the processes- simlink before running and copy after running
    # we are also traversing backwards- so we insert each set of cmds first
    for process in process_cmds_sim:
        cmds = []
        cmds.append("\t\tln -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_" + process + ".sh\n")
        cmds.append("\t\tln -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stderr\n")
        cmds.append("\t\tln -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stdout\n")
        process_cmds_sim[process] = cmds
    
    #print(process_cmds_sim)
    for process in process_cmds_cp:
        cmds = []
        cmds.append("\t\tcp -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_" + process + ".sh\n")
        cmds.append("\t\tcp -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stderr\n")
        cmds.append("\t\tcp -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stdout\n")
        process_cmds_cp[process] = cmds
    
    with open(file_path, "r+") as f:
        contents = f.readlines()
        script_idx = defaultdict(list)
        current_process = processes[-1]
        #print("current process", current_process)

        # traversing backwards through the file, so added lines don't mess up the line numbering (repeated numbers)
        for file_idx, line in reversed(list(enumerate(contents))):
            if "process" in line:
                processes.pop()
                current_process = processes[-1]
                #print("current process", current_process)
                
            if "\"\"\"" in line:
                script_idx[0].append(file_idx)
                # since we know that the script """ comes in pairs, we can add all to a list and check if it's the first in the pair

                if len(script_idx[0]) % 2 == 0:
                    # iterate through each of the unique processes
                    for idx, cmd in enumerate(process_cmds_sim[current_process]):
                        # begin inserting after the """
                        contents.insert(script_idx[0][-1] + idx + 1, cmd)
                    

                elif len(script_idx[0]) % 2 == 1:
                    for idx, cmd in enumerate(process_cmds_cp[current_process]):
                        # begin inserting before the """
                        contents.insert(script_idx[0][-1] + idx, cmd)
                

            elif "$/" in line:
                print(file_idx)
                script_idx[1].append(file_idx)
                
                if len(script_idx[1]) % 2 == 1:
                    # iterate through each of the unique processes
                    # begin inserting after the $/
                    for idx, cmd in enumerate(process_cmds_sim[current_process]):
                        contents.insert(script_idx[1][-1] + idx + 1, cmd)
                    

                elif len(script_idx[1]) % 2 == 0:
                    # begin inserting before the /$
                    for idx, cmd in enumerate(process_cmds_cp[current_process]):
                        contents.insert(script_idx[1][-1] + idx, cmd)


    with open(file_path, "w") as f:
        f.writelines(contents)

if __name__ == "__main__":
    args = parse_args()
    main_v3(args)
