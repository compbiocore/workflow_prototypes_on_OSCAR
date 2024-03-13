from collections import defaultdict
import argparse

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
                
                #print(line)


    """
        ln -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_star.sh
        ln -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_star_slurm.stderr
        ln -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_star_slurm.stdout    
    """
    #print(processes)

    # storing the commands for each of the processes
    process_cmds = []
    for process in processes:
        cmds = []
        cmds.append("\t\tcp -sf `realpath .command.sh` ${params.out_dir}/scripts/${sample_id}_" + process + ".sh\n")
        cmds.append("\t\tcp -sf `realpath .command.err` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stderr\n")
        cmds.append("\t\tcp -sf `realpath .command.log` ${params.out_dir}/logs/${sample_id}_" + process + "_slurm.stdout\n")
        process_cmds.append(cmds)


    with open(file_path, "r+") as f:
        contents = f.readlines()
        script_idx = defaultdict(list)
        
        for idx, line in enumerate(contents):
            if "\"\"\"" in line:
                script_idx[0].append(idx)
                
                # since we know that the script """ comes in pairs, we can add all to a list and check if it's the first in the pair
                if len(script_idx[0]) % 2 == 1:
                    # iterate through each of the unique processes
                    process_idx = 0
                    for idx, cmd in enumerate(process_cmds[process_idx]):
                        # begin inserting after the quotations
                        contents.insert(script_idx[0][-1] + idx + 1, cmd)

                    process_idx += 1

            elif "$/" in line:
                script_idx[1].append(idx)
                
                # similar process, but $/ is a unique pair, so we only need to consider it once
                process_idx = 0
                for idx, cmd in enumerate(process_cmds[process_idx]):
                    # begin inserting after the quotations
                    contents.insert(script_idx[1][-1] + idx + 1, cmd)

                process_idx += 1

    with open(file_path, "w") as f:
        f.writelines(contents)

if __name__ == "__main__":
    args = parse_args()
    main(args)
