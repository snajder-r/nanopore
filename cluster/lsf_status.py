#!/omics/groups/OE0540/internal/software/users/snajder/miniconda3/bin/python
import subprocess as sp
import re
import sys

def parse_status(stat_str):
    if stat_str == "DONE":
        return "success"
    elif stat_str in {"PEND", "RUN", "PSUSP", "USUSP", "SSUSP"}:
        return "running"
    elif stat_str == "EXIT":
        return "failed"
    else:
        raise ValueError


jobid = sys.argv[1]
p = sp.Popen(f"bjobs {jobid}".split(" "), stdout=sp.PIPE, stderr=sp.PIPE)
out, err = p.communicate()

if "not found" in err.decode():
    print("failed")
    sys.exit()

out = out.decode().split("\n")

for line in out:
    if line.startswith(jobid):
        line = re.split(" +", line)
        if len(line) >= 3:
            try:
                status = parse_status(line[2])
                print(status)
                sys.exit()
            except ValueError:
                pass
print("failed")
