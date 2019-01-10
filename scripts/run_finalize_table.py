#!/usr/bin/env python
import finalize_table
import sys
import os

if len(sys.argv) != 4:
    print("USAGE: out_table bed mos_mdom_fasta", file = sys.stderr)
    sys.exit(1)


finalize_table.create_table(sys.argv[1],
                            sys.argv[2],
                            sys.argv[3],
                            sys.path[0] + "/kdr_list.json",
                            "/dev/stdout")
