#!/usr/bin/env python
import finalize_table
import sys
import os

if len(sys.argv) != 5:
    print("USAGE: out_table bed mos_mdom_fasta out_file", file = sys.stderr)
    sys.exit(1)

if os.path.exists(sys.argv[4]):
    print("Error!:", file = sys.stderr)
    print(sys.argv[4] + " already exists!", file = sys.stderr)
    sys.exit(1)

finalize_table.create_table(*sys.argv[1:])
