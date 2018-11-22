import finalize_table
import sys

if len(sys.argv) != 5:
    print("USAGE: out_table bed mos_mdom_fasta out_file")
    sys.exit(1)

finalize_table.create_table(*sys.argv[1:])
