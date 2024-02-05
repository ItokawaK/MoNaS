#!/usr/bin/env python

'''
Copyright (c) 2019, Kentaro Itokawa <itokawa@nih.go.jp>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import sys
import os
import json

'''
Creating JBrowse configuration file in json.
'''

def main():
    if(len(sys.argv) != 6):
        print("USAGE: make_jbrowse_config.py ref.fa ref.bed list.txt bam_dir trackList.json")
        sys.exit(1)

    ref_fa_path = sys.argv[1]
    ref_bed = sys.argv[2]
    list_path = sys.argv[3]
    bam_dir = sys.argv[4]
    out_file = sys.argv[5]

    ref_fai_path = ref_fa_path + '.fai'

    for file in [ref_fa_path, ref_fai_path, ref_bed, list_path]:
        if not os.path.isfile(file):
            sys.stderr.write(file + ' was not found. Try again.\n')
            sys.exit(1)

    #parse list file
    sample_names = []
    with open(list_path) as f:
        sample_names = [el[0] for el in [l.split() for l in f.readlines()]]

    #sample_names = ["sample_A", "sample_B"]

    tracks = []

    tracks.append({
            "label" : "refseqs",
            "urlTemplate": ref_fa_path
            })

    tracks.append({
                    "style": {
                        "arrowheadClass": "arrowhead",
                        "className" : "feature",
                        "strandArrow" : True
                        },
                    "label": "exons",
                    "urlTemplate" : ref_bed,
                    "type": "CanvasFeatures"
                })

    for sample in sample_names:
        bam_path = bam_dir + "/" + sample + ".bam"
        if os.path.isfile(bam_path) and os.path.isfile(bam_path + ".bai"):
            dict_tmp = {}
            dict_tmp['style'] = {"height": 50}
            dict_tmp['label'] = sample
            dict_tmp['storeClass'] = "JBrowse/Store/SeqFeature/BAM"
            dict_tmp['urlTemplate'] = bam_dir + "/" + sample + ".bam"
            dict_tmp['maxFeatureScreenDensity'] = 4
            dict_tmp['metadata.category'] = "BAM"
            dict_tmp['metadata.Description'] = sample
            dict_tmp['type'] = "JBrowse/View/Track/SNPCoverage"
            tracks.append(dict_tmp)
        else:
            sys.stderr.write(bam_path + ' was not found. Skipped.\n')

    out_dict = {"refSeqs": ref_fai_path,
            "tracks": tracks}

    with open(out_file, 'w') as f:
        json.dump(out_dict,f,indent=2)
        f.write("\n")

if __name__ == '__main__':
    main()