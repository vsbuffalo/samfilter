## sf.py -- samtools filter

import sys
try:
    import pysam
except ImportError as e:
    sys.exit("Could not find pysam; if not installed, download and install from http://code.google.com/p/pysam/.")
import argparse
from os import path
import re
import csv

parser = argparse.ArgumentParser(description="filter SAM/BAM file based on input.")
parser.add_argument("file", type=str, # not file, need to determine whether binary
                    help="SAM/BAM file")
parser.add_argument('--qids', help="CSV string of query names.", type=str, default=None)
parser.add_argument('--qids-file', help="file of query names, one per line.", type=argparse.FileType("r"))
parser.add_argument('--mapped', help="query must be mapped", action="store_true", default=None)
parser.add_argument('--paired', help="query must be paired", action="store_true", default=None)
parser.add_argument('--unmapped', help="query must be unmapped", action="store_true", default=None)
parser.add_argument('--reverse', help="query must be on reverse strand", action="store_true", default=None)
parser.add_argument('--forward', help="query must be on forward strand", action="store_true", default=None)
parser.add_argument('--mapq', help="mapping quality must be greater than or equal to this value.", default=None, type=int)
parser.add_argument('--proper-pair', help="must be in proper pair.", action="store_true", default=None)
parser.add_argument('--mate-mapped', help="mate must be mapped.", action="store_true", default=None)
parser.add_argument('--output', help="output file (default stdin)",
                    default="stdin", type=str)
parser.add_argument('--output-bam', help="use binary (BAM) as output format",
                    default=False, action="store_true")

def build_sam_filters(qids=None, paired=None, mapped=None,
                      unmapped=None, reverse=None,
                      forward=None, proper_pair=None,
                      mate_mapped=None, mapq=None):
    """
    Return a list of functions that return True if the condition is
    met. Each x is read from samfile.fetch().
    """
    filters = dict()

    if qids is not None:
        filters['qid'] = lambda x: qids.get(x.qname, False)
    if paired is not None:
        filters['paired'] = lambda x: x.is_paired
    if mapped is not None:
        filters['mapped'] = lambda x: not x.is_unmapped
    if unmapped is not None:
        filters['unmapped'] = lambda x: x.is_unmapped
    if reverse is not None:
        filters['reverse'] = lambda x: x.is_reverse
    if forward is not None:
        filters['forward'] = lambda x: not x.is_reverse
    if proper_pair is not None:
        filters['proper_pair'] = lambda x: x.is_paired and x.is_proper_pair
    if mate_mapped is not None:
        filters['mate_mapped'] = lambda x: x.is_paired and not x.mate_is_unmapped
    if mapq is not None:
        filters['mapq'] = lambda x: x.mapq >= mapq

    return filters

if __name__ == "__main__":
    args = parser.parse_args()
    
    # figure out if SAM/BAM from extension
    read_mode = "rb" if path.splitext(args.file)[1] == ".bam" else "r"
    try:
       samfile = pysam.Samfile(args.file, read_mode)
    except IOError as e:
        sys.exit("Cannot open file '%s'." % args.file)

    # open output
    if args.output != "stdin":
        write_mode = "wb" if args.output_bam else "w"
        try:
            out_samfile = pysam.Samfile(args.output, write_mode, template=samfile)
        except IOError as e:
            sys.exit("Cannot open output file '%s' for writing." % args.output)
    else:
        out_samfile = pysam.Samfile("-", "w", template=samfile)


    # check if there are query ids to subset
    qids = list()
    if args.qids:
        qids.extend(re.split(" ?,", args.qids))
    if args.qids_file:
        try:
            with args.qids_file as f:
                qids.extend([q.strip() for q in f.readlines()])
        except IOError as e:
            sys.exit("Cannot open file '%s'." % args.qids_file)

    # TODO check against duplicates?
    qids_hash = dict(zip(qids, [True]*len(qids)))
    if len(qids_hash) == 0:
        qids_hash = None

    # build list of filters, all must be true
    sam_filters = build_sam_filters(qids=qids_hash, paired=args.paired,
                                    mapped=args.mapped,
                                    unmapped=args.unmapped,
                                    reverse=args.reverse,
                                    forward=args.forward,
                                    proper_pair=args.proper_pair,
                                    mate_mapped=args.mate_mapped,
                                    mapq=args.mapq)
    
    nkept = 0
    ntotal = 0
    for read in samfile:
        ntotal += 1
        if all([f(read) for f in sam_filters.values()]):
            nkept += 1
            out_samfile.write(read)


    sys.stderr.write("%d kept, %d total.\n" % (nkept, ntotal))
    str_filters = "filters used: " + ', '.join(sam_filters.keys()) + "\n"
    sys.stderr.write(str_filters)

    # try to close stdout, stderr
    try:
        sys.stdout.close()
    except:
        pass
    try:
        sys.stderr.close()
    except:
        pass    
