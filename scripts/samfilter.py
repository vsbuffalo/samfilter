## samfilter.py -- samtools filter
__version__ = 0.3
import sys
try:
    import pysam
except ImportError as e:
    sys.exit("Could not find pysam; if not installed, download and install from http://code.google.com/p/pysam/.")
import argparse
from os import path
import pdb

parser = argparse.ArgumentParser(description="""
Filter SAM/BAM file based on specified filters.

Note: These filters are not necessarily independent because in some
contexts, this would not make sense.

Filtering by --rids and --rids-files does not make sense if the read
is not mapped (and thus doesn't have a reference name), so either of
these filters also add the condition that the read is mapped (i.e., it
implicitly adds --mapped). Also, --mapq, --forward, and --reverse also
only make sense if the read is mapped, so this is an implicit filter
in these options as well. --different-rnames also only makes sense if
*both* target sequence and mate are mapped. Likewise, --xt-unique implies
--mapped.

Filtering by --proper-pair and --mate-mapped also only makes sense if
a read is paired, since a single-ended read would lead this to be
false.
""")
parser.add_argument("file", type=str, # not file, need to determine whether binary
                    help="SAM/BAM file")
parser.add_argument('--qids', help="CSV string of query names.", type=str, default=None)
parser.add_argument('--qids-file', help="file of query names, one per line.", type=argparse.FileType("r"))
parser.add_argument('--rids', help="CSV string of reference names.", type=str, default=None)
parser.add_argument('--rids-file', help="file of reference names, one per line.", type=argparse.FileType("r"))
parser.add_argument('--mapped', help="query must be mapped", action="store_true", default=None)
parser.add_argument('--paired', help="query must be paired", action="store_true", default=None)
parser.add_argument('--single', help="query must be unpaired (single-end read)", action="store_true", default=None)
parser.add_argument('--unmapped', help="query must be unmapped", action="store_true", default=None)
parser.add_argument('--reverse', help="query must be on reverse strand", action="store_true", default=None)
parser.add_argument('--xt-unique', help="query sequence mapped uniquely (BWA's XT must equal U)", action="store_true", default=None)
parser.add_argument('--forward', help="query must be on forward strand", action="store_true", default=None)
parser.add_argument('--mapq', help="mapping quality must be greater than or equal to this value.", default=None, type=int)
parser.add_argument('--different-rnames', help="the target and mate reference names differ (and both are mapped).", action="store_true", default=None)
parser.add_argument('--not-proper-pair', help="must not be in a proper pair.", action="store_true", default=None)
parser.add_argument('--proper-pair', help="must be in proper pair.", action="store_true", default=None)
parser.add_argument('--mate-mapped', help="mate must be mapped.", action="store_true", default=None)
parser.add_argument('--output', help="output file (default stdin)",
                    default="stdin", type=str)
parser.add_argument('--output-bam', help="use binary (BAM) as output format",
                    default=False, action="store_true")

def build_sam_filters(samfile, qids=None, rids=None, paired=None,
                      single=None, mapped=None,
                      xt_unique=None,
                      unmapped=None, reverse=None,
                      forward=None, proper_pair=None,
                      not_proper_pair=None, different_rnames=None,
                      mate_mapped=None, mapq=None):
    """
    Return a list of functions that return True if the condition is
    met. Each x is read from samfile.fetch().
    """
    filters = dict()

    # both of these are O(1) lookup, assuming low collisions (which
    # should be handled well, as Python's dicts use open addressing,
    # not seperate chaining).
    if qids is not None:
        filters['qid'] = lambda x: qids.get(x.qname, False)
    if rids is not None:
        filters['rid'] = lambda x: not x.is_unmapped and rids.get(samfile.getrname(x.tid), False)
    if paired is not None:
        filters['paired'] = lambda x: x.is_paired
    if single is not None:
        filters['single'] = lambda x: not x.is_paired
    if mapped is not None:
        filters['mapped'] = lambda x: not x.is_unmapped
    if xt_unique is not None:
        filters['xt_unique'] = lambda x: not x.is_unmapped and dict(x.tags)["XT"] == "U"
    if unmapped is not None:
        filters['unmapped'] = lambda x: x.is_unmapped
    if reverse is not None:
        filters['reverse'] = lambda x: not x.is_unmapped and x.is_reverse
    if forward is not None:
        filters['forward'] = lambda x: not x.is_unmapped and not x.is_reverse
    if proper_pair is not None:
        filters['proper_pair'] = lambda x: x.is_paired and x.is_proper_pair
    if not_proper_pair is not None:
        filters['not_proper_pair'] = lambda x: x.is_paired and not x.is_proper_pair
    if different_rnames is not None:
        filters['different_rnames'] = lambda x: not x.is_unmapped and not x.mate_is_unmapped and x.is_paired and x.tid != x.rnext
    if mate_mapped is not None:
        filters['mate_mapped'] = lambda x: x.is_paired and not x.mate_is_unmapped
    if mapq is not None:
        filters['mapq'] = lambda x: not x.is_unmapped and x.mapq >= mapq

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

    def build_id_hash(id_file=None, id_string=None):
        """
        For either a file of IDs or a string of CSV ids, return a hash
        of them.
        """
        ids = list()
        if id_string:
            ids.extend(re.split(" ?,", id_string))
        if id_file:
            try:
                with id_file as f:
                    ids.extend([q.strip() for q in f.readlines()])
            except IOError as e:
                sys.exit("Cannot open file '%s'." % id_file)

        # TODO check against duplicates?
        ids_hash = dict(zip(ids, [True]*len(ids)))
        if len(ids_hash) == 0:
            ids_hash = None
        return ids_hash

    qids_hash = build_id_hash(args.qids_file, args.qids)
    rids_hash = build_id_hash(args.rids_file, args.rids)

    # build list of filters, all must be true
    sam_filters = build_sam_filters(samfile, qids=qids_hash, rids=rids_hash,
                                    single=args.single, paired=args.paired,
                                    mapped=args.mapped, xt_unique=args.xt_unique, 
                                    unmapped=args.unmapped,
                                    reverse=args.reverse,
                                    forward=args.forward,
                                    proper_pair=args.proper_pair,
                                    not_proper_pair=args.not_proper_pair,
                                    different_rnames=args.different_rnames,
                                    mate_mapped=args.mate_mapped,
                                    mapq=args.mapq)

    
    nkept = 0
    ntotal = 0
    try:
        for read in samfile:
            ntotal += 1

            # If we're looking for unique, check this condition and
            # die if not true. This is because indicators of read
            # mapping uniqueness are non-orthogonal.
            if args.xt_unique is not None:
                dtags = dict(read.tags)
                if not read.is_unmapped and dtags["XT"] == "U" and dtags["X0"] > 1:
                    raise Exception, "tags X0 and XT do not agree."

            if all([f(read) for f in sam_filters.values()]):
                nkept += 1
                out_samfile.write(read)
                
        sys.stderr.write("%d kept, %d total.\n" % (nkept, ntotal))
        str_filters = "filters used: " + ', '.join(sam_filters.keys()) + "\n"
        sys.stderr.write(str_filters)
    except KeyboardInterrupt:
        sys.stderr.write("Keyboard interrupt - aborting.\n")
    finally:
        # try to close stdout, stderr
        try:
            sys.stdout.close()
        except:
            pass
        try:
            sys.stderr.close()
        except:
            pass
