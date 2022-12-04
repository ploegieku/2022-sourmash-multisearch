#! /usr/bin/env python
"""
Perform
"""
import os
import sys
import argparse
from math import log, exp
import sourmash
from sourmash import sourmash_args

from sourmash.cli.utils import (
    add_ksize_arg,
    add_moltype_args,
    add_picklist_args,
    add_scaled_arg,
    add_pattern_args,
)
from sourmash.logging import notify, error, print_results, set_quiet
from sourmash.sourmash_args import FileOutputCSV, SaveSignaturesToLocation

from sourmash.search import (
    search_databases_with_flat_query,
    search_databases_with_abund_query,
)
from scipy.special import logsumexp


# def get_prob_overlap(df, freq_a, freq_b, logged=False):
#     # import pdb; pdb.set_trace()
#     hashes = df["hashval"].values
#
#     # Not all hashes are present for some reason -> take the intersection for now while I figure that bug out
#     # Was something done at scaled=10 maybe? and then this is scaled=5 which has half the k-mers?
#     observed_hashes = freq_a.index.intersection(freq_b.index.intersection(hashes))
#     if observed_hashes.empty:
#         return 0
#
#     freq_a_subset = freq_a[observed_hashes]
#     freq_b_subset = freq_b[observed_hashes]
#     if logged:
#         # Log-probabilities to prevent underflow
#         log_probabilities = np.logaddexp(freq_a_subset, freq_b_subset)
#         prob_overlap = logsumexp(log_probabilities)
#         # import pdb; pdb.set_trace()
#     else:
#         # Raw probabilities
#         prob_overlap = freq_a_subset.multiply(freq_b_subset).sum()
#
#     return prob_overlap


def sum_hash_occurences(signature):
    return sum(v for k, v in signature.minhash.hashes.items())


def get_prob_overlap(query: sourmash.SourmashSignature, result):
    query_hashes = query.minhash.hashes
    result_hashes = result.match.minhash.hashes
    result_sum = sum_hash_occurences(result.match)
    query_sum = sum_hash_occurences(query)
    prob_overlap = sum(
        n_query / query_sum * result_hashes[hashval] / result_sum
        for hashval, n_query in query_hashes.items()
        if hashval in result_hashes
    )

    # logprob_overlap = logsumexp([
    #     log(n_query)
    #     - log(query_sum)
    #     + log(result_hashes[hashval])
    #     - log(result_sum)
    #     for hashval, n_query in query_hashes.items()
    #     if hashval in result_hashes]
    # )
    # print(f'prob_overlap: {prob_overlap}')
    # print(f'logprob_overlap: {logprob_overlap}')
    # print(f"exp(logprob_overlap): {exp(logprob_overlap)}")
    return prob_overlap


def main():
    args, p = multisearch_parser()

    set_quiet(args.quiet, args.debug)
    moltype = sourmash_args.calculate_moltype(args)
    picklist = sourmash_args.load_picklist(args)
    pattern_search = sourmash_args.load_include_exclude_db_patterns(args)

    # flatten --db and --query
    args.databases = [item for sublist in args.databases for item in sublist]
    inp_files = [item for sublist in args.query for item in sublist]
    if args.query_from_file:
        more_files = sourmash_args.load_pathlist_from_file(args.query_from_file)
        inp_files.extend(more_files)

    print(f"inp_files {inp_files}")

    # need a query to get ksize, moltype for db loading
    query = next(
        iter(
            sourmash_args.load_file_as_signatures(
                inp_files[0], ksize=args.ksize, select_moltype=moltype
            )
        )
    )

    if args.scaled:
        if not query.minhash.scaled:
            error("cannot downsample a signature not created with --scaled")
            sys.exit(-1)
        if args.scaled != query.minhash.scaled:
            notify(
                f"downsampling query from scaled={query.minhash.scaled} to {int(args.scaled)}"
            )
            with query.update() as query:
                query.minhash = query.minhash.downsample(scaled=args.scaled)

    # set up the search databases
    is_containment = args.containment or args.max_containment
    if is_containment:
        if args.containment and args.max_containment:
            notify("ERROR: cannot specify both --containment and --max-containment!")
            sys.exit(-1)

    databases = sourmash_args.load_dbs_and_sigs(
        args.databases,
        query,
        not is_containment,
        picklist=picklist,
        pattern=pattern_search,
        fail_on_empty_database=args.fail_on_empty_database,
    )

    query_background_mh = None
    db_background_mh = None

    if args.add_kmer_stats:
        if args.query_background:
            query_background_mh = background_kmer_minhash_from_sigfiles(
                args.query_background, moltype, args.ksize
            )
        else:
            notify(
                "Specified computing k-mer stats, but no --query-background file "
                "found, computing query background by merging all query signatures"
            )
            query_background_mh = background_kmer_minhash_from_sigfiles(
                inp_files, moltype, args.ksize
            )

        if args.db_background:
            # notify(f"Loading query background minhashes from {args.db_background}")
            db_background_mh = background_kmer_minhash_from_sigfiles(
                args.db_background, moltype, args.ksize
            )

        else:
            notify(
                "Specified computing k-mer stats, but no --db-background file "
                "found, computing database background by merging all db signatures"
            )
            db_background_mh = background_kmer_minhash_from_dbs(databases)

    # run SEARCH on all the queries.
    n = 0
    size_may_be_inaccurate = False
    aggregated_results = []
    for queryfile in inp_files:
        # load the query signature(s) & figure out all the things
        for query in sourmash_args.load_file_as_signatures(
            queryfile, ksize=args.ksize, select_moltype=moltype
        ):
            notify(
                f"loaded query: {str(query)[:30]}... (k={query.minhash.ksize}, "
                f"{sourmash_args.get_moltype(query)})"
            )

            # handle signatures with abundance
            if query.minhash.track_abundance:
                if args.ignore_abundance:
                    if query.minhash.track_abundance:
                        # abund sketch + ignore abundance => flatten sketch.
                        with query.update() as query:
                            query.minhash = query.minhash.flatten()
                elif args.containment or args.max_containment:
                    # abund sketch + keep abundance => no containment searches
                    notify(
                        "ERROR: cannot do containment searches on an abund signature; maybe specify --ignore-abundance?"
                    )
                    sys.exit(-1)
            else:
                # forcibly ignore abundances if query has no abundances
                args.ignore_abundance = True

            # do the actual search
            if query.minhash.track_abundance:
                try:
                    results = search_databases_with_abund_query(
                        query,
                        databases,
                        threshold=args.threshold,
                        do_containment=args.containment,
                        do_max_containment=args.max_containment,
                        best_only=args.best_only,
                        unload_data=True,
                    )
                except TypeError as exc:
                    error(f"ERROR: {str(exc)}")
                    sys.exit(-1)
            else:
                results = search_databases_with_flat_query(
                    query,
                    databases,
                    threshold=args.threshold,
                    do_containment=args.containment,
                    do_max_containment=args.max_containment,
                    best_only=args.best_only,
                    unload_data=True,
                    estimate_ani_ci=args.estimate_ani_ci,
                )

            # Add an e-value of the match based on the underlying distribution of
            # k-mers from the query signatures and databases
            for result in results:
                prob_overlap = get_prob_overlap(query, result)
                result.prob_overlap = prob_overlap
                result.weight = result.similarity / result.prob_overlap

                # Need to only do this once because the write_cols attribute is a c
                # lass attribute (set before __init__)
                if 'prob_overlap' not in result.write_cols and 'weight' not in result.write_cols:
                    result.write_cols.extend(["prob_overlap", "weight"])
                # import pdb; pdb.set_trace()
            # import pdb
            #
            # pdb.set_trace()

            n_matches = len(results)
            if args.best_only:
                args.num_results = 1

            if not args.num_results or n_matches <= args.num_results:
                print_results(
                    f"{len(results)} matches above threshold {args.threshold:0.3f}:"
                )
            else:
                print_results(
                    f"{len(results)} matches above threshold {args.threshold:0.3f}; showing first {args.num_results}:"
                )
                n_matches = args.num_results

            size_may_be_inaccurate = False
            jaccard_ani_untrustworthy = False

            # output!
            print_results("similarity    prob    weight   match")
            print_results("----------   ------   ------   -----")
            for sr in results[:n_matches]:
                pct_overlap = "{:.1f}%".format(sr.prob_overlap * 100)
                pct = "{:.1f}%".format(sr.similarity * 100)
                name = sr.match._display_name(60)
                print_results(
                    "{:>6}      {:>6}     {:.2f}        {}",
                    pct,
                    pct_overlap,
                    sr.weight,
                    name,
                )
                if sr.cmp_scaled is not None:
                    if not size_may_be_inaccurate and sr.size_may_be_inaccurate:
                        size_may_be_inaccurate = True
                    if not is_containment and sr.cmp.jaccard_ani_untrustworthy:
                        jaccard_ani_untrustworthy = True

            if args.best_only:
                notify("** reporting only one match because --best-only was set")

            output_csv = get_output_csv_name(args.output_dir, query)

            write_matches_csv(results, output_csv)

        # save matching signatures upon request
        if args.save_matches:
            notify(f'saving all matched signatures to "{args.save_matches}"')

            with SaveSignaturesToLocation(args.save_matches) as save_sig:
                for sr in results:
                    save_sig.add(sr.match)

        if picklist:
            sourmash_args.report_picklist(args, picklist)

        if size_may_be_inaccurate:
            notify(
                "WARNING: size estimation for at least one of these sketches may be inaccurate. ANI values will not be reported for these comparisons."
            )
        if jaccard_ani_untrustworthy:
            notify(
                "WARNING: Jaccard estimation for at least one of these comparisons is likely inaccurate. Could not estimate ANI for these comparisons."
            )


def get_output_csv_name(output_dir, query):
    query_filename = query.filename
    if not query_filename:
        # use md5sum if query.filename not properly set
        query_filename = query.md5sum()
    else:
        # Append first 8 chars of md5sum to query filename to deduplicate output
        query_filename += "." + query.md5sum()[:8]
    output_base = os.path.basename(query_filename)
    if output_dir:
        output_base = os.path.join(output_dir, output_base)
    output_csv = output_base + ".csv"
    return output_csv


def background_kmer_minhash_from_sigfiles(sigfiles, moltype: str, ksize: int):
    """Aggregate all signature files into one mega-signature to compute k-mer overlap statistics"""
    for i, sigfile in enumerate(sigfiles):
        notify(f"Loading query background minhashes from {sigfile}")

        for sig in sourmash.load_file_as_signatures(
            sigfile, select_moltype=moltype, ksize=ksize
        ):
            if i == 0:
                # Initialize the minhash object
                aggregated_mh = minhash_from_sig(sig)
            aggregated_mh.merge(sig.minhash)
    return aggregated_mh


def background_kmer_minhash_from_dbs(dbs):
    """Aggregate all signature files into one mega-signature to compute k-mer overlap statistics"""
    for i, db in enumerate(dbs):
        for sig in db.signatures():
            if i == 0:
                # Initialize the minhash object
                aggregated_mh = minhash_from_sig(sig)
            aggregated_mh.merge(sig.minhash)
    return aggregated_mh


def minhash_from_sig(sig):
    aggregated_mh = sourmash.MinHash(
        sig.minhash.num,
        sig.minhash.ksize,
        is_protein=sig.minhash.is_protein,
        dayhoff=sig.minhash.dayhoff,
        hp=sig.minhash.hp,
        track_abundance=sig.minhash.track_abundance,
        seed=sig.minhash.seed,
        max_hash=sig.minhash._max_hash,
    )
    return aggregated_mh


def write_matches_csv(found, output_csv):
    notify(f'saving all CSV matches to "{output_csv}"')
    w = None
    with FileOutputCSV(output_csv) as fp:
        for result in found:
            if w is None:
                w = result.init_dictwriter(fp)
            result.write(w)


def multisearch_parser():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--query", nargs="*", default=[], action="append", help="query signature"
    )
    p.add_argument(
        "--query-from-file", help="file containing list of signature files to query"
    )
    p.add_argument(
        "--databases",
        nargs="+",
        action="append",
        help="signatures/SBTs to search",
    )
    p.add_argument(
        "-q", "--quiet", action="store_true", help="suppress non-error output"
    )
    p.add_argument(
        "-d", "--debug", action="store_true", help="output debug information"
    )
    p.add_argument(
        "-t",
        "--threshold",
        metavar="T",
        default=0.08,
        type=float,
        help="minimum threshold for reporting matches; default=0.08",
    )
    p.add_argument(
        "--save-matches",
        metavar="FILE",
        help="output matching signatures to the specified file",
    )
    p.add_argument(
        "--best-only",
        action="store_true",
        help="report only the best match (with greater speed)",
    )
    p.add_argument(
        "-n",
        "--num-results",
        default=3,
        type=int,
        metavar="N",
        help="number of results to display to user; 0 to report all",
    )
    p.add_argument(
        "--containment",
        action="store_true",
        help="score based on containment rather than similarity",
    )
    p.add_argument(
        "--max-containment",
        action="store_true",
        help="score based on max containment rather than similarity",
    )
    p.add_argument(
        "--threshold-bp",
        metavar="REAL",
        type=float,
        default=5e4,
        help="threshold (in bp) for reporting results (default=50,000)",
    )
    p.add_argument(
        "--ignore-abundance",
        action="store_true",
        help="do NOT use k-mer abundances if present",
    )
    p.add_argument(
        "--estimate-ani-ci",
        action="store_true",
        help="also output confidence intervals for ANI estimates",
    )
    p.add_argument(
        "--fail-on-empty-database",
        action="store_true",
        help="stop at databases that contain no compatible signatures",
    )
    p.add_argument(
        "--no-fail-on-empty-database",
        action="store_false",
        dest="fail_on_empty_database",
        help="continue past databases that contain no compatible signatures",
    )
    p.set_defaults(fail_on_empty_database=True)
    p.add_argument(
        "--add-kmer-stats",
        action="store_true",
        help='Add statistical significance of k-mer matches based on background "'
        '"distribution of k-mers in query and database',
    )
    p.add_argument(
        "--query-background",
        nargs="*",
        help="Non-singleton (aggregated) signature of all k-mers in query. If not "
        "provided, computed by merging all query signatures and summing "
        "abundances",
    )
    p.add_argument(
        "--db-background",
        nargs="*",
        help="Non-singleton (aggregated) signature of all k-mers in databases. If not "
        "provided, computed by merging all database signatures and summing "
        "abundances",
    )
    p.add_argument(
        "--query-fasta",
        nargs="*",
        help="Non-singleton (aggregated) signature of all k-mers in query. If not "
        "provided, computed by merging all query signatures and summing "
        "abundances",
    )
    p.add_argument(
        "--db-fasta",
        nargs="*",
        help="Non-singleton (aggregated) signature of all k-mers in databases. If not "
        "provided, computed by merging all database signatures and summing "
        "abundances",
    )
    p.add_argument(
        "--output-dir",
        "--outdir",
        help="output CSV results to this directory",
    )
    add_ksize_arg(p)
    add_moltype_args(p)
    add_scaled_arg(p, 0)
    add_picklist_args(p)
    add_pattern_args(p)

    args = p.parse_args()
    return args, p


if __name__ == "__main__":
    sys.exit(main())
