#include "common.hpp"
#include "rlz_lp.hpp"

int main(int argc, char** argv) {
    // Read the command-line options.
    const Args args = parse_args(argc, argv);

    // Load the input list, reference, and suffix array.
    const LoadedData data = load_data(args);

    // Create the linear-probing RLZ parser.
    RLZLPParser<Symbol, SAType> parser(data.ref, data.sa, args.div_p);

    // Run RLZ on all input files and collect the results.
    const RunResults results = run_rlz(data, parser);

    // Print the run details and per-file results.
    print_results(args, data, results);

    // Print the final timing and cache summary.
    print_lp_summary(results, parser.cache_info());

    return 0;
}