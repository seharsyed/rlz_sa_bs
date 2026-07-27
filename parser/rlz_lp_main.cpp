#include "common.hpp"
#include "rlz_lp.hpp"

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);
    const LoadedData data = load_data(args);

    RLZLPParser<Symbol, SAType> parser(
        data.ref,
        data.sa,
        args.p
    );

    print_parser_ready(parser.cache_info());

    print_run_details(args, data);
    print_result_header();

    const RunResults results =
        run_rlz(data, parser);

    print_lp_summary(
        results,
        parser.cache_info()
    );

    return 0;
}