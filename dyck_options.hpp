#ifndef DYCK_OPTIONS_HPP
#define DYCK_OPTIONS_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <functional>
#include <utility>
#include <exception>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

namespace po = boost::program_options;

/**
 * Additional command line parser which interprets '@something' as a
 * option "config-file" with the value "something".
 */ 
std::pair<std::string, std::string> at_option_parser(std::string const& s) {
  if ('@' == s[0]) {
    return std::make_pair (std::string("response-file"), s.substr(1));
  } else { 
    return std::make_pair (std::string(), std::string()); 
  }
}

/**
 * A structure that is used to pass options to the DYCK solver. This structure
 * has default values embedded. No accessor functions are being written.
 */
struct dyck_options_t {
  /** Integer options */
  int n;
  int seed;
  int mle_seed;
  int lca_seed;
  int verbosity;
  int k;
  int lambda;
  int num_trials;
  int num_threads;
  int batch_size;

  /** Boolean options */
  bool use_random_weights;
  bool print_tree;
  bool print_path;
  bool measure_mle;
  bool measure_lca;
  bool test_confidence;
  bool dump_numbers;
  bool sample_leaves;
  bool use_stupid_tree_gen;

  /** Floating point options */
  double p;
  double l;
  double percent_lca_samples;

  /** String options */
  std::string gen_method;
  std::string dot_in;
  std::string dot_in_dir;
  std::string dot_out;
  std::string dyck_out;

  /** A parameter indicating if we need to continue or not */
  bool exit_on_return;

  /**
   * The constructor takes in all the command line parameters and parses them.
   */ 
  dyck_options_t (int argc, char** argv) : exit_on_return(false) {
    /** Set up the options that we want */
    po::options_description desc ("Dyck Path Program Options");
    desc.add_options()
  ("help", "produce a help message")
  ("n", po::value<int>(&n)->default_value(10), 
        "Number of nodes in the tree")
  ("seed", po::value<int>(&seed)->default_value(-1), "Random number seed")
  ("mle-seed", po::value<int>(&mle_seed)->default_value(-1), 
               "Random number seed for uniformly picked vertex for MLE")
  ("lca-seed", po::value<int>(&lca_seed)->default_value(-1), 
               "Random number seed for uniformly picked vertices for LCA")
  ("verbosity", po::value<int>(&verbosity)->default_value(1), 
   "Verbosity level for error messages")
  ("k", po::value<int>(&k)->default_value(4), "Number of Bernoulli trials")
  ("lambda", po::value<int>(&lambda)->default_value(10), "Mean for Poisson")
  ("num-trials", po::value<int>(&num_trials)->default_value(1), 
                              "Number of trials to run for this experiment")
  ("num-threads", po::value<int>(&num_threads)->default_value(4), 
                              "Number of threads to use for this experiment")
  ("batch-size", po::value<int>(&batch_size)->default_value(10), 
                              "Batch size for testing confidence")
  ("use-random-weights", po::value<bool>(&use_random_weights)
                         ->default_value(false),
   "Should we have random edge weights (between 0 and 1): false)")
  ("print-tree", po::value<bool>(&print_tree)
                       ->default_value(false),
   "Should we print the tree as a GraphViz graph: false)")
  ("print-path", po::value<bool>(&print_path)
                       ->default_value(false),
   "Should we print the Dyck path of the tree: false)")
  ("measure-mle", po::value<bool>(&measure_mle)
                       ->default_value(false),
   "Should we measure the MLE of a randomly picked vertex in the tree: false)")
  ("measure-lca", po::value<bool>(&measure_lca)
                       ->default_value(false),
   "Should we measure the LCA of a randomly formed tree: false)")
  ("test-confidence", po::value<bool>(&test_confidence)
                       ->default_value(false),
   "Should we measure the CI a randomly picked vertex in the tree: false)")
  ("dump-numbers", po::value<bool>(&dump_numbers)
                       ->default_value(false),
   "Should we dump mean heights and sums of squares: false)")
  ("sample-leaves", po::value<bool>(&sample_leaves)
                       ->default_value(true),
   "Should we only sample the leaf nodes and not internal nodes (LCA): true")
  ("use-stupid-tree-gen", po::value<bool>(&use_stupid_tree_gen)
                       ->default_value(false),
  "Should we use the stupid algorithm to generate unconditioned trees:false)")
  ("p", po::value<double>(&p)->default_value(0.01),
   "Probability of success for the Bernoulli trials (default: 0.01)")
  ("l", po::value<double>(&l)->default_value(0.0),
   "Probability required for the new Binary-0-2 (default: 0.01)")
  ("percent-lca-samples", 
                  po::value<double>(&percent_lca_samples)->default_value(.4), 
                   "Percentange of original nodes to sample for LCA tree")
  ("gen-method", po::value<std::string>(&gen_method)->default_value("Poisson"),
   "Distribution to use for generating the tree (default:Poisson)")
  ("dot-in", po::value<std::string>(&dot_in)->default_value(""),
   "File name to read graph from --- used with \"Graphviz\" gen-method")
  ("dot-in-dir", po::value<std::string>(&dot_in_dir)->default_value(""),
   "File listings to read graphs from --- used with \"Graphviz\" gen-method")
  ("dot-out", po::value<std::string>(&dot_out)->default_value("stdout"),
   "File name to print the tree GraphViz out")
  ("dyck-out", po::value<std::string>(&dyck_out)->default_value("stdout"),
   "File name to print the tree Dyck path out")
  ("response-file", po::value<std::string>(), 
                          "can be specified with '@name', too")
      ; /* end options */

    /** create a variable map to hold all these things */
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc)
              .extra_parser(at_option_parser).run(), vm);
 
    /** Print help and return if needed */
    if (vm.count ("help")) { 
      std::cout << desc;
      exit_on_return = true;
      return;
    } 
    
    try {
      if (vm.count("response-file")) {
        /** Load the file and tokenize it */
        std::ifstream ifs (vm["response-file"].as<std::string>().c_str());
        if (false==ifs) {
          std::cout << "Could not open response file:" 
                    << vm["response-file"].as<std::string>() << std::endl;
          exit_on_return = true;
          return;
        }
      
        /** Read the whole file into a string */
        std::stringstream ss;
        ss << ifs.rdbuf ();
        boost::char_separator<char> sep(" \n\r");
        std::string sstr(ss.str());
        boost::tokenizer<boost::char_separator<char> > tok(sstr, sep);
        std::vector<std::string> args;
        std::copy (tok.begin(), tok.end(), std::back_inserter(args));
      
        /** Parse the file and store the options */
        po::store(po::command_line_parser(args).options(desc).run(), vm);
      }
      
      po::notify(vm);
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
      exit_on_return = true;
    }

    /** Added new for Karthik */
    if (0.0 == l) l = sqrt(p/(1-p));

    /** Now, check for error conditions */
    if (0 == strcmp("Poisson", gen_method.c_str())) {
    } else if (0 == strcmp("Binomial", gen_method.c_str())) {
    } else if (0 == strcmp("Geometric", gen_method.c_str())) {
    } else if (0 == strcmp("Binary-0-2", gen_method.c_str())) {
    } else if (0 == strcmp("Binary-0-1-2", gen_method.c_str())) {
    } else if (0 == strcmp("Graphviz", gen_method.c_str())) {
      if (0==strcmp("", dot_in.c_str()) && 0==strcmp("", dot_in_dir.c_str())){
        std::cout << "Graphviz input method requires dot file" << std::endl;
      }
    } else {
      std::cout << "Graph input method not supported" << std::endl;
    }
  }

  void pretty_print () const {
    std::cout << "n = " << n << std::endl;
    std::cout << "seed = " << seed << std::endl;
    std::cout << "MLE-seed = " << mle_seed << std::endl;
    std::cout << "verbosity = " << verbosity << std::endl;
    std::cout << "k = " << k << std::endl;
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "l = " << l << std::endl;
    std::cout << "num-trials = " << num_trials << std::endl;
    std::cout << "num-threads = " << num_threads << std::endl;
    std::cout << "batch-size = " << batch_size << std::endl;
    std::cout << "gen-method = " << gen_method << std::endl;
    std::cout << "Use random weights = " << use_random_weights << std::endl;
    std::cout << "Print tree = " << print_tree << std::endl;
    std::cout << "Print Dyck path = " << print_path << std::endl;
    std::cout << "Estimate MLE measure = " << measure_mle << std::endl;
    std::cout << "Estimate confidence = " << test_confidence << std::endl;
    std::cout << "Tree file = " << dot_out << std::endl;
    std::cout << "Path file = " << dyck_out << std::endl;
  }
};

#endif // DYCK_OPTIONS_HPP
