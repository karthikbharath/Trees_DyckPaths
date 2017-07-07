/**
 * @author pkambadu
 */

#include <iostream>
#include <fstream>
#include <functional>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <utility>
#include <cstdio>
#include <map>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/filesystem.hpp>

#include "dyck_options.hpp"
#include "generate_Galton_Watson.hpp"
#include "dyck_path.hpp"

typedef generate_Galton_Watson_t::graph_type graph_type;
typedef generate_Galton_Watson_t::vertex_index_map_t vertex_index_map_t;
typedef generate_Galton_Watson_t::vertex_name_map_t vertex_name_map_t;
typedef generate_Galton_Watson_t::edge_weight_map_t edge_weight_map_t;
typedef generate_Galton_Watson_t::vertex_distance_map_t vertex_weight_map_t;
typedef boost::graph_traits<graph_type>::vertex_iterator vertex_iter_t;

static const double PI = 
    3.1415926535897932384626433832795028841971693993751058209749445923078164;

/**
 * Thin wrapper around boost::uniform<> so that it conforms to random_shuffle.
 */
template <typename EngineType>
struct random_functor_t : std::unary_function<int, int> {
  typedef EngineType engine_type;
  typedef boost::uniform_int<int> uni_int_type;
  engine_type engine;
  
  random_functor_t (engine_type engine) : engine(engine) {}

  int operator()(const int& ULIMIT) {
    uni_int_type dist (0, ULIMIT-1);
    return dist(engine);
  }
};

/**
 * Print the tree as a DOT file 
 */ 
template <typename Graph>
void print_tree (Graph& G, const dyck_options_t& options) {
#pragma omp critical
  if (options.print_tree) {
    /** Set the properties that you need to extract from the graph */
    boost::dynamic_properties dp;
    vertex_index_map_t index_map = boost::get(boost::vertex_index, G);
    vertex_name_map_t name_map = boost::get(boost::vertex_name, G);
    vertex_weight_map_t weight_map = boost::get(boost::vertex_distance, G);

    dp.property("index", index_map);
    dp.property("name", name_map);
    dp.property("weight", weight_map);

    /** Create the output */
    std::ostream* output; 
    if (0==strcmp("stdout", options.dot_out.c_str())) output = &(std::cout);
    else output = new std::ofstream(options.dot_out.c_str());

    boost::write_graphviz_dp(*output, G, dp, "index");

    if (0!=strcmp("stdout", options.dot_out.c_str())) delete output;
  }
}

/**
 * A function that lists all the DOT files in a directory.
 */
std::vector<std::string> get_dot_files (const std::string& dir_name) {
  std::vector<std::string> dot_file_list;

  if (boost::filesystem::exists(dir_name) && 
      boost::filesystem::is_directory(dir_name)) {
    boost::filesystem::directory_iterator end_iter;
    for(boost::filesystem::directory_iterator dir_iter(dir_name); 
        dir_iter != end_iter; 
        ++dir_iter) {
      if (boost::filesystem::is_regular_file(dir_iter->status())) {
        const boost::filesystem::path current_file = dir_iter->path();
        const boost::filesystem::path current_file_ext = 
                                             current_file.extension();
        if (current_file_ext.native() == ".dot") {
          dot_file_list.push_back(current_file.native());
        }
      }
    }
  }

  return dot_file_list;
}

/**
 * The one annoying problem is that OutputIterator's do not define
 * a value type, so it may not always work! Hrmp!
 */ 
int generate_seed () {
  int r_int;
#pragma omp critical
  {
    std::ifstream rf("/dev/random", std::ios::binary);
    if(rf.is_open())rf.read(reinterpret_cast<char*>(&r_int),sizeof(r_int));
    else r_int = std::time(NULL);
  }

  return r_int;
}

graph_type generate_one_graph (const dyck_options_t& options,
                               int& gen_seed) {
  graph_type G;

  gen_seed = (0>options.seed)?generate_seed():options.seed;
  if (0 == strcmp("Poisson", options.gen_method.c_str())) {
    G = generate_Galton_Watson_t::generate_Poisson(options.n,
                                                   gen_seed,
                                                   options.lambda,
                                                   options.use_random_weights,
                                                   options.verbosity,
                                                   options.use_stupid_tree_gen);
  } else if (0 == strcmp("Binomial", options.gen_method.c_str())) {
    G = generate_Galton_Watson_t::generate_Binomial (options.n,
                                                     gen_seed,
                                                     options.k,
                                                     options.p,
                                                     options.use_random_weights,
                                                     options.verbosity,
                                                 options.use_stupid_tree_gen);
  } else if (0 == strcmp("Geometric", options.gen_method.c_str())) {
    G = generate_Galton_Watson_t::generate_Geometric(options.n,
                                                     gen_seed,
                                                     options.p,
                                                     options.use_random_weights,
                                                     options.verbosity,
                                                 options.use_stupid_tree_gen);
  } else if (0 == strcmp("Binary-0-2", options.gen_method.c_str())) {
    if (false == (options.n & 0x1)) {
      std::cerr << "Binary-0-2 can only generate odd number of nodes\
 due to a shortcoming of Devroye's algorithm. Sorry" << std::endl;
      exit(-1);
    }
    G = generate_Galton_Watson_t::generate_Binary_0_2
                                              (options.n,
                                               gen_seed,
                                               options.p,
                                               options.use_random_weights,
                                               options.verbosity,
                                               options.use_stupid_tree_gen);
  } else if (0 == strcmp("Binary-0-1-2", options.gen_method.c_str())) {
    G = generate_Galton_Watson_t::generate_Binary_0_1_2
                                              (options.n,
                                               gen_seed,
                                               options.use_random_weights,
                                               options.verbosity,
                                               options.use_stupid_tree_gen);
  } else if (0 == strcmp("Graphviz", options.gen_method.c_str())) {
    /** Get the input file */
    std::ifstream input_dot_file(options.dot_in.c_str());

    /** Set the properties that you need to extract from the graph */
    boost::dynamic_properties dp;
    vertex_index_map_t index_map = boost::get(boost::vertex_index, G);
    vertex_name_map_t name_map = boost::get(boost::vertex_name, G);
    vertex_weight_map_t weight_map = boost::get(boost::vertex_distance, G);

    dp.property("index", index_map);
    dp.property("name", name_map);
    dp.property("weight", weight_map);

    /** Read the graph --- note that indices are relabeled, which is fine */
    boost::read_graphviz(input_dot_file, G, dp, "name");

    /** Label it */
    bfs_label(G);

  } else {
    std::cout << "Chosen method of graph generation not supported" 
              << std::endl;
  }
  return G;
}

template <typename NameToIndexMap, 
          typename WeightMap>
double path_length (const std::string& parent, 
                    const std::string& child,
                    const NameToIndexMap& name_map,
                    const WeightMap& weight_map,
                    const dyck_options_t& options) {

  double total_path_length = 0.0;

  if (0<options.verbosity) {
    std::cout << "Path length " << name_map.find(parent)->second << " ---> " 
              << name_map.find(child)->second << std::endl;
  }

  /** Get the first place where things mismatch */
  std::string::const_iterator position = 
    std::mismatch (parent.begin(), parent.end(), child.begin()).second;

  /** If they don't, we get to go one by one */
  while (position != child.end()) {
    const std::string current_vertex(child.begin(), position+1);
    const int current_vertex_index = name_map.find(current_vertex)->second;
    const double weight = boost::get(weight_map, current_vertex_index);
    if (0<options.verbosity) {
      std::cout << "Adding " << current_vertex_index << " (code:"
              << current_vertex << ") " 
              << " of weight " << weight << std::endl;
    }
    total_path_length += weight;
    ++position;
  }

  return total_path_length;
}

/** total path length and generate trees by cutting things off */

std::pair<int,double>
run_one_lca_experiment(dyck_options_t& options,
                       int& gen_seed,
                       int& lca_seed,
                       int& num_nodes_in_G) {
  graph_type G;
  vertex_name_map_t name_map;
  vertex_weight_map_t weight_map;

  /** Generate the seeds if needed */
  lca_seed = (0>options.lca_seed)?generate_seed():options.lca_seed;
  gen_seed = (0>options.seed)?generate_seed():options.seed;

  if (options.verbosity) std::cout << "Generating graph .... " << std::flush;

  do {
    /** Generate a graph */
    G = generate_one_graph (options, gen_seed);

    /** Get the number of nodes in the graph */
    vertex_iter_t v_begin, v_end;
    boost::tie(v_begin, v_end) = boost::vertices(G);
    num_nodes_in_G = std::distance(v_begin, v_end);

    if (false==options.dot_in.empty()) {
      options.n = num_nodes_in_G;
    }
  } while (num_nodes_in_G<options.n);
  if (options.verbosity) std::cout << "DONE" << std::endl;

  name_map = boost::get(boost::vertex_name, G);
  weight_map = boost::get(boost::vertex_distance, G);

  int num_lca_samples = num_nodes_in_G*options.percent_lca_samples;
  std::vector<int> lca_vertices(num_lca_samples);
  const std::string file_name = 
           boost::filesystem::path(options.dot_in).filename().string();

  /** 
   * We want to generate some vertices from the given set of vertices 
   * randomly and form the LCA tree. The thing wrong here is that my
   * compiler doesn't implement std::random_sample and std::iota. Hence
   * that long stupid way to sample elements from the vertex.
   *
   */
  std::vector<int> all_vertices(num_nodes_in_G-1);
  for (size_t i=0; i<all_vertices.size(); ++i) all_vertices[i]=(i+1);
  boost::mt19937 engine(lca_seed);
  random_functor_t<boost::mt19937> shuffle_prng (engine);
  std::random_shuffle (all_vertices.begin(), all_vertices.end(), shuffle_prng);

  if (options.sample_leaves) {
    /**
     * We want to only sample leaf nodes. The logic behind this is that we want
     * to only sample "original" data points. In case of Arvind's data, the 
     * intermediate nodes are formed by hierarchical clustering and hence are 
     * not original. There are two things to do:
     *
     * (1) Iterate through the shuffled vertices and 
     *    -- copy those vertices that are leaves (up to num_lca_samples)
     *    -- count the number of leaves
     * (2) If num_leaves>num_lca_samples, we are golden. However, if this is 
     *     not the case (num_lca_samples>=num_leaves), then we want to cut 
     *     down on the number of samples.
     */
    int num_leaves = 0;
    std::vector<int>::iterator lca_vertices_iter = lca_vertices.begin();
    for (size_t i=0; i<all_vertices.size(); ++i) {
      const int vertex_id = all_vertices[i];
      if (0 == boost::out_degree(vertex_id, G)) {
        if (num_leaves<num_lca_samples) { 
          *lca_vertices_iter = vertex_id;
          ++lca_vertices_iter; 
        }
        ++num_leaves;
      }
    }

    if (num_lca_samples>=num_leaves) {
      std::cerr << file_name << ": Sampling " << num_lca_samples 
                << " from " << num_leaves;
      num_lca_samples = 0.5*num_leaves;
      std::cerr << ". Changing to sample " << num_lca_samples << std::endl;
      if (0==num_lca_samples) {
        std::cerr << "Error: sampling " << num_lca_samples << std::endl;
        exit(-3);
      }
      lca_vertices.resize(num_lca_samples);
    }
  } else {
    /**
     * Approach 1:
     * We don't care whether we are sampling leaf or non-leaf nodes. As long as
     * it's a non-root node, we are golden.  All we are doing here is creating
     * an array of elements 1...(n-1), randomly shuffling it and then picking
     * the top num_lca_samples.
     */
    std::copy (all_vertices.begin(), 
               all_vertices.begin() + num_lca_samples,
               lca_vertices.begin());
  }

  /**
   * Now, build the new tree with intermediate vertices and figure out the 
   * number of edges in the new tree using the LCA algorithm to create 
   * intermediate nodes.
   */
  if (0<options.verbosity) {
    std::cout << "Sampled LCA vertices are: " << std::endl;
    for (size_t i=0; i<lca_vertices.size(); ++i) 
      std::cout << lca_vertices[i] << " (code:" 
              << boost::get(name_map, lca_vertices[i]) << ")" << std::endl;
  }

  std::map<std::string, int> all_vertices_name_map;

  /** Create a name map */
  for (size_t i=0; i<num_nodes_in_G; ++i)
    all_vertices_name_map[boost::get(name_map, i)] = i;

  /** 2. Join the vertices one by one to get the LCA tree */
  std::map<std::string,int> new_lca_vertices;
  typedef std::map<std::string, int>::const_iterator MapIterator;
  new_lca_vertices[boost::get(name_map, 0)] = 0;

  std::set<std::pair<int,int> > added_edges_set;
  double total_path_length = 0.0;
  int total_num_edges = 0; 

  for (size_t i=0; i<(lca_vertices.size()-1); ++i) {
    for (size_t j=(i+1); j<lca_vertices.size(); ++j) {
      /** Figure out the vertices we are dealing with */
      const int me = lca_vertices[i];
      const std::string my_name = boost::get(name_map, me);
      const int you = lca_vertices[j];
      const std::string your_name = boost::get(name_map, you);
      
      /** Push these two vertices into new LCA tree --- duplicate elim by map */
      new_lca_vertices[my_name] = me;
      new_lca_vertices[your_name] = you;
      
      /** Find the LCA, which also happens to be the common prefix */
      const std::string our_common_prefix = 
        std::string(my_name.begin(),
        std::mismatch(my_name.begin(), my_name.end(), your_name.begin()).first);
      
      /** Search for the common ancestor */
      const bool found_ancestor = (new_lca_vertices.end() != 
                            new_lca_vertices.find(our_common_prefix));
      
      if (false == found_ancestor) {
      const int our_common_ancestor = all_vertices_name_map[our_common_prefix];
        new_lca_vertices[our_common_prefix] = our_common_ancestor;
      }
      
      /** If the common prefix is not found, then insert it in */
      if (0<options.verbosity) {
        std::cout << "Ancestor of " << me << " (code:" << my_name << ") and " 
                  << you << " (code:" << your_name << ") is "
                  << new_lca_vertices[our_common_prefix] << " (code:" 
            << our_common_prefix << ")" << (found_ancestor?" (already present)": 
                                                    " (had to be added)")
            << std::endl;
      }
    }
  }

  if (options.sample_leaves) {
    std::cerr << file_name << ": " << num_lca_samples 
              << " lca-tree (|V|=" << new_lca_vertices.size() 
              << ", |E|=" << 2*num_lca_samples-1 << ")" << std::endl;
  }

  if (0<options.verbosity) {
    std::cout << "Sampled (and built) LCA tree has " << new_lca_vertices.size()
              << " nodes and these vertices: " << std::endl;
    for (MapIterator iter = new_lca_vertices.begin(); 
         iter != new_lca_vertices.end(); 
         ++iter) {
      std::cout << iter->second << " (code:" << iter->first << ")" << std::endl;
    }
  }

  const MapIterator range_first = new_lca_vertices.begin();
  MapIterator range_last = range_first; 
  ++range_last;
  MapIterator iter = range_last;
  while (iter != new_lca_vertices.end()) {
    /** find parent of *iter */
    MapIterator lowest_parent = range_first;
    MapIterator first = range_first;
    while (first != range_last) {
      if (0<options.verbosity) {
        std::cout << "Testing " << first->second << " for " 
                    << iter->second << std::endl;
      }
      if (first->first.end()==
          std::mismatch(first->first.begin(), 
                        first->first.end(),
                        iter->first.begin()).first) {
        if (0<options.verbosity) {
          std::cout << first->second << " is an ancestor of " 
                    << iter->second << std::endl;
        }
        lowest_parent = first;
      }
      ++first;
    }

    if (0<options.verbosity) {
      std::cout << lowest_parent->second << " is the LCA of " 
                << iter->second << std::endl;
    }

    /** compute distance from parent to *iter */
    total_path_length += path_length(lowest_parent->first, iter->first, 
                                    all_vertices_name_map, weight_map, options);
    ++total_num_edges;

    /** increment range_last */
    ++range_last;
    ++iter;
  }

  print_tree(G, options);

  return std::pair<int,double>(total_num_edges, total_path_length);
}

void run_one_experiment (const dyck_options_t& options,
                         int& gen_seed,
                         int& mle_seed,
                         int& num_nodes_in_G,
                         int& vertex_to_find,
                         double& vertex_height,
                         double& mean_height) {
  graph_type G;
  vertex_name_map_t name_map;
  vertex_weight_map_t weight_map;

  do {
    /** Generate a graph */
    G = generate_one_graph (options, gen_seed);

    /** Get the number of nodes in the graph */
    vertex_iter_t v_begin, v_end;
    boost::tie(v_begin, v_end) = boost::vertices(G);
    num_nodes_in_G = std::distance(v_begin, v_end);
  } while (num_nodes_in_G < (options.n));

  name_map = boost::get(boost::vertex_name, G);
  weight_map = boost::get(boost::vertex_distance, G);

  /** 
   * If we want to identify the height of a random (non-root) vertex in the
   * tree, then we want generate a random number between [1,n), where n is 
   * the number of nodes in the generated tree.
   */
  vertex_to_find = -1;
  mle_seed = (0>options.mle_seed)?generate_seed():options.mle_seed;
  boost::mt19937 engine (mle_seed);
  boost::uniform_int<int> dist (1, num_nodes_in_G-1);
  vertex_to_find = dist(engine);

  /** Now, compute the Dyck path using depth-first search */
  std::vector<double> y_axis;
  vertex_height = dyck_path (G, name_map, weight_map, 
                             std::back_inserter(y_axis), 
                             vertex_to_find, 
                             options.verbosity);

  /** Compute the mean height of the tree */
  assert (y_axis.size() == 2*num_nodes_in_G);
  mean_height = 0.0;
  for (size_t i=0; i<y_axis.size(); ++i) mean_height += y_axis[i];
  mean_height /= (2*pow(num_nodes_in_G,1.5));

  print_tree (G, options);

#pragma omp critical
  if (options.print_path) {
    /** Write stuff out where needed */
    /** We need to ensure that the x-axis is such that the slope is always 45*/
    std::vector<double> x_axis(y_axis.size());
    x_axis[0] = 0;
    for (size_t i=1; i<x_axis.size(); ++i)
      x_axis[i] = x_axis[i-1] + std::abs(y_axis[i] - y_axis[i-1]);

    /** Remember, DO NOT PRINT THE LAST INDEX -- this is an artifact of DFS */
    std::ostream* output; 
    if (0==strcmp("stdout", options.dyck_out.c_str())) output = &(std::cout);
    else output = new std::ofstream(options.dyck_out.c_str());
    for (size_t i=0; i<(x_axis.size()-1); ++i) {
      *output << x_axis[i] << "  " << y_axis[i] << std::endl;
    }
    if (0!=strcmp("stdout", options.dyck_out.c_str())) delete output;
  }

}

int main (int argc, char** argv) {
  dyck_options_t options (argc, argv);
  if (options.exit_on_return) { return -1; }

  if (1<options.verbosity) options.pretty_print();

  std::vector<std::string> dot_file_names;
  if (options.dot_in_dir != "") {
    dot_file_names = get_dot_files(options.dot_in_dir);
    options.num_trials = dot_file_names.size();
    options.gen_method = "Graphviz";
  } 

  std::vector<int> gen_seed_vec(options.num_trials);
  std::vector<int> mle_seed_vec(options.num_trials);
  std::vector<int> lca_seed_vec(options.num_trials);
  std::vector<int> num_nodes_vec(options.num_trials);
  std::vector<int> vertex_to_find_vec(options.num_trials);
  std::vector<double> vertex_height_vec(options.num_trials);
  std::vector<double> mean_height_vec(options.num_trials);

  if (options.measure_lca) {

    char header_string[1024];
    int h_count = sprintf (header_string, 
            "%9s %9s %2s %10s %2s %7s %7s %7s %5s %12s %12s",
                "dist-type", "f-name", "k", "p", "lm", "n",
                "n-tru", "n-smpl", "lca-edges", "lca-path-ln");
    std::cout << header_string << std::endl;

#pragma omp parallel for num_threads(options.num_threads)
    for (int i=0; i<options.num_trials; ++i) {
      /**
       * TODO: This is a hack. What we are doing is that we are checking if 
       * the input method is a DOT file. If so, we are checking if the input
       * mentioned is actually a directory listing of files. In this case, 
       * we will populate options.dot_in to be a different file name every 
       * single time.
       */
      if (false==options.dot_in_dir.empty()) {
        options.dot_in = dot_file_names[i];
      }
      std::pair<int,double> lca_result = 
        run_one_lca_experiment (options, gen_seed_vec[i], 
                                lca_seed_vec[i], num_nodes_vec[i]); 

      char value_string[1024];
      const std::string file_name = 
              boost::filesystem::path(options.dot_in).filename().string();
      int v_count = sprintf (value_string, 
               "%9s %9s %2d %.4e %2d %7d %7d %6.2f %12d %.6e",
                               options.gen_method.c_str(),
                               file_name.c_str(),
                               options.k,
                               options.p,
                               options.lambda,
                               options.n,
                               num_nodes_vec[i],
                               options.percent_lca_samples,
                               lca_result.first,
                               lca_result.second);
      std::cout << value_string << std::endl;
    }

    goto END;
  }

#pragma omp parallel for num_threads(options.num_threads)
  for (int i=0; i<options.num_trials; ++i) {
    /**
     * TODO: This is a hack. What we are doing is that we are checking if 
     * the input method is a DOT file. If so, we are checking if the input
     * mentioned is actually a directory listing of files. In this case, 
     * we will populate options.dot_in to be a different file name every 
     * single time.
     */
    if (false==options.dot_in_dir.empty()) {
      options.dot_in = dot_file_names[i];
    }
    run_one_experiment (options, gen_seed_vec[i], mle_seed_vec[i], 
                        num_nodes_vec[i], vertex_to_find_vec[i], 
                        vertex_height_vec[i], mean_height_vec[i]);
  }

  double true_var;
  if (0 == strcmp("Poisson", options.gen_method.c_str())) 
    true_var=options.lambda;
  else if (0 == strcmp("Binomial", options.gen_method.c_str())) 
    true_var=options.k*options.p*(1-options.p);
  else if (0 == strcmp("Geometric", options.gen_method.c_str()))
    true_var=(1-options.p)/(options.p*options.p);
  else if (0 == strcmp("Binary-0-2", options.gen_method.c_str()))
    true_var=1;
  else if (0 == strcmp("Binary-0-1-2", options.gen_method.c_str()))
    true_var=2./3.;

  if (options.dump_numbers) {
    for (int i=0; i<options.num_trials; ++i) {
      double x_i = (1.0/(sqrt((double)num_nodes_vec[i])))*vertex_height_vec[i];
      printf (" %-8.3f", x_i);
    }
    printf ("\n");
    for (int i=0;i<options.num_trials;++i) {
      printf(" %-8.3f", mean_height_vec[i]);
    }
    printf ("\n");
  }

  if (options.measure_mle) {
    if (options.verbosity) {
      char header_string[1024];
      int h_count = sprintf (header_string, 
            "%9s %9s %2s %10s %2s %7s %7s %7s %7s %12s %12s %12s", 
                "dist-type", "filename", "k", "p", "lm", "n",
                "n-tru", "node-#","height-raw","height-nrmzd","height-mean");
      if (1<options.verbosity)
   sprintf (&(header_string[h_count]), " %10s %10s", "gen-seed", "MLE-seed");
      std::cout << header_string << std::endl;
    }

    double sum_of_x_i_sqr = 0.0;
    double mean_of_x_i_sqr = 0.0;
    for (int i=0; i<options.num_trials; ++i) {
      double x_i = (1.0/(sqrt((double)num_nodes_vec[i])))*vertex_height_vec[i];
      sum_of_x_i_sqr += x_i*x_i;
      mean_of_x_i_sqr += x_i;

      if (options.verbosity) {
        char value_string[1024];
        const std::string file_name = 
              boost::filesystem::path(dot_file_names[i]).filename().string();
        int v_count = sprintf (value_string, 
               "%9s %9s %2d %.4e %2d %7d %7d %7d %.6e %.6e %.6e",
                               options.gen_method.c_str(),
                               file_name.c_str(),
                               options.k,
                               options.p,
                               options.lambda,
                               options.n,
                               num_nodes_vec[i],
                               vertex_to_find_vec[i],
                               vertex_height_vec[i],
                               x_i,
                               mean_height_vec[i]);
        if (1<options.verbosity)
          sprintf (&(value_string[v_count]), " %10u %10u", 
                          gen_seed_vec[i], mle_seed_vec[i]);
        std::cout << value_string << std::endl;
      }
    }
    mean_of_x_i_sqr /= options.num_trials;
    mean_of_x_i_sqr *= mean_of_x_i_sqr;

    double mle_estimate_var = (2*options.num_trials)/sum_of_x_i_sqr;
    double ano_estimate_var = PI/(2*mean_of_x_i_sqr);

    std::cout << "True variance = " << true_var
              << " MLE = " << mle_estimate_var 
              << " Another = " << ano_estimate_var << std::endl;
    std::cout << " Mean of x_i^2 =  " << mean_of_x_i_sqr << std::endl;
  }

  if (options.test_confidence) {
    int num_batches = options.num_trials / options.batch_size;
    int num_contained = 0;

    int num_nodes = num_nodes_vec[0];
    bool experimental_error = false;
    for (int i=1; i<options.num_trials; ++i) {
      if (num_nodes != num_nodes_vec[i]) {
        experimental_error = true;
        break;
      }
    }

    if (experimental_error) {
      std::cout << "Run this experiment again, please" << std::endl;
      goto END;
    }

    std::vector<double> other_variances;
    if (0 == strcmp("Poisson", options.gen_method.c_str())) {
      other_variances.push_back(1 - 1./2.);
      other_variances.push_back(1 - 1./3.);
      other_variances.push_back(1 - 1./4.);
      other_variances.push_back(1 - 1./5.);
      other_variances.push_back(2);
    } else if (0 == strcmp("Binomial", options.gen_method.c_str())) {
      other_variances.push_back(1);
      other_variances.push_back(2);
      other_variances.push_back(2.0/3.0);
    } else if (0 == strcmp("Geometric", options.gen_method.c_str())) {
      other_variances.push_back(1);
      other_variances.push_back(1 - 1./2.);
      other_variances.push_back(1 - 1./3.);
      other_variances.push_back(1 - 1./4.);
      other_variances.push_back(1 - 1./5.);
    } else if (0 == strcmp("Binary-0-2", options.gen_method.c_str())) {
      other_variances.push_back(1 - 1./2.);
      other_variances.push_back(1 - 1./4.);
      other_variances.push_back(1 - 1./5.);
      other_variances.push_back(2);
      other_variances.push_back(2.0/3.0);
    } else if (0 == strcmp("Binary-0-1-2", options.gen_method.c_str())) {
      other_variances.push_back(1 - 1./2.);
      other_variances.push_back(1 - 1./3.);
      other_variances.push_back(1 - 1./4.);
      other_variances.push_back(1 - 1./5.);
      other_variances.push_back(1);
      other_variances.push_back(2);
    }
    std::vector<double> power_against_others (other_variances.size(), 0.0);

    boost::math::chi_squared chi_dist (2*options.batch_size);
    const double a = boost::math::quantile(chi_dist, 0.025);
    const double b = boost::math::quantile(chi_dist, 0.975);
    
    for (int i=0; i<options.num_trials; i+=options.batch_size) {
      double sum_of_x_j_sqr = 0.0;
      for (int j=0; j<options.batch_size; ++j) {
        double x_j = (1.0/(sqrt((double)num_nodes_vec[i+j])))*
                                        vertex_height_vec[i+j];
        sum_of_x_j_sqr += x_j*x_j;
      }
      
      const double A = a/sum_of_x_j_sqr;
      const double B = b/sum_of_x_j_sqr;
      if (A <= true_var && B >= true_var) ++num_contained;
      for (int l=0; l<other_variances.size(); ++l)
        if (A <= other_variances[l] && B >= other_variances[l]) 
          power_against_others[l]+=1;

      if (1<options.verbosity) {
        std::cout << "Lower = " << A << " Upper = " << B  
                  << " sum_of_x_j_sqr = " << sum_of_x_j_sqr << std::endl;
      }
    }
    const double contained = 
       100*static_cast<double>(num_contained)/static_cast<double>(num_batches);

    printf ("%12s %5d %3d %7.3f%%(tru-var=%4.3f) ", options.gen_method.c_str(),
                                                    num_batches,
                                                    options.batch_size,
                                                    contained,
                                                    true_var);
    for (int l=0; l<other_variances.size(); ++l) 
      printf (" %7.3f%%(var=%4.3f)", 100*(power_against_others[l]/num_batches),
                                     other_variances[l]);
    printf ("\n");
  }

END:

  return 0;
}
