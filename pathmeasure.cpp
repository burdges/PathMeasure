
//
// usage : allpaths-montecarlo list/[graph].csv > [graph].c
// compile : g++ -W -Wall -pedantic -O3 -funroll-loops -lboost_graph -s allpaths-montecarlo.cpp -g -o allpaths-montecarlo
//
// Linux profiling :
// g++ -W -Wall -pedantic -O3 -funroll-loops -lboost_graph allpaths-montecarlo.cpp -g -pg -c -o allpaths-montecarlo.o
// g++ -g -pg -o allpaths-montecarlo allpaths-montecarlo.o
// objdump -t allpaths-montecarlo 
// ./allpaths-montecarlo 10 <list/random.csv >t
// gprof ./allpaths-montecarlo | less
//
// Mac OS X profiling :
// instruments -t /Developer/Applications/Instruments.app/Contents/Resources/templates/Time\ Profiler.tracetemplate ./allpaths-montecarlo 
//    after modifying allpaths-montecarlo.cpp to load list/random.csv and limiting the path length
//
//	recall : 
//		/opt/local/bin/g++-mp-4.4 
//  	g++ -I/opt/local/include -L/opt/local/lib 



#include <stdlib.h>
#include <getopt.h> // getlongopt
// #include <unistd.h> // getopt

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

#include <climits>
#include <ctime>
#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/filtered_graph.hpp>
// #include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
// #include <boost/graph/johnson_all_pairs_shortest.hpp>
// #include <boost/graph/floyd_warshall_shortest.hpp>

using namespace boost;  // Ain't no boost::graph  :(


// Ick Global Variables //

static int verbosity = 1;


// Bernoulli Trials Generator //

#include "randomc.h"

class Bernoulli : public CRandomMersenne {
public:
	static const unsigned int_max = INT_MAX; // 32768-1; // 2147483648-1;
	double p;
	Bernoulli(double _p,int seed = (int)time(0)) :
		CRandomMersenne(seed), p(_p) { }
	bool trial(double _p)
		{ return IRandom(0,int_max) < _p*int_max;  }
	bool trial()
		{  return this->trial(p);  }
};


// Graph Types //

struct vertex_properties { string id; };  // float p;
typedef property<edge_index_t, size_t, property<edge_weight_t, int> > edge_properties;
// struct edge_properties {  size_t edge_index;  int edge_weight;  };
typedef property<graph_name_t, string> graph_properties;

typedef adjacency_list < vecS, vecS, undirectedS,
  vertex_properties, edge_properties, graph_properties > raw_graph_t;
typedef subgraph< raw_graph_t > graph_t;

typedef graph_traits<graph_t>::edge_descriptor edge_descriptor_t;
typedef graph_traits<graph_t>::vertex_descriptor vertex_discripter_t;


// Algorithm //

typedef ublas::symmetric_matrix<double> matrix_t;

// Ideally, this template should handle a non-square symmetric_matrix
// using upper and lower.
template <typename T> inline T norm_max(ublas::symmetric_matrix<T>& M) {
	size_t n1 = M.size1();   size_t n2 = M.size2();   T m = 0;
	for (size_t i=0; i<n1; i++)  for (size_t j=i; j<n2; j++)
		m = max(M(i,j),m);
	return m;
}

template <typename T> inline T norm_max(ublas::matrix<T>& M) {
	size_t n1 = M.size1();   size_t n2 = M.size2();   T m = 0;
	for (size_t i=0; i<n1; i++)  for (size_t j=0; j<n2; j++)
		m = max(M(i,j),m);
	return m;
}

struct bernoulli_edge_predicate {
	typedef graph_traits<graph_t>::vertices_size_type vertex_size_type;
	typedef ublas::symmetric_matrix<int> bool_matrix;

	Bernoulli *bernoulli;
	graph_t *G;
	bool_matrix edge_matrix;

	bernoulli_edge_predicate() { }
	bernoulli_edge_predicate(graph_t& g, Bernoulli& b) :
		bernoulli(&b), G(&g), edge_matrix(num_vertices(g))  { }

	bool_matrix::reference
	edge(const edge_descriptor_t& e) {
		vertex_size_type i = get(vertex_index, *G, source(e,*G));
		vertex_size_type j = get(vertex_index, *G, target(e,*G));
		return edge_matrix(i,j);
	}	// cout << '(' << i << ',' << j << ')' << endl;
	bool_matrix::const_reference
	edge(const edge_descriptor_t& e) const {
		vertex_size_type i = get(vertex_index, *G, source(e,*G));
		vertex_size_type j = get(vertex_index, *G, target(e,*G));
		return edge_matrix(i,j);
	}	// cout << '(' << i << ',' << j << ')' << endl;

	void trial(double p) {
		graph_traits<graph_t>::edge_iterator e, e_end;
		for (tie(e, e_end) = edges(*G); e != e_end; ++e)
			edge(*e) = bernoulli->trial(p)? 1:0;
	}
	void trial()  { trial(bernoulli->p); }

	bool operator()(const edge_descriptor_t& e) const
		{  return (bool)edge(e);  }
};

struct run_trials_paramaters {
	int algorithm,seed;
	int trials, reports;
	double p,precision;

	run_trials_paramaters(int s = (int)time(0), double _p = -1) {
		p = _p;  seed = s;
		algorithm = 0;
		trials = 1000;
		reports = 100;
		precision = 0; // 0.00000001;
	}
};

matrix_t run_trials(graph_t& G, struct run_trials_paramaters params) {
	int n = num_vertices(G);
	matrix_t M(n);
	switch (n) {
		case 2:  M(1,1)=1;  M(0,1) = params.p;  M(1,0)=params.p;
		case 1:  M(0,0) = 1;  return M;
		case 0:  cerr << "run_trials cannot process an empty graph\n"; abort();
	}

	Bernoulli bernoulli(params.p,params.seed);
	bernoulli_edge_predicate edge_pred(G, bernoulli);
	
	ublas::symmetric_matrix<unsigned long> L = ublas::zero_matrix<unsigned long>(n);
	matrix_t P(n);  int t,nc=0;
	if (verbosity)
		{ cerr << "Round\t max. norm\n"; }
	for (t=0; t<params.trials; t++) {
		edge_pred.trial();
		filtered_graph<graph_t, bernoulli_edge_predicate> F(G, edge_pred);
		vector<int> components(n);
		nc += connected_components(F, &components[0]);
		for (int i=0; i<n; i++)  for (int j=i; j<n; j++)
			if (components[i] == components[j])
				++L(i,j);
		if (params.reports && !(t % params.reports)) {
			M = (matrix_t)L / (double)t;
			matrix_t Q = M - P;
			double norm = norm_max(Q);
			if (t > params.reports) {
				if (params.precision && norm < params.precision) return M;
				if (verbosity)
					{ cerr << '#' << t << "\t" << norm << endl; }
			}
			P = M;
		}
	}
	M = (matrix_t)L / (double)t;
	return M;
}

/*
 *		johnson_all_pairs_shortest_paths(g,d);  // sparce
 *		floyd_warshall_all_pairs_shortest_paths(g,d);  // dense
 */


// Files //

inline void add_vertex_once(vertex_discripter_t v, graph_t& g)
	{  if (! g.find_vertex(v).second)  add_vertex(v,g);  }

size_t biconnected_components(graph_t& g, vector<graph_t*>& bicomponents) {
	map< edge_descriptor_t, unsigned > bicomponent_edges;
	vector<vertex_discripter_t> articulation_points;
	articulation_points.reserve(num_vertices(g));

	int n = biconnected_components(g,
			make_assoc_property_map(bicomponent_edges),
			back_inserter(articulation_points) ).first;
	int o = bicomponents.size();
	for (int i=0; i<n; i++)
		bicomponents.push_back( &(g.create_subgraph()) );
	graph_traits<graph_t>::edge_iterator e, e_end;
	for (tie(e, e_end) = edges(g); e != e_end; ++e) {
		add_vertex_once( source(*e,g), *bicomponents[o+bicomponent_edges[*e]] );
		add_vertex_once( target(*e,g), *bicomponents[o+bicomponent_edges[*e]] );
	}
	// cout << "Graph #" << ... << " : " << get(graph_name, g);
	// cout << " |V|=" << num_vertices(g) << endl;
	// cout << " |Component|=" << bicomponents.size() << endl;

	vector<vertex_discripter_t>::iterator v,v_end;
	cout << "Art.pts.\tComponents\n";
	// cout << articulation_points.size() << endl;
	for (v=articulation_points.begin(), v_end=articulation_points.end();
			v != v_end; v++) {
		cout << get(&vertex_properties::id, g, *v) << "\t\t";
		for (int i=0; i<n; i++)
			if ((bicomponents[o+i]->find_vertex(*v)).second)
				cout << ' ' << o+i << ' ';
		cout << endl;
	}
	return n;
}

bool dot_load_istream(istream& in, vector<graph_t>& graphs, char *fn = "") {
	graph_t g(0, string(fn) );
	dynamic_properties dp;
	dp.property("id", get(&vertex_properties::id, g) );
	if (! read_graphviz(in,g,dp,"id") )  return false;
	graphs.push_back(g);
	return true;
}

void output_matrix(ostream& out, matrix_t& M, unsigned n) {
	out << "matrix_t M[" << n << "] = {\n";
	for (unsigned i=0; i<M.size1(); ++i) {
		out << "\t{ ";
		out << M(i,0);
		for (unsigned j=1; j<M.size2() ; ++j)
			out << ',' << M(i,j);
		out << " },\n";
	}
	out << "}\n";
}

// void merge_trials(vector<graph_t*>& bicomponents, vector<matrix_t>& distances) 


// Main //

void usage() {
	cerr << "Usage : allpaths-montecarlo -p=<probability> [-s=<seed>] [-t=<trials>] [-f=<reports>/-q] <graph>" << endl;
	abort();
}

int main(int argc, char **argv) {
	vector<graph_t> graphs;
	vector<graph_t*> bicomponents;
	struct run_trials_paramaters params((int)time(0));

	static struct option long_options[] = {
		{"johnson",			no_argument, &params.algorithm, 1}, // sparse 
		{"floydwarshall",	no_argument, &params.algorithm, 2}, // dense
		{"trials",		required_argument, 0, 'q'},
		{"quiet",		no_argument, 0, 'q'},
		{"help", 		no_argument, 0, 'h'},
		{"precision",	required_argument, 0, 1},
	};

	if (argc == 1)  usage();
	while (1) {
		int option_index = 0;
		int c = getopt_long (argc, argv, "p:s:t:f:qh?", long_options, &option_index);
		if (c == -1)  break;
		char* arg=optarg;
		if (*arg == '-') do { ++arg; } while(*arg && *arg != '=');
		switch (c) { 
			case 'p':  params.p = atof(arg);  break;
			case 's':  params.seed = atoi(arg);  break;
			case 't':  params.trials = atoi(arg);  break;
			case 'f':  params.reports = atoi(arg);  break;
			case 'q':  verbosity = 0;  break;
			case 'h':  case '?':  default:  usage();
			case 1:  params.precision = atof(arg);  break;
		}
	}
	if (params.p <= 0 || params.p > 1)  usage();
	if (optind == argc)
		dot_load_istream(cin,graphs,(char*)"");
	while (optind < argc) {
		ifstream in;
		in.open(argv[optind], ifstream::in);
		dot_load_istream(in,graphs,argv[optind]);
		in.close();
		++optind;
	}
	for (unsigned i=0; i<graphs.size(); i++)
		biconnected_components(graphs[i],bicomponents);
	for (unsigned i=0; i<bicomponents.size(); i++) {
		matrix_t M;
		M = run_trials(*bicomponents[i],params);
		output_matrix(cout,M,i);
	}
//	if (bicomponents.size())
//		merge_trials(bicomponents);
}

