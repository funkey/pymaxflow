#ifndef PYMAXFLOW_MAXFLOW_GRAPH_H__
#define PYMAXFLOW_MAXFLOW_GRAPH_H__

#include <lemon/list_graph.h>
#include <solver/SolverFactory.h>

template <typename CapacityType, typename FlowType>
class Graph {

public:

	typedef std::size_t NodeId;

	/**
	 * Terminal nodes.
	 */
	typedef enum {

		Source = 0,
		Sink   = 1

	} TermType;

	/**
	 * Create a new Graph for max-flow/min-cut computation for the given 
	 * estimated number of nodes and edges. More nodes and edges can be added 
	 * later, but might require re-allocation.
	 */
	Graph(std::size_t node_num_max, std::size_t edge_num_max);

	/**
	 * Add n nodes to the graph, return the index to the first one.
	 */
	NodeId add_nodes(std::size_t n);

	/**
	 * Add a bidrectional edge between u and v, with cap being the capacity from 
	 * u->v, and rev_cap from v->u.
	 */
	void add_edge(NodeId u, NodeId v, CapacityType cap, CapacityType rev_cap);

	/**
	 * Add edges source->n and n->sink with the given capacities.
	 */
	void add_tweights(NodeId n, CapacityType cap_source, CapacityType cap_sink);

	/**
	 * Compute the max-flow.
	 */
	FlowType maxflow();

	/**
	 * Return the segment to which node n belongs after computation of the 
	 * max-flow/min-cut. If there are several min-cuts and the node can belong 
	 * to either source or sink, default_segm will be returned.
	 */
	TermType what_segment(NodeId n, TermType default_segm = Source);

private:

	typedef lemon::ListDigraph GraphType;

	GraphType _graph;

	GraphType::NodeMap<CapacityType> _source_caps;
	GraphType::NodeMap<CapacityType> _sink_caps;
	GraphType::ArcMap<CapacityType>  _caps;

	std::size_t _num_nodes;

	std::unique_ptr<LinearSolverBackend> _solver;
};

template <typename CapacityType, typename FlowType>
Graph<CapacityType, FlowType>::Graph(std::size_t node_num_max, std::size_t edge_num_max) :
	_source_caps(_graph),
	_sink_caps(_graph),
	_caps(_graph),
	_num_nodes(0) {

	SolverFactory factory;
	_solver = std::unique_ptr<LinearSolverBackend>(factory.createLinearSolverBackend());
}

#endif // PYMAXFLOW_MAXFLOW_GRAPH_H__

