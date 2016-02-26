#ifndef PYMAXFLOW_MAXFLOW_GRAPH_H__
#define PYMAXFLOW_MAXFLOW_GRAPH_H__

#include <lemon/list_graph.h>
#include <solver/SolverFactory.h>
#include <util/helpers.hpp>

/**
	* Terminal nodes.
	*/
enum TerminalType {

	Source = 0,
	Sink   = 1

};

template <typename CapacityType, typename FlowType>
class Graph {

public:

	typedef std::size_t NodeId;

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
	void add_tedge(NodeId n, CapacityType cap_source, CapacityType cap_sink);

	/**
	 * Compute the max-flow.
	 */
	FlowType maxflow();

	/**
	 * Return the segment to which node n belongs after computation of the 
	 * max-flow/min-cut.
	 */
	TerminalType what_segment(NodeId n);

private:

	typedef lemon::ListDigraph GraphType;

	GraphType _graph;

	GraphType::Node _source;
	GraphType::Node _sink;

	GraphType::ArcMap<CapacityType>  _caps;

	std::size_t _num_nodes;
	std::size_t _num_arcs;

	std::unique_ptr<LinearSolverBackend> _solver;

	Solution _solution;
};

template <typename CapacityType, typename FlowType>
Graph<CapacityType, FlowType>::Graph(std::size_t node_num_max, std::size_t edge_num_max) :
	_caps(_graph),
	_num_nodes(2),
	_num_arcs(0) {

	_graph.reserveNode(node_num_max + 2); // plus source and sink
	_graph.reserveArc(2*edge_num_max);

	_source = _graph.addNode();
	_sink   = _graph.addNode();

	SolverFactory factory;
	_solver = std::unique_ptr<LinearSolverBackend>(factory.createLinearSolverBackend());
}

template <typename CapacityType, typename FlowType>
typename Graph<CapacityType, FlowType>::NodeId
Graph<CapacityType, FlowType>::add_nodes(std::size_t n) {

	if (n < 1)
		UTIL_THROW_EXCEPTION(
				UsageError,
				"at least one node has to be added with a call to add_nodes");

	_num_nodes += n;

	NodeId first = _graph.id(_graph.addNode());
	n--;

	for (std::size_t i = 0; i < n; i++)
		_graph.addNode();

	return first;
}


template <typename CapacityType, typename FlowType>
void
Graph<CapacityType, FlowType>::add_edge(NodeId u, NodeId v, CapacityType cap, CapacityType rev_cap) {

	GraphType::Arc fwd = _graph.addArc(_graph.nodeFromId(u), _graph.nodeFromId(v));
	GraphType::Arc bwd = _graph.addArc(_graph.nodeFromId(v), _graph.nodeFromId(u));

	_caps[fwd] = cap;
	_caps[bwd] = rev_cap;

	_num_arcs += 2;
}

template <typename CapacityType, typename FlowType>
void
Graph<CapacityType, FlowType>::add_tedge(NodeId n, CapacityType cap_source, CapacityType cap_sink) {

	GraphType::Arc source_arc = _graph.addArc(_source, _graph.nodeFromId(n));
	GraphType::Arc sink_arc   = _graph.addArc(_graph.nodeFromId(n), _sink);

	_caps[source_arc] = cap_source;
	_caps[sink_arc]   = cap_sink;

	_num_arcs += 2;
}

template <typename CapacityType, typename FlowType>
FlowType
Graph<CapacityType, FlowType>::maxflow() {

	std::size_t num_vars   = _num_arcs + _num_nodes;
	std::size_t source_var = _num_arcs + _graph.id(_source);
	std::size_t sink_var   = _num_arcs + _graph.id(_sink);

	LinearObjective objective(num_vars);
	LinearConstraints constraints;

	std::size_t n = 0;
	for (GraphType::ArcIt a(_graph); a != lemon::INVALID; ++a) {

		std::size_t arc_var = n;
		std::size_t node_u_var = _num_arcs + _graph.id(_graph.source(a));
		std::size_t node_v_var = _num_arcs + _graph.id(_graph.target(a));

		objective.setCoefficient(arc_var, _caps[a]);

		LinearConstraint positivity;
		positivity.setCoefficient(arc_var, 1.0);
		positivity.setRelation(GreaterEqual);
		positivity.setValue(0.0);

		LinearConstraint cut;
		cut.setCoefficient(arc_var, 1.0);
		cut.setCoefficient(node_u_var, -1.0);
		cut.setCoefficient(node_v_var,  1.0);
		cut.setRelation(GreaterEqual);
		cut.setValue(0.0);

		constraints.add(positivity);
		constraints.add(cut);

		n++;
	}

	for (GraphType::NodeIt i(_graph); i != lemon::INVALID; ++i) {

		std::size_t node_var = _num_arcs + _graph.id(i);

		LinearConstraint positivity;
		positivity.setCoefficient(node_var, 1.0);
		positivity.setRelation(GreaterEqual);
		positivity.setValue(0.0);

		constraints.add(positivity);
	}

	LinearConstraint source1sink0;
	source1sink0.setCoefficient(source_var, 1.0);
	source1sink0.setCoefficient(sink_var,  -1.0);
	source1sink0.setRelation(GreaterEqual);
	source1sink0.setValue(1);

	constraints.add(source1sink0);

	_solver->initialize(num_vars, Continuous);
	_solver->setObjective(objective);
	_solver->setConstraints(constraints);

	//std::cout << objective << std::endl;
	//for (auto c : constraints)
	//	std::cout << c << std::endl;

	std::string message;
	if (!_solver->solve(_solution, message))
		UTIL_THROW_EXCEPTION(
				LinearSolverBackendException,
				"linear program could not be solved: " << message);

	//std::cout << _solution.getVector() << std::endl;

	return _solution.getValue();
}

template <typename CapacityType, typename FlowType>
TerminalType
Graph<CapacityType, FlowType>::what_segment(NodeId n) {

	std::size_t node_var = _num_arcs + n;

	if (_solution[node_var] > 0.5)
		return Source;
	else
		return Sink;
}

#endif // PYMAXFLOW_MAXFLOW_GRAPH_H__
