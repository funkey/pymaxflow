#include <boost/python.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <util/exceptions.h>
#include <crag/Crag.h>
#include <crag/CragVolumes.h>
#include <io/Hdf5CragStore.h>
#include <io/Hdf5VolumeStore.h>
#include <io/volumes.h>
#include <inference/Costs.h>
#include <inference/RandomForest.h>
#include <inference/CragSolution.h>
#include <learning/BundleOptimizer.h>
#include <learning/Loss.h>
#include "PyOracle.h"
#include "logging.h"

template <typename Map, typename K, typename V>
const V& genericGetter(const Map& map, const K& k) { return map[k]; }
template <typename Map, typename K, typename V>
void genericSetter(Map& map, const K& k, const V& value) { map[k] = value; }

template <typename Map, typename K, typename V>
void featuresSetter(Map& map, const K& k, const V& value) { map.set(k, value); }

#if defined __clang__ && __clang_major__ < 6
// std::shared_ptr support
	template<class T> T* get_pointer(std::shared_ptr<T> p){ return p.get(); }
#endif

namespace pymaxflow {

/**
 * Translates an Exception into a python exception.
 *
 **/
void translateException(const Exception& e) {

	if (boost::get_error_info<error_message>(e))
		PyErr_SetString(PyExc_RuntimeError, boost::get_error_info<error_message>(e)->c_str());
	else
		PyErr_SetString(PyExc_RuntimeError, e.what());
}

/**
 * Defines all the python classes in the module libpymaxflow. Here we decide 
 * which functions and data members we wish to expose.
 */
BOOST_PYTHON_MODULE(pymaxflow) {

	boost::python::register_exception_translator<Exception>(&translateException);

	// Logging
	boost::python::enum_<logger::LogLevel>("LogLevel")
			.value("Quiet", logger::Quiet)
			.value("Error", logger::Error)
			.value("Debug", logger::Debug)
			.value("All", logger::All)
			.value("User", logger::User)
			;
	boost::python::def("setLogLevel", setLogLevel);
	boost::python::def("getLogLevel", getLogLevel);

	//// Crag
	//boost::python::class_<Crag, boost::noncopyable>("Crag")
			//.def("addNode", static_cast<Crag::CragNode(Crag::*)()>(&Crag::addNode))
			//.def("addNode", static_cast<Crag::CragNode(Crag::*)(Crag::NodeType)>(&Crag::addNode))
			//.def("addAdjacencyEdge", static_cast<Crag::CragEdge(Crag::*)(Crag::CragNode, Crag::CragNode)>(&Crag::addAdjacencyEdge))
			//.def("addAdjacencyEdge", static_cast<Crag::CragEdge(Crag::*)(Crag::CragNode, Crag::CragNode, Crag::EdgeType)>(&Crag::addAdjacencyEdge))
			//.def("addSubsetArc", &Crag::addSubsetArc)
			//.def("erase", static_cast<void(Crag::*)(Crag::CragNode)>(&Crag::erase))
			//.def("erase", static_cast<void(Crag::*)(Crag::CragEdge)>(&Crag::erase))
			//.def("erase", static_cast<void(Crag::*)(Crag::CragArc)>(&Crag::erase))
			//.def("id", static_cast<int(Crag::*)(Crag::CragNode) const>(&Crag::id))
			//.def("id", static_cast<int(Crag::*)(Crag::CragEdge) const>(&Crag::id))
			//.def("id", static_cast<int(Crag::*)(Crag::CragArc)  const>(&Crag::id))
			//.def("nodeFromId", &Crag::nodeFromId)
			//.def("oppositeNode", &Crag::oppositeNode)
			//.def("nodes", &Crag::nodes)
			//.def("edges", &Crag::edges)
			//.def("arcs", &Crag::arcs)
			//.def("adjEdges", &Crag::adjEdges)
			//.def("outArcs", &Crag::outArcs)
			//.def("inArcs", &Crag::inArcs)
			//.def("type", static_cast<Crag::NodeType(Crag::*)(Crag::CragNode) const>(&Crag::type))
			//.def("type", static_cast<Crag::EdgeType(Crag::*)(Crag::CragEdge) const>(&Crag::type))
			//.def("getLevel", &Crag::getLevel)
			//.def("isRootNode", &Crag::isRootNode)
			//.def("isLeafNode", &Crag::isLeafNode)
			//.def("isLeafEdge", &Crag::isLeafEdge)
			//.def("leafNodes", &Crag::leafNodes)
			//.def("leafEdges", static_cast<std::set<Crag::CragEdge>(Crag::*)(Crag::CragNode) const>(&Crag::leafEdges))
			//.def("leafEdges", static_cast<std::set<Crag::CragEdge>(Crag::*)(Crag::CragEdge) const>(&Crag::leafEdges))
			//;

}

} // namespace pymaxflow
