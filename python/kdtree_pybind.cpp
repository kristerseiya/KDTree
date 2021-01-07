#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <tuple>
#include "kdtree/kdtree.h"

namespace py = pybind11;

class PyKDTree : private KDTree<double> {
public:
    PyKDTree() : KDTree<double>() {}

    void assign_numpy(py::array_t<double, py::array::c_style | py::array::forcecast> data) {
        if (data.ndim() != 2) {
            throw std::runtime_error("Input must be 2 dimensional");
        }
        this->assign((double*)data.data(),data.shape(1),data.shape(0));
    }

    std::tuple< py::array_t<size_t>,py::array_t<double> >
    search_knn(py::array_t<double, py::array::c_style | py::array::forcecast> query,
                   int k) {
        if (data_ == nullptr) {
            throw std::runtime_error("Data is not assigned");
        }
        if ((query.ndim() != 1)||(query.shape(0) != dimension_)) {
            throw std::runtime_error("Query must be a shape of ("+ std::to_string(dimension_) +",)");
        }
        // std::vector<size_t> neighbor_idx;
        // std::vector<double> distances;
        // this->searchKNN(std::vector<double>(query.data(),query.data()+query.shape(0)),
        //                 k, neighbor_idx, distances);
        // py::array_t<size_t> py_neighbor_idx(neighbor_idx.size(), neighbor_idx.data());
        // py::array_t<double> py_distances(distances.size(), distances.data());
        // return std::make_tuple(py_neighbor_idx, py_distances);
        auto neighbor_idx = new std::vector<size_t>();
        auto distances = new std::vector<double>();
        this->searchKNN(std::vector<double>(query.data(),query.data()+query.shape(0)),
                        k, *neighbor_idx, *distances);

        auto capsule1 = py::capsule(neighbor_idx, [](void* neighbor_idx) {
            delete reinterpret_cast<std::vector<size_t>*>(neighbor_idx); });
        auto capsule2 = py::capsule(distances, [](void* distances) {
            delete reinterpret_cast<std::vector<double>*>(distances); });

        py::array_t<size_t> py_neighbor_idx(neighbor_idx->size(),
                                            neighbor_idx->data(),
                                            capsule1);

        py::array_t<double> py_distances(distances->size(),
                                         distances->data(),
                                         capsule2);
        return std::make_tuple(py_neighbor_idx, py_distances);
    }

    std::tuple< py::array_t<size_t>,py::array_t<double> >
    search_radius(py::array_t<double, py::array::c_style | py::array::forcecast> query,
               double radius) {
        if (data_ == nullptr) {
            throw std::runtime_error("Data is not assigned");
        }
        if ((query.ndim() != 1)||(query.shape(0) != dimension_)) {
            throw std::runtime_error("Query must be a shape of ("+ std::to_string(dimension_) +",)");
        }
        std::vector<size_t> neighbor_idx;
        std::vector<double> distances;
        this->searchRadius(std::vector<double>(query.data(),query.data()+query.shape(0)),
                            radius, neighbor_idx, distances);
        py::array_t<size_t> py_neighbor_idx(neighbor_idx.size(), neighbor_idx.data());
        py::array_t<double> py_distances(distances.size(), distances.data());
        return std::make_tuple(py_neighbor_idx, py_distances);
    }

    std::tuple< py::array_t<size_t>,py::array_t<double> >
    search_hybrid(py::array_t<double, py::array::c_style | py::array::forcecast> query,
                  double radius, int max_n) {
        if (data_ == nullptr) {
            throw std::runtime_error("Data is not assigned");
        }
        if ((query.ndim() != 1)||(query.shape(0) != dimension_)) {
            throw std::runtime_error("Query must be a shape of ("+ std::to_string(dimension_) +",)");
        }
        std::vector<size_t> neighbor_idx;
        std::vector<double> distances;
        this->searchHybrid(std::vector<double>(query.data(),query.data()+query.shape(0)),
                           radius, neighbor_idx, distances, max_n);
        py::array_t<size_t> py_neighbor_idx(neighbor_idx.size(), neighbor_idx.data());
        py::array_t<double> py_distances(distances.size(), distances.data());
        return std::make_tuple(py_neighbor_idx, py_distances);
    }

    int get_n_calls() {
        return visited_;
    }
};

PYBIND11_MODULE(pykdtree, m) {
    py::class_<PyKDTree> (m, "KDTree")
            .def(py::init<>())
            .def("assign", &PyKDTree::assign_numpy)
            .def("searchKNN", &PyKDTree::search_knn)
            .def("searchRadius", &PyKDTree::search_radius)
            .def("searchHybrid", &PyKDTree::search_hybrid)
            .def("get_n_calls", &PyKDTree::get_n_calls);
}
