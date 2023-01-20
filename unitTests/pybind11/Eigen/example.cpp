#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "Eigen/Dense"
#include "iostream"
class Matrix
{
    public:
        Matrix(int dim1, int dim2) : dim1(dim1), dim2(dim2)
        {
            a = Eigen::MatrixXd::Random(dim1,dim2);
        }
        ~Matrix() {}

        Eigen::MatrixXd &get_matrix() {return a;}
        const Eigen::MatrixXd &view_matrix() {return a;}

    private:
        int dim1;
        int dim2;
        Eigen::MatrixXd a;
};


namespace py = pybind11;

PYBIND11_MODULE(example, m)
{
    // bindings to Matrix class
    py::class_<Matrix>(m, "Matrix")
    .def(py::init<int, int>())
    .def("get_matrix", &Matrix::get_matrix, py::return_value_policy::reference_internal)
    .def("view_matrix", &Matrix::view_matrix, py::return_value_policy::reference_internal);
}
