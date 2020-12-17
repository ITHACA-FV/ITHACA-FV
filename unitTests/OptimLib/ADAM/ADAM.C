#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

inline
double
sphere_fn(const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)
{
    double obj_val = vals_inp.dot(vals_inp);

    if (grad_out) {
        *grad_out = 2.0*vals_inp;
    }

    return obj_val;
}

int main()
{
    const int test_dim = 6;

    Eigen::VectorXd x = Eigen::VectorXd::Ones(test_dim); // initial values (1,1,...,1)

    bool success = optim::gd(x, sphere_fn, nullptr);

    if (success) {
        std::cout << "gd: sphere test completed successfully." << "\n";
    } else {
        std::cout << "gd: sphere test completed unsuccessfully." << "\n";
    }

    std::cout << "gd: solution to sphere test:\n" << x << std::endl;

    return 0;
}
