#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

struct sphere_properties
{
    double sphere_const;
};

inline
double
sphere_fn(const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out,
          void* opt_data)
{
    sphere_properties* objfn_data = reinterpret_cast<sphere_properties*>
                                    (opt_data);
    Eigen::VectorXd x = vals_inp.array() + objfn_data->sphere_const;
    // Sphere Function
    double obj_val = x.dot(x);

    // Gradient Computation
    if (grad_out)
    {
        *grad_out = 2.0 * x;
    }

    return obj_val;
}

int main()
{
    const int test_dim = 6;
    // initial values (1,1,...,1)
    Eigen::VectorXd x = Eigen::VectorXd::Ones(test_dim);
    // Create a structure with sphere properties it can be any type of object
    sphere_properties sphere;
    // set value of the sphere property object
    sphere.sphere_const = 4;

    // settings of the algorithm see ITHACA-FV/src/thirdparty/OptimLib/misc/optim_structs.hpp file for more details
    optim::algo_settings_t settings;
    // set adam as gd method
    settings.gd_settings.method = 6;
    // Create a void pointer and cast the pointer of the sphere object in it
    void* p = reinterpret_cast<void*>(&sphere);
    // Call the optimization algorithm
    bool success = optim::gd(x, sphere_fn, p, settings);

    if (success)
    {
        std::cout << "gd: sphere test completed successfully." << "\n";
    }

    else
    {
        std::cout << "gd: sphere test completed unsuccessfully." << "\n";
    }

    std::cout << "gd: solution to sphere test:\n" << x << std::endl;
    return 0;
}
