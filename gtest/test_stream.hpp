#include "ITHACAstream.H"
#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
//include the google test dependencies
#include <gtest/gtest.h>



TEST(ITHACAstream, ReadWriteTensor)
{
    Eigen::Tensor<double, 3> output;
    Eigen::Tensor<double, 3> input;
    output.resize(3, 4, 5);
    output.setRandom();
    ITHACAstream::SaveDenseTensor(output, "./", "output");
    ITHACAstream::ReadDenseTensor(input, "./", "output");
    Eigen::Tensor<double, 0> difference = ((output - input).abs().sum());
    Eigen::Tensor<double, 0> comp;
    comp(0) = 1e-18;
    EXPECT_LT(difference(0),comp(0));
    system("rm output");
}

TEST(ITHACAstream, ReadAndWriteNPYMatrix)
{
    Eigen::MatrixXi MI_out = Eigen::MatrixXi::Random(5, 5);
    cnpy::save(MI_out, "arr_int.npy");
    Eigen::MatrixXi MI_inp;
    cnpy::load(MI_inp, "arr_int.npy");
    int difference_I = (MI_out - MI_inp).sum();
    //
    Eigen::MatrixXd MD_out = Eigen::MatrixXd::Random(5, 5);
    cnpy::save(MD_out, "arr_dou.npy");
    Eigen::MatrixXd MD_inp;
    cnpy::load(MD_inp, "arr_dou.npy");
    double difference_D = (MD_out - MD_inp).sum();
    //
    Eigen::MatrixXf MF_out = Eigen::MatrixXf::Random(5, 5);
    cnpy::save(MF_out, "arr_flo.npy");
    Eigen::MatrixXf MF_inp;
    cnpy::load(MF_inp, "arr_flo.npy");
    float difference_F = (MF_out - MF_inp).sum();
    //
    Eigen::VectorXi VI_out = Eigen::VectorXi::Random(5);
    cnpy::save(VI_out, "vec_int.npy");
    Eigen::VectorXi VI_inp;
    cnpy::load(VI_inp, "vec_int.npy");
    int difference_IV = (VI_out - VI_inp).sum();
    //
    Eigen::VectorXd VD_out = Eigen::VectorXd::Random(5);
    cnpy::save(VD_out, "vec_dou.npy");
    Eigen::VectorXd VD_inp;
    cnpy::load(VD_inp, "vec_dou.npy");
    double difference_DV = (VD_out - VD_inp).sum();
    //
    Eigen::VectorXf VF_out = Eigen::VectorXf::Random(5);
    cnpy::save(VF_out, "vec_flo.npy");
    Eigen::VectorXf VF_inp;
    cnpy::load(VF_inp, "vec_flo.npy");
    float difference_FV = (VF_out - VF_inp).sum();
    //
    Eigen::Tensor<int, 3> TI_inp(2, 3, 4);
    TI_inp.setRandom();
    cnpy::save(TI_inp, "ten_int.npy");
    Eigen::Tensor<int, 3> TI_out(2, 3, 4);
    cnpy::load(TI_out, "ten_int.npy");
    double difference_TI = ((Eigen::Tensor<int, 0>)(TI_out - TI_inp).sum())(0);
    //
    Eigen::Tensor<double, 3> TD_inp(2, 3, 4);
    TD_inp.setRandom();
    cnpy::save(TD_inp, "ten_dou.npy");
    Eigen::Tensor<double, 3> TD_out(2, 3, 4);
    cnpy::load(TD_out, "ten_dou.npy");
    double difference_TD = ((Eigen::Tensor<double, 0>)(TD_out - TD_inp).sum())(0);
    //
    Eigen::Tensor<float, 3> TF_inp(2, 3, 4);
    TF_inp.setRandom();
    cnpy::save(TF_inp, "ten_flo.npy");
    Eigen::Tensor<float, 3> TF_out(2, 3, 4);
    cnpy::load(TF_out, "ten_flo.npy");
    double difference_TF = ((Eigen::Tensor<float, 0>)(TF_out - TF_inp).sum())(0);

    EXPECT_TRUE(difference_I == 0) << "Reading and Writing Integer array failed";
    EXPECT_TRUE(difference_D == 0) << "Reading and Writing double array failed";
    EXPECT_TRUE(difference_F == 0) << "Reading and Writing float array failed";
    EXPECT_TRUE(difference_TI == 0) << "Reading and Writing Integer tensor failed";
    EXPECT_TRUE(difference_TI == 0) << "Reading and Writing double tensor failed";
    EXPECT_TRUE(difference_TF == 0) << "Reading and Writing float tensor failed";
    system("rm *.npy");
}


//TEST(ITHACAstream, cnpyTEST)
//{
    //const int Nx = 128;
    //const int Ny = 64;
    //const int Nz = 32;
    //Eigen::MatrixXf m(2, 2);
    ////create random data
    //std::complex<double>* data = new std::complex<double>[Nx * Ny * Nz];

    //for (int i = 0; i < Nx * Ny * Nz;
            //i++) data[i] = std::complex<double>(rand(), rand());

    ////save it to file
    //std::vector<size_t> shape = {Nz, Ny, Nx};
    //cnpy::npy_save("arr1.npy", data, shape, "w");
    ////load it into a new array
    //cnpy::NpyArray arr = cnpy::npy_load("arr1.npy");
    //std::complex<double>* loaded_data = arr.data<std::complex<double>>();
    ////make sure the loaded data matches the saved data
    //EXPECT_TRUE(arr.word_size == sizeof(std::complex<double>));
    //EXPECT_TRUE(arr.shape.size() == 3 && arr.shape[0] == Nz && arr.shape[1] == Ny &&
           //arr.shape[2] == Nx);

    //for (int i = 0; i < Nx * Ny * Nz; i++) EXPECT_TRUE(data[i] == loaded_data[i]);

    ////append the same data to file
    ////npy array on file now has shape (Nz+Nz,Ny,Nx)
    //cnpy::npy_save("arr1.npy", data, shape, "a");
    ////now write to an npz file
    ////non-array variables are treated as 1D arrays with 1 element
    //double myVar1 = 1.2;
    //char myVar2 = 'a';
    //std::vector<size_t> shape2 = {1};
    //cnpy::npz_save("out.npz", "myVar1", &myVar1, shape2,
                   //"w"); //"w" overwrites any existing file
    //cnpy::npz_save("out.npz", "myVar2", &myVar2, shape2,
                   //"a"); //"a" appends to the file we created above
    //cnpy::npz_save("out.npz", "arr1", data, shape,
                   //"a"); //"a" appends to the file we created above
    ////load a single var from the npz file
    //cnpy::NpyArray arr2 = cnpy::npz_load("out.npz", "arr1");
    ////load the entire npz file
    //cnpy::npz_t my_npz = cnpy::npz_load("out.npz");
    ////check that the loaded myVar1 matches myVar1
    //cnpy::NpyArray arr_mv1 = my_npz["myVar1"];
    //double* mv1 = arr_mv1.data<double>();
    //EXPECT_TRUE(arr_mv1.shape.size() == 1 && arr_mv1.shape[0] == 1);
    //EXPECT_TRUE(mv1[0] == myVar1);
    ////cleanup: note that we are responsible for deleting all loaded data
    //delete[] data;
    //delete[] loaded_data;
    //Eigen::MatrixXf pc1 = Eigen::MatrixXf::Random(5, 5);
    //cnpy::save(pc1, "arr.npy");
    //Eigen::MatrixXf pc0 = cnpy::load(pc0, "arr.npy");
    //std::cout << pc0 << std::endl;
    //cnpy::save(pc0, "test.npy");
//}


TEST(ITHACAstream, TestSparseMatrix)
{
    Eigen::SparseMatrix<double> mat;
    cnpy::load(mat,"forCnpy.npz");
    cnpy::save(mat,"forNpy.npz");
    Eigen::SparseMatrix<double> mat2;
    cnpy::load(mat2,"forNpy.npz");
    EXPECT_TRUE((mat-mat2).norm() == 0) << "Read And Write NPZ sparse matrix failed!" << std::endl;
}
