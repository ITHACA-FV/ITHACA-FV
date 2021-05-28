//Copyright (C) 2011  Carl Rogers
//Released under MIT License
//license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#include "ITHACAstream.H"
#include <complex>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <iomanip>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"

char cnpy::BigEndianTest()
{
    unsigned char x[] = {1, 0};
    short y = *(short*) x;
    return y == 1 ? '<' : '>';
}

char cnpy::map_type(const std::type_info& t)
{
    if (t == typeid(float) )
    {
        return 'f';
    }

    if (t == typeid(double) )
    {
        return 'f';
    }

    if (t == typeid(long double) )
    {
        return 'f';
    }

    if (t == typeid(int) )
    {
        return 'i';
    }

    if (t == typeid(char) )
    {
        return 'i';
    }

    if (t == typeid(short) )
    {
        return 'i';
    }

    if (t == typeid(long) )
    {
        return 'i';
    }

    if (t == typeid(long long) )
    {
        return 'i';
    }

    if (t == typeid(unsigned char) )
    {
        return 'u';
    }

    if (t == typeid(unsigned short) )
    {
        return 'u';
    }

    if (t == typeid(unsigned long) )
    {
        return 'u';
    }

    if (t == typeid(unsigned long long) )
    {
        return 'u';
    }

    if (t == typeid(unsigned int) )
    {
        return 'u';
    }

    if (t == typeid(bool) )
    {
        return 'b';
    }

    if (t == typeid(std::complex<float>) )
    {
        return 'c';
    }

    if (t == typeid(std::complex<double>) )
    {
        return 'c';
    }

    if (t == typeid(std::complex<long double>) )
    {
        return 'c';
    }
    else
    {
        return '?';
    }
}

template<> std::vector<char>& cnpy::operator+=(std::vector<char>& lhs,
        const std::string rhs)
{
    lhs.insert(lhs.end(), rhs.begin(), rhs.end());
    return lhs;
}

template<> std::vector<char>& cnpy::operator+=(std::vector<char>& lhs,
        const char* rhs)
{
    //write in little endian
    size_t len = strlen(rhs);
    lhs.reserve(len);

    for (size_t byte = 0; byte < len; byte++)
    {
        lhs.push_back(rhs[byte]);
    }

    return lhs;
}

void cnpy::parse_npy_header(FILE* fp, unsigned int& word_size,
                            unsigned int*& shape, unsigned int& ndims, bool& fortran_order)
{
    char buffer[256];
    size_t res = fread(buffer, sizeof(char), 11, fp);

    if (res != 11)
    {
        throw std::runtime_error("parse_npy_header: failed fread");
    }

    std::string header = fgets(buffer, 256, fp);
    assert(header[header.size() - 1] == '\n');
    int loc1, loc2;
    //fortran order
    loc1 = header.find("fortran_order") + 16;
    fortran_order = (header.substr(loc1, 5) == "True" ? true : false);
    //shape
    loc1 = header.find("(");
    loc2 = header.find(")");
    std::string str_shape = header.substr(loc1 + 1, loc2 - loc1 - 1);

    if (str_shape[str_shape.size() - 1] == ',')
    {
        ndims = 1;
    }
    else
    {
        ndims = std::count(str_shape.begin(), str_shape.end(), ',') + 1;
    }

    shape = new unsigned int[ndims];

    for (unsigned int i = 0; i < ndims; i++)
    {
        loc1 = str_shape.find(",");
        shape[i] = atoi(str_shape.substr(0, loc1).c_str());
        str_shape = str_shape.substr(loc1 + 1);
    }

    //endian, word size, data type
    //byte order code | stands for not applicable.
    //not sure when this applies except for byte array
    loc1 = header.find("descr") + 9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    assert(littleEndian);
    //char type = header[loc1+1];
    //assert(type == map_type(T));
    std::string str_ws = header.substr(loc1 + 2);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(0, loc2).c_str());
}

void cnpy::parse_zip_footer(FILE* fp, unsigned short& nrecs,
                            unsigned int& global_header_size, unsigned int& global_header_offset)
{
    std::vector<char> footer(22);
    fseek(fp, -22, SEEK_END);
    size_t res = fread(&footer[0], sizeof(char), 22, fp);

    if (res != 22)
    {
        throw std::runtime_error("parse_zip_footer: failed fread");
    }

    unsigned short disk_no, disk_start, nrecs_on_disk, comment_len;
    disk_no = *(unsigned short*) &footer[4];
    disk_start = *(unsigned short*) &footer[6];
    nrecs_on_disk = *(unsigned short*) &footer[8];
    nrecs = *(unsigned short*) &footer[10];
    global_header_size = *(unsigned int*) &footer[12];
    global_header_offset = *(unsigned int*) &footer[16];
    comment_len = *(unsigned short*) &footer[20];
    assert(disk_no == 0);
    assert(disk_start == 0);
    assert(nrecs_on_disk == nrecs);
    assert(comment_len == 0);
}

cnpy::NpyArray load_the_npy_file(FILE* fp)
{
    unsigned int* shape;
    unsigned int ndims, word_size;
    bool fortran_order;
    cnpy::parse_npy_header(fp, word_size, shape, ndims, fortran_order);
    unsigned long long size =
        1; //long long so no overflow when multiplying by word_size

    for (unsigned int i = 0; i < ndims; i++)
    {
        size *= shape[i];
    }

    cnpy::NpyArray arr;
    arr.word_size = word_size;
    arr.shape = std::vector<unsigned int>(shape, shape + ndims);
    delete[] shape;
    arr.data = new char[size * word_size];
    arr.fortran_order = fortran_order;
    size_t nread = fread(arr.data, word_size, size, fp);

    if (nread != size)
    {
        throw std::runtime_error("load_the_npy_file: failed fread");
    }

    return arr;
}

cnpy::npz_t cnpy::npz_load(std::string fname)
{
    FILE* fp = fopen(fname.c_str(), "rb");

    if (!fp)
    {
        printf("npz_load: Error! Unable to open file %s!\n", fname.c_str());
    }

    assert(fp);
    cnpy::npz_t arrays;

    while (1)
    {
        std::vector<char> local_header(30);
        size_t headerres = fread(&local_header[0], sizeof(char), 30, fp);

        if (headerres != 30)
        {
            throw std::runtime_error("npz_load: failed fread");
        }

        //if we've reached the global header, stop reading
        if (local_header[2] != 0x03 || local_header[3] != 0x04)
        {
            break;
        }

        //read in the variable name
        unsigned short name_len = *(unsigned short*) &local_header[26];
        std::string varname(name_len, ' ');
        size_t vname_res = fread(&varname[0], sizeof(char), name_len, fp);

        if (vname_res != name_len)
        {
            throw std::runtime_error("npz_load: failed fread");
        }

        //erase the lagging .npy
        varname.erase(varname.end() - 4, varname.end());
        //read in the extra field
        unsigned short extra_field_len = *(unsigned short*) &local_header[28];

        if (extra_field_len > 0)
        {
            std::vector<char> buff(extra_field_len);
            size_t efield_res = fread(&buff[0], sizeof(char), extra_field_len, fp);

            if (efield_res != extra_field_len)
            {
                throw std::runtime_error("npz_load: failed fread");
            }
        }

        arrays[varname] = load_the_npy_file(fp);
    }

    fclose(fp);
    return arrays;
}

cnpy::NpyArray cnpy::npz_load(std::string fname,
                              std::string varname)
{
    FILE* fp = fopen(fname.c_str(), "rb");

    if (!fp)
    {
        printf("npz_load: Error! Unable to open file %s!\n", fname.c_str());
        abort();
    }

    while (1)
    {
        std::vector<char> local_header(30);
        size_t header_res = fread(&local_header[0], sizeof(char), 30, fp);

        if (header_res != 30)
        {
            throw std::runtime_error("npz_load: failed fread");
        }

        //if we've reached the global header, stop reading
        if (local_header[2] != 0x03 || local_header[3] != 0x04)
        {
            break;
        }

        //read in the variable name
        unsigned short name_len = *(unsigned short*) &local_header[26];
        std::string vname(name_len, ' ');
        size_t vname_res = fread(&vname[0], sizeof(char), name_len, fp);

        if (vname_res != name_len)
        {
            throw std::runtime_error("npz_load: failed fread");
        }

        vname.erase(vname.end() - 4, vname.end()); //erase the lagging .npy
        //read in the extra field
        unsigned short extra_field_len = *(unsigned short*) &local_header[28];
        fseek(fp, extra_field_len, SEEK_CUR); //skip past the extra field

        if (vname == varname)
        {
            NpyArray array = load_the_npy_file(fp);
            fclose(fp);
            return array;
        }
        else
        {
            //skip past the data
            unsigned int size = *(unsigned int*) &local_header[22];
            fseek(fp, size, SEEK_CUR);
        }
    }

    fclose(fp);
    printf("npz_load: Error! Variable name %s not found in %s!\n", varname.c_str(),
           fname.c_str());
    abort();
}

cnpy::NpyArray cnpy::npy_load(std::string fname)
{
    FILE* fp = fopen(fname.c_str(), "rb");

    if (!fp)
    {
        printf("npy_load: Error! Unable to open file %s!\n", fname.c_str());
        abort();
    }

    NpyArray arr = load_the_npy_file(fp);
    fclose(fp);
    return arr;
}

template<class typeNumber, int dim>
void cnpy::save(const Eigen::Matrix<typeNumber, Eigen::Dynamic, dim>&
                mat, const std::string fname)
{
    std::vector<typeNumber> matvec(mat.rows() * mat.cols());
    unsigned int shape[] = {(unsigned int) mat.rows(), (unsigned int)mat.cols()};

    for (int i = 0; i < mat.rows(); ++i)
    {
        for (int j = 0; j < mat.cols(); ++j)
        {
            matvec[i * mat.cols() + j] = mat(i, j);
        }
    }

    npy_save(fname, matvec.data(), shape, 2);
}

template<class typeNumber>
void cnpy::save(const Eigen::Tensor<typeNumber, 3>&
                tens, const std::string fname)
{
    typename Eigen::Tensor<typeNumber, 3>::Dimensions dim = tens.dimensions();
    int tot = 1;

    for (int k = 0; k < dim.size(); k++)
    {
        tot *= dim[k];
    }

    std::vector<typeNumber> matvec(tot);
    unsigned int shape[] = {(unsigned int) dim[0], (unsigned int) dim[1], (unsigned int) dim[2]};

    for (int i = 0; i < dim[0]; ++i)
    {
        for (int j = 0; j < dim[1]; ++j)
        {
            for (int k = 0; k < dim[2]; ++k)
            {
                matvec[i * dim[2]*dim[1] + j * dim[2] + k] = tens(i, j, k);
            }
        }
    }

    npy_save(fname, matvec.data(), shape, 3);
}

template<class typeNumber, int dim>
Eigen::Matrix<typeNumber, Eigen::Dynamic, dim> cnpy::load(
    Eigen::Matrix<typeNumber, Eigen::Dynamic, dim>& mat,
    const std::string fname, std::string order)
{
    M_Assert(order == "rowMajor" ||
             order == "colMajor", "Order can be only rowMajor or colMajor");
    NpyArray arr = npy_load(fname);
    assert(arr.shape.size() == 2);
    mat.resize(arr.shape[0], arr.shape[1]);
    typeNumber* data = reinterpret_cast<typeNumber*>(arr.data);

    if (order == "rowMajor")
    {
        for (size_t i = 0; i < arr.shape[0]; ++i)
        {
            for (size_t j = 0; j < arr.shape[1]; ++j)
            {
                mat(i, j) = (typeNumber) data[arr.shape[1] * i + j];
            }
        }
    }
    else if (order == "colMajor")
    {
        for (size_t i = 0; i < arr.shape[0]; ++i)
        {
            for (size_t j = 0; j < arr.shape[1]; ++j)
            {
                mat(i, j) = (typeNumber) data[arr.shape[0] * j + i];
            }
        }
    }

    delete[] arr.data;
    return mat;
}

template<typename typeNumber>
Eigen::Tensor<typeNumber, 3> cnpy::load(Eigen::Tensor<typeNumber, 3>& tens,
                                        const std::string fname, std::string order)
{
    M_Assert(order == "rowMajor" ||
             order == "colMajor", "Order can be only rowMajor or colMajor");
    NpyArray arr = npy_load(fname);
    assert(arr.shape.size() == 3);
    tens.resize(arr.shape[0], arr.shape[1], arr.shape[2]);
    typeNumber* data = reinterpret_cast<typeNumber*>(arr.data);

    if (order == "rowMajor")
    {
        for (size_t i = 0; i < arr.shape[0]; ++i)
        {
            for (size_t j = 0; j < arr.shape[1]; ++j)
            {
                for (size_t k = 0; k < arr.shape[2]; ++k)
                {
                    tens(i, j, k) = (typeNumber) data[arr.shape[1] * arr.shape[2] * i + j *
                                                                   arr.shape[2] + k];
                }
            }
        }
    }
    else if (order == "colMajor")
    {
        for (size_t i = 0; i < arr.shape[0]; ++i)
        {
            for (size_t j = 0; j < arr.shape[1]; ++j)
            {
                for (size_t k = 0; k < arr.shape[2]; ++k)
                {
                    tens(i, j, k) = (typeNumber) data[arr.shape[0] * arr.shape[1] * k + j *
                                                                   arr.shape[0] + i];
                }
            }
        }
    }

    delete[] arr.data;
    return tens;
}

template<typename T>
Eigen::SparseMatrix<T> cnpy::load(Eigen::SparseMatrix<T>& smatrix,
                                  const std::string fname)
{
    auto d1 = cnpy::npz_load(fname, "data");
    auto d2 = cnpy::npz_load(fname, "indices");
    auto d3 = cnpy::npz_load(fname, "indptr");
    auto d4 = cnpy::npz_load(fname, "shape");
    double* data = reinterpret_cast<double*>(d1.data);
    int32_t* indices = reinterpret_cast<int32_t*>(d2.data);
    int32_t* indptr = reinterpret_cast<int32_t*>(d3.data);
    int* shape = reinterpret_cast<int*>(d4.data);
    int rows, cols;
    M_Assert(d4.shape[0] == 2, "Method works only with 2-D matrices");
    rows = (int) shape[0];
    cols = (int) d3.shape[0] - 1;
    int nel = d1.shape[0];
    smatrix.resize(rows, cols);
    smatrix.reserve(nel);
    typedef Eigen::Triplet<T> Trip;
    std::vector<Trip> tripletList;

    for (int i = 0; i < cols; ++i)
    {
        for (int j = indptr[i]; j < indptr[i + 1]; j++)
        {
            tripletList.push_back(Trip(indices[j], i, (T) data[j]));
        }
    }

    smatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return smatrix;
}

template<typename T>
void cnpy::save(const Eigen::SparseMatrix<T>& mat, const std::string fname)
{
    unsigned int shape1[] = {mat.nonZeros()};
    unsigned int shape2[] = {mat.outerSize() + 1};
    unsigned int shape3[] = {2};
    unsigned int shape4[] = {256};
    cnpy::npz_save(fname, "data", mat.valuePtr(), shape1, 1, "w");
    cnpy::npz_save(fname, "indices", mat.innerIndexPtr(), shape1, 1, "a");
    cnpy::npz_save(fname, "indptr", mat.outerIndexPtr(), shape2, 1, "a");
    int64_t t1 = mat.rows();
    int64_t t2 = mat.cols();
    int64_t* shape = new int64_t[2];
    shape[0] = mat.rows();
    shape[1] = mat.cols();
    cnpy::npz_save(fname, "shape", shape, shape3, 1, "a");
    char myVar2 = 'abc';
    cnpy::npz_save(fname, "format", &myVar2, shape4, 1, "a");
}

template void cnpy::save(const Eigen::MatrixXi& mat, const std::string fname);
template Eigen::MatrixXi cnpy::load(Eigen::MatrixXi& mat,
                                    const std::string fname, std::string order);
template void cnpy::save(const Eigen::MatrixXd& mat, const std::string fname);
template Eigen::MatrixXd cnpy::load(Eigen::MatrixXd& mat,
                                    const std::string fname, std::string order);
template void cnpy::save(const Eigen::MatrixXf& mat, const std::string fname);
template Eigen::MatrixXf cnpy::load(Eigen::MatrixXf& mat,
                                    const std::string fname, std::string order);
template void cnpy::save(const Eigen::MatrixXcd& mat, const std::string fname);
template Eigen::MatrixXcd cnpy::load(Eigen::MatrixXcd& mat,
                                     const std::string fname, std::string order);
template void cnpy::save(const Eigen::VectorXi& mat, const std::string fname);
template Eigen::VectorXi cnpy::load(Eigen::VectorXi& mat,
                                    const std::string fname, std::string order);
template void cnpy::save(const Eigen::VectorXd& mat, const std::string fname);
template Eigen::VectorXd cnpy::load(Eigen::VectorXd& mat,
                                    const std::string fname, std::string order);
template void cnpy::save(const Eigen::VectorXf& mat, const std::string fname);
template Eigen::VectorXf cnpy::load(Eigen::VectorXf& mat,
                                    const std::string fname, std::string order);
template void cnpy::save(const Eigen::VectorXcd& mat, const std::string fname);
template Eigen::VectorXcd cnpy::load(Eigen::VectorXcd& mat,
                                     const std::string fname, std::string order);
template void cnpy::save(const Eigen::SparseMatrix<double>& mat,
                         const std::string fname);
template Eigen::SparseMatrix<double> cnpy::load(Eigen::SparseMatrix<double>&
        smatrix, const std::string fname);

template void cnpy::save(const Eigen::Tensor<int, 3>& mat,
                         const std::string fname);
template Eigen::Tensor<int, 3> cnpy::load(Eigen::Tensor<int, 3>& tens,
        const std::string fname, std::string order);
template void cnpy::save(const Eigen::Tensor<double, 3>& mat,
                         const std::string fname);
template Eigen::Tensor<double, 3> cnpy::load(Eigen::Tensor<double, 3>& tens,
        const std::string fname, std::string order);
template void cnpy::save(const Eigen::Tensor<float, 3>& mat,
                         const std::string fname);
template Eigen::Tensor<float, 3> cnpy::load(Eigen::Tensor<float, 3>& tens,
        const std::string fname, std::string order);
template void cnpy::save(const Eigen::Tensor<std::complex<double>, 3>& mat,
                         const std::string fname);
template Eigen::Tensor<std::complex<double>, 3> cnpy::load(
    Eigen::Tensor<std::complex<double>, 3>& tens,
    const std::string fname, std::string order);
#pragma GCC diagnostic pop
