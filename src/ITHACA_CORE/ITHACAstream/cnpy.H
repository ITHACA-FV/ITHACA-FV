//Copyright (C) 2011  Carl Rogers
//Copyright (C) 2017  Matheus Gadelha
//Copyright (C) 2019  Giovanni Stabile
//Released under MIT License
//license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#ifndef LIBCNPY_H_
#define LIBCNPY_H_

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>
#include<string>
#include<stdexcept>
#include<sstream>
#include<vector>
#include<cstdio>
#include<typeinfo>
#include<iostream>
#include<cassert>
#include<zlib.h>
#include<map>
#include<memory>
#include<stdint.h>
#include<numeric>
#include <regex>


namespace cnpy
{

struct NpyArray
{
    NpyArray(const std::vector<size_t>& _shape, size_t _word_size,
             bool _fortran_order, std::string _number_type) :
        shape(_shape), word_size(_word_size), fortran_order(_fortran_order),
        number_type(_number_type)
    {
        num_vals = 1;

        for (size_t i = 0; i < shape.size(); i++)
        {
            num_vals *= shape[i];
        }

        data_holder = std::shared_ptr<std::vector<char>>(
                          new std::vector<char>(num_vals * word_size));
    }

    NpyArray() : shape(0), word_size(0), fortran_order(0), num_vals(0) { }

    template<typename T>
    T* data()
    {
        return reinterpret_cast<T*>(&(*data_holder)[0]);
    }

    template<typename T>
    const T* data() const
    {
        return reinterpret_cast<T*>(&(*data_holder)[0]);
    }

    template<typename T>
    std::vector<T> as_vec()
    {
        if (number_type == "f4")
        {
            float* p = data<float>();
            return std::vector<T>(static_cast<float*>(p),
                                  static_cast<float*>(p) + num_vals);
        }

        if (number_type == "f8")
        {
            double* p = data<double>();
            return std::vector<T>(static_cast<double*>(p),
                                  static_cast<double*>(p) + num_vals);
        }

        if (number_type == "f16")
        {
            long double* p = data<long double>();
            return std::vector<T>(static_cast<long double*>(p),
                                  static_cast<long double*>(p) + num_vals);
        }

        if (number_type == "i2")
        {
            short* p = data<short>();
            return std::vector<T>(static_cast<short*>(p),
                                  static_cast<short*>(p) + num_vals);
        }

        if (number_type == "i4")
        {
            int* p = data<int>();
            return std::vector<T>(static_cast<int*>(p), static_cast<int*>(p) + num_vals);
        }

        if (number_type == "i8")
        {
            long* p = data<long>();
            return std::vector<T>(static_cast<long*>(p), static_cast<long*>(p) + num_vals);
        }

        // if (number_type == "c8")
        // {
        //     _Complex float* p = data<_Complex float>();
        //     return std::vector<T>(static_cast<_Complex float*>(p), static_cast<_Complex float*>(p) + num_vals);
        // }
        // if (number_type == "c16")
        // {
        //     std::complex<double>* p = data<std::complex<double>>();
        //     return std::vector<T>(static_cast<std::complex<double>*>(p), static_cast<std::complex<double>*>(p) + num_vals);
        // }
        // if (number_type == "c32")
        // {
        //     std::complex<long double>* p = data<std::complex<long double>>();
        //     return std::vector<T>(static_cast<std::complex<long double>*>(p), static_cast<std::complex<long double>*>(p) + num_vals);
        // }
        if (number_type == "u2")
        {
            unsigned short* p = data<unsigned short>();
            return std::vector<T>(static_cast<unsigned short*>(p),
                                  static_cast<unsigned short*>(p) + num_vals);
        }

        if (number_type == "u4")
        {
            unsigned int* p = data<unsigned int>();
            return std::vector<T>(static_cast<unsigned int*>(p),
                                  static_cast<unsigned int*>(p) + num_vals);
        }

        if (number_type == "u8")
        {
            unsigned long* p = data<unsigned long>();
            return std::vector<T>(static_cast<unsigned long*>(p),
                                  static_cast<unsigned long*>(p) + num_vals);
        }

        if (number_type == "f8")
        {
            double* p = data<double>();
            return std::vector<T>(static_cast<double*>(p),
                                  static_cast<double*>(p) + num_vals);
        }
    }

    size_t num_bytes() const
    {
        return data_holder->size();
    }

    std::shared_ptr<std::vector<char>> data_holder;
    std::vector<size_t> shape;
    size_t word_size;
    bool fortran_order;
    std::string number_type;
    size_t num_vals;
};

using npz_t = std::map<std::string, NpyArray>;

char BigEndianTest();
char map_type(const std::type_info& t);
template<typename T> std::vector<char> create_npy_header(
    const std::vector<size_t>& shape);
void parse_npy_header(FILE* fp, size_t& word_size, std::vector<size_t>& shape,
                      bool& fortran_order, std::string& number_type);
void parse_npy_header(unsigned char* buffer, size_t& word_size,
                      std::vector<size_t>& shape, bool& fortran_order, std::string& number_type);
void parse_zip_footer(FILE* fp, uint16_t& nrecs, size_t& global_header_size,
                      size_t& global_header_offset);
npz_t npz_load(std::string fname);
NpyArray npz_load(std::string fname, std::string varname);
NpyArray npy_load(std::string fname);

template<typename T> std::vector<char>& operator+=(std::vector<char>& lhs,
        const T rhs)
{
    //write in little endian
    for (size_t byte = 0; byte < sizeof(T); byte++)
    {
        char val = *((char*)&rhs + byte);
        lhs.push_back(val);
    }

    return lhs;
}

template<> std::vector<char>& operator+=(std::vector<char>& lhs,
        const std::string rhs);
template<> std::vector<char>& operator+=(std::vector<char>& lhs,
        const char* rhs);


template<typename T> void npy_save(std::string fname, const T* data,
                                   const std::vector<size_t> shape, std::string mode = "w")
{
    FILE* fp = NULL;
    std::vector<size_t>
    true_data_shape; //if appending, the shape of existing + new data

    if (mode == "a")
    {
        fp = fopen(fname.c_str(), "r+b");
    }

    if (fp)
    {
        //file exists. we need to append to it. read the header, modify the array size
        size_t word_size;
        bool fortran_order;
        std::string number_type;
        parse_npy_header(fp, word_size, true_data_shape, fortran_order, number_type);
        assert(!fortran_order);

        if (word_size != sizeof(T))
        {
            std::cout << "libnpy error: " << fname << " has word size " << word_size <<
                      " but npy_save appending data sized " << sizeof(T) << "\n";
            assert( word_size == sizeof(T) );
        }

        if (true_data_shape.size() != shape.size())
        {
            std::cout <<
                      "libnpy error: npy_save attempting to append misdimensioned data to " << fname
                      << "\n";
            assert(true_data_shape.size() != shape.size());
        }

        for (size_t i = 1; i < shape.size(); i++)
        {
            if (shape[i] != true_data_shape[i])
            {
                std::cout << "libnpy error: npy_save attempting to append misshaped data to " <<
                          fname << "\n";
                assert(shape[i] == true_data_shape[i]);
            }
        }

        true_data_shape[0] += shape[0];
    }
    else
    {
        fp = fopen(fname.c_str(), "wb");
        true_data_shape = shape;
    }

    std::vector<char> header = create_npy_header<T>(true_data_shape);
    size_t nels = std::accumulate(shape.begin(), shape.end(), 1,
                                  std::multiplies<size_t>());
    fseek(fp, 0, SEEK_SET);
    fwrite(&header[0], sizeof(char), header.size(), fp);
    fseek(fp, 0, SEEK_END);
    fwrite(data, sizeof(T), nels, fp);
    fclose(fp);
}

template<typename T> void npz_save(std::string zipname, std::string fname,
                                   const T* data, const std::vector<size_t>& shape, std::string mode = "w")
{
    //first, append a .npy to the fname
    fname += ".npy";
    //now, on with the show
    FILE* fp = NULL;
    uint16_t nrecs = 0;
    size_t global_header_offset = 0;
    std::vector<char> global_header;

    if (mode == "a")
    {
        fp = fopen(zipname.c_str(), "r+b");
    }

    if (fp)
    {
        //zip file exists. we need to add a new npy file to it.
        //first read the footer. this gives us the offset and size of the global header
        //then read and store the global header.
        //below, we will write the the new data at the start of the global header then append the global header and footer below it
        size_t global_header_size;
        parse_zip_footer(fp, nrecs, global_header_size, global_header_offset);
        fseek(fp, global_header_offset, SEEK_SET);
        global_header.resize(global_header_size);
        size_t res = fread(&global_header[0], sizeof(char), global_header_size, fp);

        if (res != global_header_size)
        {
            throw std::runtime_error("npz_save: header read error while adding to existing zip");
        }

        fseek(fp, global_header_offset, SEEK_SET);
    }
    else
    {
        fp = fopen(zipname.c_str(), "wb");
    }

    std::vector<char> npy_header = create_npy_header<T>(shape);
    size_t nels = std::accumulate(shape.begin(), shape.end(), 1,
                                  std::multiplies<size_t>());
    size_t nbytes = nels * sizeof(T) + npy_header.size();
    //get the CRC of the data to be added
    uint32_t crc = crc32(0L, (uint8_t*)&npy_header[0], npy_header.size());
    crc = crc32(crc, (uint8_t*)data, nels * sizeof(T));
    //build the local header
    std::vector<char> local_header;
    local_header += "PK"; //first part of sig
    local_header += (uint16_t) 0x0403; //second part of sig
    local_header += (uint16_t) 20; //min version to extract
    local_header += (uint16_t) 0; //general purpose bit flag
    local_header += (uint16_t) 0; //compression method
    local_header += (uint16_t) 0; //file last mod time
    local_header += (uint16_t) 0;     //file last mod date
    local_header += (uint32_t) crc; //crc
    local_header += (uint32_t) nbytes; //compressed size
    local_header += (uint32_t) nbytes; //uncompressed size
    local_header += (uint16_t) fname.size(); //fname length
    local_header += (uint16_t) 0; //extra field length
    local_header += fname;
    //build global header
    global_header += "PK"; //first part of sig
    global_header += (uint16_t) 0x0201; //second part of sig
    global_header += (uint16_t) 20; //version made by
    global_header.insert(global_header.end(), local_header.begin() + 4,
                         local_header.begin() + 30);
    global_header += (uint16_t) 0; //file comment length
    global_header += (uint16_t) 0; //disk number where file starts
    global_header += (uint16_t) 0; //internal file attributes
    global_header += (uint32_t) 0; //external file attributes
    global_header += (uint32_t)
                     global_header_offset; //relative offset of local file header, since it begins where the global header used to begin
    global_header += fname;
    //build footer
    std::vector<char> footer;
    footer += "PK"; //first part of sig
    footer += (uint16_t) 0x0605; //second part of sig
    footer += (uint16_t) 0; //number of this disk
    footer += (uint16_t) 0; //disk where footer starts
    footer += (uint16_t) (nrecs + 1); //number of records on this disk
    footer += (uint16_t) (nrecs + 1); //total number of records
    footer += (uint32_t) global_header.size(); //nbytes of global headers
    footer += (uint32_t) (global_header_offset + nbytes +
                          local_header.size()); //offset of start of global headers, since global header now starts after newly written array
    footer += (uint16_t) 0; //zip file comment length
    //write everything
    fwrite(&local_header[0], sizeof(char), local_header.size(), fp);
    fwrite(&npy_header[0], sizeof(char), npy_header.size(), fp);
    fwrite(data, sizeof(T), nels, fp);
    fwrite(&global_header[0], sizeof(char), global_header.size(), fp);
    fwrite(&footer[0], sizeof(char), footer.size(), fp);
    fclose(fp);
}

template<typename T> void npy_save(std::string fname, const std::vector<T> data,
                                   std::string mode = "w")
{
    std::vector<size_t> shape;
    shape.push_back(data.size());
    npy_save(fname, &data[0], shape, mode);
}

template<typename T> void npz_save(std::string zipname, std::string fname,
                                   const std::vector<T> data, std::string mode = "w")
{
    std::vector<size_t> shape;
    shape.push_back(data.size());
    npz_save(zipname, fname, &data[0], shape, mode);
}

template<typename T> std::vector<char> create_npy_header(
    const std::vector<size_t>& shape)
{
    std::vector<char> dict;
    dict += "{'descr': '";
    dict += BigEndianTest();
    dict += map_type(typeid(T));
    dict += std::to_string(sizeof(T));
    dict += "', 'fortran_order': False, 'shape': (";
    dict += std::to_string(shape[0]);

    for (size_t i = 1; i < shape.size(); i++)
    {
        dict += ", ";
        dict += std::to_string(shape[i]);
    }

    if (shape.size() == 1)
    {
        dict += ",";
    }

    dict += "), }";
    //pad with spaces so that preamble+dict is modulo 16 bytes. preamble is 10 bytes. dict needs to end with \n
    int remainder = 16 - (10 + dict.size()) % 16;
    dict.insert(dict.end(), remainder, ' ');
    dict.back() = '\n';
    std::vector<char> header;
    header += (char) 0x93;
    header += "NUMPY";
    header += (char) 0x01; //major version of numpy format
    header += (char) 0x00; //minor version of numpy format
    header += (uint16_t) dict.size();
    header.insert(header.end(), dict.begin(), dict.end());
    return header;
}

template<typename T, int dim>
Eigen::Matrix < T, -1, dim > load(Eigen::Matrix < T, -1, dim > & mat,
                                  const std::string fname, std::string order = "rowMajor");

template<typename T>
Eigen::SparseMatrix<T> load(Eigen::SparseMatrix<T>& mat,
                            const std::string fname);

template<typename T>
Eigen::Tensor<T, 3> load(Eigen::Tensor<T, 3>& tens,
                         const std::string fname, std::string order = "rowMajor");

template<typename T, int dim>
void save(const Eigen::Matrix < T, -1, dim > & mat, const std::string fname);

template<typename T>
void save(const Eigen::Tensor<T, 3>& tensor, const std::string fname);

template<typename T>
void save(const Eigen::SparseMatrix<T>& mat, const std::string fname);

cnpy::NpyArray load_the_npy_file(FILE* fp);
cnpy::NpyArray load_the_npy_file(std::string fname);

}

#endif
