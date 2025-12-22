// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "core/io/fieldWriter.h"
#include <H5Cpp.h>
#include <iostream>

void write_fieldarray_hdf5(const std::string& path, const FieldArray& fld) {
    if (!fld.is_valid()) {
        throw std::runtime_error("Invalid FieldArray: cannot write to HDF5");
    }

    H5::H5File file(path, H5F_ACC_TRUNC);

    auto write_1d = [&](const std::string& name, const std::vector<double>& v) {
        hsize_t dims[1] = { v.size() };
        H5::DataSpace space(1, dims);
        H5::DataSet dset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
        dset.write(v.data(), H5::PredType::NATIVE_DOUBLE);
    };

    write_1d("x", fld.xs);
    write_1d("y", fld.ys);
    write_1d("z", fld.zs);

    hsize_t dims3[3] = { fld.nx, fld.ny, fld.nz };
    H5::DataSpace space3(3, dims3);

    H5::DataSet dEx = file.createDataSet("Ex", H5::PredType::NATIVE_DOUBLE, space3);
    H5::DataSet dEy = file.createDataSet("Ey", H5::PredType::NATIVE_DOUBLE, space3);
    H5::DataSet dEz = file.createDataSet("Ez", H5::PredType::NATIVE_DOUBLE, space3);
    H5::DataSet dPhi = file.createDataSet("phi", H5::PredType::NATIVE_DOUBLE, space3);

    dEx.write(fld.Ex.data(), H5::PredType::NATIVE_DOUBLE);
    dEy.write(fld.Ey.data(), H5::PredType::NATIVE_DOUBLE);
    dEz.write(fld.Ez.data(), H5::PredType::NATIVE_DOUBLE);
    dPhi.write(fld.phi.data(), H5::PredType::NATIVE_DOUBLE);

    // attributes
    {
        H5::Attribute a_nx = file.createAttribute("nx", H5::PredType::NATIVE_ULONG, H5::DataSpace(H5S_SCALAR));
        unsigned long nx = fld.nx;
        a_nx.write(H5::PredType::NATIVE_ULONG, &nx);
    }
    {
        H5::Attribute a_ny = file.createAttribute("ny", H5::PredType::NATIVE_ULONG, H5::DataSpace(H5S_SCALAR));
        unsigned long ny = fld.ny;
        a_ny.write(H5::PredType::NATIVE_ULONG, &ny);
    }
    {
        H5::Attribute a_nz = file.createAttribute("nz", H5::PredType::NATIVE_ULONG, H5::DataSpace(H5S_SCALAR));
        unsigned long nz = fld.nz;
        a_nz.write(H5::PredType::NATIVE_ULONG, &nz);
    }

    file.close();
    std::cout << "Wrote field HDF5: " << path << std::endl;
}
