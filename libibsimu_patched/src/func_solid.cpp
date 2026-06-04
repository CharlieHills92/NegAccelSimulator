#include <iostream>
#include "func_solid.hpp"
#include "ibsimu.hpp"
#include "error.hpp"

// Modified by C. Poggi on 10/02/2025 to include use of Lambda functions, for creating classes of geometries

bool FuncSolid::inside(const Vec3D &x) const
{
    if (!_func)
        throw Error(ERROR_LOCATION, "solid function not defined");

    // Transform 3D -> 3D
    Vec3D y = _T.transform_point(x);
    return _func(y[0], y[1], y[2]);
}

FuncSolid::FuncSolid(std::istream &is)
    : _func(nullptr)  // Initialize with nullptr
{
    ibsimu.message(MSG_WARNING, 1) << "Warning: loading of FuncSolid not implemented\n";
    ibsimu.flush();
}

void FuncSolid::debug_print(std::ostream &os) const
{
    os << "**FuncSolid\n";
    os << "func = " << (_func ? "Defined" : "NULL") << "\n";
}

void FuncSolid::save(std::ostream &os) const
{
    write_int32(os, FILEID_FUNCSOLID);
    ibsimu.message(MSG_WARNING, 1) << "Warning: saving of FuncSolid not implemented\n";
    ibsimu.flush();
}
