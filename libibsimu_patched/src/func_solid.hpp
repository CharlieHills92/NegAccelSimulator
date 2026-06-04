#ifndef FUNC_SOLID_HPP
#define FUNC_SOLID_HPP 1

#include <iostream>
#include <functional>  // Include for std::function
#include "solid.hpp"

// Modified by C. Poggi on 10/02/2025 to include use of Lambda functions, for creating classes of geometries

/*! \brief Function solid class.
 *
 *  FuncSolid class holds the definition for one solid defining
 *  C-function. This solid implementation suffers from the inability
 *  of saving to file. If FuncSolid is constructed from stream the
 *  function pointer inside is set to NULL and error is thrown if
 *  function is evaluated using inside().
 */
class FuncSolid : public Solid {

    std::function<bool(double, double, double)> _func;  // Use std::function instead of function pointer

public:

    /*! \brief Constructor.
     */
    FuncSolid(std::function<bool(double, double, double)> func) : _func(func) {}

    /*! \brief Constructor for function pointers (backward compatibility).
     */
    FuncSolid(bool (*func)(double, double, double)) : _func(func) {}

    /*! \brief Constructor for loading solid data from stream \a is.
     */
    FuncSolid(std::istream &is);

    /*! \brief Destructor.
     */
    ~FuncSolid() {}

    /*! \brief Return if point x is inside funcsolid.
     */
    bool inside(const Vec3D &x) const;

    /*! \brief Print debugging information to \a os.
     */
    void debug_print(std::ostream &os) const;

    /*! \brief Saves solid data to stream \a os.
     */
    void save(std::ostream &os) const;
};

#endif
