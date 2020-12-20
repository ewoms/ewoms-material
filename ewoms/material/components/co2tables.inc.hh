/* Tables for CO2 fluid properties calculated according to Span and
 * Wagner (1996).
 *
 * THIS AN AUTO-GENERATED FILE! DO NOT EDIT IT!
 *
 * Temperature range: 280.000 K to 400.000 K, using 200 sampling points
 * Pressure range: 0.100 MPa to 100.000 MPa, using 500 sampling points
 *
 * Generated using:
 *
 * ./extractproperties 280.0 400.0 200 1e5 100e6 500
 */
struct TabulatedDensityTraits {
    typedef double Scalar;
    static constexpr char name[] = "density";
    static constexpr int numX = 200;
    static constexpr Scalar xMin = 2.800000000000000e+02;
    static constexpr Scalar xMax = 4.000000000000000e+02;
    static constexpr int    numY = 500;
    static constexpr Scalar yMin = 1.000000000000000e+05;
    static constexpr Scalar yMax = 1.000000000000000e+08;

    static const double vals[numX][numY];
};

struct TabulatedEnthalpyTraits {
    typedef double Scalar;
    static constexpr char name[] = "enthalpy";
    static constexpr int numX = 200;
    static constexpr Scalar xMin = 2.800000000000000e+02;
    static constexpr Scalar xMax = 4.000000000000000e+02;
    static constexpr int    numY = 500;
    static constexpr Scalar yMin = 1.000000000000000e+05;
    static constexpr Scalar yMax = 1.000000000000000e+08;
    static const double vals[numX][numY];
};

typedef Ewoms::UniformTabulated2DFunction< double > TabulatedFunction;

// this class collects all the tabulated quantities in one convenient place
struct CO2Tables {
    static TabulatedFunction tabulatedEnthalpy;
    static TabulatedFunction tabulatedDensity;
    static constexpr double brineSalinity = 1.000000000000000e-01;
};
