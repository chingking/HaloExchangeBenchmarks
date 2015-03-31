#pragma once

#include <fstream>
#include <cassert>
#include "ComparingModule.h"


/**
 *  @class GnuplotLogger generates a text gnuplotfile that draws 2D maps of boolean grids.
 *  The gnuplot file plots red cells where the grid point did not pass the validation.
 */
class GnuplotLogger
{
public:

    //default constructor, destructor
    GnuplotLogger() {}
    ~GnuplotLogger() {}

    /**
     * Method that initializes the GnuplotLogger module
     * @param saveDigits precision in number of digits
     * @param calculationDomain size of the calculation domain to run the validation
     * @param filename name of the file to write gnuplot
     * @param kmin minimum K level to produce gnu plots
     * @param kmax maximum K level to produce gnu plots
     * @param numRowPlots number of rows of plots for every gnuplot page
     * @param numColumnPlots number of columns of plots for every gnuplot page
     */
    void Init(int saveDigits, int iSize, int jSize, int kSize, std::string filename, int kmin, int kmax, int numRowPlots, int numColumnPlots);

    /**
     * @return description of the module
     */
    std::string Name();

    /**
     * Setups the module
     */
    void SetUp();

    /* Compares two values at the Index specified by i, j and k in the datafield.
     * @param int i - the Index of the i-counter
     * @param int j - the Index of the j-counter
     * @param int k - the Index of the k-counter
     * @param double demanded - a value that acts as reference value
     * @param double actual - a value that will be compared to the reference value
     */
    void Compare(int i, int j, int k, double* demanded, double* actual);

    /*
     * @return true if the module passes validation
     */
    bool Passed();

    void ConvertValidatedMatrix(double* valid);

    /*
     * It generates the gnuplot with the resulting grid map of booleans
     * @return string indicating the location of the gnuplot file
     */
    std::string getResult();

    void SetRelativePrecision(double relativePrecision){ relativePrecision_ = relativePrecision; }

private:

    bool BasicRelativeCompare(double demanded, double actual);
    bool BasicAbsoluteCompare(double demanded, double actual);
    bool BasicCompare(double demanded, double actual);


    int numLevels_;
    int iSize_;
    int jSize_;
    int kSize_;
    int kmin_;
    int kmax_;
    int numRowPlots_;
    int numColumnPlots_;
    bool passed_;
    std::string filename_;
    int* pValidatedGrid_;
    std::ofstream outFile_;
    bool passedValidation_;

    double relativePrecision_;

};

