#include <cmath>
#include <climits>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>


#include "GnuplotLogger.h"

void GnuplotLogger::Init(int saveDigits, int iSize, int jSize, int kSize, std::string filename, int kmin = 0, int kmax = 0, int numRowPlots = 2,
        int numColumnPlots = 3)
{

    assert(numRowPlots > 0);
    assert(numColumnPlots > 0);
    SetRelativePrecision(pow(10.0, -saveDigits));
    passed_ = true;

    assert(kmin >= 0 && kmin < kSize);
    assert(kmax >= 0 && kmax >= kmin && kmax <= kSize);

    iSize_ = iSize;
    jSize_ = jSize;
    kSize_ = kSize;
    pValidatedGrid_ = (int*)malloc( iSize_*jSize_*kSize*sizeof(int) );
    filename_ = filename+".gpi";
    kmin_ = kmin;
    kmax_ = kmax;
    numRowPlots_ = numRowPlots;
    numColumnPlots_ = numColumnPlots;

    numLevels_ = kmax_ - kmin_;
    if (kmin_ == 0)
    {
        numLevels_++;
    }

    if( numLevels_ < numRowPlots_ * numColumnPlots_)
    {
        bool inLoop = true;
        while(inLoop)
        {
            inLoop = false;
            if (numRowPlots_ > numColumnPlots_)
            {
                if (numLevels_ < (numRowPlots_ - 1) * numColumnPlots_)
                {
                    numRowPlots_--;
                    inLoop = true;
                }
            }
            else
            {
                if(numLevels_ < (numRowPlots_ * (numColumnPlots_-1)) )
                {
                    numColumnPlots_--;
                    inLoop = true;
                }
            }
        }
    }

    SetUp();
}

std::string GnuplotLogger::Name()
{
    std::stringstream buffer;
    buffer << "GnuplotLogger, Epsilon = " << relativePrecision_;
    return buffer.str();
}

void GnuplotLogger::SetUp()
{
    for (int i = 0; i < iSize_ ; ++i)
    {
        for (int j = 0; j < jSize_ ; ++j)
        {
            for (int k = 0; k < kSize_ ; ++k)
            {
                pValidatedGrid_[i+iSize_*(j+k*jSize_) ]  = -1;
            }
        }
    }

}

void GnuplotLogger::ConvertValidatedMatrix(double* valid)
{
  int asize = iSize_*jSize_*kSize_;
  for(int l=0; l < asize; ++l) 
  {
// value == 0 : correct inner domain
    if((int) valid[l] == 0 ) 
      pValidatedGrid_[l] = 0;
// value == 1 : correct halo
    else if((int) valid[l] == 1)
      pValidatedGrid_[l] = 1;
// value == 5 : wrong inner domain
    else if((int) valid[l] == 5)
      pValidatedGrid_[l] = 10;
// value == 6 : wrong halo
    else if((int) valid[l] == 6)
      pValidatedGrid_[l] = 6.5;
// value == 5 : uninitialized
    else if((int) valid[l] < 0 )
      pValidatedGrid_[l] = 5;

    if( (int) valid[l] != 0 && (int) valid[l] != 1)
      passed_ = false;
   
  }
}

void GnuplotLogger::Compare(int i, int j, int k, double* demanded, double* actual)
{

    if (BasicCompare(demanded[ i+iSize_*(j+k*jSize_)  ], actual[ i+iSize_*(j+k*jSize_)  ]) == false)
    {
        passed_ = false;
        pValidatedGrid_[i+iSize_*(j+k*jSize_) ] = 3;
    }
    else
    {
        pValidatedGrid_[i+iSize_*(j+k*jSize_) ] = 1;
    }
}

bool GnuplotLogger::Passed()
{
    return passed_;
}

std::string GnuplotLogger::getResult()
{

  if( ! passed_ ) {

    outFile_.open(filename_.c_str(), std::ios::out);

    if (!outFile_)
    {
        throw std::runtime_error("MatlabLogger::Init(...) Can't open file\n");
    }

    if (numLevels_ > 1)
    {
        outFile_ << "set xrange [-1:1]" << std::endl;
        outFile_ << "set yrange [-1:1]" << std::endl;
        outFile_ << "set size 1,1" << std::endl;
        outFile_ << "set origin 0,0" << std::endl;
        outFile_ << "set object 2 rect from  screen 0.87, 0.94 to  screen 0.88, 0.96 fc ls 6  fs bo 1" << std::endl;
        outFile_ << "set label \"Uninitialized\" at screen 0.885, 0.952 font \",9\"" << std::endl;
        outFile_ << "set object 3 rect from  screen 0.87, 0.9 to  screen 0.88, 0.92 fc ls 3  fs bo 1" << std::endl;
        outFile_ << "set label \"Correct\" at screen 0.885, 0.912 font \",9\"" << std::endl;
        outFile_ << "set object 4 rect from  screen 0.87, 0.86 to  screen 0.88, 0.88 fc ls 1  fs bo 1" << std::endl;
        outFile_ << "set label \"Incorrect\" at screen 0.885, 0.872 font \",9\"" << std::endl;
        outFile_ << "set multiplot" << std::endl;
    }
    int kInPage = 0;
    int kInRow = 0;
    int kInColumn = 0;
    float xOrigin = 0;
    float yOrigin = 0;
    float xSize = 0;
    float ySize = 0;
    std::stringstream ss;
    for (int k = kmin_; k < kmax_ + 1; ++k)
    {
        kInPage = (k - kmin_);

        kInPage = kInPage % (numRowPlots_ * numColumnPlots_);
        kInPage++;
        if (kInPage == 1 && k != kmin_)
        {
            outFile_ << "unset multiplot" << std::endl;
            outFile_ << "pause -1 \"Hit return to continue\"" << std::endl;
            outFile_ << "set xrange [-1:1]" << std::endl;
            outFile_ << "set size 1,1" << std::endl;
            outFile_ << "set origin 0,0" << std::endl;
            outFile_ << "set multiplot" << std::endl;
        }

        kInRow = (kInPage - 1) / numColumnPlots_ + 1;
        kInColumn = (kInPage - 1) % numColumnPlots_ + 1;

        float ratio = iSize_ / (float)(jSize_) ;

        yOrigin = 1 / (float) (numRowPlots_) * (numRowPlots_ - kInRow);
        xOrigin = 1 / (float) (numColumnPlots_) * (kInColumn - 1);
        xSize = 1 / (float) (numRowPlots_);
        ySize = 1 / (float) (numColumnPlots_);

        if (ySize / xSize < ratio)
        {
            xSize = ySize * 1 / ratio;

        }
        else
        {
            ySize = xSize * ratio;

        }

        if (numRowPlots_ * numColumnPlots_ > 5 && numRowPlots_ * numColumnPlots_ < 8)
        {
            xSize = xSize * 1.15;

            ySize = ySize * 1.15;
        }
        if (numRowPlots_ * numColumnPlots_ > 7)
        {
            xSize = xSize * 1.25;
            ySize = ySize * 1.25;
        }
        outFile_ << "set title \"Level " << k << "\"" << std::endl;
        outFile_ << "set size " << xSize << "," << ySize << std::endl;
        outFile_ << "set origin " << xOrigin << "," << yOrigin << std::endl;

        outFile_ << "unset key" << std::endl;
        outFile_ << "set tic scale 0" << std::endl;
        outFile_ << "set border 3 front" << std::endl;
        outFile_ << "unset colorbox" << std::endl;
        outFile_ << "set cbrange [-2:7]" << std::endl;
        outFile_ << "set cbtics 0,1,5" << std::endl;
        outFile_ << "set size ratio " << ratio << std::endl;
        float xrangemin = -0.5;
        float xrangemax = jSize_ - 0.5;
        outFile_ << "set xrange ["<<xrangemin<< ":" << xrangemax << "]" << std::endl;
        float yrangemin = -0.5;
        float yrangemax = iSize_ - 0.5;
        outFile_ << "set yrange ["<<yrangemin<< ":" << yrangemax << "] reverse" << std::endl;
        outFile_ << "unset ytics" << std::endl;
        outFile_ << "unset xtics" << std::endl;

        outFile_ << "set view map" << std::endl;

        outFile_ << "set title \"Level " << k << "\"" << std::endl;
        outFile_ << "set xtics" << std::endl;
        outFile_ << "set ytics" << std::endl;
        outFile_ << "set xlabel 'j'" << std::endl;
        outFile_ << "set ylabel 'i'" << std::endl;

        outFile_ << "set palette defined (0 \"blue\", 14 \"yellow\", 20 \"red\" )" << std::endl;
        outFile_ << "splot '-' matrix with image failsafe" << std::endl;
        for (int i = 0; i < iSize_ ; ++i)
        {
            for (int j = 0; j < jSize_ ; ++j)
            {
                //            for(int k=0; k < calculationDomain_.kSize() + 1; ++k)
                outFile_ << pValidatedGrid_[i+iSize_*(j+k*jSize_) ] << " ";
            }
            outFile_ << std::endl;
        }
        outFile_ << "e" << std::endl;
        outFile_ << "e" << std::endl;

    }
    if (numLevels_ > 1)
    {
        outFile_ << "unset multiplot" << std::endl;
    }
    outFile_ << "pause -1 \"Hit return to continue\"" << std::endl;

    outFile_.close();
  }
  

  if ( ! passed_ )
    return std::string("GnuplotLogger: failures detected written in gnuplot file " + filename_);
  else
    return std::string("GnuplotLogger: passed");
}


bool GnuplotLogger::BasicRelativeCompare(double demanded, double actual)
{
    return fabs(demanded - actual) < fabs(demanded*relativePrecision_);
}

bool GnuplotLogger::BasicAbsoluteCompare(double demanded, double actual)
{
    return fabs(demanded - actual) < relativePrecision_ ;
}

bool GnuplotLogger::BasicCompare(double demanded, double actual)
{
    const double absolutLimit = 1;
    if(fabs(demanded) < absolutLimit && fabs(actual) < absolutLimit)
    {
        return BasicAbsoluteCompare(demanded, actual);
    }
    else
    {
        return BasicRelativeCompare(demanded, actual);
    }
}

