#include <limits>
#include <cmath>
#include "ComparingModule.h"

/*void ComparingModule::Init(std::string name, double relativePrecision)
{
    name_ = name;
    relativePrecision_ = relativePrecision;
} */

bool ComparingModule::BasicRelativeCompare(double demanded, double actual)
{
    return fabs(demanded - actual) < fabs(demanded*relativePrecision_);
}

bool ComparingModule::BasicAbsoluteCompare(double demanded, double actual)
{
    return fabs(demanded - actual) < relativePrecision_ ;
}
        
bool ComparingModule::BasicCompare(double demanded, double actual)
{
    /*
    const double epsilon = fabs(demanded * 1e-11) + 1e-26;
    if(fabs(demanded) > 1e-15)
    {
        //double relativeError = fabs(actual-demanded) / fabs(demanded);
        //if(relativeError == 1.0)
        //{
        //    relativeError = 0.0;
        //}
        if(fabs(actual-demanded) <= 1e-8 * fabs(demanded))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return fabs(actual-demanded) < 1e-14;
    } 
*/
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
        
  
