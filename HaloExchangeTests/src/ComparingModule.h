#pragma once

#include <iostream>

/** 
* @file ComparingModule.h
* Abstract class that defines which methodes a Comparator must provide.
* @author Florian Scheidegger
* @date 01/02/2011
*/

class ComparingModule
{

    public:
        
        /**
        * Default constructor 
        */
        ComparingModule(){ relativePrecision_  = -1.0; }

        /** 
        * function that checks wheater the relativePrecision is setted for this module
        * @return bool - true if the relativePrecision was setted
        */
        bool IsRelativePrecisionSetted(){ return relativePrecision_ != -1.0; }

        /**
        * Sets the relativepPrecision for this module.
        * @param double relativePrecision - the relativePrecision on which this modules basic compare functions works
        */
        void SetRelativePrecision(double relativePrecision){ relativePrecision_ = relativePrecision; }

        /**
        * Function that compares two double values according to the relative error. If the relative error is smaller as relativePrecision of this module, the both arguments will be considered as equal.
        * @param double demanded - a value that acts as reference value
        * @param double actual - a value that will be compared to the reference value
        * @return bool - true if the values are equal or nearly-equal, false otherwise
        */
        bool BasicRelativeCompare(double demanded, double actual);

        /**
        * Function that compares two double values according to the absolute error. If the absolute error is smaller as relativePrecision of this module, the both arguments will be considered as equal.
        * @param double demanded - a value that acts as reference value
        * @param double actual - a value that will be compared to the reference value
        * @return bool - true if the values are equal or nearly-equal, false otherwise
        */
        bool BasicAbsoluteCompare(double demanded, double actual);

        /**
        * It compares two values, it in most cases it returns the same result as BasicRelativeCompare(...) expect that values near zero are compared with BasicAbsoluteCompare(...).
        * @param double demanded - a value that acts as reference value
        * @param double actual - a value that will be compared to the reference value
        * @return bool - true if the values are equal or nearly-equal, false otherwise
        */
        bool BasicCompare(double demanded, double actual);

        /**
        * @return std::string - returns the selfgenerated name of that module, of this string consists of a constant name and some inital settings.
        */
        virtual std::string Name() = 0;

        /**
        * This function is called before the module is used, giving the derived module a chance to setup its inner state.
        */
        virtual void SetUp() = 0;

        /**
        * Compares two values at the Index specified by i, j and k in the datafield.
        * @param int i - the Index of the i-counter
        * @param int j - the Index of the j-counter
        * @param int k - the Index of the k-counter
        * @param double demanded - a value that acts as reference value
        * @param double actual - a value that will be compared to the reference value
        */
        virtual void Compare(int i, int j, int k, double* demanded, double* actual) = 0;

        /**
        * After all compares are done, the tester asks the module wheater it overall passed or not.
        * @return bool - true if all compares were okey, false otherwise.
        */
        virtual bool Passed() = 0;

        /**
        * @return std::string - a message that specifies the computed result, this methode will be called after the test. If passed() returns true, in some mostly it makes sense that this function will return a empty string in this case. Otherwise it returns a explication why the modul has faild and/or some additional information. 
        */
        virtual std::string getResult() = 0;


    protected:

        double relativePrecision_;

};
  
