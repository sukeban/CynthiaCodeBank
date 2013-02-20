//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	ComplexNumber.h
//
//		a useful complex number class
//          from http://developer.apple.com/library/mac/#samplecode/AUPinkNoise/Listings/Utility_ComplexNumber_h.html
//
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef __CoreAudio_ComplexNumber
#define __CoreAudio_ComplexNumber

#include <math.h>
#include <stdio.h>

class Complex
{
public:
	Complex()
		: mReal(0.0), mImag(0.0) {};
	
	Complex(float inReal, float inImag )
		: mReal(inReal), mImag(inImag) {};

	Complex(float inReal)				// construct complex from real
		: mReal(inReal), mImag(0) {};
		


	inline float			GetReal() const {return mReal;};
	inline float			GetImag() const {return mImag;};
	
	void			SetReal(float inReal) {mReal = inReal;};
	void			SetImag(float inImag) {mImag = inImag;};
	
	float			Phase() const {return atan2(mImag, mReal);};
	float			GetPhase() const {return atan2(mImag, mReal);};
	float			Magnitude() const {return sqrt(mImag*mImag + mReal*mReal);};
	float			GetMagnitude() const {return sqrt(mImag*mImag + mReal*mReal);};
		
	void			SetMagnitudePhase(float inMagnitude, float inPhase)
	{
		mReal = inMagnitude * cos(inPhase);
		mImag = inMagnitude * sin(inPhase);
	};

	
	Complex			Pow(float inPower)
	{
		float mag = GetMagnitude();
		float phase = GetPhase();
		
		Complex result;
		result.SetMagnitudePhase(pow(mag, inPower), phase*inPower );
		
		return result;
	};
	
	Complex			GetConjugate() const {return Complex(mReal, -mImag);};
	
	
	Complex			inline operator += (const Complex &a);
	Complex			inline operator -= (const Complex &a);


	void			Print() {printf("(%f,%f)", mReal, mImag ); };
	void			PrintMagnitudePhase() {printf("(%f,%f)\n", GetMagnitude(), GetPhase() ); };
	
	
	float			mReal;
	float			mImag;
};

Complex			inline operator+ (const Complex &a, const Complex &b )
	{return Complex(a.GetReal() + b.GetReal(), a.GetImag() + b.GetImag() ); };

Complex			inline operator - (const Complex &a, const Complex &b )
	{return Complex(a.GetReal() - b.GetReal(), a.GetImag() - b.GetImag() ); };
	
Complex			inline operator * (const Complex &a, const Complex &b )
	{return Complex(	a.GetReal()*b.GetReal() - a.GetImag()*b.GetImag(),
						a.GetReal()*b.GetImag() + a.GetImag()*b.GetReal() ); };
	
Complex			inline operator * (const Complex &a, float b)
	{return Complex(a.GetReal()*b, a.GetImag()*b );};
	
Complex			inline operator * (float b, const Complex &a )
	{return Complex(a.GetReal()*b, a.GetImag()*b );};
	
Complex			inline operator/(const Complex& a, const Complex& b)
{
	float mag1 = a.GetMagnitude();
	float mag2 = b.GetMagnitude();
	
	float phase1 = a.GetPhase();
	float phase2 = b.GetPhase();
	
	Complex c;
	c.SetMagnitudePhase(mag1/mag2, phase1 - phase2 );
	
	return c;
}

Complex			inline Complex::operator += (const Complex &a)
{
	*this = *this + a;
	return *this;
};

Complex			inline Complex::operator -= (const Complex &a)
{
	*this = *this - a;
	return *this;
};

bool			inline	operator == (const Complex &a, const Complex &b )
{
	return a.GetReal() == b.GetReal() && a.GetImag() == b.GetImag();
}

inline Complex		UnitCircle(float mag, float phase)
{
	return Complex(mag * cos(phase), mag * sin(phase) );
}

#endif // __ComplexNumber
