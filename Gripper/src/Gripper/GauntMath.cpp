#include <math/gripper_math.hpp>

#ifdef _WIN32
int round(double x)
{
	return int(x > 0.0 ? x + 0.5 : x - 0.5);
}
#endif


unsigned int nearestInferiorPowerOf2( unsigned int n )
{
        return (unsigned int) pow( 2, floor( log( (float)n ) / log( 2.0f ) ) );
}


int greatestCommonFactor(int a, int b){

    if(a > b)
    {
        if(a % b == 0) return b;
        else
        {
            b = a % b;
            greatestCommonFactor(a, b);
        }
    }
    else
    {
        if(b % a == 0) return a;
        else
        {
            a = b % a;
            greatestCommonFactor(a, b);
        }
    }
    return 0;
}