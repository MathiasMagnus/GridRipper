#pragma once

// Gripper includes
#include <Gripper/stl/stlConfig.hpp>
#include <Gripper/stl/stlMathConstants.hpp>

// Standard C++ includes
#include <cmath>
#include <complex>


namespace math
{
	template <typename Floating>
    std::complex<Floating> Y_lms(Floating theta, Floating phi, int l, int m, int s)
	{
		using F = Floating;

		constexpr auto e = math::constants::e<F>;
		constexpr auto pi = math::constants::pi<F>;
		constexpr auto i = std::complex<F>{ static_cast<F>(0),
		                                    static_cast<F>(1) };

		constexpr int max = 256; // max index value
		switch(l*max*max+m*max+s)
		{
		case 0:	// {0,0,0}
			return (F)1/((F)2*std::sqrt(pi));
		case 65280:	// {1,-1,0}
			return (std::sqrt((F)3/((F)2*pi))*std::sin(theta))/((F)2*std::exp(i*phi));
		case 65281:	// {1,-1,1}
			return -(std::sqrt((F)3/pi)*std::pow(std::cos(theta/(F)2),(F)2))/((F)2*std::exp(i*phi));
		case 65536:	// {1,0,0}
			return (std::sqrt((F)3/pi)*std::cos(theta))/(F)2;
		case 65537:	// {1,0,1}
			return (std::sqrt((F)3/((F)2*pi))*std::sin(theta))/(F)2;
		case 65792:	// {1,1,0}
			return -(std::exp(i*phi)*std::sqrt((F)3/((F)2*pi))*std::sin(theta))/(F)2;
		case 65793:	// {1,1,1}
			return -(std::exp(i*phi)*std::sqrt((F)3/pi)*std::pow(std::sin(theta/(F)2),(F)2))/(F)2;
		case 130560:	// {2,-2,0}
			return (std::sqrt((F)15/((F)2*pi))*std::pow(std::sin(theta),(F)2))/((F)4*std::exp((F)2*i*phi));
		case 130561:	// {2,-2,1}
			return -((std::sqrt((F)5/pi)*std::pow(std::cos(theta/(F)2),(F)3)*std::sin(theta/(F)2))/std::exp((F)2*i*phi));
		case 130562:	// {2,-2,2}
			return (std::sqrt((F)5/pi)*std::pow(std::cos(theta/(F)2),(F)4))/((F)2*std::exp((F)2*i*phi));
		case 130816:	// {2,-1,0}
			return (std::sqrt((F)15/((F)2*pi))*std::sin((F)2*theta))/((F)4*std::exp(i*phi));
		case 130817:	// {2,-1,1}
			return -(std::sqrt((F)5/pi)*std::pow(std::cos(theta/(F)2),(F)2)*(-(F)1 + (F)2*std::cos(theta)))/((F)2*std::exp(i*phi));
		case 130818:	// {2,-1,2}
			return -((std::sqrt((F)5/pi)*std::pow(std::cos(theta/(F)2),(F)3)*std::sin(theta/(F)2))/std::exp(i*phi));
		case 131072:	// {2,0,0}
			return (std::sqrt((F)5/pi)*((F)1 + (F)3*std::cos((F)2*theta)))/(F)8;
		case 131073:	// {2,0,1}
			return (std::sqrt((F)15/((F)2*pi))*std::cos(theta)*std::sin(theta))/(F)2;
		case 131074:	// {2,0,2}
			return (std::sqrt((F)15/((F)2*pi))*std::pow(std::sin(theta),(F)2))/(F)4;
		case 131328:	// {2,1,0}
			return -(std::exp(i*phi)*std::sqrt((F)15/((F)2*pi))*std::cos(theta)*std::sin(theta))/(F)2;
		case 131329:	// {2,1,1}
			return -(std::exp(i*phi)*std::sqrt((F)5/pi)*((F)1 + (F)2*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)2))/(F)2;
		case 131330:	// {2,1,2}
			return -(std::exp(i*phi)*std::sqrt((F)5/pi)*std::pow(std::sin(theta/(F)2),(F)2)*std::sin(theta))/(F)2;
		case 131584:	// {2,2,0}
			return (std::exp((F)2*i*phi)*std::sqrt((F)15/((F)2*pi))*std::pow(std::sin(theta),(F)2))/(F)4;
		case 131585:	// {2,2,1}
			return (std::exp((F)2*i*phi)*std::sqrt((F)5/pi)*std::pow(std::sin(theta/(F)2),(F)2)*std::sin(theta))/(F)2;
		case 131586:	// {2,2,2}
			return (std::exp((F)2*i*phi)*std::sqrt((F)5/pi)*std::pow(std::sin(theta/(F)2),(F)4))/(F)2;
		case 195840:	// {3,-3,0}
			return (std::sqrt((F)35/pi)*std::pow(std::sin(theta),(F)3))/((F)8*std::exp((F)3*i*phi));
		case 195841:	// {3,-3,1}
			return (std::sqrt((F)105/pi)*(-(F)1 + std::cos(theta))*std::pow((F)1 + std::cos(theta),(F)2))/((F)16*std::exp((F)3*i*phi));
		case 195842:	// {3,-3,2}
			return (std::sqrt((F)21/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*std::sin(theta/(F)2))/std::exp((F)3*i*phi);
		case 195843:	// {3,-3,3}
			return -(std::sqrt((F)7/pi)*std::pow(std::cos(theta/(F)2),(F)6))/((F)2*std::exp((F)3*i*phi));
		case 196096:	// {3,-2,0}
			return (std::sqrt((F)105/((F)2*pi))*std::cos(theta)*std::pow(std::sin(theta),(F)2))/((F)4*std::exp((F)2*i*phi));
		case 196097:	// {3,-2,1}
			return -(std::sqrt((F)35/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*(-(F)1 + (F)3*std::cos(theta))*std::sin(theta/(F)2))/((F)2*std::exp((F)2*i*phi));
		case 196098:	// {3,-2,2}
			return (std::sqrt((F)7/pi)*std::pow(std::cos(theta/(F)2),(F)4)*(-(F)2 + (F)3*std::cos(theta)))/((F)2*std::exp((F)2*i*phi));
		case 196099:	// {3,-2,3}
			return (std::sqrt((F)21/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*std::sin(theta/(F)2))/std::exp((F)2*i*phi);
		case 196352:	// {3,-1,0}
			return (std::sqrt((F)21/pi)*((F)3 + (F)5*std::cos((F)2*theta))*std::sin(theta))/((F)16*std::exp(i*phi));
		case 196353:	// {3,-1,1}
			return -(std::sqrt((F)7/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)13 - (F)20*std::cos(theta) + (F)15*std::cos((F)2*theta)))/((F)16*std::exp(i*phi));
		case 196354:	// {3,-1,2}
			return -(std::sqrt((F)35/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*(-(F)1 + (F)3*std::cos(theta))*std::sin(theta/(F)2))/((F)2*std::exp(i*phi));
		case 196355:	// {3,-1,3}
			return (std::sqrt((F)105/pi)*(-(F)1 + std::cos(theta))*std::pow((F)1 + std::cos(theta),(F)2))/((F)16*std::exp(i*phi));
		case 196608:	// {3,0,0}
			return (std::sqrt((F)7/pi)*((F)3*std::cos(theta) + (F)5*std::cos((F)3*theta)))/(F)16;
		case 196609:	// {3,0,1}
			return (std::sqrt((F)21/pi)*((F)3 + (F)5*std::cos((F)2*theta))*std::sin(theta))/(F)16;
		case 196610:	// {3,0,2}
			return (std::sqrt((F)105/((F)2*pi))*std::cos(theta)*std::pow(std::sin(theta),(F)2))/(F)4;
		case 196611:	// {3,0,3}
			return (std::sqrt((F)35/pi)*std::pow(std::sin(theta),(F)3))/(F)8;
		case 196864:	// {3,1,0}
			return -(std::exp(i*phi)*std::sqrt((F)21/pi)*((F)3 + (F)5*std::cos((F)2*theta))*std::sin(theta))/(F)16;
		case 196865:	// {3,1,1}
			return -(std::exp(i*phi)*std::sqrt((F)7/pi)*((F)13 + (F)20*std::cos(theta) + (F)15*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)2))/(F)16;
		case 196866:	// {3,1,2}
			return -(std::exp(i*phi)*std::sqrt((F)35/((F)2*pi))*std::cos(theta/(F)2)*((F)1 + (F)3*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)3))/(F)2;
		case 196867:	// {3,1,3}
			return -(std::exp(i*phi)*std::sqrt((F)105/pi)*std::pow(std::sin(theta/(F)2),(F)2)*std::pow(std::sin(theta),(F)2))/(F)8;
		case 197120:	// {3,2,0}
			return (std::exp((F)2*i*phi)*std::sqrt((F)105/((F)2*pi))*std::cos(theta)*std::pow(std::sin(theta),(F)2))/(F)4;
		case 197121:	// {3,2,1}
			return (std::exp((F)2*i*phi)*std::sqrt((F)35/((F)2*pi))*std::cos(theta/(F)2)*((F)1 + (F)3*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)3))/(F)2;
		case 197122:	// {3,2,2}
			return (std::exp((F)2*i*phi)*std::sqrt((F)7/pi)*((F)2 + (F)3*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)2;
		case 197123:	// {3,2,3}
			return (std::exp((F)2*i*phi)*std::sqrt((F)21/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)4)*std::sin(theta))/(F)2;
		case 197376:	// {3,3,0}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)35/pi)*std::pow(std::sin(theta),(F)3))/(F)8;
		case 197377:	// {3,3,1}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)105/pi)*std::pow(std::sin(theta/(F)2),(F)2)*std::pow(std::sin(theta),(F)2))/(F)8;
		case 197378:	// {3,3,2}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)21/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)4)*std::sin(theta))/(F)2;
		case 197379:	// {3,3,3}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)7/pi)*std::pow(std::sin(theta/(F)2),(F)6))/(F)2;
		case 261120:	// {4,-4,0}
			return ((F)3*std::sqrt((F)35/((F)2*pi))*std::pow(std::sin(theta),(F)4))/((F)16*std::exp((F)4*i*phi));
		case 261121:	// {4,-4,1}
			return (-(F)3*std::sqrt((F)14/pi)*std::pow(std::cos(theta/(F)2),(F)5)*std::pow(std::sin(theta/(F)2),(F)3))/std::exp((F)4*i*phi);
		case 261122:	// {4,-4,2}
			return (-(F)3*std::sqrt((F)7/pi)*(-(F)1 + std::cos(theta))*std::pow((F)1 + std::cos(theta),(F)3))/((F)16*std::exp((F)4*i*phi));
		case 261123:	// {4,-4,3}
			return (-(F)3*std::pow((F)1 + std::cos(theta),(F)35)*std::sqrt(std::pow(std::sin(theta/(F)2),(F)2)))/((F)8*std::exp((F)4*i*phi)*std::sqrt(pi));
		case 261124:	// {4,-4,4}
			return ((F)3*std::pow(std::cos(theta/(F)2),(F)8))/((F)2*std::exp((F)4*i*phi)*std::sqrt(pi));
		case 261376:	// {4,-3,0}
			return ((F)3*std::sqrt((F)35/pi)*std::cos(theta)*std::pow(std::sin(theta),(F)3))/((F)8*std::exp((F)3*i*phi));
		case 261377:	// {4,-3,1}
			return (-(F)3*std::sqrt((F)7/pi)*std::pow(std::cos(theta/(F)2),(F)4)*(-(F)1 + (F)4*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)2*std::exp((F)3*i*phi));
		case 261378:	// {4,-3,2}
			return ((F)3*std::sqrt((F)7/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*(-(F)1 + (F)2*std::cos(theta))*std::sin(theta/(F)2))/std::exp((F)3*i*phi);
		case 261379:	// {4,-3,3}
			return (-(F)3*std::pow(std::cos(theta/(F)2),(F)6)*(-(F)3 + (F)4*std::cos(theta)))/((F)2*std::exp((F)3*i*phi)*std::sqrt(pi));
		case 261380:	// {4,-3,4}
			return (-(F)3*std::pow((F)1 + std::cos(theta),(F)35)*std::sqrt(std::pow(std::sin(theta/(F)2),(F)2)))/((F)8*std::exp((F)3*i*phi)*std::sqrt(pi));
		case 261632:	// {4,-2,0}
			return ((F)3*std::sqrt((F)5/((F)2*pi))*((F)5 + (F)7*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)2))/((F)16*std::exp((F)2*i*phi));
		case 261633:	// {4,-2,1}
			return (-(F)3*std::pow(std::cos(theta/(F)2),(F)3)*((F)6 - (F)7*std::cos(theta) + (F)7*std::cos((F)2*theta))*std::sin(theta/(F)2))/((F)2*std::exp((F)2*i*phi)*std::sqrt((F)2*pi));
		case 261634:	// {4,-2,2}
			return ((F)3*std::pow(std::cos(theta/(F)2),(F)4)*((F)9 - (F)14*std::cos(theta) + (F)7*std::cos((F)2*theta)))/((F)4*std::exp((F)2*i*phi)*std::sqrt(pi));
		case 261635:	// {4,-2,3}
			return ((F)3*std::sqrt((F)7/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*(-(F)1 + (F)2*std::cos(theta))*std::sin(theta/(F)2))/std::exp((F)2*i*phi);
		case 261636:	// {4,-2,4}
			return (-(F)3*std::sqrt((F)7/pi)*(-(F)1 + std::cos(theta))*std::pow((F)1 + std::cos(theta),(F)3))/((F)16*std::exp((F)2*i*phi));
		case 261888:	// {4,-1,0}
			return ((F)3*std::sqrt((F)5/pi)*((F)9*std::cos(theta) + (F)7*std::cos((F)3*theta))*std::sin(theta))/((F)32*std::exp(i*phi));
		case 261889:	// {4,-1,1}
			return (-(F)3*std::pow(std::cos(theta/(F)2),(F)2)*(-(F)15 + (F)30*std::cos(theta) - (F)21*std::cos((F)2*theta) + (F)14*std::cos((F)3*theta)))/((F)16*std::exp(i*phi)*std::sqrt(pi));
		case 261890:	// {4,-1,2}
			return (-(F)3*std::pow(std::cos(theta/(F)2),(F)3)*((F)6 - (F)7*std::cos(theta) + (F)7*std::cos((F)2*theta))*std::sin(theta/(F)2))/((F)2*std::exp(i*phi)*std::sqrt((F)2*pi));
		case 261891:	// {4,-1,3}
			return (-(F)3*std::sqrt((F)7/pi)*std::pow(std::cos(theta/(F)2),(F)4)*(-(F)1 + (F)4*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)2*std::exp(i*phi));
		case 261892:	// {4,-1,4}
			return (-(F)3*std::sqrt((F)14/pi)*std::pow(std::cos(theta/(F)2),(F)5)*std::pow(std::sin(theta/(F)2),(F)3))/std::exp(i*phi);
		case 262144:	// {4,0,0}
			return ((F)3*((F)9 + (F)20*std::cos((F)2*theta) + (F)35*std::cos((F)4*theta)))/((F)128*std::sqrt(pi));
		case 262145:	// {4,0,1}
			return ((F)3*std::sqrt((F)5/pi)*((F)9*std::cos(theta) + (F)7*std::cos((F)3*theta))*std::sin(theta))/(F)32;
		case 262146:	// {4,0,2}
			return ((F)3*std::sqrt((F)5/((F)2*pi))*((F)5 + (F)7*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)2))/(F)16;
		case 262147:	// {4,0,3}
			return ((F)3*std::sqrt((F)35/pi)*std::cos(theta)*std::pow(std::sin(theta),(F)3))/(F)8;
		case 262148:	// {4,0,4}
			return ((F)3*std::sqrt((F)35/((F)2*pi))*std::pow(std::sin(theta),(F)4))/(F)16;
		case 262400:	// {4,1,0}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)5/pi)*((F)9*std::cos(theta) + (F)7*std::cos((F)3*theta))*std::sin(theta))/(F)32;
		case 262401:	// {4,1,1}
			return (-(F)3*std::exp(i*phi)*((F)15 + (F)30*std::cos(theta) + (F)21*std::cos((F)2*theta) + (F)14*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)16*std::sqrt(pi));
		case 262402:	// {4,1,2}
			return (-(F)3*std::exp(i*phi)*std::cos(theta/(F)2)*((F)6 + (F)7*std::cos(theta) + (F)7*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)2*std::sqrt((F)2*pi));
		case 262403:	// {4,1,3}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)7/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)1 + (F)4*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)2;
		case 262404:	// {4,1,4}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)7/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)2)*std::pow(std::sin(theta),(F)3))/(F)4;
		case 262656:	// {4,2,0}
			return ((F)3*std::exp((F)2*i*phi)*std::sqrt((F)5/((F)2*pi))*((F)5 + (F)7*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)2))/(F)16;
		case 262657:	// {4,2,1}
			return ((F)3*std::exp((F)2*i*phi)*std::cos(theta/(F)2)*((F)6 + (F)7*std::cos(theta) + (F)7*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)2*std::sqrt((F)2*pi));
		case 262658:	// {4,2,2}
			return ((F)3*std::exp((F)2*i*phi)*((F)9 + (F)14*std::cos(theta) + (F)7*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)4))/((F)4*std::sqrt(pi));
		case 262659:	// {4,2,3}
			return (F)3*std::exp((F)2*i*phi)*std::sqrt((F)7/((F)2*pi))*std::cos(theta/(F)2)*((F)1 + (F)2*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)5);
		case 262660:	// {4,2,4}
			return ((F)3*std::exp((F)2*i*phi)*std::sqrt((F)7/pi)*std::pow(std::sin(theta/(F)2),(F)4)*std::pow(std::sin(theta),(F)2))/(F)4;
		case 262912:	// {4,3,0}
			return (-(F)3*std::exp((F)3*i*phi)*std::sqrt((F)35/pi)*std::cos(theta)*std::pow(std::sin(theta),(F)3))/(F)8;
		case 262913:	// {4,3,1}
			return (-(F)3*std::exp((F)3*i*phi)*std::sqrt((F)7/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)1 + (F)4*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)2;
		case 262914:	// {4,3,2}
			return -(F)3*std::exp((F)3*i*phi)*std::sqrt((F)7/((F)2*pi))*std::cos(theta/(F)2)*((F)1 + (F)2*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)5);
		case 262915:	// {4,3,3}
			return (-(F)3*std::exp((F)3*i*phi)*((F)3 + (F)4*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)6))/((F)2*std::sqrt(pi));
		case 262916:	// {4,3,4}
			return (-(F)3*std::exp((F)3*i*phi)*std::pow(std::sin(theta/(F)2),(F)6)*std::sin(theta))/std::sqrt((F)2*pi);
		case 263168:	// {4,4,0}
			return ((F)3*std::exp((F)4*i*phi)*std::sqrt((F)35/((F)2*pi))*std::pow(std::sin(theta),(F)4))/(F)16;
		case 263169:	// {4,4,1}
			return ((F)3*std::exp((F)4*i*phi)*std::sqrt((F)7/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)2)*std::pow(std::sin(theta),(F)3))/(F)4;
		case 263170:	// {4,4,2}
			return ((F)3*std::exp((F)4*i*phi)*std::sqrt((F)7/pi)*std::pow(std::sin(theta/(F)2),(F)4)*std::pow(std::sin(theta),(F)2))/(F)4;
		case 263171:	// {4,4,3}
			return ((F)3*std::exp((F)4*i*phi)*std::pow(std::sin(theta/(F)2),(F)6)*std::sin(theta))/std::sqrt((F)2*pi);
		case 263172:	// {4,4,4}
			return ((F)3*std::exp((F)4*i*phi)*std::pow(std::sin(theta/(F)2),(F)8))/((F)2*std::sqrt(pi));
		case 326400:	// {5,-5,0}
			return ((F)3*std::sqrt((F)77/pi)*std::pow(std::sin(theta),(F)5))/((F)32*std::exp((F)5*i*phi));
		case 326401:	// {5,-5,1}
			return -((std::sqrt((F)1155/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)6)*std::pow(std::sin(theta/(F)2),(F)4))/std::exp((F)5*i*phi));
		case 326402:	// {5,-5,2}
			return (std::sqrt((F)330/pi)*std::pow(std::cos(theta/(F)2),(F)7)*std::pow(std::sin(theta/(F)2),(F)3))/std::exp((F)5*i*phi);
		case 326403:	// {5,-5,3}
			return ((F)3*std::sqrt((F)55/pi)*(-(F)1 + std::cos(theta))*std::pow((F)1 + std::cos(theta),(F)4))/((F)64*std::exp((F)5*i*phi));
		case 326404:	// {5,-5,4}
			return (std::sqrt((F)55/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)9)*std::sin(theta/(F)2))/std::exp((F)5*i*phi);
		case 326405:	// {5,-5,5}
			return -(std::sqrt((F)11/pi)*std::pow(std::cos(theta/(F)2),(F)10))/((F)2*std::exp((F)5*i*phi));
		case 326656:	// {5,-4,0}
			return ((F)3*std::sqrt((F)385/((F)2*pi))*std::cos(theta)*std::pow(std::sin(theta),(F)4))/((F)16*std::exp((F)4*i*phi));
		case 326657:	// {5,-4,1}
			return (std::sqrt((F)231/pi)*std::sqrt((F)1 - std::cos(theta))*std::pow((F)1 + std::cos(theta),(F)25)*((F)1 - (F)6*std::cos(theta) + (F)5*std::pow(std::cos(theta),(F)2)))/((F)32*std::exp((F)4*i*phi));
		case 326658:	// {5,-4,2}
			return (std::sqrt((F)33/pi)*std::pow(std::cos(theta/(F)2),(F)6)*(-(F)2 + (F)5*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)2))/std::exp((F)4*i*phi);
		case 326659:	// {5,-4,3}
			return (-(F)3*std::sqrt((F)11/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)7)*(-(F)3 + (F)5*std::cos(theta))*std::sin(theta/(F)2))/((F)2*std::exp((F)4*i*phi));
		case 326660:	// {5,-4,4}
			return (std::sqrt((F)11/pi)*std::pow(std::cos(theta/(F)2),(F)8)*(-(F)4 + (F)5*std::cos(theta)))/((F)2*std::exp((F)4*i*phi));
		case 326661:	// {5,-4,5}
			return (std::sqrt((F)55/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)9)*std::sin(theta/(F)2))/std::exp((F)4*i*phi);
		case 326912:	// {5,-3,0}
			return (std::sqrt((F)385/pi)*((F)7 + (F)9*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)3))/((F)64*std::exp((F)3*i*phi));
		case 326913:	// {5,-3,1}
			return -(std::sqrt((F)231/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)4)*((F)13 - (F)12*std::cos(theta) + (F)15*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)8*std::exp((F)3*i*phi));
		case 326914:	// {5,-3,2}
			return (std::sqrt((F)33/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*((F)17 - (F)24*std::cos(theta) + (F)15*std::cos((F)2*theta))*std::sin(theta/(F)2))/((F)4*std::exp((F)3*i*phi));
		case 326915:	// {5,-3,3}
			return -(std::sqrt((F)11/pi)*std::pow(std::cos(theta/(F)2),(F)6)*((F)71 - (F)108*std::cos(theta) + (F)45*std::cos((F)2*theta)))/((F)16*std::exp((F)3*i*phi));
		case 326916:	// {5,-3,4}
			return (-(F)3*std::sqrt((F)11/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)7)*(-(F)3 + (F)5*std::cos(theta))*std::sin(theta/(F)2))/((F)2*std::exp((F)3*i*phi));
		case 326917:	// {5,-3,5}
			return ((F)3*std::sqrt((F)55/pi)*(-(F)1 + std::cos(theta))*std::pow((F)1 + std::cos(theta),(F)4))/((F)64*std::exp((F)3*i*phi));
		case 327168:	// {5,-2,0}
			return (std::sqrt((F)1155/((F)2*pi))*((F)5*std::cos(theta) + (F)3*std::cos((F)3*theta))*std::pow(std::sin(theta),(F)2))/((F)32*std::exp((F)2*i*phi));
		case 327169:	// {5,-2,1}
			return -(std::sqrt((F)77/pi)*std::pow(std::cos(theta/(F)2),(F)3)*(-(F)14 + (F)33*std::cos(theta) - (F)18*std::cos((F)2*theta) + (F)15*std::cos((F)3*theta))*std::sin(theta/(F)2))/((F)16*std::exp((F)2*i*phi));
		case 327170:	// {5,-2,2}
			return (std::sqrt((F)11/pi)*std::pow(std::cos(theta/(F)2),(F)4)*(-(F)32 + (F)57*std::cos(theta) - (F)36*std::cos((F)2*theta) + (F)15*std::cos((F)3*theta)))/((F)8*std::exp((F)2*i*phi));
		case 327171:	// {5,-2,3}
			return (std::sqrt((F)33/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*((F)17 - (F)24*std::cos(theta) + (F)15*std::cos((F)2*theta))*std::sin(theta/(F)2))/((F)4*std::exp((F)2*i*phi));
		case 327172:	// {5,-2,4}
			return (std::sqrt((F)33/pi)*std::pow(std::cos(theta/(F)2),(F)6)*(-(F)2 + (F)5*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)2))/std::exp((F)2*i*phi);
		case 327173:	// {5,-2,5}
			return (std::sqrt((F)330/pi)*std::pow(std::cos(theta/(F)2),(F)7)*std::pow(std::sin(theta/(F)2),(F)3))/std::exp((F)2*i*phi);
		case 327424:	// {5,-1,0}
			return (std::sqrt((F)165/((F)2*pi))*((F)15 + (F)28*std::cos((F)2*theta) + (F)21*std::cos((F)4*theta))*std::sin(theta))/((F)128*std::exp(i*phi));
		case 327425:	// {5,-1,1}
			return -(std::sqrt((F)11/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)155 - (F)280*std::cos(theta) + (F)252*std::cos((F)2*theta) - (F)168*std::cos((F)3*theta) + (F)105*std::cos((F)4*theta)))/((F)128*std::exp(i*phi));
		case 327426:	// {5,-1,2}
			return -(std::sqrt((F)77/pi)*std::pow(std::cos(theta/(F)2),(F)3)*(-(F)14 + (F)33*std::cos(theta) - (F)18*std::cos((F)2*theta) + (F)15*std::cos((F)3*theta))*std::sin(theta/(F)2))/((F)16*std::exp(i*phi));
		case 327427:	// {5,-1,3}
			return -(std::sqrt((F)231/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)4)*((F)13 - (F)12*std::cos(theta) + (F)15*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)8*std::exp(i*phi));
		case 327428:	// {5,-1,4}
			return (std::sqrt((F)231/pi)*std::sqrt((F)1 - std::cos(theta))*std::pow((F)1 + std::cos(theta),(F)25)*((F)1 - (F)6*std::cos(theta) + (F)5*std::pow(std::cos(theta),(F)2)))/((F)32*std::exp(i*phi));
		case 327429:	// {5,-1,5}
			return -((std::sqrt((F)1155/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)6)*std::pow(std::sin(theta/(F)2),(F)4))/std::exp(i*phi));
		case 327680:	// {5,0,0}
			return (std::sqrt((F)11/pi)*((F)30*std::cos(theta) + (F)35*std::cos((F)3*theta) + (F)63*std::cos((F)5*theta)))/(F)256;
		case 327681:	// {5,0,1}
			return (std::sqrt((F)165/((F)2*pi))*((F)15 + (F)28*std::cos((F)2*theta) + (F)21*std::cos((F)4*theta))*std::sin(theta))/(F)128;
		case 327682:	// {5,0,2}
			return (std::sqrt((F)1155/((F)2*pi))*((F)5*std::cos(theta) + (F)3*std::cos((F)3*theta))*std::pow(std::sin(theta),(F)2))/(F)32;
		case 327683:	// {5,0,3}
			return (std::sqrt((F)385/pi)*((F)7 + (F)9*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)3))/(F)64;
		case 327684:	// {5,0,4}
			return ((F)3*std::sqrt((F)385/((F)2*pi))*std::cos(theta)*std::pow(std::sin(theta),(F)4))/(F)16;
		case 327685:	// {5,0,5}
			return ((F)3*std::sqrt((F)77/pi)*std::pow(std::sin(theta),(F)5))/(F)32;
		case 327936:	// {5,1,0}
			return -(std::exp(i*phi)*std::sqrt((F)165/((F)2*pi))*((F)15 + (F)28*std::cos((F)2*theta) + (F)21*std::cos((F)4*theta))*std::sin(theta))/(F)128;
		case 327937:	// {5,1,1}
			return -(std::exp(i*phi)*std::sqrt((F)11/pi)*((F)155 + (F)280*std::cos(theta) + (F)252*std::cos((F)2*theta) + (F)168*std::cos((F)3*theta) + (F)105*std::cos((F)4*theta))*std::pow(std::sin(theta/(F)2),(F)2))/(F)128;
		case 327938:	// {5,1,2}
			return -(std::exp(i*phi)*std::sqrt((F)77/pi)*std::cos(theta/(F)2)*((F)14 + (F)33*std::cos(theta) + (F)18*std::cos((F)2*theta) + (F)15*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)3))/(F)16;
		case 327939:	// {5,1,3}
			return -(std::exp(i*phi)*std::sqrt((F)231/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)13 + (F)12*std::cos(theta) + (F)15*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)8;
		case 327940:	// {5,1,4}
			return -(std::exp(i*phi)*std::sqrt((F)231/pi)*std::pow(std::cos(theta/(F)2),(F)3)*((F)1 + (F)5*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)2;
		case 327941:	// {5,1,5}
			return -(std::exp(i*phi)*std::sqrt((F)1155/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)4)*std::pow(std::sin(theta/(F)2),(F)6));
		case 328192:	// {5,2,0}
			return (std::exp((F)2*i*phi)*std::sqrt((F)1155/((F)2*pi))*((F)5*std::cos(theta) + (F)3*std::cos((F)3*theta))*std::pow(std::sin(theta),(F)2))/(F)32;
		case 328193:	// {5,2,1}
			return (std::exp((F)2*i*phi)*std::sqrt((F)77/pi)*std::cos(theta/(F)2)*((F)14 + (F)33*std::cos(theta) + (F)18*std::cos((F)2*theta) + (F)15*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)3))/(F)16;
		case 328194:	// {5,2,2}
			return (std::exp((F)2*i*phi)*std::sqrt((F)11/pi)*((F)32 + (F)57*std::cos(theta) + (F)36*std::cos((F)2*theta) + (F)15*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)8;
		case 328195:	// {5,2,3}
			return (std::exp((F)2*i*phi)*std::sqrt((F)33/((F)2*pi))*std::cos(theta/(F)2)*((F)17 + (F)24*std::cos(theta) + (F)15*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)4;
		case 328196:	// {5,2,4}
			return std::exp((F)2*i*phi)*std::sqrt((F)33/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)2 + (F)5*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)6);
		case 328197:	// {5,2,5}
			return std::exp((F)2*i*phi)*std::sqrt((F)330/pi)*std::pow(std::cos(theta/(F)2),(F)3)*std::pow(std::sin(theta/(F)2),(F)7);
		case 328448:	// {5,3,0}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)385/pi)*((F)7 + (F)9*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)3))/(F)64;
		case 328449:	// {5,3,1}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)231/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)13 + (F)12*std::cos(theta) + (F)15*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)8;
		case 328450:	// {5,3,2}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)33/((F)2*pi))*std::cos(theta/(F)2)*((F)17 + (F)24*std::cos(theta) + (F)15*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)4;
		case 328451:	// {5,3,3}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)11/pi)*((F)71 + (F)108*std::cos(theta) + (F)45*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)6))/(F)16;
		case 328452:	// {5,3,4}
			return (-(F)3*std::exp((F)3*i*phi)*std::sqrt((F)11/((F)2*pi))*std::cos(theta/(F)2)*((F)3 + (F)5*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)2;
		case 328453:	// {5,3,5}
			return (-(F)3*std::exp((F)3*i*phi)*std::sqrt((F)55/pi)*std::pow(std::sin(theta/(F)2),(F)6)*std::pow(std::sin(theta),(F)2))/(F)8;
		case 328704:	// {5,4,0}
			return ((F)3*std::exp((F)4*i*phi)*std::sqrt((F)385/((F)2*pi))*std::cos(theta)*std::pow(std::sin(theta),(F)4))/(F)16;
		case 328705:	// {5,4,1}
			return (std::exp((F)4*i*phi)*std::sqrt((F)231/pi)*std::pow(std::cos(theta/(F)2),(F)3)*((F)1 + (F)5*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)2;
		case 328706:	// {5,4,2}
			return std::exp((F)4*i*phi)*std::sqrt((F)33/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)2 + (F)5*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)6);
		case 328707:	// {5,4,3}
			return ((F)3*std::exp((F)4*i*phi)*std::sqrt((F)11/((F)2*pi))*std::cos(theta/(F)2)*((F)3 + (F)5*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)2;
		case 328708:	// {5,4,4}
			return (std::exp((F)4*i*phi)*std::sqrt((F)11/pi)*((F)4 + (F)5*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)8))/(F)2;
		case 328709:	// {5,4,5}
			return (std::exp((F)4*i*phi)*std::sqrt((F)55/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)8)*std::sin(theta))/(F)2;
		case 328960:	// {5,5,0}
			return (-(F)3*std::exp((F)5*i*phi)*std::sqrt((F)77/pi)*std::pow(std::sin(theta),(F)5))/(F)32;
		case 328961:	// {5,5,1}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)1155/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)4)*std::pow(std::sin(theta/(F)2),(F)6));
		case 328962:	// {5,5,2}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)330/pi)*std::pow(std::cos(theta/(F)2),(F)3)*std::pow(std::sin(theta/(F)2),(F)7));
		case 328963:	// {5,5,3}
			return (-(F)3*std::exp((F)5*i*phi)*std::sqrt((F)55/pi)*std::pow(std::sin(theta/(F)2),(F)6)*std::pow(std::sin(theta),(F)2))/(F)8;
		case 328964:	// {5,5,4}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)55/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)8)*std::sin(theta))/(F)2;
		case 328965:	// {5,5,5}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)11/pi)*std::pow(std::sin(theta/(F)2),(F)10))/(F)2;
		case 391680:	// {6,-6,0}
			return (std::sqrt((F)3003/pi)*std::pow(std::sin(theta),(F)6))/((F)64*std::exp((F)6*i*phi));
		case 391681:	// {6,-6,1}
			return (-(F)3*std::sqrt((F)286/pi)*std::pow(std::cos(theta/(F)2),(F)7)*std::pow(std::sin(theta/(F)2),(F)5))/std::exp((F)6*i*phi);
		case 391682:	// {6,-6,2}
			return ((F)3*std::sqrt((F)715/pi)*std::pow(-(F)1 + std::cos(theta),(F)2)*std::pow((F)1 + std::cos(theta),(F)4))/((F)128*std::exp((F)6*i*phi));
		case 391683:	// {6,-6,3}
			return -((std::sqrt((F)715/pi)*std::pow(std::cos(theta/(F)2),(F)9)*std::pow(std::sin(theta/(F)2),(F)3))/std::exp((F)6*i*phi));
		case 391684:	// {6,-6,4}
			return (std::sqrt((F)429/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)10)*std::pow(std::sin(theta/(F)2),(F)2))/std::exp((F)6*i*phi);
		case 391685:	// {6,-6,5}
			return -((std::sqrt((F)39/pi)*std::pow(std::cos(theta/(F)2),(F)11)*std::sin(theta/(F)2))/std::exp((F)6*i*phi));
		case 391686:	// {6,-6,6}
			return (std::sqrt((F)13/pi)*std::pow(std::cos(theta/(F)2),(F)12))/((F)2*std::exp((F)6*i*phi));
		case 391936:	// {6,-5,0}
			return ((F)3*std::sqrt((F)1001/pi)*std::cos(theta)*std::pow(std::sin(theta),(F)5))/((F)32*std::exp((F)5*i*phi));
		case 391937:	// {6,-5,1}
			return -((std::sqrt((F)429/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)6)*(-(F)1 + (F)6*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)4))/std::exp((F)5*i*phi));
		case 391938:	// {6,-5,2}
			return (std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)7)*(-(F)1 + (F)3*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)2*std::exp((F)5*i*phi));
		case 391939:	// {6,-5,3}
			return -(std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)8)*(-(F)1 + (F)2*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)2*std::exp((F)5*i*phi));
		case 391940:	// {6,-5,4}
			return (std::sqrt((F)143/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)9)*(-(F)2 + (F)3*std::cos(theta))*std::sin(theta/(F)2))/std::exp((F)5*i*phi);
		case 391941:	// {6,-5,5}
			return -(std::sqrt((F)13/pi)*std::pow(std::cos(theta/(F)2),(F)10)*(-(F)5 + (F)6*std::cos(theta)))/((F)2*std::exp((F)5*i*phi));
		case 391942:	// {6,-5,6}
			return -((std::sqrt((F)39/pi)*std::pow(std::cos(theta/(F)2),(F)11)*std::sin(theta/(F)2))/std::exp((F)5*i*phi));
		case 392192:	// {6,-4,0}
			return ((F)3*std::sqrt((F)91/((F)2*pi))*((F)9 + (F)11*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)4))/((F)64*std::exp((F)4*i*phi));
		case 392193:	// {6,-4,1}
			return -(std::sqrt((F)39/pi)*std::pow(std::cos(theta/(F)2),(F)5)*((F)29 - (F)22*std::cos(theta) + (F)33*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)4*std::exp((F)4*i*phi));
		case 392194:	// {6,-4,2}
			return (std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)6)*((F)35 - (F)44*std::cos(theta) + (F)33*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)8*std::exp((F)4*i*phi));
		case 392195:	// {6,-4,3}
			return -(std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)7)*((F)15 - (F)22*std::cos(theta) + (F)11*std::cos((F)2*theta))*std::sin(theta/(F)2))/((F)4*std::exp((F)4*i*phi));
		case 392196:	// {6,-4,4}
			return (std::sqrt((F)13/pi)*std::pow(std::cos(theta/(F)2),(F)8)*((F)59 - (F)88*std::cos(theta) + (F)33*std::cos((F)2*theta)))/((F)8*std::exp((F)4*i*phi));
		case 392197:	// {6,-4,5}
			return (std::sqrt((F)143/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)9)*(-(F)2 + (F)3*std::cos(theta))*std::sin(theta/(F)2))/std::exp((F)4*i*phi);
		case 392198:	// {6,-4,6}
			return (std::sqrt((F)429/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)10)*std::pow(std::sin(theta/(F)2),(F)2))/std::exp((F)4*i*phi);
		case 392448:	// {6,-3,0}
			return (std::sqrt((F)1365/pi)*((F)21*std::cos(theta) + (F)11*std::cos((F)3*theta))*std::pow(std::sin(theta),(F)3))/((F)128*std::exp((F)3*i*phi));
		case 392449:	// {6,-3,1}
			return (-(F)3*std::sqrt((F)65/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)4)*(-(F)9 + (F)25*std::cos(theta) - (F)11*std::cos((F)2*theta) + (F)11*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)8*std::exp((F)3*i*phi));
		case 392450:	// {6,-3,2}
			return ((F)3*std::sqrt((F)13/pi)*std::pow(std::cos(theta/(F)2),(F)5)*(-(F)98 + (F)185*std::cos(theta) - (F)110*std::cos((F)2*theta) + (F)55*std::cos((F)3*theta))*std::sin(theta/(F)2))/((F)32*std::exp((F)3*i*phi));
		case 392451:	// {6,-3,3}
			return -(std::sqrt((F)13/pi)*std::pow(std::cos(theta/(F)2),(F)6)*(-(F)167 + (F)285*std::cos(theta) - (F)165*std::cos((F)2*theta) + (F)55*std::cos((F)3*theta)))/((F)16*std::exp((F)3*i*phi));
		case 392452:	// {6,-3,4}
			return -(std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)7)*((F)15 - (F)22*std::cos(theta) + (F)11*std::cos((F)2*theta))*std::sin(theta/(F)2))/((F)4*std::exp((F)3*i*phi));
		case 392453:	// {6,-3,5}
			return -(std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)8)*(-(F)1 + (F)2*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)2*std::exp((F)3*i*phi));
		case 392454:	// {6,-3,6}
			return -((std::sqrt((F)715/pi)*std::pow(std::cos(theta/(F)2),(F)9)*std::pow(std::sin(theta/(F)2),(F)3))/std::exp((F)3*i*phi));
		case 392704:	// {6,-2,0}
			return (std::sqrt((F)1365/pi)*((F)35 + (F)60*std::cos((F)2*theta) + (F)33*std::cos((F)4*theta))*std::pow(std::sin(theta),(F)2))/((F)512*std::exp((F)2*i*phi));
		case 392705:	// {6,-2,1}
			return -(std::sqrt((F)65/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*((F)161 - (F)252*std::cos(theta) + (F)252*std::cos((F)2*theta) - (F)132*std::cos((F)3*theta) + (F)99*std::cos((F)4*theta))*std::sin(theta/(F)2))/((F)64*std::exp((F)2*i*phi));
		case 392706:	// {6,-2,2}
			return (std::sqrt((F)13/pi)*std::pow(std::cos(theta/(F)2),(F)4)*((F)1709 - (F)3096*std::cos(theta) + (F)2340*std::cos((F)2*theta) - (F)1320*std::cos((F)3*theta) + (F)495*std::cos((F)4*theta)))/((F)256*std::exp((F)2*i*phi));
		case 392707:	// {6,-2,3}
			return ((F)3*std::sqrt((F)13/pi)*std::pow(std::cos(theta/(F)2),(F)5)*(-(F)98 + (F)185*std::cos(theta) - (F)110*std::cos((F)2*theta) + (F)55*std::cos((F)3*theta))*std::sin(theta/(F)2))/((F)32*std::exp((F)2*i*phi));
		case 392708:	// {6,-2,4}
			return (std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)6)*((F)35 - (F)44*std::cos(theta) + (F)33*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)8*std::exp((F)2*i*phi));
		case 392709:	// {6,-2,5}
			return (std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)7)*(-(F)1 + (F)3*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)2*std::exp((F)2*i*phi));
		case 392710:	// {6,-2,6}
			return ((F)3*std::sqrt((F)715/pi)*std::pow(-(F)1 + std::cos(theta),(F)2)*std::pow((F)1 + std::cos(theta),(F)4))/((F)128*std::exp((F)2*i*phi));
		case 392960:	// {6,-1,0}
			return (std::sqrt((F)273/((F)2*pi))*((F)50*std::cos(theta) + (F)45*std::cos((F)3*theta) + (F)33*std::cos((F)5*theta))*std::sin(theta))/((F)256*std::exp(i*phi));
		case 392961:	// {6,-1,1}
			return -(std::sqrt((F)13/pi)*std::pow(std::cos(theta/(F)2),(F)2)*(-(F)175 + (F)350*std::cos(theta) - (F)300*std::cos((F)2*theta) + (F)255*std::cos((F)3*theta) - (F)165*std::cos((F)4*theta) + (F)99*std::cos((F)5*theta)))/((F)128*std::exp(i*phi));
		case 392962:	// {6,-1,2}
			return -(std::sqrt((F)65/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*((F)161 - (F)252*std::cos(theta) + (F)252*std::cos((F)2*theta) - (F)132*std::cos((F)3*theta) + (F)99*std::cos((F)4*theta))*std::sin(theta/(F)2))/((F)64*std::exp(i*phi));
		case 392963:	// {6,-1,3}
			return (-(F)3*std::sqrt((F)65/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)4)*(-(F)9 + (F)25*std::cos(theta) - (F)11*std::cos((F)2*theta) + (F)11*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)8*std::exp(i*phi));
		case 392964:	// {6,-1,4}
			return -(std::sqrt((F)39/pi)*std::pow(std::cos(theta/(F)2),(F)5)*((F)29 - (F)22*std::cos(theta) + (F)33*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)4*std::exp(i*phi));
		case 392965:	// {6,-1,5}
			return -((std::sqrt((F)429/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)6)*(-(F)1 + (F)6*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)4))/std::exp(i*phi));
		case 392966:	// {6,-1,6}
			return (-(F)3*std::sqrt((F)286/pi)*std::pow(std::cos(theta/(F)2),(F)7)*std::pow(std::sin(theta/(F)2),(F)5))/std::exp(i*phi);
		case 393216:	// {6,0,0}
			return (std::sqrt((F)13/pi)*((F)50 + (F)105*std::cos((F)2*theta) + (F)126*std::cos((F)4*theta) + (F)231*std::cos((F)6*theta)))/(F)1024;
		case 393217:	// {6,0,1}
			return (std::sqrt((F)273/((F)2*pi))*((F)50*std::cos(theta) + (F)45*std::cos((F)3*theta) + (F)33*std::cos((F)5*theta))*std::sin(theta))/(F)256;
		case 393218:	// {6,0,2}
			return (std::sqrt((F)1365/pi)*((F)35 + (F)60*std::cos((F)2*theta) + (F)33*std::cos((F)4*theta))*std::pow(std::sin(theta),(F)2))/(F)512;
		case 393219:	// {6,0,3}
			return (std::sqrt((F)1365/pi)*((F)21*std::cos(theta) + (F)11*std::cos((F)3*theta))*std::pow(std::sin(theta),(F)3))/(F)128;
		case 393220:	// {6,0,4}
			return ((F)3*std::sqrt((F)91/((F)2*pi))*((F)9 + (F)11*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)4))/(F)64;
		case 393221:	// {6,0,5}
			return ((F)3*std::sqrt((F)1001/pi)*std::cos(theta)*std::pow(std::sin(theta),(F)5))/(F)32;
		case 393222:	// {6,0,6}
			return (std::sqrt((F)3003/pi)*std::pow(std::sin(theta),(F)6))/(F)64;
		case 393472:	// {6,1,0}
			return -(std::exp(i*phi)*std::sqrt((F)273/((F)2*pi))*((F)50*std::cos(theta) + (F)45*std::cos((F)3*theta) + (F)33*std::cos((F)5*theta))*std::sin(theta))/(F)256;
		case 393473:	// {6,1,1}
			return -(std::exp(i*phi)*std::sqrt((F)13/pi)*((F)175 + (F)350*std::cos(theta) + (F)300*std::cos((F)2*theta) + (F)255*std::cos((F)3*theta) + (F)165*std::cos((F)4*theta) + (F)99*std::cos((F)5*theta))*std::pow(std::sin(theta/(F)2),(F)2))/(F)128;
		case 393474:	// {6,1,2}
			return -(std::exp(i*phi)*std::sqrt((F)65/((F)2*pi))*std::cos(theta/(F)2)*((F)161 + (F)252*std::cos(theta) + (F)252*std::cos((F)2*theta) + (F)132*std::cos((F)3*theta) + (F)99*std::cos((F)4*theta))*std::pow(std::sin(theta/(F)2),(F)3))/(F)64;
		case 393475:	// {6,1,3}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)65/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)9 + (F)25*std::cos(theta) + (F)11*std::cos((F)2*theta) + (F)11*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)8;
		case 393476:	// {6,1,4}
			return -(std::exp(i*phi)*std::sqrt((F)39/pi)*std::pow(std::cos(theta/(F)2),(F)3)*((F)29 + (F)22*std::cos(theta) + (F)33*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)4;
		case 393477:	// {6,1,5}
			return -(std::exp(i*phi)*std::sqrt((F)429/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)4)*((F)1 + (F)6*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)6));
		case 393478:	// {6,1,6}
			return -(F)3*std::exp(i*phi)*std::sqrt((F)286/pi)*std::pow(std::cos(theta/(F)2),(F)5)*std::pow(std::sin(theta/(F)2),(F)7);
		case 393728:	// {6,2,0}
			return (std::exp((F)2*i*phi)*std::sqrt((F)1365/pi)*((F)35 + (F)60*std::cos((F)2*theta) + (F)33*std::cos((F)4*theta))*std::pow(std::sin(theta),(F)2))/(F)512;
		case 393729:	// {6,2,1}
			return (std::exp((F)2*i*phi)*std::sqrt((F)65/((F)2*pi))*std::cos(theta/(F)2)*((F)161 + (F)252*std::cos(theta) + (F)252*std::cos((F)2*theta) + (F)132*std::cos((F)3*theta) + (F)99*std::cos((F)4*theta))*std::pow(std::sin(theta/(F)2),(F)3))/(F)64;
		case 393730:	// {6,2,2}
			return (std::exp((F)2*i*phi)*std::sqrt((F)13/pi)*((F)1709 + (F)3096*std::cos(theta) + (F)2340*std::cos((F)2*theta) + (F)1320*std::cos((F)3*theta) + (F)495*std::cos((F)4*theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)256;
		case 393731:	// {6,2,3}
			return ((F)3*std::exp((F)2*i*phi)*std::sqrt((F)13/pi)*std::cos(theta/(F)2)*((F)98 + (F)185*std::cos(theta) + (F)110*std::cos((F)2*theta) + (F)55*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)32;
		case 393732:	// {6,2,4}
			return (std::exp((F)2*i*phi)*std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)35 + (F)44*std::cos(theta) + (F)33*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)6))/(F)8;
		case 393733:	// {6,2,5}
			return (std::exp((F)2*i*phi)*std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)3)*((F)1 + (F)3*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)2;
		case 393734:	// {6,2,6}
			return ((F)3*std::exp((F)2*i*phi)*std::sqrt((F)715/pi)*std::pow(std::sin(theta/(F)2),(F)4)*std::pow(std::sin(theta),(F)4))/(F)32;
		case 393984:	// {6,3,0}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)1365/pi)*((F)21*std::cos(theta) + (F)11*std::cos((F)3*theta))*std::pow(std::sin(theta),(F)3))/(F)128;
		case 393985:	// {6,3,1}
			return (-(F)3*std::exp((F)3*i*phi)*std::sqrt((F)65/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)9 + (F)25*std::cos(theta) + (F)11*std::cos((F)2*theta) + (F)11*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)8;
		case 393986:	// {6,3,2}
			return (-(F)3*std::exp((F)3*i*phi)*std::sqrt((F)13/pi)*std::cos(theta/(F)2)*((F)98 + (F)185*std::cos(theta) + (F)110*std::cos((F)2*theta) + (F)55*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)32;
		case 393987:	// {6,3,3}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)13/pi)*((F)167 + (F)285*std::cos(theta) + (F)165*std::cos((F)2*theta) + (F)55*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)6))/(F)16;
		case 393988:	// {6,3,4}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)195/((F)2*pi))*std::cos(theta/(F)2)*((F)15 + (F)22*std::cos(theta) + (F)11*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)4;
		case 393989:	// {6,3,5}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)2145/pi)*std::pow(std::sin(theta/(F)2),(F)6)*std::sin(theta)*(std::sin(theta) + std::sin((F)2*theta)))/(F)8;
		case 393990:	// {6,3,6}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)715/pi)*std::pow(std::sin(theta/(F)2),(F)6)*std::pow(std::sin(theta),(F)3))/(F)8;
		case 394240:	// {6,4,0}
			return ((F)3*std::exp((F)4*i*phi)*std::sqrt((F)91/((F)2*pi))*((F)9 + (F)11*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)4))/(F)64;
		case 394241:	// {6,4,1}
			return (std::exp((F)4*i*phi)*std::sqrt((F)39/pi)*std::pow(std::cos(theta/(F)2),(F)3)*((F)29 + (F)22*std::cos(theta) + (F)33*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)4;
		case 394242:	// {6,4,2}
			return (std::exp((F)4*i*phi)*std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)35 + (F)44*std::cos(theta) + (F)33*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)6))/(F)8;
		case 394243:	// {6,4,3}
			return (std::exp((F)4*i*phi)*std::sqrt((F)195/((F)2*pi))*std::cos(theta/(F)2)*((F)15 + (F)22*std::cos(theta) + (F)11*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)4;
		case 394244:	// {6,4,4}
			return (std::exp((F)4*i*phi)*std::sqrt((F)13/pi)*((F)59 + (F)88*std::cos(theta) + (F)33*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)8))/(F)8;
		case 394245:	// {6,4,5}
			return std::exp((F)4*i*phi)*std::sqrt((F)143/((F)2*pi))*std::cos(theta/(F)2)*((F)2 + (F)3*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)9);
		case 394246:	// {6,4,6}
			return (std::exp((F)4*i*phi)*std::sqrt((F)429/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)8)*std::pow(std::sin(theta),(F)2))/(F)4;
		case 394496:	// {6,5,0}
			return (-(F)3*std::exp((F)5*i*phi)*std::sqrt((F)1001/pi)*std::cos(theta)*std::pow(std::sin(theta),(F)5))/(F)32;
		case 394497:	// {6,5,1}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)429/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)4)*((F)1 + (F)6*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)6));
		case 394498:	// {6,5,2}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)3)*((F)1 + (F)3*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)2;
		case 394499:	// {6,5,3}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)2145/pi)*std::pow(std::sin(theta/(F)2),(F)6)*std::sin(theta)*(std::sin(theta) + std::sin((F)2*theta)))/(F)8;
		case 394500:	// {6,5,4}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)143/((F)2*pi))*std::cos(theta/(F)2)*((F)2 + (F)3*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)9));
		case 394501:	// {6,5,5}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)13/pi)*((F)5 + (F)6*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)10))/(F)2;
		case 394502:	// {6,5,6}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)39/pi)*std::pow(std::sin(theta/(F)2),(F)10)*std::sin(theta))/(F)2;
		case 394752:	// {6,6,0}
			return (std::exp((F)6*i*phi)*std::sqrt((F)3003/pi)*std::pow(std::sin(theta),(F)6))/(F)64;
		case 394753:	// {6,6,1}
			return (F)3*std::exp((F)6*i*phi)*std::sqrt((F)286/pi)*std::pow(std::cos(theta/(F)2),(F)5)*std::pow(std::sin(theta/(F)2),(F)7);
		case 394754:	// {6,6,2}
			return ((F)3*std::exp((F)6*i*phi)*std::sqrt((F)715/pi)*std::pow(std::sin(theta/(F)2),(F)4)*std::pow(std::sin(theta),(F)4))/(F)32;
		case 394755:	// {6,6,3}
			return (std::exp((F)6*i*phi)*std::sqrt((F)715/pi)*std::pow(std::sin(theta/(F)2),(F)6)*std::pow(std::sin(theta),(F)3))/(F)8;
		case 394756:	// {6,6,4}
			return (std::exp((F)6*i*phi)*std::sqrt((F)429/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)8)*std::pow(std::sin(theta),(F)2))/(F)4;
		case 394757:	// {6,6,5}
			return (std::exp((F)6*i*phi)*std::sqrt((F)39/pi)*std::pow(std::sin(theta/(F)2),(F)10)*std::sin(theta))/(F)2;
		case 394758:	// {6,6,6}
			return (std::exp((F)6*i*phi)*std::sqrt((F)13/pi)*std::pow(std::sin(theta/(F)2),(F)12))/(F)2;
		case 456960:	// {7,-7,0}
			return ((F)3*std::sqrt((F)715/((F)2*pi))*std::pow(std::sin(theta),(F)7))/((F)64*std::exp((F)7*i*phi));
		case 456961:	// {7,-7,1}
			return ((F)3*std::sqrt((F)5005/pi)*std::pow(-(F)1 + std::cos(theta),(F)3)*std::pow((F)1 + std::cos(theta),(F)4))/((F)256*std::exp((F)7*i*phi));
		case 456962:	// {7,-7,2}
			return (std::sqrt((F)15015/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)9)*std::pow(std::sin(theta/(F)2),(F)5))/std::exp((F)7*i*phi);
		case 456963:	// {7,-7,3}
			return -(std::sqrt((F)15015/pi)*std::pow(-(F)1 + std::cos(theta),(F)2)*std::pow((F)1 + std::cos(theta),(F)5))/((F)256*std::exp((F)7*i*phi));
		case 456964:	// {7,-7,4}
			return (std::sqrt((F)1365/pi)*std::pow(std::cos(theta/(F)2),(F)11)*std::pow(std::sin(theta/(F)2),(F)3))/std::exp((F)7*i*phi);
		case 456965:	// {7,-7,5}
			return -(std::sqrt((F)1365/pi)*std::pow(std::cos(theta/(F)2),(F)12)*std::pow(std::sin(theta/(F)2),(F)2))/((F)2*std::exp((F)7*i*phi));
		case 456966:	// {7,-7,6}
			return (std::sqrt((F)105/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)13)*std::sin(theta/(F)2))/std::exp((F)7*i*phi);
		case 456967:	// {7,-7,7}
			return -(std::sqrt((F)15/pi)*std::pow(std::cos(theta/(F)2),(F)14))/((F)2*std::exp((F)7*i*phi));
		case 457216:	// {7,-6,0}
			return ((F)3*std::sqrt((F)5005/pi)*std::cos(theta)*std::pow(std::sin(theta),(F)6))/((F)64*std::exp((F)6*i*phi));
		case 457217:	// {7,-6,1}
			return (-(F)3*std::sqrt((F)715/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)7)*(-(F)1 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)5))/((F)2*std::exp((F)6*i*phi));
		case 457218:	// {7,-6,2}
			return (std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)8)*(-(F)2 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)4))/((F)2*std::exp((F)6*i*phi));
		case 457219:	// {7,-6,3}
			return -(std::sqrt((F)2145/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)9)*(-(F)3 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)2*std::exp((F)6*i*phi));
		case 457220:	// {7,-6,4}
			return (std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)10)*(-(F)4 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)2))/std::exp((F)6*i*phi);
		case 457221:	// {7,-6,5}
			return -(std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)11)*(-(F)5 + (F)7*std::cos(theta))*std::sin(theta/(F)2))/((F)2*std::exp((F)6*i*phi));
		case 457222:	// {7,-6,6}
			return (std::sqrt((F)15/pi)*std::pow(std::cos(theta/(F)2),(F)12)*(-(F)6 + (F)7*std::cos(theta)))/((F)2*std::exp((F)6*i*phi));
		case 457223:	// {7,-6,7}
			return (std::sqrt((F)105/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)13)*std::sin(theta/(F)2))/std::exp((F)6*i*phi);
		case 457472:	// {7,-5,0}
			return ((F)3*std::sqrt((F)385/((F)2*pi))*((F)11 + (F)13*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)5))/((F)128*std::exp((F)5*i*phi));
		case 457473:	// {7,-5,1}
			return (-(F)3*std::sqrt((F)55/pi)*std::pow(std::cos(theta/(F)2),(F)6)*((F)81 - (F)52*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)4))/((F)16*std::exp((F)5*i*phi));
		case 457474:	// {7,-5,2}
			return (std::sqrt((F)165/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)7)*((F)93 - (F)104*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)8*std::exp((F)5*i*phi));
		case 457475:	// {7,-5,3}
			return -(std::sqrt((F)165/pi)*std::pow(std::cos(theta/(F)2),(F)8)*((F)113 - (F)156*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)16*std::exp((F)5*i*phi));
		case 457476:	// {7,-5,4}
			return (std::sqrt((F)15/pi)*std::pow(std::cos(theta/(F)2),(F)9)*((F)141 - (F)208*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::sin(theta/(F)2))/((F)8*std::exp((F)5*i*phi));
		case 457477:	// {7,-5,5}
			return -(std::sqrt((F)15/pi)*std::pow(std::cos(theta/(F)2),(F)10)*((F)177 - (F)260*std::cos(theta) + (F)91*std::cos((F)2*theta)))/((F)16*std::exp((F)5*i*phi));
		case 457478:	// {7,-5,6}
			return -(std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)11)*(-(F)5 + (F)7*std::cos(theta))*std::sin(theta/(F)2))/((F)2*std::exp((F)5*i*phi));
		case 457479:	// {7,-5,7}
			return -(std::sqrt((F)1365/pi)*std::pow(std::cos(theta/(F)2),(F)12)*std::pow(std::sin(theta/(F)2),(F)2))/((F)2*std::exp((F)5*i*phi));
		case 457728:	// {7,-4,0}
			return ((F)3*std::sqrt((F)385/((F)2*pi))*((F)27*std::cos(theta) + (F)13*std::cos((F)3*theta))*std::pow(std::sin(theta),(F)4))/((F)128*std::exp((F)4*i*phi));
		case 457729:	// {7,-4,1}
			return (-(F)3*std::sqrt((F)55/pi)*std::pow(std::cos(theta/(F)2),(F)5)*(-(F)66 + (F)213*std::cos(theta) - (F)78*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)32*std::exp((F)4*i*phi));
		case 457730:	// {7,-4,2}
			return (std::sqrt((F)165/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)6)*(-(F)140 + (F)285*std::cos(theta) - (F)156*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)16*std::exp((F)4*i*phi));
		case 457731:	// {7,-4,3}
			return -(std::sqrt((F)165/pi)*std::pow(std::cos(theta/(F)2),(F)7)*(-(F)230 + (F)405*std::cos(theta) - (F)234*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::sin(theta/(F)2))/((F)32*std::exp((F)4*i*phi));
		case 457732:	// {7,-4,4}
			return (std::sqrt((F)15/pi)*std::pow(std::cos(theta/(F)2),(F)8)*(-(F)344 + (F)573*std::cos(theta) - (F)312*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta)))/((F)16*std::exp((F)4*i*phi));
		case 457733:	// {7,-4,5}
			return (std::sqrt((F)15/pi)*std::pow(std::cos(theta/(F)2),(F)9)*((F)141 - (F)208*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::sin(theta/(F)2))/((F)8*std::exp((F)4*i*phi));
		case 457734:	// {7,-4,6}
			return (std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)10)*(-(F)4 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)2))/std::exp((F)4*i*phi);
		case 457735:	// {7,-4,7}
			return (std::sqrt((F)1365/pi)*std::pow(std::cos(theta/(F)2),(F)11)*std::pow(std::sin(theta/(F)2),(F)3))/std::exp((F)4*i*phi);
		case 457984:	// {7,-3,0}
			return ((F)3*std::sqrt((F)35/((F)2*pi))*((F)189 + (F)308*std::cos((F)2*theta) + (F)143*std::cos((F)4*theta))*std::pow(std::sin(theta),(F)3))/((F)512*std::exp((F)3*i*phi));
		case 457985:	// {7,-3,1}
			return (-(F)3*std::sqrt((F)5/pi)*std::pow(std::cos(theta/(F)2),(F)4)*(-(F)2376*std::cos(theta) + (F)2684*std::cos((F)2*theta) + (F)13*((F)135 - (F)88*std::cos((F)3*theta) + (F)77*std::cos((F)4*theta)))*std::pow(std::sin(theta/(F)2),(F)2))/((F)256*std::exp((F)3*i*phi));
		case 457986:	// {7,-3,2}
			return (std::sqrt((F)15/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*((F)3115 - (F)5456*std::cos(theta) + (F)4268*std::cos((F)2*theta) - (F)2288*std::cos((F)3*theta) + (F)1001*std::cos((F)4*theta))*std::sin(theta/(F)2))/((F)128*std::exp((F)3*i*phi));
		case 457987:	// {7,-3,3}
			return -(std::sqrt((F)15/pi)*std::pow(std::cos(theta/(F)2),(F)6)*((F)5595 - (F)9944*std::cos(theta) + (F)6908*std::cos((F)2*theta) - (F)3432*std::cos((F)3*theta) + (F)1001*std::cos((F)4*theta)))/((F)256*std::exp((F)3*i*phi));
		case 457988:	// {7,-3,4}
			return -(std::sqrt((F)165/pi)*std::pow(std::cos(theta/(F)2),(F)7)*(-(F)230 + (F)405*std::cos(theta) - (F)234*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::sin(theta/(F)2))/((F)32*std::exp((F)3*i*phi));
		case 457989:	// {7,-3,5}
			return -(std::sqrt((F)165/pi)*std::pow(std::cos(theta/(F)2),(F)8)*((F)113 - (F)156*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)16*std::exp((F)3*i*phi));
		case 457990:	// {7,-3,6}
			return -(std::sqrt((F)2145/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)9)*(-(F)3 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)2*std::exp((F)3*i*phi));
		case 457991:	// {7,-3,7}
			return -(std::sqrt((F)15015/pi)*std::pow(-(F)1 + std::cos(theta),(F)2)*std::pow((F)1 + std::cos(theta),(F)5))/((F)256*std::exp((F)3*i*phi));
		case 458240:	// {7,-2,0}
			return ((F)3*std::sqrt((F)35/pi)*((F)350*std::cos(theta) + (F)275*std::cos((F)3*theta) + (F)143*std::cos((F)5*theta))*std::pow(std::sin(theta),(F)2))/((F)1024*std::exp((F)2*i*phi));
		case 458241:	// {7,-2,1}
			return (-(F)3*std::sqrt((F)5/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*(-(F)1890 + (F)4130*std::cos(theta) - (F)3080*std::cos((F)2*theta) + (F)2805*std::cos((F)3*theta) - (F)1430*std::cos((F)4*theta) + (F)1001*std::cos((F)5*theta))*std::sin(theta/(F)2))/((F)512*std::exp((F)2*i*phi));
		case 458242:	// {7,-2,2}
			return (std::sqrt((F)15/pi)*std::pow(std::cos(theta/(F)2),(F)4)*(-(F)5220 + (F)9810*std::cos(theta) - (F)7920*std::cos((F)2*theta) + (F)5445*std::cos((F)3*theta) - (F)2860*std::cos((F)4*theta) + (F)1001*std::cos((F)5*theta)))/((F)512*std::exp((F)2*i*phi));
		case 458243:	// {7,-2,3}
			return (std::sqrt((F)15/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*((F)3115 - (F)5456*std::cos(theta) + (F)4268*std::cos((F)2*theta) - (F)2288*std::cos((F)3*theta) + (F)1001*std::cos((F)4*theta))*std::sin(theta/(F)2))/((F)128*std::exp((F)2*i*phi));
		case 458244:	// {7,-2,4}
			return (std::sqrt((F)165/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)6)*(-(F)140 + (F)285*std::cos(theta) - (F)156*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)2))/((F)16*std::exp((F)2*i*phi));
		case 458245:	// {7,-2,5}
			return (std::sqrt((F)165/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)7)*((F)93 - (F)104*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)8*std::exp((F)2*i*phi));
		case 458246:	// {7,-2,6}
			return (std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)8)*(-(F)2 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)4))/((F)2*std::exp((F)2*i*phi));
		case 458247:	// {7,-2,7}
			return (std::sqrt((F)15015/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)9)*std::pow(std::sin(theta/(F)2),(F)5))/std::exp((F)2*i*phi);
		case 458496:	// {7,-1,0}
			return (std::sqrt((F)105/((F)2*pi))*((F)350 + (F)675*std::cos((F)2*theta) + (F)594*std::cos((F)4*theta) + (F)429*std::cos((F)6*theta))*std::sin(theta))/((F)2048*std::exp(i*phi));
		case 458497:	// {7,-1,1}
			return -(std::sqrt((F)15/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)6650 - (F)12600*std::cos(theta) + (F)11925*std::cos((F)2*theta) - (F)9900*std::cos((F)3*theta) + (F)8118*std::cos((F)4*theta) - (F)5148*std::cos((F)5*theta) + (F)3003*std::cos((F)6*theta)))/((F)4096*std::exp(i*phi));
		case 458498:	// {7,-1,2}
			return (-(F)3*std::sqrt((F)5/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*(-(F)1890 + (F)4130*std::cos(theta) - (F)3080*std::cos((F)2*theta) + (F)2805*std::cos((F)3*theta) - (F)1430*std::cos((F)4*theta) + (F)1001*std::cos((F)5*theta))*std::sin(theta/(F)2))/((F)512*std::exp(i*phi));
		case 458499:	// {7,-1,3}
			return (-(F)3*std::sqrt((F)5/pi)*std::pow(std::cos(theta/(F)2),(F)4)*(-(F)2376*std::cos(theta) + (F)2684*std::cos((F)2*theta) + (F)13*((F)135 - (F)88*std::cos((F)3*theta) + (F)77*std::cos((F)4*theta)))*std::pow(std::sin(theta/(F)2),(F)2))/((F)256*std::exp(i*phi));
		case 458500:	// {7,-1,4}
			return (-(F)3*std::sqrt((F)55/pi)*std::pow(std::cos(theta/(F)2),(F)5)*(-(F)66 + (F)213*std::cos(theta) - (F)78*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)3))/((F)32*std::exp(i*phi));
		case 458501:	// {7,-1,5}
			return (-(F)3*std::sqrt((F)55/pi)*std::pow(std::cos(theta/(F)2),(F)6)*((F)81 - (F)52*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)4))/((F)16*std::exp(i*phi));
		case 458502:	// {7,-1,6}
			return (-(F)3*std::sqrt((F)715/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)7)*(-(F)1 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)5))/((F)2*std::exp(i*phi));
		case 458503:	// {7,-1,7}
			return ((F)3*std::sqrt((F)5005/pi)*std::pow(-(F)1 + std::cos(theta),(F)3)*std::pow((F)1 + std::cos(theta),(F)4))/((F)256*std::exp(i*phi));
		case 458752:	// {7,0,0}
			return (std::sqrt((F)15/pi)*((F)175*std::cos(theta) + (F)189*std::cos((F)3*theta) + (F)231*std::cos((F)5*theta) + (F)429*std::cos((F)7*theta)))/(F)2048;
		case 458753:	// {7,0,1}
			return (std::sqrt((F)105/((F)2*pi))*((F)350 + (F)675*std::cos((F)2*theta) + (F)594*std::cos((F)4*theta) + (F)429*std::cos((F)6*theta))*std::sin(theta))/(F)2048;
		case 458754:	// {7,0,2}
			return ((F)3*std::sqrt((F)35/pi)*((F)350*std::cos(theta) + (F)275*std::cos((F)3*theta) + (F)143*std::cos((F)5*theta))*std::pow(std::sin(theta),(F)2))/(F)1024;
		case 458755:	// {7,0,3}
			return ((F)3*std::sqrt((F)35/((F)2*pi))*((F)189 + (F)308*std::cos((F)2*theta) + (F)143*std::cos((F)4*theta))*std::pow(std::sin(theta),(F)3))/(F)512;
		case 458756:	// {7,0,4}
			return ((F)3*std::sqrt((F)385/((F)2*pi))*((F)27*std::cos(theta) + (F)13*std::cos((F)3*theta))*std::pow(std::sin(theta),(F)4))/(F)128;
		case 458757:	// {7,0,5}
			return ((F)3*std::sqrt((F)385/((F)2*pi))*((F)11 + (F)13*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)5))/(F)128;
		case 458758:	// {7,0,6}
			return ((F)3*std::sqrt((F)5005/pi)*std::cos(theta)*std::pow(std::sin(theta),(F)6))/(F)64;
		case 458759:	// {7,0,7}
			return ((F)3*std::sqrt((F)715/((F)2*pi))*std::pow(std::sin(theta),(F)7))/(F)64;
		case 459008:	// {7,1,0}
			return -(std::exp(i*phi)*std::sqrt((F)105/((F)2*pi))*((F)350 + (F)675*std::cos((F)2*theta) + (F)594*std::cos((F)4*theta) + (F)429*std::cos((F)6*theta))*std::sin(theta))/(F)2048;
		case 459009:	// {7,1,1}
			return -(std::exp(i*phi)*std::sqrt((F)15/pi)*((F)6650 + (F)12600*std::cos(theta) + (F)11925*std::cos((F)2*theta) + (F)9900*std::cos((F)3*theta) + (F)8118*std::cos((F)4*theta) + (F)5148*std::cos((F)5*theta) + (F)3003*std::cos((F)6*theta))*std::pow(std::sin(theta/(F)2),(F)2))/(F)4096;
		case 459010:	// {7,1,2}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)5/((F)2*pi))*std::cos(theta/(F)2)*((F)1890 + (F)4130*std::cos(theta) + (F)3080*std::cos((F)2*theta) + (F)2805*std::cos((F)3*theta) + (F)1430*std::cos((F)4*theta) + (F)1001*std::cos((F)5*theta))*std::pow(std::sin(theta/(F)2),(F)3))/(F)512;
		case 459011:	// {7,1,3}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)5/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)2376*std::cos(theta) + (F)2684*std::cos((F)2*theta) + (F)13*((F)135 + (F)88*std::cos((F)3*theta) + (F)77*std::cos((F)4*theta)))*std::pow(std::sin(theta/(F)2),(F)4))/(F)256;
		case 459012:	// {7,1,4}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)55/pi)*std::pow(std::cos(theta/(F)2),(F)3)*((F)66 + (F)213*std::cos(theta) + (F)78*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)32;
		case 459013:	// {7,1,5}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)55/pi)*std::pow(std::cos(theta/(F)2),(F)4)*((F)81 + (F)52*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)6))/(F)16;
		case 459014:	// {7,1,6}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)715/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*((F)1 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)2;
		case 459015:	// {7,1,7}
			return (-(F)3*std::exp(i*phi)*std::sqrt((F)5005/pi)*std::pow(std::sin(theta/(F)2),(F)2)*std::pow(std::sin(theta),(F)6))/(F)128;
		case 459264:	// {7,2,0}
			return ((F)3*std::exp((F)2*i*phi)*std::sqrt((F)35/pi)*((F)350*std::cos(theta) + (F)275*std::cos((F)3*theta) + (F)143*std::cos((F)5*theta))*std::pow(std::sin(theta),(F)2))/(F)1024;
		case 459265:	// {7,2,1}
			return ((F)3*std::exp((F)2*i*phi)*std::sqrt((F)5/((F)2*pi))*std::cos(theta/(F)2)*((F)1890 + (F)4130*std::cos(theta) + (F)3080*std::cos((F)2*theta) + (F)2805*std::cos((F)3*theta) + (F)1430*std::cos((F)4*theta) + (F)1001*std::cos((F)5*theta))*std::pow(std::sin(theta/(F)2),(F)3))/(F)512;
		case 459266:	// {7,2,2}
			return (std::exp((F)2*i*phi)*std::sqrt((F)15/pi)*((F)5220 + (F)9810*std::cos(theta) + (F)7920*std::cos((F)2*theta) + (F)5445*std::cos((F)3*theta) + (F)2860*std::cos((F)4*theta) + (F)1001*std::cos((F)5*theta))*std::pow(std::sin(theta/(F)2),(F)4))/(F)512;
		case 459267:	// {7,2,3}
			return (std::exp((F)2*i*phi)*std::sqrt((F)15/((F)2*pi))*std::cos(theta/(F)2)*((F)3115 + (F)5456*std::cos(theta) + (F)4268*std::cos((F)2*theta) + (F)2288*std::cos((F)3*theta) + (F)1001*std::cos((F)4*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)128;
		case 459268:	// {7,2,4}
			return (std::exp((F)2*i*phi)*std::sqrt((F)165/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)140 + (F)285*std::cos(theta) + (F)156*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)6))/(F)16;
		case 459269:	// {7,2,5}
			return (std::exp((F)2*i*phi)*std::sqrt((F)165/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*((F)93 + (F)104*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)8;
		case 459270:	// {7,2,6}
			return (std::exp((F)2*i*phi)*std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)4)*((F)2 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)8))/(F)2;
		case 459271:	// {7,2,7}
			return std::exp((F)2*i*phi)*std::sqrt((F)15015/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*std::pow(std::sin(theta/(F)2),(F)9);
		case 459520:	// {7,3,0}
			return (-(F)3*std::exp((F)3*i*phi)*std::sqrt((F)35/((F)2*pi))*((F)189 + (F)308*std::cos((F)2*theta) + (F)143*std::cos((F)4*theta))*std::pow(std::sin(theta),(F)3))/(F)512;
		case 459521:	// {7,3,1}
			return (-(F)3*std::exp((F)3*i*phi)*std::sqrt((F)5/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)2376*std::cos(theta) + (F)2684*std::cos((F)2*theta) + (F)13*((F)135 + (F)88*std::cos((F)3*theta) + (F)77*std::cos((F)4*theta)))*std::pow(std::sin(theta/(F)2),(F)4))/(F)256;
		case 459522:	// {7,3,2}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)15/((F)2*pi))*std::cos(theta/(F)2)*((F)3115 + (F)5456*std::cos(theta) + (F)4268*std::cos((F)2*theta) + (F)2288*std::cos((F)3*theta) + (F)1001*std::cos((F)4*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)128;
		case 459523:	// {7,3,3}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)15/pi)*((F)5595 + (F)9944*std::cos(theta) + (F)6908*std::cos((F)2*theta) + (F)3432*std::cos((F)3*theta) + (F)1001*std::cos((F)4*theta))*std::pow(std::sin(theta/(F)2),(F)6))/(F)256;
		case 459524:	// {7,3,4}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)165/pi)*std::cos(theta/(F)2)*((F)230 + (F)405*std::cos(theta) + (F)234*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)32;
		case 459525:	// {7,3,5}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)165/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)113 + (F)156*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)8))/(F)16;
		case 459526:	// {7,3,6}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)2145/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*((F)3 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)9))/(F)2;
		case 459527:	// {7,3,7}
			return -(std::exp((F)3*i*phi)*std::sqrt((F)15015/pi)*std::pow(std::sin(theta/(F)2),(F)6)*std::pow(std::sin(theta),(F)4))/(F)32;
		case 459776:	// {7,4,0}
			return ((F)3*std::exp((F)4*i*phi)*std::sqrt((F)385/((F)2*pi))*((F)27*std::cos(theta) + (F)13*std::cos((F)3*theta))*std::pow(std::sin(theta),(F)4))/(F)128;
		case 459777:	// {7,4,1}
			return ((F)3*std::exp((F)4*i*phi)*std::sqrt((F)55/pi)*std::pow(std::cos(theta/(F)2),(F)3)*((F)66 + (F)213*std::cos(theta) + (F)78*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)5))/(F)32;
		case 459778:	// {7,4,2}
			return (std::exp((F)4*i*phi)*std::sqrt((F)165/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)140 + (F)285*std::cos(theta) + (F)156*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)6))/(F)16;
		case 459779:	// {7,4,3}
			return (std::exp((F)4*i*phi)*std::sqrt((F)165/pi)*std::cos(theta/(F)2)*((F)230 + (F)405*std::cos(theta) + (F)234*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)32;
		case 459780:	// {7,4,4}
			return (std::exp((F)4*i*phi)*std::sqrt((F)15/pi)*((F)344 + (F)573*std::cos(theta) + (F)312*std::cos((F)2*theta) + (F)91*std::cos((F)3*theta))*std::pow(std::sin(theta/(F)2),(F)8))/(F)16;
		case 459781:	// {7,4,5}
			return (std::exp((F)4*i*phi)*std::sqrt((F)15/pi)*std::cos(theta/(F)2)*((F)141 + (F)208*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)9))/(F)8;
		case 459782:	// {7,4,6}
			return std::exp((F)4*i*phi)*std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)4 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)10);
		case 459783:	// {7,4,7}
			return (std::exp((F)4*i*phi)*std::sqrt((F)1365/pi)*std::pow(std::sin(theta/(F)2),(F)8)*std::pow(std::sin(theta),(F)3))/(F)8;
		case 460032:	// {7,5,0}
			return (-(F)3*std::exp((F)5*i*phi)*std::sqrt((F)385/((F)2*pi))*((F)11 + (F)13*std::cos((F)2*theta))*std::pow(std::sin(theta),(F)5))/(F)128;
		case 460033:	// {7,5,1}
			return (-(F)3*std::exp((F)5*i*phi)*std::sqrt((F)55/pi)*std::pow(std::cos(theta/(F)2),(F)4)*((F)81 + (F)52*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)6))/(F)16;
		case 460034:	// {7,5,2}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)165/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*((F)93 + (F)104*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)8;
		case 460035:	// {7,5,3}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)165/pi)*std::pow(std::cos(theta/(F)2),(F)2)*((F)113 + (F)156*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)8))/(F)16;
		case 460036:	// {7,5,4}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)15/pi)*std::cos(theta/(F)2)*((F)141 + (F)208*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)9))/(F)8;
		case 460037:	// {7,5,5}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)15/pi)*((F)177 + (F)260*std::cos(theta) + (F)91*std::cos((F)2*theta))*std::pow(std::sin(theta/(F)2),(F)10))/(F)16;
		case 460038:	// {7,5,6}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)195/((F)2*pi))*std::cos(theta/(F)2)*((F)5 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)11))/(F)2;
		case 460039:	// {7,5,7}
			return -(std::exp((F)5*i*phi)*std::sqrt((F)1365/pi)*std::pow(std::sin(theta/(F)2),(F)10)*std::pow(std::sin(theta),(F)2))/(F)8;
		case 460288:	// {7,6,0}
			return ((F)3*std::exp((F)6*i*phi)*std::sqrt((F)5005/pi)*std::cos(theta)*std::pow(std::sin(theta),(F)6))/(F)64;
		case 460289:	// {7,6,1}
			return ((F)3*std::exp((F)6*i*phi)*std::sqrt((F)715/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*((F)1 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)7))/(F)2;
		case 460290:	// {7,6,2}
			return (std::exp((F)6*i*phi)*std::sqrt((F)2145/pi)*std::pow(std::cos(theta/(F)2),(F)4)*((F)2 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)8))/(F)2;
		case 460291:	// {7,6,3}
			return (std::exp((F)6*i*phi)*std::sqrt((F)2145/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)3)*((F)3 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)9))/(F)2;
		case 460292:	// {7,6,4}
			return std::exp((F)6*i*phi)*std::sqrt((F)195/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)2)*((F)4 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)10);
		case 460293:	// {7,6,5}
			return (std::exp((F)6*i*phi)*std::sqrt((F)195/((F)2*pi))*std::cos(theta/(F)2)*((F)5 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)11))/(F)2;
		case 460294:	// {7,6,6}
			return (std::exp((F)6*i*phi)*std::sqrt((F)15/pi)*((F)6 + (F)7*std::cos(theta))*std::pow(std::sin(theta/(F)2),(F)12))/(F)2;
		case 460295:	// {7,6,7}
			return (std::exp((F)6*i*phi)*std::sqrt((F)105/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)12)*std::sin(theta))/(F)2;
		case 460544:	// {7,7,0}
			return (-(F)3*std::exp((F)7*i*phi)*std::sqrt((F)715/((F)2*pi))*std::pow(std::sin(theta),(F)7))/(F)64;
		case 460545:	// {7,7,1}
			return (-(F)3*std::exp((F)7*i*phi)*std::sqrt((F)5005/pi)*std::pow(std::sin(theta/(F)2),(F)2)*std::pow(std::sin(theta),(F)6))/(F)128;
		case 460546:	// {7,7,2}
			return -(std::exp((F)7*i*phi)*std::sqrt((F)15015/((F)2*pi))*std::pow(std::cos(theta/(F)2),(F)5)*std::pow(std::sin(theta/(F)2),(F)9));
		case 460547:	// {7,7,3}
			return -(std::exp((F)7*i*phi)*std::sqrt((F)15015/pi)*std::pow(std::sin(theta/(F)2),(F)6)*std::pow(std::sin(theta),(F)4))/(F)32;
		case 460548:	// {7,7,4}
			return -(std::exp((F)7*i*phi)*std::sqrt((F)1365/pi)*std::pow(std::sin(theta/(F)2),(F)8)*std::pow(std::sin(theta),(F)3))/(F)8;
		case 460549:	// {7,7,5}
			return -(std::exp((F)7*i*phi)*std::sqrt((F)1365/pi)*std::pow(std::sin(theta/(F)2),(F)10)*std::pow(std::sin(theta),(F)2))/(F)8;
		case 460550:	// {7,7,6}
			return -(std::exp((F)7*i*phi)*std::sqrt((F)105/((F)2*pi))*std::pow(std::sin(theta/(F)2),(F)12)*std::sin(theta))/(F)2;
		case 460551:	// {7,7,7}
			return -(std::exp((F)7*i*phi)*std::sqrt((F)15/pi)*std::pow(std::sin(theta/(F)2),(F)14))/(F)2;
		default:
			throw std::domain_error{ std::string{"unknown {l,m,s} values for Y_lms with {"} +
			                         std::to_string(l) + "," +
			                         std::to_string(m) + "," +
			                         std::to_string(s) + "}"};
			return static_cast<std::complex<F>>(0);
			break;
		}
	}
} // namespace math
