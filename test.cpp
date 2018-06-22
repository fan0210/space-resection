#include "spaceresection.h"
#include <vector>

int main()
{
	SpaceResection sr;

	sr.setInitialValue({ 5195.7049,5195.7049 ,0,0,0,0,0,0,4272,2848,5000,0,-1200,25.0,0,0 });

	if (!sr.importImagePoints("像点.txt") || !sr.importControlPoints("控制点.txt"))
	{
		std::cout << "file path error." << std::endl;

		system("pause");
		return -1;
	}

	std::vector<double> elems = sr.slove("解算结果.txt");

	std::cout << "{fx, fy, x0, y0, k1, k2, p1, p2, X, Y, Z, phi, omega, kappa}" << std::endl;
	for (auto i = elems.cbegin(); i != elems.cend(); ++i)
		std::cout << *i << "  ";
	std::cout << std::endl;

	system("pause");
	return 0;
}