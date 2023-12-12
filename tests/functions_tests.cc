#include <iostream>
#include "functions/functions.h"
#include <complex>

using namespace polynom;
using namespace std;

int main() {
	int e = 10;
	Polynomial <int> a = (5);
	cout << a << endl;
	int* c = new int [3] {2, 4, 6};
	Polynomial <int> b(c, 3);
	for (int i = 0; i < b.size(); i++) {
		cout << i << ':' << b[i] << endl;
	}
	cout << b << endl;

	cout << b[1] << endl;

	cout << b * e << endl;

	b.set(12, 1);
	cout << b << endl;
	
	Polynomial <int> d(c, 3);
	cout << b + d << endl;

	int* f = new int [4] {2, 4, 6, 8};
	Polynomial <int> h(f, 4);
	cout << b + h << endl;

	cout << d.calculation_Polynomial_x(2) << endl;



	return 0;
}