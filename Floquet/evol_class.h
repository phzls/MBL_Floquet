#include <Eigen/Core>
#include <cmath>

using namespace std;

/*
A base class defines the evolution of a quantum system. The size of the system needs to
be passed in. In this case the local dimension at each site is assumed to be 2, and the 
total dimension is calculated. The local dimension can also be specified when constructing
the class.
*/

class EvolMatrix
{
	protected:
		const int size_; // Size of the system
		const int dim_; // Dimension of the space given size
		const int local_dim_; // Local dimension of the Hilber space
		virtual void Matrix_Init() = 0; // Initial corresponding matrices/vectors
	public:
		EvolMatrix(int size): 
			size_(size), local_dim_(2), dim_(1 << size){}
		EvolMatrix(int size, int local_dim):
			size_(size), local_dim_(local_dim), dim_(int(pow(double(local_dim), double(size)))){}
		virtual ~EvolMatrix();
};