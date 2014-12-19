#include <cmath>
#include <string>

using namespace std;

/*
A pure base structure of parameters which gives data used for the construction of time evolution
operator.
*/
template <class T>
struct ModelConstructData
{
	T data_set; // The data set that will be used
	virtual void Initialize() = 0; // Initialize the data
};

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

	public:
		// When local dimension is not given
		EvolMatrix(int size): 
			size_(size), local_dim_(2), dim_(1 << size){}

		// When local dimension is explicitly given
		EvolMatrix(int size, int local_dim):
			size_(size), local_dim_(local_dim), dim_(int(pow(double(local_dim), double(size)))){}

		// Constructing time evolution matrix
		virtual void Evol_Construct(const ModelConstructData<class T>&) = 0; 
		
		// Diagnolize time evolution matrix with eigenvectors kept
		virtual void Evol_Diag() = 0;

		// Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
		// False is not kept; True is kept
		virtual void Evol_Diag(bool keep) = 0;

		// Return the string format of representation string stream 
		virtual string Repr() const = 0; 

		virtual ~EvolMatrix();
};