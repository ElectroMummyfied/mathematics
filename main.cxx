#include <math.h>
#include <numeric>
#include <vector>
#include <iostream>

#pragma region class-real
	namespace numeric {
		class real;
		std::ostream& operator<<(std::ostream& os, const numeric::real& obj);
	};
	class numeric::real {
		#pragma region aliases
			public:
				typedef double float_t;
		#pragma endregion aliases
		#pragma region variables
			private:
				float_t m_real;
		#pragma endregion variables
		#pragma region initializers
			public:
				real();
				real(const float_t&& val);
		#pragma endregion initializers

		#pragma region peripheries
			public:
				real H(const real& rhs, const float_t& h_param);
		#pragma endregion peripheries

		#pragma region operator-overloading
			public:
				real operator+(const real& rhs) const;
				real operator-(const real& rhs) const;
				real operator*(const real& rhs) const;
				real operator/(const real& rhs) const;
				real operator^(const real& rhs) const;
				real operator&(const real& rhs) const;
				void operator=(const float_t val);
		#pragma endregion operator-overloading
		#pragma region type-casting
			public:
				operator float_t() const;
		#pragma endregion type-casting
		#pragma region friends
				friend std::ostream& operator<<(std::ostream& os, const real& obj);
		#pragma endregion friends
	};
	#pragma region class-real-implementation
		#pragma region initializers
			numeric::real::real() : m_real(0.0) {}
			numeric::real::real(const float_t&& val) : m_real(val) {}
		#pragma endregion initializers
		#pragma region peripheries
			numeric::real numeric::real::H (
				const real& rhs,
				const float_t& h_param
			) {
				bool sign;
				float_t abs_h;
				real result;

				sign = std::signbit(h_param);
				abs_h = std::abs(h_param);

				if((0.0 <= abs_h) && (abs_h < 1.0)) {
					result = rhs + real(1.0);
				} else if((1.0 <= abs_h) && (abs_h < 2.0)) {
					result = (*this) + rhs;
				} else if((2.0 <= abs_h) && (abs_h < 3.0)) {
					result = (*this) * rhs;
				} else if((3.0 <= abs_h) && (abs_h < 4.0)) {
					result = (*this) ^ rhs;
				} else if((4.0 <= abs_h) && (abs_h < 5.0)) {
					result = (*this) & rhs;
				} else {
					result.m_real = std::nan("");
				}

				return result;
			}
		#pragma endregion peripheries
		#pragma region operator-overloading
			#pragma region operator+
				numeric::real
				numeric::real::operator+ (
					const real& rhs
				) const {
					real result;
					result.m_real = this->m_real + rhs.m_real;
					return result;
				}
			#pragma endregion operator+
			#pragma region operator-
				numeric::real
				numeric::real::operator- (
					const real& rhs
				) const {
					real result;
					result.m_real = this->m_real - rhs.m_real;
					return result;
				}
			#pragma endregion operator-
			#pragma region operator*
				numeric::real
				numeric::real::operator* (
					const real& rhs
				) const {
					real result;
					result.m_real = this->m_real * rhs.m_real;
					return result;
				}
			#pragma endregion operator*
			#pragma region ooperator/
				numeric::real
				numeric::real::operator/ (
					const real& rhs
				) const {
					real result;
					result.m_real = this->m_real / rhs.m_real;
					return result;
				}
			#pragma endregion operator/
			#pragma region operator^
				numeric::real
				numeric::real::operator^(const real& rhs) const {
					real result;
					result.m_real = std::pow(this->m_real, rhs.m_real);
					return result;
				}
			#pragma endregion operator^
			#pragma region operator&
				numeric::real numeric::real::operator&(const real& rhs) const {
					real result;
					float_t _x = std::floor(rhs.m_real);
					if(_x < 0) {
						result.m_real = std::nan("");
					} else {
						result = (*this);
					}
					for(long long it = 1; it < _x; it++) {
						result = (*this) ^ result;
					}
					return result;
				}
			#pragma endregion operator&
			#pragma region operator=
				void
				numeric::real::operator= (
					const float_t val
				){
					m_real = val;
				}
			#pragma endregion operator=
		#pragma endregion operator-overloading
		#pragma region type-cast
			numeric::real::operator float_t() const {
				return m_real;
			}
		#pragma endregion type-cast
		#pragma region friends
			#pragma region operator<<
				std::ostream&
				numeric::operator<< (
					std::ostream& os,
					const numeric::real& obj
				) {
					os << obj.m_real;
					return os;
				}
			#pragma endregion operator<<
		#pragma endregion friends
	#pragma endregion class-real-implementation
#pragma endregion class-real
#pragma region class-dense_matrix
	namespace linear_algebra {
		class dense_tensor;
	};
	class linear_algebra::dense_tensor {
		#pragma region alises
			public:
				typedef numeric::real float_t;
				typedef std::uint64_t u64_t;
				typedef std::vector<float_t> data_t;
				typedef std::vector<u64_t> u64vec_t;
		#pragma endregion alises
		#pragma region varibles
			private:
				data_t m_data;
				u64vec_t m_dims;
				u64vec_t m_strides;
		#pragma endregion varibles
		#pragma region helpers
			private:
				u64vec_t strides();
				u64_t id_nochecks (const u64vec_t& _addr);
				u64vec_t addr_nochecks (u64_t _id);
			public:
				u64_t id (const u64vec_t& _addr);
				u64vec_t addr (u64_t _id);
		#pragma endregion helpers
		#pragma region initializers
			public:
				void set_size (const u64vec_t& dims, const float_t& default_val = 0);
		#pragma endregion initializers
		#pragma region getters
			u64_t rank();
			u64_t n_elem();
		#pragma endregion getters

		#pragma region peripheries
			public:
				float_t& at (const u64_t& id);
				float_t& at (const u64vec_t& id);
		#pragma endregion peripheries
	};
	#pragma region class-dense_tensor-implementation
		#pragma region helpers	
			#pragma region strides
				linear_algebra::dense_tensor::u64vec_t
				linear_algebra::dense_tensor::strides() {
					u64vec_t result(rank(), 1);

					for(u64_t it = 1; it < rank(); it++) {
						result.at(it) = m_dims.at(it - 1) * result.at(it - 1);	
					}

					return result;
				}
			#pragma endregion strides
			#pragma region id_nochecks
				linear_algebra::dense_tensor::u64_t
				linear_algebra::dense_tensor::id_nochecks (
					const	u64vec_t& _addr 
				) {
					u64_t result = 0;

					for(u64_t it = 0; it < rank(); it++) {
						result = m_strides.at(it) * _addr.at(it);
					}
					return result;
				}	
			#pragma endregion id_nochecks
			#pragma region id
				linear_algebra::dense_tensor::u64_t
				linear_algebra::dense_tensor::id (
					const	u64vec_t& _addr 
				) {
					return id_nochecks(_addr);
				}	
			#pragma endregion id
			#pragma region addr_nochecks
				linear_algebra::dense_tensor::u64vec_t
				linear_algebra::dense_tensor::addr_nochecks (
					u64_t _id
				) {
					u64vec_t result;
					u64vec_t::reverse_iterator addr_it;
					
					result.assign(m_dims.size(), 0);	
					addr_it = result.rbegin();
					for (
						u64vec_t::const_reverse_iterator it = m_strides.crbegin();
						it != m_strides.crend();
						it++
					) {
						*addr_it = _id/(*it);	
						_id %= (*it);
						addr_it++;
					}

					return result;	
				}
			#pragma endregion addr_nochecks
			#pragma region addr
				linear_algebra::dense_tensor::u64vec_t
				linear_algebra::dense_tensor::addr (
					u64_t _id
				) {
					return addr_nochecks(_id);	
				}
			#pragma endregion addr
		#pragma endregion helpers	
		#pragma region initializers
			#pragma region set_size
				void
				linear_algebra::dense_tensor::set_size (
					const u64vec_t& dims,
					const float_t& default_val
				)
				{
					u64_t _n_elem = 1;
					m_dims = dims;
					m_strides = strides();

					for(u64_t it = 0; it < dims.size(); it++) {
						_n_elem *= dims.at(it);
					}	

					m_data.assign(_n_elem, default_val);
				}
			#pragma endregion set_size
		#pragma endregion initializers
		#pragma region getters
			#pragma region rank
				linear_algebra::dense_tensor::u64_t
				linear_algebra::dense_tensor::rank ()
				{
					return m_dims.size();
				}
			#pragma endregion rank
			#pragma region n_elem
				linear_algebra::dense_tensor::u64_t
				linear_algebra::dense_tensor::n_elem()
				{
					return m_data.size();
				}
			#pragma endregion n_elem
		#pragma endregion getters
		#pragma region peripheries
				#pragma region at
					linear_algebra::dense_tensor::float_t&
					linear_algebra::dense_tensor::at (
						const u64_t& _id
					) {
						return m_data.at(_id);
					}
					linear_algebra::dense_tensor::float_t&
					linear_algebra::dense_tensor::at (
						const u64vec_t& _addr
					)
					{
						// if(_addr.size() != rank()) {
						// 	throw;
						// }
						return m_data.at(id_nochecks(_addr));
					}
				#pragma endregion at
		#pragma endregion peripheries
	#pragma endregion class-dense_tensor-implementation
#pragma endregion class-dense_matrix

int main(int argc, char* argv[]) {
	linear_algebra::dense_tensor A, B;
	auto vec = A.addr(0);

	A.set_size({{3, 4, 2}}, 1.0);

	std::cout << "\n";
	for(int it = 0; it < A.n_elem(); it++) {
		vec = A.addr(it);
		for(int it2 = 0; it2 < vec.size(); it2++) { std::cout << vec.at(it2); }
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "rank: " << A.rank() << "\n";
	std::cout << "# entries: " << A.n_elem() << "\n";

	#pragma region output
		// std::cout << "\nclass real: " << sizeof(double) << " bytes" << "\n\n";
		// std::cout << "\ta = " << a << "\n";
		// std::cout << "\tb = " << b << "\n\n";

		// c = a + b;
		// std::cout << "\tc = a + b : " << c << "\n";
		// c = a - b;
		// std::cout << "\tc = a - b : " << c << "\n";
		// c = a * b;
		// std::cout << "\tc = a * b : " << c << "\n";
		// c = a / b;
		// std::cout << "\tc = a / b : " << c << "\n";
		// c = a ^ b;
		// std::cout << "\tc = a ^ b : " << c << "\n";
		// c = a & b;
		// std::cout << "\tc = a & b : " << c << "\n";

		// c = a.H(b, 4.5);
		// std::cout << "\tc = a.H(b) : " << c << "\n";

		// std::cout << std::endl;
	#pragma endregion output
	return 0;
}