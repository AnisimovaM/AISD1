#pragma once
#include <iostream>
#include <cmath>
#include <complex>
#include <stdexcept>
# define M_PI           3.14159265358979323846  /* pi */
using namespace std;

namespace polynom {
    template<typename T>
    class Polynomial {
    private:
        T* _data;
        size_t _size;
        inline static const double EPSILION = 0.001;
    public:
        Polynomial(size_t size) : _size(size), _data(new T[size]()) {}
        Polynomial(T* data, size_t size) : _data(new T[size]), _size(size) {
            for (size_t i = 0; i < _size; i++) {
                _data[i] = data[i];
            }
        }
        
        size_t size() const {
            return _size;
        }
        ~Polynomial() {
            delete[] _data;
            _data = nullptr;
            _size = 0;
        }
        Polynomial(const Polynomial<T>& a) : _data(new T[a._size]), _size(a._size) {
            for (size_t i = 0; i < _size; i++) {
                _data[i] = a[i];
            }
        }

        T operator[](size_t index) const {
            if (_size <= index) {
                return 0;
            }
            return _data[index];
        }

        void set(T data, size_t index) {
            if (index > _size) {
                expand(index + 1);
            }
            _data[index] = data;
        }

        size_t get_size() const {
            return _size;
        }

        T* get_data() const {
            return _data;
        }

        void expand(size_t size) {
            if (size < _size) {
                throw std::out_of_range("this power is alredy exsists.");
            }
            auto temp = (new T[size]());
            for (size_t i = 0; i < _size; i++) {
                temp[i] = _data[i];
            }
            delete[] _data;
            _data = temp;
            _size = size;
        }

        void swap(Polynomial<T>& a) {
            std::swap(_data, a._data);
            std::swap(_size, a._size);
        }

        

        Polynomial<T>& operator=(Polynomial<T> a) {
            swap(a);
            return *this;
        }

        Polynomial<T>& operator+= (const Polynomial<T>& a) {
            if (a.size() > _size) {

                expand(a.size());
            }
            for (size_t i = 0; i < a._size; i++) {
                _data[i] += a[i];
            }
            return *this;
        }

        Polynomial<T>& operator-= (const Polynomial<T>& a) {
            if (a.size() > _size) {
                expand(a.size());
            }
            for (size_t i = 0; i < a._size; i++) {
                _data[i] -= a[i];
            }
            return *this;
        }

        Polynomial<T> operator- (Polynomial<T> a) const {
            return a -= *this;
        }

        Polynomial<T> operator+ (Polynomial<T> a)  const {
            return a += *this;
        }

        Polynomial<T>& operator*(const T& a) {
            for (size_t i = 0; i < _size; i++) {
                _data[i] *= a;
            }
            return *this;
        }

        bool operator== (Polynomial<T> a) const {
            a.shrink_to_fit();
            Polynomial<T> copy(*this);
            copy.shrink_to_fit();
            if (copy._size != a._size) { return false; }
            for (size_t i = 0; i < a._size; i++) {
                if (std::abs(_data[i] - a._data[i]) > EPSILION) {
                    return false;
                }
            }
            return true;
        }

        bool operator!= (const Polynomial<T>& a) const {
            return !(a == *this);
        }

        friend Polynomial<T> operator* (const T& a, const Polynomial<T>& pol) {
            Polynomial<T> res = pol;
            res = res * a;
            return res;
        }
        T calculation_Polynomial_x(const T& x) {
            T sum = 0;
            for (size_t i = 0; i < _size; i++) {
                sum += _data[i] * pow(x, i);
            }
            return sum;
        }

        void shrink_to_fit() {
            T zero = T(0);
            for (size_t i = _size - 1; i > 0; i--) {
                if (_data[i] != zero) {
                    i++;
                    auto temp = (new T[i]());
                    for (size_t j = 0; j < i; j++) {
                        temp[j] = _data[j];
                    }
                    delete[] _data;
                    _data = temp;
                    _size = i;
                    break;
                }
            }
        }

    };
    template<typename T>
    void findCubicRoots(Polynomial<T> polynomial) {


        // Приведение к упрощенной форме
        T* mass = polynomial.get_data();
        T a = mass[3];
        T b = mass[2];
        T c = mass[1];
        T d = mass[0];

        T p = c / a - b * b / (3 * a * a);
        T q = 2 * b * b * b / (27 * a * a * a) - b * c / (3 * a * a) + d / a;

        T Q = p / 3;
        T R = q / 2;

        T D = Q * Q * Q + R * R; // Дискриминант


        if (D > 0) {
            // Один действительный корень
            T S = std::cbrt(-R + std::sqrt(D));
            T t = std::cbrt(-R - std::sqrt(D));
            T x1 = (S + t - b) / (3 * a);
            std::cout << "One real root: " << x1 << std::endl;
        }
        else {
            // Три действительных корня (D <= 0)
            T theta = std::acos(-R / std::sqrt(-Q * Q * Q));
            T x1 = 2 * std::sqrt(-Q) * std::cos(theta / 3) - b / (3 * a);
            T x2 = 2 * std::sqrt(-Q) * std::cos((theta + 2 * M_PI) / 3) - b / (3 * a);
            T x3 = 2 * std::sqrt(-Q) * std::cos((theta + 4 * M_PI) / 3) - b / (3 * a);

            std::cout << "Three real roots: " << x1 << ", " << x2 << ", " << x3 << std::endl;
        }
    }

   

    template<typename T>
    std::ostream& operator << (std::ostream& stream, const Polynomial<T>& a) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (i == 0) {
                stream << a[i];
            }
            else if (i == 1) {
                stream << " + " << a[i] << "*x";
            }
            else {
                stream << " + " << a[i] << "*x^" << i;
            }
        }
        return stream;
    }
};
