#ifndef AUTODIFF_HPP
#define AUTODIFF_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

namespace ASC_ode {
    template<typename T = double>
    class AutoDiffDynamic {
    private:
        T m_val;
        std::vector<T> m_deriv; // Dynamic storage

    public:
        // Constructor for constants (no derivatives yet)
        AutoDiffDynamic(T v) : m_val(v), m_deriv() {}

        AutoDiffDynamic() : m_val(0), m_deriv() {}

        // Constructor for variables (value + known dimension of system)
        AutoDiffDynamic(T v, size_t dim) : m_val(v), m_deriv(dim, 0.0) {}

        T value() const { return m_val; }

        // Mutable access to derivatives
        std::vector<T> &deriv() { return m_deriv; }

        // Const access to derivatives
        const std::vector<T> &deriv() const { return m_deriv; }

        // Helper to check dimension
        size_t dim() const { return m_deriv.size(); }
    };

    // --- Output Operator ---
    template<typename T>
    std::ostream &operator<<(std::ostream &os, const AutoDiffDynamic<T> &ad) {
      os << "Value: " << ad.value() << ", Deriv: [";
      for (size_t i = 0; i < ad.deriv().size(); i++) {
        os << ad.deriv()[i];
        if (i < ad.deriv().size() - 1) os << ", ";
      }
      os << "]";
      return os;
    }

    // --- Arithmetic Operators ---

    template<typename T>
    size_t get_dim(const AutoDiffDynamic<T> &a, const AutoDiffDynamic<T> &b) {
      size_t da = a.dim();
      size_t db = b.dim();
      if (da > 0 && db > 0) {
        assert(da == db && "AutoDiffDynamic dimension mismatch!");
        return da;
      }
      return std::max(da, db);
    }

    // ADDITION
    template<typename T>
    AutoDiffDynamic<T> operator+(const AutoDiffDynamic<T> &a, const AutoDiffDynamic<T> &b) {
      size_t n = get_dim(a, b);
      AutoDiffDynamic<T> result(a.value() + b.value(), n);

      if (a.dim() > 0 && b.dim() > 0) {
        for (size_t i = 0; i < n; i++) result.deriv()[i] = a.deriv()[i] + b.deriv()[i];
      } else if (a.dim() > 0) {
        for (size_t i = 0; i < n; i++) result.deriv()[i] = a.deriv()[i];
      } else if (b.dim() > 0) {
        for (size_t i = 0; i < n; i++) result.deriv()[i] = b.deriv()[i];
      }
      return result;
    }

    template<typename T>
    AutoDiffDynamic<T> operator+(T a, const AutoDiffDynamic<T> &b) {
      AutoDiffDynamic<T> result(a + b.value(), b.dim());
      for (size_t i = 0; i < b.dim(); ++i) result.deriv()[i] = b.deriv()[i];
      return result;
    }

    template<typename T>
    AutoDiffDynamic<T> operator+(const AutoDiffDynamic<T> &a, T b) { return b + a; }

    // SUBTRACTION
    template<typename T>
    AutoDiffDynamic<T> operator-(const AutoDiffDynamic<T> &a, const AutoDiffDynamic<T> &b) {
      size_t n = get_dim(a, b);
      AutoDiffDynamic<T> result(a.value() - b.value(), n);

      if (a.dim() > 0 && b.dim() > 0) {
        for (size_t i = 0; i < n; i++) result.deriv()[i] = a.deriv()[i] - b.deriv()[i];
      } else if (a.dim() > 0) {
        // b is constant
        for (size_t i = 0; i < n; i++) result.deriv()[i] = a.deriv()[i];
      } else if (b.dim() > 0) {
        // a is constant
        for (size_t i = 0; i < n; i++) result.deriv()[i] = -b.deriv()[i];
      }
      return result;
    }

    template<typename T>
    AutoDiffDynamic<T> operator-(T a, const AutoDiffDynamic<T> &b) {
      AutoDiffDynamic<T> result(a - b.value(), b.dim());
      for (size_t i = 0; i < b.dim(); ++i) result.deriv()[i] = -b.deriv()[i];
      return result;
    }

    template<typename T>
    AutoDiffDynamic<T> operator-(const AutoDiffDynamic<T> &a, T b) {
      AutoDiffDynamic<T> result(a.value() - b, a.dim());
      for (size_t i = 0; i < a.dim(); ++i) result.deriv()[i] = a.deriv()[i];
      return result;
    }

    // MULTIPLICATION
    template<typename T>
    AutoDiffDynamic<T> operator*(const AutoDiffDynamic<T> &a, const AutoDiffDynamic<T> &b) {
      size_t n = get_dim(a, b);
      AutoDiffDynamic<T> result(a.value() * b.value(), n);

      if (a.dim() > 0 && b.dim() > 0) {
        for (size_t i = 0; i < n; i++)
          result.deriv()[i] = a.deriv()[i] * b.value() + a.value() * b.deriv()[i];
      } else if (a.dim() > 0) {
        // b is constant
        for (size_t i = 0; i < n; i++) result.deriv()[i] = a.deriv()[i] * b.value();
      } else if (b.dim() > 0) {
        // a is constant
        for (size_t i = 0; i < n; i++) result.deriv()[i] = a.value() * b.deriv()[i];
      }
      return result;
    }

    template<typename T>
    AutoDiffDynamic<T> operator*(T a, const AutoDiffDynamic<T> &b) {
      AutoDiffDynamic<T> result(a * b.value(), b.dim());
      for (size_t i = 0; i < b.dim(); ++i) result.deriv()[i] = a * b.deriv()[i];
      return result;
    }

    template<typename T>
    AutoDiffDynamic<T> operator*(const AutoDiffDynamic<T> &b, T a) {
      return a * b;
    }

    // DIVISION
    template<typename T>
    AutoDiffDynamic<T> operator/(const AutoDiffDynamic<T> &a, const AutoDiffDynamic<T> &b) {
      size_t n = get_dim(a, b);
      AutoDiffDynamic<T> result(a.value() / b.value(), n);

      if (a.dim() > 0 && b.dim() > 0) {
        T b2 = b.value() * b.value();
        for (size_t i = 0; i < n; i++)
          result.deriv()[i] = (a.deriv()[i] * b.value() - a.value() * b.deriv()[i]) / b2;
      } else if (a.dim() > 0) {
        // b is constant: a' / b
        T inv_b = 1.0 / b.value();
        for (size_t i = 0; i < n; i++) result.deriv()[i] = a.deriv()[i] * inv_b;
      } else if (b.dim() > 0) {
        // a is constant: -a * b' / b^2
        T factor = -a.value() / (b.value() * b.value());
        for (size_t i = 0; i < n; i++) result.deriv()[i] = factor * b.deriv()[i];
      }
      return result;
    }


    template<typename T>
    void operator+=(AutoDiffDynamic<T> &a,
                    const AutoDiffDynamic<T> &b
    ) {
      a = a + b;
    }

    template<typename T>
    void operator-=(AutoDiffDynamic<T> &a, const AutoDiffDynamic<T> &b) {
      a = a - b;
    }

    template<typename T>
    void operator+=(AutoDiffDynamic<T> &a, const T b) {
      a = a + b;
    }

    template<typename T>
    void operator-=(AutoDiffDynamic<T> &a, const T b) {
      a = a - b;
    }

    template<typename T>
    AutoDiffDynamic<T> operator/(const AutoDiffDynamic<T> &a, double b) {
      // Check for division by zero
      if (b == 0.0) {
        throw std::runtime_error("Division by zero in AutoDiffDynamic / double");
      }

      // 1. Create result with new value: a.value() / b
      // 2. Initialize gradient vector size to match 'a'
      AutoDiffDynamic<T> result(a.value() / b, a.dim());

      // 3. Compute derivatives: a.deriv()[i] / b
      // Multiplying by reciprocal (1/b) is usually slightly faster than repeated division
      double inv_b = 1.0 / b;
      for (size_t i = 0; i < a.dim(); i++) {
        result.deriv()[i] = a.deriv()[i] * inv_b;
      }

      return result;
    }
// MATH FUNCTIONS
    using std::sin;
    using std::cos;
    using std::sqrt;

    template<typename T>
    AutoDiffDynamic<T> sqrt(const AutoDiffDynamic<T> &a) {
      AutoDiffDynamic<T> result(std::sqrt(a.value()), a.dim());
      for (size_t i = 0; i < a.dim(); i++)
        result.deriv()[i] = (1.0 / (2.0 * result.value())) * a.deriv()[i];
      return result;
    }

    template<size_t N, typename T = double>
    AutoDiffDynamic<T> sin(const AutoDiffDynamic<T> &a) {
      AutoDiffDynamic<T> result(sin(a.value()));
      for (size_t i = 0; i < N; i++)
        result.deriv()[i] = cos(a.value()) * a.deriv()[i];
      return result;
    }

    template<size_t N, typename T = double>
    AutoDiffDynamic<T> cos(const AutoDiffDynamic<T> &a) {
      AutoDiffDynamic<T> result(cos(a.value()));
      for (size_t i = 0; i < N; i++)
        result.deriv()[i] = -sin(a.value()) * a.deriv()[i];
      return result;
    }

    template<size_t N, typename T = double>
    AutoDiffDynamic<T> exp(const AutoDiffDynamic<T> &a) {
      AutoDiffDynamic<T> result(exp(a.value()));
      for (size_t i = 0; i < N; i++)
        result.deriv()[i] = exp(a.value()) * a.deriv()[i];
      return result;
    }

    template<size_t N, typename T = double>
    AutoDiffDynamic<T> log(const AutoDiffDynamic<T> &a) {
      AutoDiffDynamic<T> result(log(a.value()));
      for (size_t i = 0; i < N; i++)
        result.deriv()[i] = (1 / a.value()) * a.deriv()[i];
      return result;
    }
// (Add sin/cos/exp similarly using a.dim() in loops)

// Put this in ASC_ode namespace, ideally in autodiff_dynamic.hpp
    template<typename T>
    nanoblas::Vector<AutoDiffDynamic<T>>
    operator*(const nanoblas::VecExpr<AutoDiffDynamic<T>> &v, const AutoDiffDynamic<T> &s) {
      nanoblas::Vector<AutoDiffDynamic<T>> result(v.size());
      for (int i = 0; i < v.size(); i++) {
        result(i) = s * v(i); // Multiply scalar s with each component v(i)
      }
      return result;
    }

// Put this in ASC_ode namespace, ideally in autodiff_dynamic.hpp
    template<typename T>
    nanoblas::Vector<AutoDiffDynamic<T>>
    operator*(const AutoDiffDynamic<T> &s, const nanoblas::VecExpr<AutoDiffDynamic<T>> &v) {
      nanoblas::Vector<AutoDiffDynamic<T>> result(v.size());
      for (int i = 0; i < v.size(); i++) {
        result(i) = s * v(i); // Multiply scalar s with each component v(i)
      }
      return result;
    }

    // Put this in ASC_ode namespace, ideally in autodiff_dynamic.hpp
    template<unsigned long int D, typename T>
    nanoblas::Vec<D, AutoDiffDynamic<T>>
    operator*(const AutoDiffDynamic<T> &s, const nanoblas::Vec<D, AutoDiffDynamic<T>> &v) {
      nanoblas::Vec<D, AutoDiffDynamic<T>> result;
      for (int i = 0; i < D; i++) {
        result(i) = s * v(i); // Multiply scalar s with each component v(i)
      }
      return result;
    }

    template<typename T>
    AutoDiffDynamic<T> operator/(double a, const AutoDiffDynamic<T> &b) {
      AutoDiffDynamic<T> result(a / b.value(), b.dim());
      double b_sq = b.value() * b.value();
      for (size_t i = 0; i < b.dim(); i++)
        result.deriv()[i] = -a * b.deriv()[i] / b_sq;

      return result;
    }


    template <typename T>
    AutoDiffDynamic<T> norm2 (AutoDiffDynamic<T> x) { return x*x; }
}
#endif
