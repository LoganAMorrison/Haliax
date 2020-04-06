//
// Created by Logan Morrison on 2/21/20.
//

#ifndef LANRE_AUTODIFF_DUAL_HPP
#define LANRE_AUTODIFF_DUAL_HPP

#include "lanre/special_functions/lambertw.hpp"
#include "lanre/special_functions/besselk.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <type_traits>
#include <Eigen/Core>

namespace lanre {
namespace autodiff {

/**
 * Class to implement Dual numbers in c++
 *
 * A Dual number `z` is defined as: z = a + b eps, where `a` and `b`
 * are real numbers and `eps` is defined such that eps * eps = 0.
 * These are useful because, if we evaluate a function at a Dual
 * number, we get: f(z) = f(a) + b * f'(a) * eps. That is, we can
 * obtain the derivative of the function by grabbing the eps
 * component of the result of f(z). We can also get higher-order
 * derivative by nesting Dual numbers. For example, if x = z + w * eps
 * and z = x0 + eps and w = 1 + 0 * eps, then
 * f(z) = ((f(x0), f'(x0)), (f'(x0), f''(x0))). By templating the class,
 * we can nest Dual numbers.
 *
 * @tparam T
 */
template<typename T>
class Dual {
public:
    /**
     * The real component of the Dual number.
     */
    T val;
    /**
     * The infinitesimal component of the Dual number.
     * This component is define such that eps * eps = 0.
     */
    T eps;

    /**
     * Default construct for a Dual number. Both `val` and `eps` are set
     * to zero.
     */
    Dual() : val(static_cast<T>(0)), eps(static_cast<T>(0)) {
    }

    /**
     * Constructor to set real component of the Dual number. `eps` is
     * set to zero.
     * @param val real component of the Dual number.
     */
    template<class U>
    explicit Dual(U val) : val(static_cast<T>(val)), eps(static_cast<T>(0)) {}

    explicit Dual(T val) : val(val), eps(static_cast<T>(0)) {}

    /**
     * Full constructor for a Dual number.
     * @param val real component of the Dual number.
     * @param eps infinitesimal component of the Dual number.
     */
    template<class U1, class U2>
    Dual(U1 val, U2 eps) : val(static_cast<T>(val)), eps(static_cast<T>(eps)) {}

    /**
     * Copy constructor for a Dual number.
     * @param d the Dual number to copy.
     */
    template<class U>
    explicit Dual(const Dual<U> &d) : val(static_cast<T>(d.val)), eps(static_cast<T>(d.eps)) {}

    Dual(const Dual<T> &d) : val(d.val), eps(d.eps) {}

    /**
     * Default destructor.
     */
    ~Dual() = default;

    Dual<T> &operator=(const Dual<T> &z) {
        if (this == &z)
            return *this;
        val = z.val;
        eps = z.eps;
        return *this;
    }

    template<class U>
    Dual<T> &operator=(const U &z) {
        auto res = static_cast<Dual<T>>(z);
        val = res.val;
        eps = res.eps;
        return *this;
    }

    /**
     * Overload of the stream operator to print Dual number.
     * @param os stream.
     * @param z Dual number.
     * @return stream.
     */
    friend std::ostream &operator<<(std::ostream &os, const Dual<T> &z) {
        os << "(" << z.val << ", " << z.eps << ")";
        return os;
    }

    /**
     * Comparison operator for two Dual numbers.
     *
     * Two duals are equal if their real and infinitesimal components
     * are equal.
     *
     * @param z first Dual number.
     * @param w second Dual number.
     * @return true if the duals are equal, false if not.
     */
    friend bool operator==(const Dual<T> &z, const Dual<T> &w) {
        return z.val == w.val; // && z.eps == w.eps;
    }

    /**
     * Comparison operator for a Dual number and a real number.
     *
     * A Dual number is equal to a real number if its real part
     * is equal to the number we are comparing to and its
     * infinitesimal component is zero.
     *
     * @param z Dual number.
     * @param x real number.
     * @return true if equal, false if not.
     */
    template<class U>
    friend bool operator==(const Dual<T> &z, const U &x) {
        return static_cast<Dual<T>>(x) == z;
    }

    /**
     * Comparison operator for a Dual number and a real number.
     *
     * A Dual number is equal to a real number if its real part
     * is equal to the number we are comparing to and its
     * infinitesimal component is zero.
     *
     * @param z Dual number.
     * @param x real number.
     * @return true if equal, false if not.
     */
    template<class U>
    friend bool operator==(const U &x, const Dual<T> &z) {
        return static_cast<Dual<T>>(x) == z;
    }

    /**
    * Comparison operator for two Dual numbers.
    *
    * Dual numbers differ if their values differ or if
     * their infinitesimal components differ.
    *
    * @param z first Dual number.
    * @param w second Dual number.
    * @return true if not equal, false if equal.
    */
    friend bool operator!=(const Dual<T> &z, const Dual<T> &w) {
        return z.val != w.val; //|| z.eps != w.eps;
    }

    /**
    * Comparison operator for a Dual number and a real.
    *
    * A Dual numbers differs from a real if its values differs
    * from the real or if its infinitesimal component is non-zero.
    *
    * @param z first Dual number.
    * @param x real number.
    * @return true if not equal, false if equal.
    */
    template<class U>
    friend bool operator!=(const Dual<T> &z, const U &x) {
        return static_cast<Dual<T>>(x) != z;
    }

    /**
    * Comparison operator for a Dual number and a real.
    *
    * A Dual numbers differs from a real if its values differs
    * from the real or if its infinitesimal component is non-zero.
    *
    * @param z first Dual number.
    * @param x real number.
    * @return true if not equal, false if equal.
    */
    template<class U>
    friend bool operator!=(const U &x, const Dual<T> &z) {
        return static_cast<Dual<T>>(x) != z;
    }

    /**
    * Comparison less-than operator for two Dual numbers.
    *
    * Dual numbers are compared based on their real components.
    *
    * @param z first Dual number.
    * @param w second Dual number.
    * @return true if z.val < w.val.
    */
    friend bool operator<(const Dual<T> &z, const Dual<T> &w) {
        return z.val < w.val;
    }

    /**
    * Comparison less-than operator for a Dual number and a real.
    *
    * A Dual and a real are compared based on the real components.
    *
    * @param z Dual number.
    * @param x real number.
    * @return true if z.val < x.
    */
    template<class U>
    friend bool operator<(const Dual<T> &z, const U &x) {
        return z < static_cast<Dual<T>>(x);
    }

    /**
    * Comparison less-than operator for a Dual number and a real.
    *
    * A Dual and a real are compared based on the real components.
    *
    * @param z Dual number.
    * @param x real number.
    * @return true if z.val > x.
    */
    template<class U>
    friend bool operator<(const U &x, const Dual<T> &z) {
        return static_cast<Dual<T>>(x) < z;
    }

    /**
    * Comparison greater-than operator for two Dual numbers.
    *
    * Dual numbers are compared based on the real components.
    *
    * @param z first Dual number.
    * @param w second Dual number.
    * @return true if z.val > w.val.
    */
    friend bool operator>(const Dual<T> &z, const Dual<T> &w) {
        return z.val > w.val;
    }

    /**
    * Comparison greater-than operator for a Dual and a real.
    *
    * A Dual and real are compared based on the real components.
    *
    * @param z Dual number.
    * @param x real number.
    * @return true if z.val > x.
    */
    template<class U>
    friend bool operator>(const Dual<T> &z, const U &x) {
        return z > static_cast<Dual<T>>(x);
    }

    /**
    * Comparison greater-than operator for a Dual and a real.
    *
    * A Dual and real are compared based on the real components.
    *
    * @param z Dual number.
    * @param x real number.
    * @return true if z.val < x.
    */
    template<class U>
    friend bool operator>(const U &x, const Dual<T> &z) {
        return static_cast<Dual<T>>(x) > z;
    }

    /**
    * Comparison greater-than or equal operator for two duals.
    *
    * Duals are compared based on the real components.
    *
    * @param z first Dual number.
    * @param w second Dual number.
    * @return true if z.val >= w.val.
    */
    friend bool operator>=(const Dual<T> &z, const Dual<T> &w) {
        return z > w || z == w;
    }

    /**
    * Comparison greater-than or equal operator for a Dual and real.
    *
    * A Dual and a real are compared based on the real components.
    *
    * @param z Dual number.
    * @param x real number.
    * @return true if z.val >= x.
    */
    template<class U>
    friend bool operator>=(const Dual<T> &z, const U &x) {
        return z >= static_cast<Dual<T>>(x);
    }

    /**
    * Comparison greater-than or equal operator for a Dual and real.
    *
    * A Dual and a real are compared based on the real components.
    *
    * @param x real number.
    * @param z Dual number.
    * @return true if z.val < x.
    */
    template<class U>
    friend bool operator>=(const U &x, const Dual<T> &z) {
        return static_cast<Dual<T>>(x) >= z;
    }

    /**
    * Comparison less-than or equal operator for two duals.
    *
    * Duals are compared based on the real components.
    *
    * @param z first Dual number.
    * @param w second Dual number.
    * @return true if z.val <= w.val.
    */
    friend bool operator<=(const Dual<T> &z, const Dual<T> &w) {
        return z < w || z == w;
    }

    /**
    * Comparison less-than or equal operator for a Dual and real.
    *
    * A Dual and a real are compared based on the real components.
    *
    * @param z Dual number.
    * @param x real number.
    * @return true if z.val <= x.
    */
    template<class U>
    friend bool operator<=(const Dual<T> &z, const U &x) {
        return z <= static_cast<Dual<T>>(x);
    }

    /**
    * Comparison less-than or equal operator for a Dual and real.
    *
    * A Dual and a real are compared based on the real components.
    *
    * @param x real number.
    * @param z Dual number.
    * @return true if z.val > x.
    */
    template<class U>
    friend bool operator<=(const U &x, const Dual<T> &z) {
        return static_cast<Dual<T>>(x) <= z;
    }

    /**
     * Increment overload for duals.
     *
     * Duals are added component-wise
     *
     * @param z first Dual
     * @param w second Dual
     */
    friend void operator+=(Dual<T> &z, const Dual<T> &w) {
        z.val += w.val;
        z.eps += w.eps;
    }

    /**
     * Increment overload for duals.
     *
     * Duals are subtracted component-wise
     *
     * @param z first Dual
     * @param w second Dual
     */
    friend void operator-=(Dual<T> &z, const Dual<T> &w) {
        z.val -= w.val;
        z.eps -= w.eps;
    }

    /**
     * Increment overload for Dual and real.
     *
     * Addes to real compoenent of Dual.
     *
     * @param z first Dual
     * @param x real number.
     */
    template<class U>
    friend void operator+=(Dual<T> &z, const U &x) {
        z.val += static_cast<Dual<T>>(x).val;
    }

    /**
     * Decrement overload for Dual and real.
     *
     * Subtracts to real compoenent of Dual.
     *
     * @param z first Dual
     * @param x real number.
     */
    template<class U>
    friend void operator-=(Dual<T> &z, const U &x) {
        z.val -= static_cast<Dual<T>>(x).val;
    }

    /**
     * Addition overload for duals.
     *
     * Adds component-wise
     *
     * @param z first Dual
     * @param w second Dual.
     */
    friend Dual<T> operator+(const Dual<T> &z, const Dual<T> &w) {
        return Dual<T>(w.val + z.val, w.eps + z.eps);
    }

    /**
     * Addition overload for Dual and real.
     *
     * Adds real components
     *
     * @param z first Dual
     * @param w second Dual.
     */
    template<class U>
    friend Dual<T> operator+(const Dual<T> &z, const U &x) {
        return static_cast<Dual<T>>(x) + z;
    }

    /**
     * Addition overload for Dual and real.
     *
     * Adds real components
     *
     * @param w second Dual.
     * @param z first Dual
     */
    template<class U>
    friend Dual<T> operator+(const U &x, const Dual<T> &z) {
        return static_cast<Dual<T>>(x) + z;
    }

    /**
     * Subtraction overload for duals.
     *
     * Subtracts component-wise
     *
     * @param z first Dual
     * @param w second Dual.
     */
    friend Dual<T> operator-(const Dual<T> &z, const Dual<T> &w) {
        return Dual<T>(-w.val + z.val, -w.eps + z.eps);
    }

    /**
     * Subtraction overload for Dual and real.
     *
     * Subtracts real components
     *
     * @param w second Dual.
     * @param x real
     */
    template<class U>
    friend Dual<T> operator-(const Dual<T> &z, const U &x) {
        return z - static_cast<Dual<T>>(x);
    }

    /**
     * Subtraction overload for Dual and real.
     *
     * Subtracts real components
     *
     * @param x real
     * @param w second Dual.
     */
    template<class U>
    friend Dual<T> operator-(const U &x, const Dual<T> &z) {
        return static_cast<Dual<T>>(x) - z;
    }

    /**
     * Negation operator for duals.
     *
     * Change signs of both components of Dual.
     *
     * @param z Dual number
     * @return negated Dual.
     */
    friend Dual<T> operator-(const Dual<T> &z) {
        return Dual<T>(-z.val, -z.eps);
    }

    /**
     * Multiplication of duals
     *
     * Multiplication follows product rule.
     *
     * @param z first Dual
     * @param w second Dual
     * @return multiplied Dual
     */
    friend Dual<T> operator*(const Dual<T> &z, const Dual<T> &w) {
        return Dual<T>(w.val * z.val, w.eps * z.val + w.val * z.eps);
    }

    /**
     * Overload for *=
     * @param z Dual number
     * @param w Dual number
     */
    friend void operator*=(Dual<T> &z, const Dual<T> &w) {
        z.val = w.val * z.val;
        z.eps = w.eps * z.val + w.val * z.eps;
    }

    /**
     * Multiplication of Dual and real
     *
     * Multiplication is done component wise.
     *
     * @param z first Dual
     * @param x real
     * @return scaled Dual
     */
    template<class U>
    friend Dual<T> operator*(const Dual<T> &z, const U &x) {
        return static_cast<Dual<T>>(x) * z;
    }

    /**
     * Multiplication of Dual and real
     *
     * Multiplication is done component wise.
     *
     * @param z first Dual
     * @param x real
     * @return scaled Dual
     */
    template<class U>
    friend Dual<T> operator*(const U &x, const Dual<T> &z) {
        return static_cast<Dual<T>>(x) * z;
    }

    /**
     * Overload for *= with Dual and real.
     * @tparam U type of real
     * @param z Dual number
     * @param x real
     */
    template<class U>
    friend void operator*=(Dual<T> &z, const U &x) {
        z = z * x;
    }

    /**
     * Division of two duals.
     *
     * Division follows the quotient rule.
     *
     * @param z first Dual
     * @param w second Dual
     * @return quotent of duals.
     */
    friend Dual<T> operator/(const Dual<T> &z, const Dual<T> &w) {
        if (w.val == 0) {
            std::cout << "division by zero with two duals: " << std::endl;
            std::cout << z << std::endl;
            std::cout << w << std::endl;
            throw std::runtime_error("division by zero with two duals.");
        }
        return Dual<T>(z.val / w.val, -w.eps * z.val / pow(w.val, 2) + z.eps / w.val);
    }

    /**
     * Overload for /= with two duals
     * @param z Dual
     * @param w Dual
     */
    friend void operator/=(Dual<T> &z, const Dual<T> &w) {
        if (w.val == 0) {
            std::cout << "division by zero with Dual and Dual: " << std::endl;
            std::cout << z << std::endl;
            std::cout << w << std::endl;
            throw std::runtime_error("division by zero with two duals.");
        }
        z = z / w;
    }

    /**
     * Divide Dual by real.
     *
     * Division component-wise
     *
     * @param z first Dual
     * @param x real
     * @return scaled Dual
     */
    template<class U>
    friend Dual<T> operator/(const Dual<T> &z, const U &x) {
        if (x == 0) {
            std::cout << "division by zero with Dual and other: " << std::endl;
            std::cout << z << std::endl;
            std::cout << x << std::endl;
            throw std::runtime_error("division by zero with Dual.");
        }
        return z / static_cast<Dual<T>>(x);
    }

    template<class U>
    friend void operator/=(Dual<T> &z, const U &x) {
        if (x == 0) {
            std::cout << "division by zero with Dual and other: " << std::endl;
            std::cout << z << std::endl;
            std::cout << x << std::endl;
            throw std::runtime_error("division by zero with Dual.");
        }
        z = z / static_cast<Dual<T>>(x);
    }

    /**
     * Divide real by Dual.
     *
     * Division follow quotent rule.
     *
     * @param x real
     * @param z Dual
     * @return quotient
     */
    template<class U>
    friend Dual<T> operator/(const U &x, const Dual<T> &z) {
        if (z.val == 0) {
            std::cout << "division by zero with other and Dual: " << std::endl;
            std::cout << z << std::endl;
            std::cout << x << std::endl;
            throw std::runtime_error("division by zero with Dual.");
        }
        return static_cast<Dual<T>>(x) / z;
    }

    /**
     * Raise Dual to a Dual power.
     *
     * @param z first Dual
     * @param w second Dual
     * @return first Dual to power of second Dual
     */
    friend Dual<T> pow(const Dual<T> &z, const Dual<T> &w) {
        if (z.val == 0)
            return Dual<T>{0};
        if (z.val < 0) {
            throw std::runtime_error("error in pow(const Dual<T> &z, const Dual<T> &w).");
        }
        return Dual<T>(pow(z.val, w.val),
                       w.val * pow(z.val, w.val - 1) * z.eps +
                               w.eps * pow(z.val, w.val) * log(z.val));
    }

    /**
     * Raise Dual to real power.
     *
     * @param z first Dual
     * @param x real
     * @return first Dual to power of real
     */
    template<class U>
    friend Dual<T> pow(const Dual<T> &z, const U &x) {
        return pow(z, static_cast<Dual<T>>(x));
    }

    friend Dual<T> pow(const Dual<T> &z, double x) {
        if (z.val == 0)
            return Dual<T>{0};
        return Dual<T>(pow(z.val, x), x * pow(z.val, x - 1) * z.eps);
    }

    /**
     * Raise Dual to integrer power.
     *
     * @param z first Dual
     * @param x integer
     * @return first Dual to power of integer
     */
    friend Dual<T> pow(const Dual<T> &z, int n) {
        /*
         * Logic:
         * - If n == 0, we return 1
         * - If n == 1, we return z
         * - If n > 1, then call function recursively
         *   using z^n = z * z^n-1. We should eventually hit
         *   z^n = z * z * ... * z * z^1.
         * - If n == -1, then just return 1 / z;
         * - If n < -1, call recursively using
         *   z^(-n) = z^(1-n) / z. Should eventual hit
         *   z^(-n) = 1 / z * 1 /z * ... * z^-1.
         */
        if (n == 0) {
            return Dual<T>(1);
        } else if (n == 1) {
            return z;
        } else if (n > 1) {
            return z * pow(z, n - 1);
        } else if (n == -1) {
            return 1 / z;
        } else if (n < -1) {
            return pow(z, n + 1) / z;
        }
    }

    /**
     * Raise real to Dual power.
     *
     * @param z first Dual
     * @param x real
     * @return real to power of Dual
     */
    template<class U>
    friend Dual<T> pow(const U &x, const Dual<T> &z) {
        return pow(static_cast<Dual<T>>(x), z);
    }

    /**
     * Square root of a Dual number.
     *
     * Derivative of sqrt(x) is 1 / 2sqrt(x)
     *
     * @param z Dual number
     * @return Dual with (sqrt(x), 1 / 2sqrt(x))
     */
    friend Dual<T> sqrt(const Dual<T> &z) {
        if (z.val <= static_cast<T>(0))
            throw std::runtime_error("Invalid sqrt argument with Dual.");
        return Dual<T>(sqrt(z.val), z.eps / (2 * sqrt(z.val)));
    }

    /**
     * Absolute value of a Dual.
     * @param z Dual number
     * @return abs(z.val)
     */
    friend Dual<T> fabs(const Dual<T> &z) {
        if (z.val < 0)
            return Dual<T>(-z);
        return Dual<T>(z);
    }

    /**
     * Absolute value of a Dual.
     * @param z Dual number
     * @return abs(z.val)
     */
    friend Dual<T> abs(const Dual<T> &z) {
        if (z.val < static_cast<T>(0))
            return Dual<T>(-z);
        return Dual<T>(z);
    }

    /**
     * Sine of a Dual number.
     *
     * If z = a + b eps, then sin(z) = sin(a) + b cos(a) eps
     *
     * @param z Dual number
     * @return Dual w with w.val = sin(z.val) and
     * w.eps = z.eps * cos(z.val)
     */
    friend Dual<T> sin(const Dual<T> &z) {
        return Dual<T>(sin(z.val), z.eps * cos(z.val));
    }

    /**
     * Cosine of a Dual number.
     *
     * If z = a + b eps, then cos(z) = cos(a) - b sin(a) eps
     *
     * @param z Dual number
     * @return Dual w with w.val = cos(z.val) and
     * w.eps = -z.eps * sin(z.val)
     */
    friend Dual<T> cos(const Dual<T> &z) {
        return Dual<T>(cos(z.val), -(z.eps * sin(z.val)));
    }

    /**
     * Tangent of a Dual number.
     *
     * If z = a + b eps, then tan(z) = tan(a) + b sec(a)^2 eps
     *
     * @param z Dual number
     * @return Dual w with w.val = tan(z.val) and
     * w.eps = z.eps * sec(z.val)^2
     */
    friend Dual<T> tan(const Dual<T> &z) {
        return Dual<T>(tan(z.val), z.eps + z.eps * pow(tan(z.val), 2));
    }

    /**
     * Exponential of a Dual number.
     *
     * If z = a + b eps, then exp(z) = exp(a) + b exp(a) eps
     *
     * @param z Dual number
     * @return Dual w with w.val = exp(z.val) and
     * w.eps = z.eps * exp(z.val)
     */
    friend Dual<T> exp(const Dual<T> &z) {
        return Dual<T>(exp(z.val), exp(z.val) * z.eps);
    }

    /**
     * Natural log of a Dual number.
     *
     * If z = a + b eps, then log(z) = log(a) + (b / a) eps
     *
     * @param z Dual number
     * @return Dual w with w.val = log(z.val) and
     * w.eps = z.eps / z.val
     */
    friend Dual<T> log(const Dual<T> &z) {
        if (z.val <= 0)
            throw std::runtime_error("Invalid log argument with Dual.");
        return Dual<T>(log(z.val), z.eps / z.val);
    }

    friend Dual<T> atan(const Dual<T> &z) {
        return Dual<T>(atan(z.val), z.eps / (1.0 + z.val * z.val));
    }

    friend Dual<T> atan2(const Dual<T> &z, const Dual<T> &w) {
        return Dual<T>(atan2(z.val, w.val), (-w.val * z.eps + z.val * w.eps) / (z.val * z.val + w.val * w.val));
    }

    friend std::complex<Dual<T>> lambertw(const Dual<T> &z, int k = 0) {
        using lanre::special_functions::lambertw;
        auto lw = lambertw(z.val, k);
        auto valpart = lw;
        auto epspart = (z.eps * lw) / (z.val * (1.0 + lw));
        return std::complex<Dual<T>>{
                Dual<T>(valpart.real(), epspart.real()),
                Dual<T>(valpart.imag(), epspart.imag())
        };
    }

    friend std::complex<Dual<T>> lambertw(const std::complex<Dual<T>> &z, int k) {
        using lanre::special_functions::lambertw;
        auto lw = lambertw(z.val, k);
        auto valpart = lw;
        auto epspart = (z.eps * lw) / (z.val * (1.0 + lw));
        return std::complex<Dual<T>>{
                Dual<T>(valpart.real(), epspart.real()),
                Dual<T>(valpart.imag(), epspart.imag())
        };
    }

    friend Dual<T> besselk0(const Dual<T> &z) {
        using special_functions::besselk1;
        using special_functions::besselk0;
        return autodiff::Dual<T>(besselk0(z.val), -z.eps * besselk1(z.val));
    }

    friend Dual<T> besselk0e(const Dual<T> &z) {
        using special_functions::besselk1e;
        using special_functions::besselk0e;
        const double x = z.val;
        return autodiff::Dual<T>(besselk0e(x), besselk0e(x) - besselk1e(x));
    }

    friend Dual<T> besselk1(const Dual<T> &z) {
        using special_functions::besselk1;
        using special_functions::besselk0;
        const double x = z.val;
        return autodiff::Dual<T>(
                besselk1(x),
                -z.eps / x * (x * besselk0(x) + besselk1(x))
        );
    }

    friend Dual<T> besselk1e(const Dual<T> &z) {
        using special_functions::besselk1e;
        using special_functions::besselk0e;
        const double x = z.val;
        return autodiff::Dual<T>(
                besselk1e(x),
                z.eps * (-besselk0e(x) + (x - 1.0) / x * besselk1e(x))
        );
    }

    friend Dual<T> besselkn(int n, const Dual<T> &z) {
        using special_functions::besselkn;
        const double x = z.val;
        return autodiff::Dual<T>(
                besselkn(n, x),
                -z.eps / x * (x * besselkn(n - 1, x) + n * besselkn(n, x))
        );
    }

    friend Dual<T> besselkne(int n, const Dual<T> &z) {
        using special_functions::besselkne;
        const double x = z.val;
        return autodiff::Dual<T>(
                besselkne(n, x),
                z.eps * ((n + x) * besselkne(n, x) / x - besselkne(1 + n, x))
        );
    }
};

inline const Dual<double> &conj(const Dual<double> &x) { return x; }

inline const Dual<double> &real(const Dual<double> &x) { return x; }

inline Dual<double> imag(const Dual<double> &) { return Dual<double>(0.0); }

inline Dual<double> abs2(const Dual<double> &x) { return x * x; }

} // namespace autodiff
} // namespace lanre

namespace Eigen {

template<>
struct NumTraits<lanre::autodiff::Dual<double>>
        : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
    typedef lanre::autodiff::Dual<double> Real;
    typedef lanre::autodiff::Dual<double> NonInteger;
    typedef lanre::autodiff::Dual<double> Nested;
    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 3,
        MulCost = 3
    };
};
} // namespace Eigen

#endif //LANRE_AUTODIFF_DUAL_HPP
