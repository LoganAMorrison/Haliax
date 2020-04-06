//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_SPECIAL_FUNCTIONS_BESSEL_HPP
#define LANRE_SPECIAL_FUNCTIONS_BESSEL_HPP

#include "lanre/special_functions/chbevl.hpp"
#include "lanre/special_functions/besseli.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <iostream>

namespace lanre {
namespace special_functions {

static const double EUL = 0.57721566490153286061;


static std::vector<double> kBESSELK0_A = {
        1.37446543561352307156E-16,
        4.25981614279661018399E-14,
        1.03496952576338420167E-11,
        1.90451637722020886025E-9,
        2.53479107902614945675E-7,
        2.28621210311945178607E-5,
        1.26461541144692592338E-3,
        3.59799365153615016266E-2,
        3.44289899924628486886E-1,
        -5.35327393233902768720E-1
};

static std::vector<double> kBESSELK0_B = {
        5.30043377268626276149E-18,
        -1.64758043015242134646E-17,
        5.21039150503902756861E-17,
        -1.67823109680541210385E-16,
        5.51205597852431940784E-16,
        -1.84859337734377901440E-15,
        6.34007647740507060557E-15,
        -2.22751332699166985548E-14,
        8.03289077536357521100E-14,
        -2.98009692317273043925E-13,
        1.14034058820847496303E-12,
        -4.51459788337394416547E-12,
        1.85594911495471785253E-11,
        -7.95748924447710747776E-11,
        3.57739728140030116597E-10,
        -1.69753450938905987466E-9,
        8.57403401741422608519E-9,
        -4.66048989768794782956E-8,
        2.76681363944501510342E-7,
        -1.83175552271911948767E-6,
        1.39498137188764993662E-5,
        -1.28495495816278026384E-4,
        1.56988388573005337491E-3,
        -3.14481013119645005427E-2,
        2.44030308206595545468E0
};

static std::vector<double> kBESSELK1_A = {
        -7.02386347938628759343E-18,
        -2.42744985051936593393E-15,
        -6.66690169419932900609E-13,
        -1.41148839263352776110E-10,
        -2.21338763073472585583E-8,
        -2.43340614156596823496E-6,
        -1.73028895751305206302E-4,
        -6.97572385963986435018E-3,
        -1.22611180822657148235E-1,
        -3.53155960776544875667E-1,
        1.52530022733894777053E0
};


static std::vector<double> kBESSELK1_B = {
        -5.75674448366501715755E-18,
        1.79405087314755922667E-17,
        -5.68946255844285935196E-17,
        1.83809354436663880070E-16,
        -6.05704724837331885336E-16,
        2.03870316562433424052E-15,
        -7.01983709041831346144E-15,
        2.47715442448130437068E-14,
        -8.97670518232499435011E-14,
        3.34841966607842919884E-13,
        -1.28917396095102890680E-12,
        5.13963967348173025100E-12,
        -2.12996783842756842877E-11,
        9.21831518760500529508E-11,
        -4.19035475934189648750E-10,
        2.01504975519703286596E-9,
        -1.03457624656780970260E-8,
        5.74108412545004946722E-8,
        -3.50196060308781257119E-7,
        2.40648494783721712015E-6,
        -1.93619797416608296024E-5,
        1.95215518471351631108E-4,
        -2.85781685962277938680E-3,
        1.03923736576817238437E-1,
        2.72062619048444266945E0
};

/**
 * Compute the modified bessel function of the second kind of order 0.
 */
double besselk0(double x) {
    double y, z;

    if (x == 0.0) {
        std::cout << "Warning: Evaluating besselk0 at x=0, which is a singular point." << std::endl;
        return HUGE_VAL;
    } else if (x < 0.0) {
        std::cout << "Warning: Evaluating besselk0 for x<0, which is not implemented. "
                     "Use besselkn for negative arguments" << std::endl;
        return NAN;
    }

    if (x <= 2.0) {
        y = x * x - 2.0;
        y = chbevl(y, kBESSELK0_A.data(), 10) - log(0.5 * x) * besseli0(x);
        return (y);
    }
    z = 8.0 / x - 2.0;
    y = exp(-x) * chbevl(z, kBESSELK0_B.data(), 25) / sqrt(x);
    return (y);
}

/**
 * Compute the modified bessel function of the second kind of order 0 scaled by
 * exp(x).
 */
double besselk0e(double x) {
    double y;

    if (x == 0.0) {
        std::cout << "Warning: Evaluating besselk0e at x=0, which is a singular point." << std::endl;
        return HUGE_VAL;
    } else if (x < 0.0) {
        std::cout << "Warning: Evaluating besselk0e for x<0, which is not implemented. "
                     "Use besselkne for negative arguments" << std::endl;
        return NAN;
    }

    if (x <= 2.0) {
        y = x * x - 2.0;
        y = chbevl(y, kBESSELK0_A.data(), 10) - log(0.5 * x) * besseli0(x);
        return (y * exp(x));
    }

    y = chbevl(8.0 / x - 2.0, kBESSELK0_B.data(), 25) / sqrt(x);
    return (y);
}

/**
 * Compute the modified bessel function of the second kind of order 1.
 */
double besselk1(double x) {
    double y, z;

    if (x == 0.0) {
        std::cout << "Warning: Evaluating besselk1 at x=0, which is a singular point." << std::endl;
        return HUGE_VAL;
    } else if (x < 0.0) {
        std::cout << "Warning: Evaluating besselk1 for x<0, which is not implemented. "
                     "Use besselkn for negative arguments" << std::endl;
        return NAN;
    }
    z = 0.5 * x;

    if (x <= 2.0) {
        y = x * x - 2.0;
        y = log(z) * besseli1(x) + chbevl(y, kBESSELK1_A.data(), 11) / x;
        return (y);
    }

    return (exp(-x) * chbevl(8.0 / x - 2.0, kBESSELK1_B.data(), 25) / sqrt(x));
}

/**
 * Compute the modified bessel function of the second kind of order 1 scaled by
 * exp(x).
 */
double besselk1e(double x) {
    double y;

    if (x == 0.0) {
        std::cout << "Warning: Evaluating besselk1 at x=0, which is a singular point." << std::endl;
        return HUGE_VAL;
    } else if (x < 0.0) {
        std::cout << "Warning: Evaluating besselk1 for x<0, which is not implemented. "
                     "Use besselkn for negative arguments" << std::endl;
        return NAN;
    }

    if (x <= 2.0) {
        y = x * x - 2.0;
        y = log(0.5 * x) * besseli1(x) + chbevl(y, kBESSELK1_A.data(), 11) / x;
        return (y * exp(x));
    }

    return (chbevl(8.0 / x - 2.0, kBESSELK1_B.data(), 25) / sqrt(x));
}

/**
 * Compute the modified bessel function of the second kind of order n.
 */
double besselkn(int nn, double x) {
    static const double EPS = std::numeric_limits<double>::epsilon();
    static const double DBL_MAX = std::numeric_limits<double>::max();
    static const double MAXLOG = 7.09782712893383996732e2;

    double k, kf, nk1f, nkf, zn, t, s, z0, z;
    double ans, fn, pn, pk, zmn, tlg, tox;
    int i, n;

    if (nn < 0)
        n = -nn;
    else
        n = nn;

    if (n > 31) {
        overf:
        //sf_error("kn", SF_ERROR_OVERFLOW, NULL);
        return HUGE_VAL;
    }

    if (x <= 0.0) {
        if (x < 0.0) {
            //sf_error("kn", SF_ERROR_DOMAIN, NULL);
            return NAN;
        } else {
            //sf_error("kn", SF_ERROR_SINGULAR, NULL);
            return HUGE_VAL;
        }
    }


    if (x > 9.55)
        goto asymp;

    ans = 0.0;
    z0 = 0.25 * x * x;
    fn = 1.0;
    pn = 0.0;
    zmn = 1.0;
    tox = 2.0 / x;

    if (n > 0) {
        /* compute factorial of n and psi(n) */
        pn = -EUL;
        k = 1.0;
        for (i = 1; i < n; i++) {
            pn += 1.0 / k;
            k += 1.0;
            fn *= k;
        }

        zmn = tox;

        if (n == 1) {
            ans = 1.0 / x;
        } else {
            nk1f = fn / n;
            kf = 1.0;
            s = nk1f;
            z = -z0;
            zn = 1.0;
            for (i = 1; i < n; i++) {
                nk1f = nk1f / (n - i);
                kf = kf * i;
                zn *= z;
                t = nk1f * zn / kf;
                s += t;
                if ((DBL_MAX - fabs(t)) < fabs(s))
                    goto overf;
                if ((tox > 1.0) && ((DBL_MAX / tox) < zmn))
                    goto overf;
                zmn *= tox;
            }
            s *= 0.5;
            t = fabs(s);
            if ((zmn > 1.0) && ((DBL_MAX / zmn) < t))
                goto overf;
            if ((t > 1.0) && ((DBL_MAX / t) < zmn))
                goto overf;
            ans = s * zmn;
        }
    }


    tlg = 2.0 * log(0.5 * x);
    pk = -EUL;
    if (n == 0) {
        pn = pk;
        t = 1.0;
    } else {
        pn = pn + 1.0 / n;
        t = 1.0 / fn;
    }
    s = (pk + pn - tlg) * t;
    k = 1.0;
    do {
        t *= z0 / (k * (k + n));
        pk += 1.0 / k;
        pn += 1.0 / (k + n);
        s += (pk + pn - tlg) * t;
        k += 1.0;
    } while (fabs(t / s) > EPS);

    s = 0.5 * s / zmn;
    if (n & 1)
        s = -s;
    ans += s;

    return (ans);



    /* Asymptotic expansion for Kn(x) */
    /* Converges to 1.4e-17 for x > 18.4 */

    asymp:

    if (x > MAXLOG) {
        //sf_error("kn", SF_ERROR_UNDERFLOW, NULL);
        return (0.0);
    }
    k = n;
    pn = 4.0 * k * k;
    pk = 1.0;
    z0 = 8.0 * x;
    fn = 1.0;
    t = 1.0;
    s = t;
    nkf = HUGE_VAL;
    i = 0;
    do {
        z = pn - pk * pk;
        t = t * z / (fn * z0);
        nk1f = fabs(t);
        if ((i >= n) && (nk1f > nkf)) {
            goto adone;
        }
        nkf = nk1f;
        s += t;
        fn += 1.0;
        pk += 2.0;
        i += 1;
    } while (fabs(t / s) > EPS);

    adone:
    ans = exp(-x) * sqrt(M_PI / (2.0 * x)) * s;
    return (ans);
}

/**
 * Compute the modified bessel function of the second kind of order n scaled by
 * exp(x).
 */
double besselkne(int nn, double x) {
    static const double EPS = std::numeric_limits<double>::epsilon();
    static const double DBL_MAX = std::numeric_limits<double>::max();
    static const double MAXLOG = 7.09782712893383996732e2;

    double k, kf, nk1f, nkf, zn, t, s, z0, z;
    double ans, fn, pn, pk, zmn, tlg, tox;
    int i, n;

    if (nn < 0)
        n = -nn;
    else
        n = nn;

    if (n > 31) {
        overf:
        //sf_error("kn", SF_ERROR_OVERFLOW, NULL);
        return HUGE_VAL;
    }

    if (x <= 0.0) {
        if (x < 0.0) {
            //sf_error("kn", SF_ERROR_DOMAIN, NULL);
            return NAN;
        } else {
            //sf_error("kn", SF_ERROR_SINGULAR, NULL);
            return HUGE_VAL;
        }
    }


    if (x > 9.55)
        goto asymp;

    ans = 0.0;
    z0 = 0.25 * x * x;
    fn = 1.0;
    pn = 0.0;
    zmn = 1.0;
    tox = 2.0 / x;

    if (n > 0) {
        /* compute factorial of n and psi(n) */
        pn = -EUL;
        k = 1.0;
        for (i = 1; i < n; i++) {
            pn += 1.0 / k;
            k += 1.0;
            fn *= k;
        }

        zmn = tox;

        if (n == 1) {
            ans = 1.0 / x;
        } else {
            nk1f = fn / n;
            kf = 1.0;
            s = nk1f;
            z = -z0;
            zn = 1.0;
            for (i = 1; i < n; i++) {
                nk1f = nk1f / (n - i);
                kf = kf * i;
                zn *= z;
                t = nk1f * zn / kf;
                s += t;
                if ((DBL_MAX - fabs(t)) < fabs(s))
                    goto overf;
                if ((tox > 1.0) && ((DBL_MAX / tox) < zmn))
                    goto overf;
                zmn *= tox;
            }
            s *= 0.5;
            t = fabs(s);
            if ((zmn > 1.0) && ((DBL_MAX / zmn) < t))
                goto overf;
            if ((t > 1.0) && ((DBL_MAX / t) < zmn))
                goto overf;
            ans = s * zmn;
        }
    }


    tlg = 2.0 * log(0.5 * x);
    pk = -EUL;
    if (n == 0) {
        pn = pk;
        t = 1.0;
    } else {
        pn = pn + 1.0 / n;
        t = 1.0 / fn;
    }
    s = (pk + pn - tlg) * t;
    k = 1.0;
    do {
        t *= z0 / (k * (k + n));
        pk += 1.0 / k;
        pn += 1.0 / (k + n);
        s += (pk + pn - tlg) * t;
        k += 1.0;
    } while (fabs(t / s) > EPS);

    s = 0.5 * s / zmn;
    if (n & 1)
        s = -s;
    ans += s;

    return ans * exp(x);



    /* Asymptotic expansion for Kn(x) */
    /* Converges to 1.4e-17 for x > 18.4 */

    asymp:

    if (x > MAXLOG) {
        //sf_error("kn", SF_ERROR_UNDERFLOW, NULL);
        return (0.0);
    }
    k = n;
    pn = 4.0 * k * k;
    pk = 1.0;
    z0 = 8.0 * x;
    fn = 1.0;
    t = 1.0;
    s = t;
    nkf = HUGE_VAL;
    i = 0;
    do {
        z = pn - pk * pk;
        t = t * z / (fn * z0);
        nk1f = fabs(t);
        if ((i >= n) && (nk1f > nkf)) {
            goto adone;
        }
        nkf = nk1f;
        s += t;
        fn += 1.0;
        pk += 2.0;
        i += 1;
    } while (fabs(t / s) > EPS);

    adone:
    ans = sqrt(M_PI / (2.0 * x)) * s;
    return (ans);
}



}
}

#endif // LANRE_SPECIAL_FUNCTIONS_BESSEL_HPP
