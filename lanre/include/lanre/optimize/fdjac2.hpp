//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_OPTIMIZE_FDJAC2_HPP
#define LANRE_OPTIMIZE_FDJAC2_HPP

#include <vector>

namespace lanre {
namespace optimize {

template<class Function>
void fdjac2(
        Function fcn,
        Vector<double> &x,
        Vector<double> &fvec,
        Matrix<double> &fjac,
        int &iflag,
        double epsfcn,
        Vector<double> &wa
) {

    const int m = fvec.size();
    const int n = x.size();

    static const double epsmch = std::numeric_limits<double>::epsilon();

    const double eps = sqrt(std::max(epsfcn, epsmch));

    for (int j = 0; j < n; j++) {
        double temp = x(j);
        double h;
        if (temp == 0.0) {
            h = eps;
        } else {
            h = eps * std::abs(temp);
        }
        x(j) = temp + h;
        fcn(x, wa, iflag);
        if (iflag < 0) return;
        x(j) = temp;
        for (int i = 0; i < m; i++) {
            fjac(i, j) = (wa(i) - fvec(i)) / h;
        }
    }
}

}
}

#endif //LANRE_OPTIMIZE_FDJAC2_HPP
