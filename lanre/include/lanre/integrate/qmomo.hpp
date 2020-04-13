//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_QMOMO_HPP
#define LANRE_INTEGRATE_QMOMO_HPP

#include <array>

namespace lanre {
namespace integrate {

void qmomo(
        double alfa,
        double beta,
        std::array<double, 25> &ri,
        std::array<double, 25> &rj,
        std::array<double, 25> &rg,
        std::array<double, 25> &rh,
        int wgtfunc
) {
    double alfp1 = alfa + 1.0;
    double betp1 = beta + 1.0;
    double alfp2 = alfa + 2.0;
    double betp2 = beta + 2.0;
    double ralf = pow(2.0, alfp1);
    double rbet = pow(2.0, betp1);

    // compute ri, rj using a forward recurrence relation.
    ri[0] = ralf / alfp1;
    rj[0] = rbet / betp1;
    ri[1] = ri[0] * alfa / alfp2;
    rj[1] = rj[0] * beta / betp2;
    double an = 2.0;
    double anm1 = 1.0;

    for (int i = 2; i < 25; i++) {
        ri[i] = -(ralf + an * (an - alfp2) * ri[i - 1]) / (anm1 * (an + alfp1));
        rj[i] = -(rbet + an * (an - betp2) * rj[i - 1]) / (anm1 * (an + betp1));
        anm1 = an;
        an = an + 1.0;
    }

    if (wgtfunc != 1) {
        if (wgtfunc != 3) {
            // compute rg using a forward recurrence relation.
            rg[0] = -ri[0] / alfp1;
            rg[1] = -(ralf + ralf) / (alfp2 * alfp2) - rg[0];
            an = 2.0;
            anm1 = 1.0;
            int im1 = 2;
            for (int i = 3; i <= 25; i++) {
                rg[i - 1] = -(an * (an - alfp2) * rg[im1 - 1] - an * ri[im1 - 1] + anm1 * ri[i - 1]) /
                        (anm1 * (an + alfp1));
                anm1 = an;
                an = an + 1.0;
                im1 = i;
            }
        }
        if (wgtfunc != 2) {
            // compute rh using a forward recurrence relation.
            rh[0] = -rj[0] / betp1;
            rh[1] = -(rbet + rbet) / (betp2 * betp2) - rh[0];
            an = 2.0;
            anm1 = 1.0;
            int im1 = 2;

            for (int i = 3; i <= 25; i++) {
                rh[i - 1] = -(an * (an - betp2) * rh[im1 - 1] - an * rj[im1 - 1] +
                        anm1 * rj[i - 1]) / (anm1 * (an + betp1));
                anm1 = an;
                an = an + 1.0;
                im1 = i;
            }
            for (int i = 1; i < 25; i += 2) {
                rh[i] *= -1;
            }
        }
    }

    // 70
    for (int i = 1; i < 25; i += 2) {
        rj[i] *= -1;
    }
}

}
}

#endif //LANRE_INTEGRATE_QMOMO_HPP
