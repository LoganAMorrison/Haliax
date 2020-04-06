//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_WGT_HPP
#define LANRE_INTEGRATE_WGT_HPP

namespace lanre {
namespace integrate {

double dqwgtc(double x, double c, double p2, double p3, double p4, int kp) {
    return 1.0 / (x - c);
}

double dqwgto(double x, double omega, double p2, double p3, double p4, int wgtfunc) {
    double omx;

    omx = omega * x;
    if (wgtfunc == 1)
        return cos(omx);
    else
        return sin(omx);
}

double dqwgts(double x, double a, double b, double alpha, double beta, int wgtfunc) {
    double bmx, xma, result;

    xma = x - a;
    bmx = b - x;
    result = pow(xma, alpha) * pow(bmx, beta);
    switch (wgtfunc) {
        case 1:
            return result;
        case 2:
            return result * log(xma);
        case 3:
            return result * log(bmx);
        case 4:
            return result * log(xma) * log(bmx);
        default:
            return result;
    }
}

}
}

#endif //LANRE_INTEGRATE_WGT_HPP
