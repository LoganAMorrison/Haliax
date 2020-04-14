//
// Created by Logan Morrison on 4/2/20.
//

#ifndef LANRE_INTEGRATE_QK_HPP
#define LANRE_INTEGRATE_QK_HPP

#include "lanre/integrate/base.hpp"
#include "lanre/integrate/wgt.hpp"
#include "lanre/integrate/cheb.hpp"
#include <cmath>
#include <array>

namespace lanre {
namespace integrate {


template<typename Function>
double qk15(
        Function f,
        double a,
        double b,
        double *abserr,
        double *resabs,
        double *resasc
) {
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.

    static double XGK15[8] = {
            0.99145537112081263921,
            0.94910791234275852453,
            0.86486442335976907279,
            0.74153118559939443986,
            0.58608723546769113029,
            0.40584515137739716691,
            0.20778495500789846760,
            0.00000000000000000000};
    static double WGK15[8] = {
            0.02293532201052922496,
            0.06309209262997855329,
            0.10479001032225018384,
            0.14065325971552591875,
            0.16900472663926790283,
            0.19035057806478540991,
            0.20443294007529889241,
            0.20948214108472782801};
    static double WG7[4] = {
            0.12948496616886969327,
            0.27970539148927666790,
            0.38183005050511894495,
            0.41795918367346938776};
    double fv1[7], fv2[7];
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc = f(centr);
    resg = fc * WG7[3];
    resk = fc * WGK15[7];
    *resabs = fabs(resk);
    for (j = 0; j < 3; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK15[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG7[j] * fsum;
        resk += WGK15[jtw] * fsum;
        *resabs = *resabs + WGK15[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 4; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK15[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK15[jtwm1] * fsum;
        *resabs = (*resabs) + WGK15[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK15[7] * fabs(fc - reskh);
    for (j = 0; j < 7; j++)
        *resasc = (*resasc) + WGK15[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * std::min(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = std::max(kEPMACH * 50.0 * (*resabs), (*abserr));
    return result;
}

template<typename Function>
double qk21(
        Function f,
        double a,
        double b,
        double *abserr,
        double *resabs,
        double *resasc
) {
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.

    static double XGK21[11] = {
            0.99565716302580808074,
            0.97390652851717172008,
            0.93015749135570822600,
            0.86506336668898451073,
            0.78081772658641689706,
            0.67940956829902440623,
            0.56275713466860468334,
            0.43339539412924719080,
            0.29439286270146019813,
            0.14887433898163121088,
            0.00000000000000000000};
    static double WGK21[11] = {
            0.01169463886737187428,
            0.03255816230796472748,
            0.05475589657435199603,
            0.07503967481091995277,
            0.09312545458369760554,
            0.10938715880229764190,
            0.12349197626206585108,
            0.13470921731147332593,
            0.14277593857706008080,
            0.14773910490133849137,
            0.14944555400291690566};
    static double WG10[5] = {
            0.06667134430868813759,
            0.14945134915058059315,
            0.21908636251598204400,
            0.26926671930999635509,
            0.29552422471475287017};
    double fv1[10], fv2[10];
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc = f(centr);
    resk = fc * WGK21[10];
    *resabs = fabs(resk);
    for (j = 0; j < 5; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK21[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG10[j] * fsum;
        resk += WGK21[jtw] * fsum;
        *resabs = *resabs + WGK21[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 5; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK21[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK21[jtwm1] * fsum;
        *resabs = (*resabs) + WGK21[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK21[10] * fabs(fc - reskh);
    for (j = 0; j < 10; j++)
        *resasc = (*resasc) + WGK21[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * std::min(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = std::max(kEPMACH * 50.0 * (*resabs), (*abserr));
    return result;
}

template<typename Function>
double qk31(
        Function f,
        double a,
        double b,
        double *abserr,
        double *resabs,
        double *resasc
) {
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.

    static double XGK31[16] = {
            0.99800229869339706029,
            0.98799251802048542849,
            0.96773907567913913426,
            0.93727339240070590431,
            0.89726453234408190088,
            0.84820658341042721620,
            0.79041850144246593297,
            0.72441773136017004742,
            0.65099674129741697053,
            0.57097217260853884754,
            0.48508186364023968069,
            0.39415134707756336990,
            0.29918000715316881217,
            0.20119409399743452230,
            0.10114206691871749903,
            0.00000000000000000000};
    static double WGK31[16] = {
            0.00537747987292334899,
            0.01500794732931612254,
            0.02546084732671532019,
            0.03534636079137584622,
            0.04458975132476487661,
            0.05348152469092808727,
            0.06200956780067064029,
            0.06985412131872825871,
            0.07684968075772037889,
            0.08308050282313302104,
            0.08856444305621177065,
            0.09312659817082532123,
            0.09664272698362367851,
            0.09917359872179195933,
            0.10076984552387559504,
            0.10133000701479154902};
    static double WG15[8] = {
            0.03075324199611726835,
            0.07036604748810812471,
            0.10715922046717193501,
            0.13957067792615431445,
            0.16626920581699393355,
            0.18616100001556221103,
            0.19843148532711157646,
            0.20257824192556127288};

    double fv1[15], fv2[15];
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc = f(centr);
    resg = fc * WG15[7];
    resk = fc * WGK31[15];
    *resabs = fabs(resk);
    for (j = 0; j < 7; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK31[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG15[j] * fsum;
        resk += WGK31[jtw] * fsum;
        *resabs = *resabs + WGK31[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 8; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK31[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK31[jtwm1] * fsum;
        *resabs = (*resabs) + WGK31[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK31[15] * fabs(fc - reskh);
    for (j = 0; j < 15; j++)
        *resasc = (*resasc) + WGK31[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * std::min(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = std::max(kEPMACH * 50.0 * (*resabs), (*abserr));
    return result;
}


template<typename Function>
double qk41(
        Function f,
        double a,
        double b,
        double *abserr,
        double *resabs,
        double *resasc
) {
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.

    static double XGK41[21] = {
            0.99885903158827766384,
            0.99312859918509492479,
            0.98150787745025025919,
            0.96397192727791379127,
            0.94082263383175475352,
            0.91223442825132590587,
            0.87827681125228197608,
            0.83911697182221882339,
            0.79504142883755119835,
            0.74633190646015079261,
            0.69323765633475138481,
            0.63605368072651502545,
            0.57514044681971031534,
            0.51086700195082709800,
            0.44359317523872510320,
            0.37370608871541956067,
            0.30162786811491300432,
            0.22778585114164507808,
            0.15260546524092267551,
            0.07652652113349733375,
            0.00000000000000000000};
    static double WGK41[21] = {
            0.00307358371852053150,
            0.00860026985564294220,
            0.01462616925697125298,
            0.02038837346126652360,
            0.02588213360495115883,
            0.03128730677703279896,
            0.03660016975820079803,
            0.04166887332797368626,
            0.04643482186749767472,
            0.05094457392372869193,
            0.05519510534828599474,
            0.05911140088063957237,
            0.06265323755478116803,
            0.06583459713361842211,
            0.06864867292852161935,
            0.07105442355344406831,
            0.07303069033278666750,
            0.07458287540049918899,
            0.07570449768455667466,
            0.07637786767208073671,
            0.07660071191799965645};
    static double WG20[10] = {
            0.01761400713915211831,
            0.04060142980038694133,
            0.06267204833410906357,
            0.08327674157670474872,
            0.10193011981724043504,
            0.11819453196151841731,
            0.13168863844917662690,
            0.14209610931838205133,
            0.14917298647260374679,
            0.15275338713072585070};
    double fv1[20], fv2[20];
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc = f(centr);
    resk = fc * WGK41[20];
    *resabs = fabs(resk);
    for (j = 0; j < 10; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK41[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG20[j] * fsum;
        resk += WGK41[jtw] * fsum;
        *resabs = *resabs + WGK41[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 10; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK41[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK41[jtwm1] * fsum;
        *resabs = (*resabs) + WGK41[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK41[20] * fabs(fc - reskh);
    for (j = 0; j < 20; j++)
        *resasc = (*resasc) + WGK41[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * std::min(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = std::max(kEPMACH * 50.0 * (*resabs), (*abserr));
    return result;
}


template<typename Function>
double qk51(
        Function f,
        double a,
        double b,
        double *abserr,
        double *resabs,
        double *resasc
) {
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.

    static double XGK51[26] = {
            0.99926210499260983419,
            0.99555696979049809791,
            0.98803579453407724764,
            0.97666392145951751150,
            0.96161498642584251242,
            0.94297457122897433941,
            0.92074711528170156175,
            0.89499199787827536885,
            0.86584706529327559545,
            0.83344262876083400142,
            0.79787379799850005941,
            0.75925926303735763058,
            0.71776640681308438819,
            0.67356636847346836449,
            0.62681009901031741279,
            0.57766293024122296772,
            0.52632528433471918260,
            0.47300273144571496052,
            0.41788538219303774885,
            0.36117230580938783774,
            0.30308953893110783017,
            0.24386688372098843205,
            0.18371893942104889202,
            0.12286469261071039639,
            0.06154448300568507889,
            0.00000000000000000000};
    static double WGK51[26] = {
            0.00198738389233031593,
            0.00556193213535671376,
            0.00947397338617415161,
            0.01323622919557167481,
            0.01684781770912829823,
            0.02043537114588283546,
            0.02400994560695321622,
            0.02747531758785173780,
            0.03079230016738748889,
            0.03400213027432933784,
            0.03711627148341554356,
            0.04008382550403238207,
            0.04287284502017004948,
            0.04550291304992178891,
            0.04798253713883671391,
            0.05027767908071567196,
            0.05236288580640747586,
            0.05425112988854549014,
            0.05595081122041231731,
            0.05743711636156783285,
            0.05868968002239420796,
            0.05972034032417405998,
            0.06053945537604586295,
            0.06112850971705304831,
            0.06147118987142531666,
            0.06158081806783293508};
    static double WG25[13] = {
            0.01139379850102628795,
            0.02635498661503213726,
            0.04093915670130631266,
            0.05490469597583519193,
            0.06803833381235691721,
            0.08014070033500101801,
            0.09102826198296364981,
            0.10053594906705064420,
            0.10851962447426365312,
            0.11485825914571164834,
            0.11945576353578477223,
            0.12224244299031004169,
            0.12317605372671545120};

    double fv1[25], fv2[25];
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc = f(centr);
    resg = fc * WG25[12];
    resk = fc * WGK51[25];
    *resabs = fabs(resk);
    for (j = 0; j < 12; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK51[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG25[j] * fsum;
        resk += WGK51[jtw] * fsum;
        *resabs = *resabs + WGK51[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 13; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK51[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK51[jtwm1] * fsum;
        *resabs = (*resabs) + WGK51[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK51[25] * fabs(fc - reskh);
    for (j = 0; j < 25; j++)
        *resasc = (*resasc) + WGK51[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * std::min(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = std::max(kEPMACH * 50.0 * (*resabs), (*abserr));
    return result;
}


template<typename Function>
double qk61(
        Function f,
        double a,
        double b,
        double *abserr,
        double *resabs,
        double *resasc
) {
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.

    static double XGK61[31] = {
            0.99948441005049063757,
            0.99689348407464954027,
            0.99163099687040459486,
            0.98366812327974720997,
            0.97311632250112626837,
            0.96002186496830751222,
            0.94437444474855997942,
            0.92620004742927432588,
            0.90557330769990779855,
            0.88256053579205268154,
            0.85720523354606109896,
            0.82956576238276839744,
            0.79972783582183908301,
            0.76777743210482619492,
            0.73379006245322680473,
            0.69785049479331579693,
            0.66006106412662696137,
            0.62052618298924286114,
            0.57934523582636169176,
            0.53662414814201989926,
            0.49248046786177857499,
            0.44703376953808917678,
            0.40040125483039439254,
            0.35270472553087811347,
            0.30407320227362507737,
            0.25463692616788984644,
            0.20452511668230989144,
            0.15386991360858354696,
            0.10280693796673703015,
            0.05147184255531769583,
            0.00000000000000000000};
    static double WGK61[31] = {
            0.00138901369867700762,
            0.00389046112709988405,
            0.00663070391593129217,
            0.00927327965951776343,
            0.01182301525349634174,
            0.01436972950704580481,
            0.01692088918905327263,
            0.01941414119394238117,
            0.02182803582160919230,
            0.02419116207808060137,
            0.02650995488233310161,
            0.02875404876504129284,
            0.03090725756238776247,
            0.03298144705748372603,
            0.03497933802806002414,
            0.03688236465182122922,
            0.03867894562472759295,
            0.04037453895153595911,
            0.04196981021516424615,
            0.04345253970135606932,
            0.04481480013316266319,
            0.04605923827100698812,
            0.04718554656929915395,
            0.04818586175708712914,
            0.04905543455502977889,
            0.04979568342707420636,
            0.05040592140278234684,
            0.05088179589874960649,
            0.05122154784925877217,
            0.05142612853745902593,
            0.05149472942945156756};
    static double WG30[15] = {
            0.00796819249616660562,
            0.01846646831109095914,
            0.02878470788332336935,
            0.03879919256962704960,
            0.04840267283059405290,
            0.05749315621761906648,
            0.06597422988218049513,
            0.07375597473770520627,
            0.08075589522942021535,
            0.08689978720108297980,
            0.09212252223778612872,
            0.09636873717464425964,
            0.09959342058679526706,
            0.10176238974840550460,
            0.10285265289355884034};

    double fv1[30], fv2[30];
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc = f(centr);
    resk = fc * WGK61[30];
    *resabs = fabs(resk);
    for (j = 0; j < 15; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK61[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG30[j] * fsum;
        resk += WGK61[jtw] * fsum;
        *resabs = *resabs + WGK61[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 15; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK61[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk += WGK61[jtwm1] * fsum;
        *resabs = *resabs + WGK61[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK61[30] * fabs(fc - reskh);
    for (j = 0; j < 30; j++)
        *resasc = (*resasc) + WGK61[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * std::min(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = std::max(kEPMACH * 50.0 * (*resabs), (*abserr));
    return result;
}


template<typename Function>
double qk15i(
        Function f,
        double boun,
        int inf,
        double a,
        double b,
        double *abserr,
        double *resabs,
        double *resasc
) {
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.

    static double XGK15[8] = {
            0.99145537112081263921,
            0.94910791234275852453,
            0.86486442335976907279,
            0.74153118559939443986,
            0.58608723546769113029,
            0.40584515137739716691,
            0.20778495500789846760,
            0.00000000000000000000};
    static double WGK15[8] = {
            0.02293532201052922496,
            0.06309209262997855329,
            0.10479001032225018384,
            0.14065325971552591875,
            0.16900472663926790283,
            0.19035057806478540991,
            0.20443294007529889241,
            0.20948214108472782801};
    static double WG7[4] = {
            0.12948496616886969327,
            0.27970539148927666790,
            0.38183005050511894495,
            0.41795918367346938776};
    double fv1[8], fv2[8];
    double absc, absc1, absc2, centr, dinf;
    double fc, fsum, fval1, fval2, hlgth, resg, resk;
    double reskh, result, tabsc1, tabsc2;
    int j;

    dinf = std::min((double) (1.0), (double) inf);
    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    tabsc1 = boun + dinf * (1.0 - centr) / centr;
    fval1 = f(tabsc1);
    if (inf == 2)
        fval1 += f(-tabsc1);
    fc = (fval1 / centr) / centr;
    resg = fc * WG7[3];
    resk = fc * WGK15[7];
    *resabs = fabs(resk);
    for (j = 0; j < 7; j++) {
        absc = hlgth * XGK15[j];
        absc1 = centr - absc;
        absc2 = centr + absc;
        tabsc1 = boun + dinf * (1.0 - absc1) / absc1;
        tabsc2 = boun + dinf * (1.0 - absc2) / absc2;
        fval1 = f(tabsc1);
        fval2 = f(tabsc2);
        if (inf == 2) {
            fval1 += f(-tabsc1);
            fval2 += f(-tabsc2);
        }
        fval1 = (fval1 / absc1) / absc1;
        fval2 = (fval2 / absc2) / absc2;
        fv1[j] = fval1;
        fv2[j] = fval2;
        fsum = fval1 + fval2;
        if (j & 1) resg += WG7[j / 2] * fsum; /* odd 'j's are truncated */
        resk += WGK15[j] * fsum;
        *resabs = (*resabs) + WGK15[j] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK15[7] * fabs(fc - reskh);
    for (j = 0; j < 7; j++)
        *resasc = (*resasc) + WGK15[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * hlgth;
    *resasc = (*resasc) * hlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * std::min(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = std::max(kEPMACH * 50.0 * (*resabs), (*abserr));
    return result;
}

template<typename Integrand, typename WeightFunc>
double qk15w(
        Integrand f,
        WeightFunc w,
        double p1,
        double p2,
        double p3,
        double p4,
        int kp,
        double a,
        double b,
        double *abserr,
        double *resabs,
        double *resasc
) {
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.

    static double XGK15[8] = {
            0.99145537112081263921,
            0.94910791234275852453,
            0.86486442335976907279,
            0.74153118559939443986,
            0.58608723546769113029,
            0.40584515137739716691,
            0.20778495500789846760,
            0.00000000000000000000};
    static double WGK15[8] = {
            0.02293532201052922496,
            0.06309209262997855329,
            0.10479001032225018384,
            0.14065325971552591875,
            0.16900472663926790283,
            0.19035057806478540991,
            0.20443294007529889241,
            0.20948214108472782801};
    static double WG7[4] = {
            0.12948496616886969327,
            0.27970539148927666790,
            0.38183005050511894495,
            0.41795918367346938776};
    double fv1[7], fv2[7];
    double absc, absc1, absc2, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc = f(centr) * (*w)(centr, p1, p2, p3, p4, kp);
    resg = fc * WG7[3];
    resk = fc * WGK15[7];
    *resabs = fabs(resk);
    for (j = 0; j < 3; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK15[jtw];
        absc1 = centr - absc;
        absc2 = centr + absc;
        fval1 = f(absc1) * (*w)(absc1, p1, p2, p3, p4, kp);
        fval2 = f(absc2) * (*w)(absc2, p1, p2, p3, p4, kp);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG7[j] * fsum;
        resk += WGK15[jtw] * fsum;
        *resabs = *resabs + WGK15[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 4; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK15[jtwm1];
        absc1 = centr - absc;
        absc2 = centr + absc;
        fval1 = f(absc1) * (*w)(absc1, p1, p2, p3, p4, kp);
        fval2 = f(absc2) * (*w)(absc2, p1, p2, p3, p4, kp);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK15[jtwm1] * fsum;
        *resabs = (*resabs) + WGK15[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK15[7] * fabs(fc - reskh);
    for (j = 0; j < 7; j++)
        *resasc = (*resasc) + WGK15[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * std::min(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = std::max(kEPMACH * 50.0 * (*resabs), (*abserr));
    return result;
}

template<typename Integrand>
double qc25c(
        Integrand f,
        double a,
        double b,
        double c,
        double *abserr,
        int *krul,
        int *neval
) {

    static double x[11] = {
            0.99144486137381041114,
            0.96592582628906828675,
            0.92387953251128675613,
            0.86602540378443864676,
            0.79335334029123516458,
            0.70710678118654752440,
            0.60876142900872063942,
            0.50000000000000000000,
            0.38268343236508977173,
            0.25881904510252076235,
            0.13052619222005159155};
    double ak22, amom0, amom1, amom2, cc, centr;
    double cheb12[13], cheb24[25], fval[25];
    double hlgth, resabs, resasc, res12, res24, u, result;
    int i, isym, k;
    int unitialized_value = 0xCCCCCCCC;
    int kp = unitialized_value;
    double p2 = unitialized_value, p3 = unitialized_value, p4 = unitialized_value;

    cc = (2.0 * c - b - a) / (b - a);
    if (fabs(cc) < 1.1) goto _10;

/*  Apply the 15-point Gauss-Kronrod scheme.    */
    (*krul)--;
    result = qk15w(f, qwgtc, c, p2, p3, p4, kp, a, b, abserr, &resabs, &resasc);
    *neval = 15;
    if (resasc == *abserr) (*krul)++;
    goto _50;

/*  Use the generalized Clenshaw-Curtis method. */
    _10:
    hlgth = 0.5 * (b - a);
    centr = 0.5 * (b + a);
    *neval = 25;
    fval[0] = 0.5 * f(hlgth + centr);
    fval[12] = f(centr);
    fval[24] = 0.5 * f(centr - hlgth);
    for (i = 1; i < 12; i++) {
        u = hlgth * x[i - 1];
        isym = 24 - i;
        fval[i] = f(u + centr);
        fval[isym] = f(centr - u);
    }

    /*  Compute the Chebyshev series expansion. */
    qcheb(x, fval, cheb12, cheb24);

    /*  The modified Chebyshev moments are computed by forward
     *  recursion, using amom0 and amom1 as starting values.
     */
    amom0 = log(fabs((1.0 - cc) / (1.0 + cc)));
    amom1 = 2.0 + cc * amom0;
    res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
    res24 = cheb24[0] * amom0 + cheb24[1] * amom1;
    for (k = 2; k < 13; k++) {
        amom2 = 2.0 * cc * amom1 - amom0;
        ak22 = (k - 1) * (k - 1);
        if ((k / 2) * 2 != k) amom2 -= (4.0 / (ak22 - 1.0));
        res12 += (cheb12[k] * amom2);
        res24 += (cheb24[k] * amom2);
        amom0 = amom1;
        amom1 = amom2;
    }
    for (k = 13; k < 25; k++) {
        amom2 = 2.0 * cc * amom1 - amom0;
        ak22 = (k - 1) * (k - 1);
        if ((k / 2) * 2 != k) amom2 -= (4.0 / (ak22 - 1.0));
        res24 += (cheb24[k] * amom2);
        amom0 = amom1;
        amom1 = amom2;
    }
    result = res24;
    *abserr = fabs(res24 - res12);
    _50:
    return result;
}


template<typename Integrand>
double qc25o(
        Integrand f,
        double a,
        double b,
        double omega,
        int sincos,
        int nrmom,
        int maxp1,
        int ksave,
        double *abserr,
        int *neval,
        double *resabs,
        double *resasc,
        int *momcom,
        double **chebmo
) {

    static double x[11] = {
            0.99144486137381041114,
            0.96592582628906828675,
            0.92387953251128675613,
            0.86602540378443864676,
            0.79335334029123516458,
            0.70710678118654752440,
            0.60876142900872063942,
            0.50000000000000000000,
            0.38268343236508977173,
            0.25881904510252076235,
            0.13052619222005159155};
    double ac, an, an2, as, asap, ass, centr, conc, cons, cospar;
    double estc, ests, hlgth, parint, par2, par22;
    double resc12, resc24, ress12, ress24, result, sinpar;
    double cheb12[13], cheb24[25], d[28], d1[28], d2[28];
    double d3[28], fval[25], v[28];
    int unitialized_value = 0xCCCCCCCC;
    double p2 = unitialized_value, p3 = unitialized_value, p4 = unitialized_value;
    int i, isym, j, k, m, noequ, noeq1, mm1;

    const int COSINE = 1;
    const int SINE = 2;

    centr = 0.5 * (b + a);
    hlgth = 0.5 * (b - a);
    parint = omega * hlgth;

/* Compute the integral using the 15-point Gauss-Kronrod formula
 * if the value of the parameter in the integrand is small or
 * is less than (bb-aa)/2^(maxp1-2), where (aa,bb) is the original
 * integration interval.
 */
    if (fabs(parint) > 2.0) goto _10;
    result = qk15w(f, qwgto, omega, p2, p3, p4, sincos, a, b,
                   abserr, resabs, resasc);
    *neval = 15;
    goto _190;

    /* Compute the integral using the generalized Clenshaw-Curtis method. */
    _10:
    conc = hlgth * cos(centr * omega);
    cons = hlgth * sin(centr * omega);
    *resasc = kOFLOW;
    *neval = 25;

/* Check whether the Chebyshev moments for this interval have
 * already been computed.
 */
    if ((nrmom < *momcom) || (ksave == 1)) goto _140;

/* Compute a new set of Chebyshev moments. */
    m = *momcom + 1;
/*** Add variable mm1 to ease transliteration from FORTRAN array
 *** indexing to C indexing.
 ***/
    mm1 = m - 1;
    par2 = parint * parint;
    par22 = par2 + 2.0;
    sinpar = sin(parint);
    cospar = cos(parint);

/* Compute the Chebyshev moments with respect to cosine. */
    v[0] = 2.0 * sinpar / parint;
    v[1] = (8.0 * cospar + (par2 + par2 - 8.0) * sinpar / parint) / par2;
    v[2] = (32.0 * (par2 - 12.0) * cospar + (2.0 * ((par2 - 80.0) *
            par2 + 192.0) * sinpar) / parint) / (par2 * par2);
    ac = 8.0 * cospar;
    as = 24.0 * parint * sinpar;
    if (fabs(parint) > 24.0) goto _70;

    /* Compute the Chebyshev moments as the solutions of a boundary
     * value problem with 1 initial value (v[2]) and 1 end value
     * (computed using an asymptotic formula).
     */

    noequ = 27 - 3;
    noeq1 = noequ - 1;
    an = 6.0;
    for (k = 0; k <= noeq1; k++) {
        an2 = an * an;
        d[k] = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
        d2[k] = (an - 1.0) * (an - 2.0) * par2;
        d1[k] = (an + 3.0) * (an + 4.0) * par2;
        v[k + 3] = as - (an2 - 4.0) * ac;
        an += 2.0;
    }
    an2 = an * an;
    d[noequ] = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
    v[noequ + 3] = as - (an2 - 4.0) * ac;
    v[3] -= (56.0 * par2 * v[2]);
    ass = parint * sinpar;
    asap = (((((210.0 * par2 - 1.0) * cospar - (105.0 * par2 - 63.0) *
            ass) / an2 - (1.0 - 15.0 * par2) * cospar + 15.0 * ass) /
            an2 - cospar + 3.0 * ass) / an2 - cospar) / an2;
    v[noequ + 3] -= (2.0 * asap * par2 * (an - 1.0) * (an - 2.0));
/* Solve the tridiagonal system by means of Gaussian elimination
 * with partial pivoting.
 */
    for (i = 0; i <= noequ; i++)
        d3[i] = 0.0;
    d2[noequ] = 0.0;
    for (i = 0; i <= noeq1; i++) {
        if (fabs(d1[i]) <= fabs(d[i])) goto _40;
        an = d1[i];
        d1[i] = d[i];
        d[i] = an;
        an = d2[i];
        d2[i] = d[i + 1];
        d[i + 1] = an;
        d3[i] = d2[i + 1];
        d2[i + 1] = 0.0;
        an = v[i + 4];
        v[i + 4] = v[i + 3];
        v[i + 3] = an;
        _40:
        d[i + 1] -= (d2[i] * d1[i] / d[i]);
        d2[i + 1] -= (d3[i] * d1[i] / d[i]);
        v[i + 4] -= (v[i + 3] * d1[i] / d[i]);
    }
    v[noequ + 3] /= d[noequ];
    v[noequ + 2] = (v[noequ + 2] - d2[noeq1] * v[noequ + 3]) / d[noeq1];
    for (i = 1; i <= noeq1; i++) {
        k = noequ - i - 1;
        v[k + 3] = (v[k + 3] - d3[k] * v[k + 5] - d2[k] * v[k + 4]) / d[k];
    }
    goto _90;

/* Compute the Chebyshev moments by means of forward recursion. */
    _70:
    an = 4.0;
    for (i = 3; i < 13; i++) {
        an2 = an * an;
        v[i] = ((an2 - 4.0) * (2.0 * (par22 - an2 - an2) *
                v[i - 1] - ac) + as - par2 * (an + 1.0) *
                (an + 2.0) * v[i - 2]) / (par2 * (an - 1.0) *
                (an - 2.0));
        an += 2.0;
    }
    _90:
    for (j = 0; j < 13; j++)
        chebmo[mm1][2 * j] = v[j];

/* Compute the Chebyshev moments with respect to sine. */
    v[0] = 2.0 * (sinpar - parint * cospar) / par2;
    v[1] = (18.0 - 48.0 / par2) * sinpar / par2 + (-2.0 + 48.0 / par2) *
            cospar / parint;
    ac = -24.0 * parint * cospar;
    as = -8.0 * sinpar;
    chebmo[mm1][1] = v[0];
    chebmo[mm1][3] = v[1];
    if (fabs(parint) > 24.0) goto _120;
    for (k = 2; k < 12; k++) {
        an = k + 1;
        chebmo[mm1][2 * k + 1] = -sinpar / (an * (2.0 * an - 2.0)) -
                0.25 * parint * (v[k + 1] / an - v[k] / (an - 1.0));
    }
    goto _140;

/* Compute the Chebyshev moments by means of forward recursion. */
    _120:
    an = 3.0;
    for (i = 2; i < 12; i++) {
        an2 = an * an;
        v[i] = ((an2 - 4.0) * (2.0 * (par22 - an2 - an2) *
                v[i - 1] + as) + ac - par2 * (an + 1.0) *
                (an + 2.0) * v[i - 2]) / (par2 * (an - 1.0) *
                (an - 2.0));
        an += 2.0;
        chebmo[mm1][2 * i + 1] = v[i];
    }
    _140:
    if (nrmom < *momcom) {
        m = nrmom + 1;
        mm1 = m - 1;
    }
    if ((*momcom < (maxp1 - 1)) && (nrmom >= (*momcom)))
        (*momcom)++;

/* Compute the coefficients of the Chebyshev expansions of degrees
 * 12 and 24 of the function f.
 */
    fval[0] = 0.5 * f(centr + hlgth);
    fval[12] = f(centr);
    fval[24] = 0.5 * f(centr - hlgth);
    for (i = 1; i < 12; i++) {
        isym = 24 - i;
        fval[i] = f(hlgth * x[i - 1] + centr);
        fval[isym] = f(centr - hlgth * x[i - 1]);
    }

    qcheb(x, fval, cheb12, cheb24);

/* Compute the integral and error estimates. */
    resc12 = cheb12[12] * chebmo[mm1][12];
    ress12 = 0.0;
    estc = fabs(cheb24[24] * chebmo[mm1][24]) + fabs((cheb12[12] -
            cheb24[12]) * chebmo[mm1][12]);
    ests = 0.0;
    k = 10;
    for (j = 0; j < 6; j++) {
        resc12 += (cheb12[k] * chebmo[mm1][k]);
        ress12 += (cheb12[k + 1] * chebmo[mm1][k + 1]);
        estc += fabs((cheb12[k] - cheb24[k]) * chebmo[mm1][k]);
        ests += fabs((cheb12[k + 1] - cheb24[k + 1]) * chebmo[mm1][k + 1]);
        k -= 2;
    }
    resc24 = cheb24[24] * chebmo[mm1][24];
    ress24 = 0.0;
    *resabs = fabs(cheb24[24]);
    k = 22;
    for (j = 0; j < 12; j++) {
        resc24 += (cheb24[k] * chebmo[mm1][k]);
        ress24 += (cheb24[k + 1] * chebmo[mm1][k + 1]);
        *resabs += (fabs(cheb24[k]) + fabs(cheb24[k + 1]));
        if (j <= 4) {
            estc += (fabs(cheb24[k] * chebmo[mm1][k]));
            ests += (fabs(cheb24[k + 1] * chebmo[mm1][k + 1]));
        }
        k -= 2;
    }
    *resabs *= fabs(hlgth);
    if (sincos == SINE)
        goto _180;
    result = conc * resc24 - cons * ress24;
    *abserr = fabs(conc * estc) + fabs(cons * ests);
    goto _190;
    _180:
    result = conc * ress24 + cons * resc24;
    *abserr = fabs(conc * ests) + fabs(cons * estc);
    _190:
    return result;
}

template<typename Integrand>
double qc25s(
        Integrand f,
        double a,
        double b,
        double bl,
        double br,
        double alfa,
        double beta,
        const double ri[],
        const double rj[],
        const double rg[],
        const double rh[],
        double *abserr,
        double *resasc,
        int wgtfunc,
        int *nev
) {
    static double x[11] = {
            0.99144486137381041114,
            0.96592582628906828675,
            0.92387953251128675613,
            0.86602540378443864676,
            0.79335334029123516458,
            0.70710678118654752440,
            0.60876142900872063942,
            0.50000000000000000000,
            0.38268343236508977173,
            0.25881904510252076235,
            0.13052619222005159155};
    double centr, dc, factor, fix, hlgth, resabs, res12, res24, u, result;
    double cheb12[13], cheb24[25], fval[25];
    int i, isym;

    *nev = 25;
    if ((bl == a) && ((alfa != 0.0) || (wgtfunc == 2) || (wgtfunc == 4)))
        goto _10;
    if ((br == b) && ((beta != 0.0) || (wgtfunc == 3) || (wgtfunc == 4)))
        goto _140;

/*  If a>bl and b<br, apply the 15-point Gauss-Kronrod scheme. */
    result = qk15w(f, qwgts, a, b, alfa, beta, wgtfunc, bl, br, abserr,
                   &resabs, resasc);
    *nev = 15;
    goto _270;

/*  This part is only executed if a = bl.
 *  Compute the Chebyshev series expansion of the following function:
 *  f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)^beta*
 *          f(0.5*(br-a)*x+0.5*(br+a))
 */
    _10:
    hlgth = 0.5 * (br - bl);
    centr = 0.5 * (br + bl);
    fix = b - centr;
    fval[0] = 0.5 * f(hlgth + centr) * pow(fix - hlgth, beta);
    fval[12] = f(centr) * pow(fix, beta);
    fval[24] = 0.5 * f(centr - hlgth) * pow(fix + hlgth, beta);
    for (i = 1; i < 12; i++) {
        u = hlgth * x[i - 1];
        isym = 24 - i;
        fval[i] = f(u + centr) * pow(fix - u, beta);
        fval[isym] = f(centr - u) * pow(fix + u, beta);
    }
    factor = pow(hlgth, alfa + 1.0);
    result = 0.0;
    *abserr = 0.0;
    res12 = 0.0;
    res24 = 0.0;
    if (wgtfunc > 2) goto _70;
    qcheb(x, fval, cheb12, cheb24);

/*  wgtfunc = 1  (or 2) */
    for (i = 0; i < 13; i++) {
        res12 += (cheb12[i] * ri[i]);
        res24 += (cheb24[i] * ri[i]);
    }
    for (i = 13; i < 25; i++) {
        res24 += (cheb24[i] * ri[i]);
    }
    if (wgtfunc == 1) goto _130;

/*  wgtfunc = 2 */
    dc = log(br - bl);
    result = res24 * dc;
    res12 = 0.0;
    res24 = 0.0;
    for (i = 0; i < 13; i++) {
        res12 += (cheb12[i] * rg[i]);
        res24 += (cheb24[i] * rg[i]);
    }
    for (i = 13; i < 25; i++) {
        res24 += (cheb24[i] * rg[i]);
    }
    goto _130;

/*  Compute the Chebyshev series expansion of the following function:
 *      f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
 */
    _70:
    fval[0] *= log(fix - hlgth);
    fval[12] *= log(fix);
    fval[24] *= log(fix + hlgth);
    for (i = 1; i < 12; i++) {
        u = hlgth * x[i - 1];
        isym = 24 - i;
        fval[i] *= log(fix - u);
        fval[isym] *= log(fix + u);
    }
    qcheb(x, fval, cheb12, cheb24);

/*  wgtfunc = 3  (or 4) */
    for (i = 0; i < 13; i++) {
        res12 += (cheb12[i] * ri[i]);
        res24 += (cheb24[i] * ri[i]);
    }
    for (i = 13; i < 25; i++) {
        res24 += (cheb24[i] * ri[i]);
    }
    if (wgtfunc == 3) goto _130;

/*  wgtfunc = 4 */
    dc = log(br - bl);
    result = res24 * dc;
    *abserr = fabs((res24 - res12) * dc);
    res12 = 0.0;
    res24 = 0.0;
    for (i = 0; i < 13; i++) {
        res12 += (cheb12[i] * rg[i]);
        res24 += (cheb24[i] * rg[i]);
    }
    for (i = 0; i < 13; i++) {
        res24 += (cheb24[i] * rg[i]);
    }
    _130:
    result = (result + res24) * factor;
    *abserr = (*abserr + fabs(res24 - res12)) * factor;
    goto _270;

/*  This part is executed only if b = br
 *
 *  Compute the Chebyshev series expansion of the following function:
 *
 *  f2 = (0.5 *(b+bl-a-a)+0.5*(b-bl)*x)^alfa *
 *      f(0.5*(b-bl)*x+0.5*(b+bl))
 */
    _140:
    hlgth = 0.5 * (b - bl);
    centr = 0.5 * (br + bl);
    fix = centr - a;
    fval[0] = 0.5 * f(hlgth + centr) * pow(fix + hlgth, alfa);
    fval[12] = f(centr) * pow(fix, alfa);
    fval[24] = 0.5 * f(centr - hlgth) * pow(fix - hlgth, alfa);
    for (i = 1; i < 12; i++) {
        u = hlgth * x[i - 1];
        isym = 24 - i;
        fval[i] = f(u + centr) * pow(fix + u, alfa);
        fval[isym] = f(centr - u) * pow(fix - u, alfa);
    }
    factor = pow(hlgth, beta + 1.0);
    result = 0.0;
    *abserr = 0.0;
    res12 = 0.0;
    res24 = 0.0;
    if ((wgtfunc == 2) || (wgtfunc == 4)) goto _200;

    /*  wgtfunc = 1  (or 3)  */
    qcheb(x, fval, cheb12, cheb24);
    for (i = 0; i < 13; i++) {
        res12 += (cheb12[i] * rj[i]);
        res24 += (cheb24[i] * rj[i]);
    }
    for (i = 13; i < 25; i++) {
        res24 += (cheb24[i] * rj[i]);
    }
    if (wgtfunc == 1) goto _260;

/*  wgtfunc = 3  */
    dc = log(br - bl);
    result = res24 * dc;
    *abserr = fabs((res24 - res12) * dc);
    res12 = 0.0;
    res24 = 0.0;
    for (i = 0; i < 13; i++) {
        res12 += (cheb12[i] * rh[i]);
        res24 += (cheb24[i] * rh[i]);
    }
    for (i = 13; i < 25; i++) {
        res24 += (cheb24[i] * rh[i]);
    }
    _190:
    goto _260;

    /*  Compute the Chebyshev series expansion of the following function:
     *
     *      f3 = f2 * log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
     */
    _200:
    fval[0] *= log(hlgth + fix);
    fval[12] *= log(fix);
    fval[24] *= log(fix - hlgth);
    for (i = 1; i < 12; i++) {
        u = hlgth * x[i - 1];
        isym = 24 - i;
        fval[i] *= log(u + fix);
        fval[isym] *= log(fix - u);
    }
    qcheb(x, fval, cheb12, cheb24);

/*  wgtfunc = 2  (or 4)  */
    for (i = 0; i < 13; i++) {
        res12 += (cheb12[i] * rj[i]);
        res24 += (cheb24[i] * rj[i]);
    }
    for (i = 13; i < 25; i++) {
        res24 += (cheb24[i] * rj[i]);
    }
    if (wgtfunc == 2) goto _260;
    dc = log(br - bl);
    result = res24 * dc;
    *abserr = fabs((res24 - res12) * dc);
    res12 = 0.0;
    res24 = 0.0;

/*  wgtfunc == 4  */
    for (i = 0; i < 13; i++) {
        res12 += (cheb12[i] * rh[i]);
        res24 += (cheb24[i] * rh[i]);
    }
    for (i = 13; i < 25; i++) {
        res24 += (cheb24[i] * rh[i]);
    }
    _260:
    result = (result + res24) * factor;
    *abserr = (*abserr + fabs(res24 - res12)) * factor;
    _270:
    return result;
}

}
}

#endif //LANRE_INTEGRATE_QK_HPP
