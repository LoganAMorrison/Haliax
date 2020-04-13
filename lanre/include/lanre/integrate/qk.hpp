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

    const static std::array<double, 4> wg = {
            0.129484966168869693270611432679082,
            0.279705391489276667901467771423780,
            0.381830050505118944950369775488975,
            0.417959183673469387755102040816327
    };
    const static std::array<double, 8> wgk = {
            0.022935322010529224963732008058970,
            0.063092092629978553290700663189204,
            0.104790010322250183839876322541518,
            0.140653259715525918745189590510238,
            0.169004726639267902826583426598550,
            0.190350578064785409913256402421014,
            0.204432940075298892414161999234649,
            0.209482141084727828012999174891714
    };
    const static std::array<double, 8> xgk = {
            0.991455371120812639206854697526329,
            0.949107912342758524526189684047851,
            0.864864423359769072789712788640926,
            0.741531185599394439863864773280788,
            0.586087235467691130294144838258730,
            0.405845151377397166906606412076961,
            0.207784955007898467600689403773245,
            0.000000000000000000000000000000000
    };

    std::array<double, 7> fv1 = {0}, fv2 = {0};
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;


    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc = f(centr);
    resg = fc * wg[3];
    resk = fc * wgk[7];
    *resabs = fabs(resk);
    for (j = 0; j < 3; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * xgk[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j] * fsum;
        resk += wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 4; j++) {
        jtwm1 = j * 2;
        absc = hlgth * xgk[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = (*resabs) + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = wgk[7] * fabs(fc - reskh);
    for (j = 0; j < 7; j++) {
        *resasc = (*resasc) + wgk[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    }
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = fmax(kEPMACH * 50.0 * (*resabs), (*abserr));
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

    static const std::array<double, 5> wg = {
            0.066671344308688137593568809893332,
            0.149451349150580593145776339657697,
            0.219086362515982043995534934228163,
            0.269266719309996355091226921569469,
            0.295524224714752870173892994651338
    };
    static const std::array<double, 11> wgk = {
            0.011694638867371874278064396062192,
            0.032558162307964727478818972459390,
            0.054755896574351996031381300244580,
            0.075039674810919952767043140916190,
            0.093125454583697605535065465083366,
            0.109387158802297641899210590325805,
            0.123491976262065851077958109831074,
            0.134709217311473325928054001771707,
            0.142775938577060080797094273138717,
            0.147739104901338491374841515972068,
            0.149445554002916905664936468389821
    };
    static const std::array<double, 11> xgk = {
            0.995657163025808080735527280689003,
            0.973906528517171720077964012084452,
            0.930157491355708226001207180059508,
            0.865063366688984510732096688423493,
            0.780817726586416897063717578345042,
            0.679409568299024406234327365114874,
            0.562757134668604683339000099272694,
            0.433395394129247190799265943165784,
            0.294392862701460198131126603103866,
            0.148874338981631210884826001129720,
            0.000000000000000000000000000000000
    };

    std::array<double, 10> fv1 = {0}, fv2 = {0};
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc = f(centr);
    resk = fc * wgk[10];
    *resabs = fabs(resk);
    for (j = 0; j < 5; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * xgk[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j] * fsum;
        resk += wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 5; j++) {
        jtwm1 = j * 2;
        absc = hlgth * xgk[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = (*resabs) + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = wgk[10] * fabs(fc - reskh);
    for (j = 0; j < 10; j++)
        *resasc = (*resasc) + wgk[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = fmax(kEPMACH * 50.0 * (*resabs), (*abserr));
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

    static const std::array<double, 8> wg = {
            0.030753241996117268354628393577204,
            0.070366047488108124709267416450667,
            0.107159220467171935011869546685869,
            0.139570677926154314447804794511028,
            0.166269205816993933553200860481209,
            0.186161000015562211026800561866423,
            0.198431485327111576456118326443839,
            0.202578241925561272880620199967519
    };
    static const std::array<double, 16> wgk = {
            0.005377479872923348987792051430128,
            0.015007947329316122538374763075807,
            0.025460847326715320186874001019653,
            0.035346360791375846222037948478360,
            0.044589751324764876608227299373280,
            0.053481524690928087265343147239430,
            0.062009567800670640285139230960803,
            0.069854121318728258709520077099147,
            0.076849680757720378894432777482659,
            0.083080502823133021038289247286104,
            0.088564443056211770647275443693774,
            0.093126598170825321225486872747346,
            0.096642726983623678505179907627589,
            0.099173598721791959332393173484603,
            0.100769845523875595044946662617570,
            0.101330007014791549017374792767493
    };
    static const std::array<double, 16> xgk = {
            0.998002298693397060285172840152271,
            0.987992518020485428489565718586613,
            0.967739075679139134257347978784337,
            0.937273392400705904307758947710209,
            0.897264532344081900882509656454496,
            0.848206583410427216200648320774217,
            0.790418501442465932967649294817947,
            0.724417731360170047416186054613938,
            0.650996741297416970533735895313275,
            0.570972172608538847537226737253911,
            0.485081863640239680693655740232351,
            0.394151347077563369897207370981045,
            0.299180007153168812166780024266389,
            0.201194093997434522300628303394596,
            0.101142066918717499027074231447392,
            0.000000000000000000000000000000000
    };

    std::array<double, 15> fv1 = {0}, fv2 = {0};
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc = f(centr);
    resg = fc * wg[7];
    resk = fc * wgk[15];
    *resabs = fabs(resk);
    for (j = 0; j < 7; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * xgk[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j] * fsum;
        resk += wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 8; j++) {
        jtwm1 = j * 2;
        absc = hlgth * xgk[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = (*resabs) + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = wgk[15] * fabs(fc - reskh);
    for (j = 0; j < 15; j++)
        *resasc = (*resasc) + wgk[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = fmax(kEPMACH * 50.0 * (*resabs), (*abserr));
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

    static const std::array<double, 11> wg = {
            0.017614007139152118311861962351853,
            0.040601429800386941331039952274932,
            0.062672048334109063569506535187042,
            0.083276741576704748724758143222046,
            0.101930119817240435036750135480350,
            0.118194531961518417312377377711382,
            0.131688638449176626898494499748163,
            0.142096109318382051329298325067165,
            0.149172986472603746787828737001969,
            0.152753387130725850698084331955098
    };
    static const std::array<double, 21> wgk = {
            0.003073583718520531501218293246031,
            0.008600269855642942198661787950102,
            0.014626169256971252983787960308868,
            0.020388373461266523598010231432755,
            0.025882133604951158834505067096153,
            0.031287306777032798958543119323801,
            0.036600169758200798030557240707211,
            0.041668873327973686263788305936895,
            0.046434821867497674720231880926108,
            0.050944573923728691932707670050345,
            0.055195105348285994744832372419777,
            0.059111400880639572374967220648594,
            0.062653237554781168025870122174255,
            0.065834597133618422111563556969398,
            0.068648672928521619345623411885368,
            0.071054423553444068305790361723210,
            0.073030690332786667495189417658913,
            0.074582875400499188986581418362488,
            0.075704497684556674659542775376617,
            0.076377867672080736705502835038061,
            0.076600711917999656445049901530102
    };
    static const std::array<double, 21> xgk = {
            0.998859031588277663838315576545863,
            0.993128599185094924786122388471320,
            0.981507877450250259193342994720217,
            0.963971927277913791267666131197277,
            0.940822633831754753519982722212443,
            0.912234428251325905867752441203298,
            0.878276811252281976077442995113078,
            0.839116971822218823394529061701521,
            0.795041428837551198350638833272788,
            0.746331906460150792614305070355642,
            0.693237656334751384805490711845932,
            0.636053680726515025452836696226286,
            0.575140446819710315342946036586425,
            0.510867001950827098004364050955251,
            0.443593175238725103199992213492640,
            0.373706088715419560672548177024927,
            0.301627868114913004320555356858592,
            0.227785851141645078080496195368575,
            0.152605465240922675505220241022678,
            0.076526521133497333754640409398838,
            0.000000000000000000000000000000000
    };

    std::array<double, 20> fv1 = {0}, fv2 = {0};
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc = f(centr);
    resk = fc * wgk[20];
    *resabs = fabs(resk);
    for (j = 0; j < 10; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * xgk[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j] * fsum;
        resk += wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 10; j++) {
        jtwm1 = j * 2;
        absc = hlgth * xgk[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = (*resabs) + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = wgk[20] * fabs(fc - reskh);
    for (j = 0; j < 20; j++)
        *resasc = (*resasc) + wgk[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = fmax(kEPMACH * 50.0 * (*resabs), (*abserr));
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

    static const std::array<double, 13> wg = {
            0.011393798501026287947902964113235,
            0.026354986615032137261901815295299,
            0.040939156701306312655623487711646,
            0.054904695975835191925936891540473,
            0.068038333812356917207187185656708,
            0.080140700335001018013234959669111,
            0.091028261982963649811497220702892,
            0.100535949067050644202206890392686,
            0.108519624474263653116093957050117,
            0.114858259145711648339325545869556,
            0.119455763535784772228178126512901,
            0.122242442990310041688959518945852,
            0.123176053726715451203902873079050
    };
    static const std::array<double, 26> wgk = {
            0.001987383892330315926507851882843,
            0.005561932135356713758040236901066,
            0.009473973386174151607207710523655,
            0.013236229195571674813656405846976,
            0.016847817709128298231516667536336,
            0.020435371145882835456568292235939,
            0.024009945606953216220092489164881,
            0.027475317587851737802948455517811,
            0.030792300167387488891109020215229,
            0.034002130274329337836748795229551,
            0.037116271483415543560330625367620,
            0.040083825504032382074839284467076,
            0.042872845020170049476895792439495,
            0.045502913049921788909870584752660,
            0.047982537138836713906392255756915,
            0.050277679080715671963325259433440,
            0.052362885806407475864366712137873,
            0.054251129888545490144543370459876,
            0.055950811220412317308240686382747,
            0.057437116361567832853582693939506,
            0.058689680022394207961974175856788,
            0.059720340324174059979099291932562,
            0.060539455376045862945360267517565,
            0.061128509717053048305859030416293,
            0.061471189871425316661544131965264,
            // note: wgk[25] was calculated from the values of wgk(1..25)
            0.061580818067832935078759824240066
    };
    static const std::array<double, 26> xgk = {
            0.999262104992609834193457486540341,
            0.995556969790498097908784946893902,
            0.988035794534077247637331014577406,
            0.976663921459517511498315386479594,
            0.961614986425842512418130033660167,
            0.942974571228974339414011169658471,
            0.920747115281701561746346084546331,
            0.894991997878275368851042006782805,
            0.865847065293275595448996969588340,
            0.833442628760834001421021108693570,
            0.797873797998500059410410904994307,
            0.759259263037357630577282865204361,
            0.717766406813084388186654079773298,
            0.673566368473468364485120633247622,
            0.626810099010317412788122681624518,
            0.577662930241222967723689841612654,
            0.526325284334719182599623778158010,
            0.473002731445714960522182115009192,
            0.417885382193037748851814394594572,
            0.361172305809387837735821730127641,
            0.303089538931107830167478909980339,
            0.243866883720988432045190362797452,
            0.183718939421048892015969888759528,
            0.122864692610710396387359818808037,
            0.061544483005685078886546392366797,
            0.000000000000000000000000000000000
    };

    std::array<double, 25> fv1 = {0}, fv2 = {0};
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc = f(centr);
    resg = fc * wg[12];
    resk = fc * wgk[25];
    *resabs = fabs(resk);
    for (j = 0; j < 12; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * xgk[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j] * fsum;
        resk += wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 13; j++) {
        jtwm1 = j * 2;
        absc = hlgth * xgk[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = (*resabs) + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = wgk[25] * fabs(fc - reskh);
    for (j = 0; j < 25; j++)
        *resasc = (*resasc) + wgk[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = fmax(kEPMACH * 50.0 * (*resabs), (*abserr));
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

    static const std::array<double, 15> wg = {
            0.007968192496166605615465883474674,
            0.018466468311090959142302131912047,
            0.028784707883323369349719179611292,
            0.038799192569627049596801936446348,
            0.048402672830594052902938140422808,
            0.057493156217619066481721689402056,
            0.065974229882180495128128515115962,
            0.073755974737705206268243850022191,
            0.080755895229420215354694938460530,
            0.086899787201082979802387530715126,
            0.092122522237786128717632707087619,
            0.096368737174644259639468626351810,
            0.099593420586795267062780282103569,
            0.101762389748405504596428952168554,
            0.102852652893558840341285636705415
    };
    static const std::array<double, 31> wgk = {
            0.001389013698677007624551591226760,
            0.003890461127099884051267201844516,
            0.006630703915931292173319826369750,
            0.009273279659517763428441146892024,
            0.011823015253496341742232898853251,
            0.014369729507045804812451432443580,
            0.016920889189053272627572289420322,
            0.019414141193942381173408951050128,
            0.021828035821609192297167485738339,
            0.024191162078080601365686370725232,
            0.026509954882333101610601709335075,
            0.028754048765041292843978785354334,
            0.030907257562387762472884252943092,
            0.032981447057483726031814191016854,
            0.034979338028060024137499670731468,
            0.036882364651821229223911065617136,
            0.038678945624727592950348651532281,
            0.040374538951535959111995279752468,
            0.041969810215164246147147541285970,
            0.043452539701356069316831728117073,
            0.044814800133162663192355551616723,
            0.046059238271006988116271735559374,
            0.047185546569299153945261478181099,
            0.048185861757087129140779492298305,
            0.049055434555029778887528165367238,
            0.049795683427074206357811569379942,
            0.050405921402782346840893085653585,
            0.050881795898749606492297473049805,
            0.051221547849258772170656282604944,
            0.051426128537459025933862879215781,
            0.051494729429451567558340433647099
    };
    static const std::array<double, 31> xgk = {
            0.999484410050490637571325895705811,
            0.996893484074649540271630050918695,
            0.991630996870404594858628366109486,
            0.983668123279747209970032581605663,
            0.973116322501126268374693868423707,
            0.960021864968307512216871025581798,
            0.944374444748559979415831324037439,
            0.926200047429274325879324277080474,
            0.905573307699907798546522558925958,
            0.882560535792052681543116462530226,
            0.857205233546061098958658510658944,
            0.829565762382768397442898119732502,
            0.799727835821839083013668942322683,
            0.767777432104826194917977340974503,
            0.733790062453226804726171131369528,
            0.697850494793315796932292388026640,
            0.660061064126626961370053668149271,
            0.620526182989242861140477556431189,
            0.579345235826361691756024932172540,
            0.536624148142019899264169793311073,
            0.492480467861778574993693061207709,
            0.447033769538089176780609900322854,
            0.400401254830394392535476211542661,
            0.352704725530878113471037207089374,
            0.304073202273625077372677107199257,
            0.254636926167889846439805129817805,
            0.204525116682309891438957671002025,
            0.153869913608583546963794672743256,
            0.102806937966737030147096751318001,
            0.051471842555317695833025213166723,
            0.000000000000000000000000000000000
    };

    std::array<double, 30> fv1 = {0}, fv2 = {0};
    double absc, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc = f(centr);
    resk = fc * wgk[30];
    *resabs = fabs(resk);
    for (j = 0; j < 15; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * xgk[jtw];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j] * fsum;
        resk += wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 15; j++) {
        jtwm1 = j * 2;
        absc = hlgth * xgk[jtwm1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk += wgk[jtwm1] * fsum;
        *resabs = *resabs + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = wgk[30] * fabs(fc - reskh);
    for (j = 0; j < 30; j++)
        *resasc = (*resasc) + wgk[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = fmax(kEPMACH * 50.0 * (*resabs), (*abserr));
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

    static const std::array<double, 8> wg = {
            0.0,
            0.129484966168869693270611432679082,
            0.0,
            0.279705391489276667901467771423780,
            0.0,
            0.381830050505118944950369775488975,
            0.0,
            0.417959183673469387755102040816327
    };
    static const std::array<double, 8> wgk = {
            0.022935322010529224963732008058970,
            0.063092092629978553290700663189204,
            0.104790010322250183839876322541518,
            0.140653259715525918745189590510238,
            0.169004726639267902826583426598550,
            0.190350578064785409913256402421014,
            0.204432940075298892414161999234649,
            0.209482141084727828012999174891714
    };
    static const std::array<double, 8> xgk = {
            0.991455371120812639206854697526329,
            0.949107912342758524526189684047851,
            0.864864423359769072789712788640926,
            0.741531185599394439863864773280788,
            0.586087235467691130294144838258730,
            0.405845151377397166906606412076961,
            0.207784955007898467600689403773245,
            0.000000000000000000000000000000000
    };
    static const double epmach = std::numeric_limits<double>::epsilon();
    static const double uflow = std::numeric_limits<double>::min();

    std::array<double, 7> fv1 = {0}, fv2 = {0};

    const int dinf = std::min(1, inf);

    double centr = 0.5 * (a + b);
    double hlgth = 0.5 * (b - a);
    double tabsc1 = boun + dinf * (1.0 - centr) / centr;
    double fval1 = f(tabsc1);
    if (inf == 2) fval1 = fval1 + f(-tabsc1);
    double fc = (fval1 / centr) / centr;

    double resg = wg[7] * fc;
    double resk = wgk[7] * fc;
    *resabs = std::abs(resk);
    for (int j = 0; j < 7; j++) {
        double absc = hlgth * xgk[j];
        double absc1 = centr - absc;
        double absc2 = centr + absc;
        tabsc1 = boun + dinf * (1.0 - absc1) / absc1;
        double tabsc2 = boun + dinf * (1.0 - absc2) / absc2;
        fval1 = f(tabsc1);
        double fval2 = f(tabsc2);
        if (inf == 2) fval1 = fval1 + f(-tabsc1);
        if (inf == 2) fval2 = fval2 + f(-tabsc2);
        fval1 = (fval1 / absc1) / absc1;
        fval2 = (fval2 / absc2) / absc2;
        fv1[j] = fval1;
        fv2[j] = fval2;
        double fsum = fval1 + fval2;
        resg += wg[j] * fsum;
        resk += wgk[j] * fsum;
        *resabs += wgk[j] * (std::abs(fval1) + std::abs(fval2));
    }
    double reskh = resk * 0.5;
    *resasc = wgk[7] * std::abs(fc - reskh);

    for (int j = 0; j < 7; j++) {
        *resasc += wgk[j] * (std::abs(fv1[j] - reskh) + std::abs(fv2[j] - reskh));
    }

    double result = resk * hlgth;
    *resasc *= hlgth;
    *resabs *= hlgth;
    *abserr = std::abs((resk - resg) * hlgth);

    if (*resasc != 0.0 && *abserr != 0.0) *abserr = *resasc * std::min(1.0, pow(200.0 * (*abserr) / (*resasc), 1.5));
    if (*resabs > uflow / (50.0 * epmach)) *abserr = std::max((epmach * 50.0) * (*resabs), (*abserr));

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

    static const std::array<double, 4> wg = {
            0.1294849661688697e+00,
            0.2797053914892767e+00,
            0.3818300505051889e+00,
            0.4179591836734694e+00
    };
    static const std::array<double, 8> wgk = {
            0.2293532201052922e-01,
            0.6309209262997855e-01,
            0.1047900103222502e+00,
            0.1406532597155259e+00,
            0.1690047266392679e+00,
            0.1903505780647854e+00,
            0.2044329400752989e+00,
            0.2094821410847278e+00
    };
    static const std::array<double, 8> xgk = {
            0.9914553711208126e+00,
            0.9491079123427585e+00,
            0.8648644233597691e+00,
            0.7415311855993944e+00,
            0.5860872354676911e+00,
            0.4058451513773972e+00,
            0.2077849550078985e+00,
            0.0000000000000000e+00
    };

    std::array<double, 7> fv1 = {0}, fv2 = {0};
    double absc, absc1, absc2, centr, dhlgth;
    double fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh, result;
    int j, jtw, jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc = f(centr) * w(centr, p1, p2, p3, p4, kp);
    resg = fc * wg[3];
    resk = fc * wgk[7];
    *resabs = fabs(resk);
    for (j = 0; j < 3; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * xgk[jtw];
        absc1 = centr - absc;
        absc2 = centr + absc;
        fval1 = f(absc1) * w(absc1, p1, p2, p3, p4, kp);
        fval2 = f(absc2) * w(absc2, p1, p2, p3, p4, kp);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j] * fsum;
        resk += wgk[jtw] * fsum;
        *resabs = *resabs + wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 4; j++) {
        jtwm1 = j * 2;
        absc = hlgth * xgk[jtwm1];
        absc1 = centr - absc;
        absc2 = centr + absc;
        fval1 = f(absc1) * w(absc1, p1, p2, p3, p4, kp);
        fval2 = f(absc2) * w(absc2, p1, p2, p3, p4, kp);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + wgk[jtwm1] * fsum;
        *resabs = (*resabs) + wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = wgk[7] * fabs(fc - reskh);
    for (j = 0; j < 7; j++)
        *resasc = (*resasc) + wgk[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * fmin(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
    if (*resabs > kUFLOW / (50.0 * kEPMACH))
        *abserr = fmax(kEPMACH * 50.0 * (*resabs), (*abserr));
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

    static const std::array<double, 11> x = {
            0.991444861373810411144557526928563,
            0.965925826289068286749743199728897,
            0.923879532511286756128183189396788,
            0.866025403784438646763723170752936,
            0.793353340291235164579776961501299,
            0.707106781186547524400844362104849,
            0.608761429008720639416097542898164,
            0.500000000000000000000000000000000,
            0.382683432365089771728459984030399,
            0.258819045102520762348898837624048,
            0.130526192220051591548406227895489
    };

    double ak22, amom0, amom1, amom2, cc, centr;
    std::array<double, 13> cheb12 = {0};
    std::array<double, 25> cheb24 = {0};
    std::array<double, 25> fval = {0};
    double hlgth, resabs, resasc, res12, res24, u, result;
    int i, isym, k;
    int unitialized_value = 0xCCCCCCCC;
    int kp = unitialized_value;
    double p2 = unitialized_value, p3 = unitialized_value, p4 = unitialized_value;

    cc = (2.0 * c - b - a) / (b - a);
    if (fabs(cc) < 1.1) goto _10;

    /*  Apply the 15-point Gauss-Kronrod scheme.    */
    (*krul)--;
    result = qk15w(f, dqwgtc, c, p2, p3, p4, kp, a, b, abserr, &resabs, &resasc);
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
    dqcheb(x, fval, cheb12, cheb24);

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
double qc25f(
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

    static const std::array<double, 11> x = {
            0.991444861373810411144557526928563,
            0.965925826289068286749743199728897,
            0.923879532511286756128183189396788,
            0.866025403784438646763723170752936,
            0.793353340291235164579776961501299,
            0.707106781186547524400844362104849,
            0.608761429008720639416097542898164,
            0.500000000000000000000000000000000,
            0.382683432365089771728459984030399,
            0.258819045102520762348898837624048,
            0.130526192220051591548406227895489
    };

    double ac, an, an2, as, asap, ass, centr, conc, cons, cospar;
    double estc, ests, hlgth, parint, par2, par22;
    double resc12, resc24, ress12, ress24, result, sinpar;
    std::array<double, 13> cheb12 = {0};
    std::array<double, 25> cheb24 = {0};
    std::array<double, 28> d = {0};
    std::array<double, 28> d1 = {0};
    std::array<double, 28> d2 = {0};
    std::array<double, 28> d3 = {0};
    std::array<double, 25> fval = {0};
    std::array<double, 28> v = {0};
    int unitialized_value = 0xCCCCCCCC;
    double p2 = unitialized_value, p3 = unitialized_value, p4 = unitialized_value;
    int i, isym, j, k, m, noequ, noeq1, mm1;

    centr = 0.5 * (b + a);
    hlgth = 0.5 * (b - a);
    parint = omega * hlgth;

    /* Compute the integral using the 15-point Gauss-Kronrod formula
     * if the value of the parameter in the integrand is small or
     * is less than (bb-aa)/2^(maxp1-2), where (aa,bb) is the original
     * integration interval.
     */
    if (fabs(parint) > 2.0) goto _10;
    result = qk15w(f, dqwgto, omega, p2, p3, p4, sincos, a, b,
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

    dqcheb(x, fval, cheb12, cheb24);

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
    if (sincos == kSINE)
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
        double ri[],
        double rj[],
        double rg[],
        double rh[],
        double *abserr,
        double *resasc,
        int wgtfunc,
        int *nev
) {

    static const std::array<double, 11> x = {
            0.991444861373810411144557526928563,
            0.965925826289068286749743199728897,
            0.923879532511286756128183189396788,
            0.866025403784438646763723170752936,
            0.793353340291235164579776961501299,
            0.707106781186547524400844362104849,
            0.608761429008720639416097542898164,
            0.500000000000000000000000000000000,
            0.382683432365089771728459984030399,
            0.258819045102520762348898837624048,
            0.130526192220051591548406227895489
    };

    double centr, dc, factor, fix, hlgth, resabs, res12, res24, u, result;
    std::array<double, 13> cheb12 = {0};
    std::array<double, 25> cheb24 = {0};
    std::array<double, 25> fval = {0};
    int i, isym;

    *nev = 25;
    if ((bl == a) && ((alfa != 0.0) || (wgtfunc == 2) || (wgtfunc == 4)))
        goto _10;
    if ((br == b) && ((beta != 0.0) || (wgtfunc == 3) || (wgtfunc == 4)))
        goto _140;

    /*  If a>bl and b<br, apply the 15-point Gauss-Kronrod scheme. */
    result = qk15w(f, dqwgts, a, b, alfa, beta, wgtfunc, bl, br, abserr,
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
    dqcheb(x, fval, cheb12, cheb24);

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
    dqcheb(x, fval, cheb12, cheb24);

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
    dqcheb(x, fval, cheb12, cheb24);
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
    dqcheb(x, fval, cheb12, cheb24);

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
