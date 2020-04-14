//
// Created by Logan Morrison on 3/19/20.
//

#ifndef LANRE_DM_MODELS_DARKSUN_HPP
#define LANRE_DM_MODELS_DARKSUN_HPP

#include "lanre/dm_models/thermally_decoupled_model.hpp"
#include "lanre/cosmology/thermodynamic_particle.hpp"
#include "lanre/constants.hpp"
#include "lanre/special_functions/besselk.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/diffeq/rodas.hpp"
#include "lanre/diffeq/radau.hpp"
#include "lanre/diffeq/integrator.hpp"
#include "lanre/diffeq/problem.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <string>
#include <memory>
#include <utility>

namespace lanre {
namespace dm_models {

using namespace lanre::diffeq;

static constexpr double _log_xmin = 0.602059991327962; // log10(4)
static constexpr double _log_xmax = 2.69897000433602; // log10(500)
static constexpr double _log_xstep = 0.006989700043360188;
static constexpr double _eta_cs_slope = 14.02653358183077;
static constexpr double _eta_cs_intercept = -6.1933023533564535;

static const std::vector<double> _log_eta_cs{
        -16.637655236501224, -3.2073778028451474, -2.064304369976662, -1.3546809466722889, -0.8281658359797516,
        -0.40253938967491654, -0.039502564053162006, 0.28364160682026235,
        0.5772239884528931, 0.838111266851666, 1.0868009202155158, 1.3165324941873415, 1.5326218146462145,
        1.7425200505064233, 1.9370828345768174, 2.125085632208347, 2.3069051266553204,
        2.478546038197779, 2.65639789484337, 2.8203069536323957, 2.978730419114667, 3.1338320589598774,
        3.288998184793619, 3.4376303031062525, 3.58417626478382, 3.734959349971237,
        3.8721140443084954, 4.009342444091389, 4.144021960780492, 4.282415582197282, 4.417730984233178,
        4.54843204980097, 4.679145110145505, 4.8046999385083815, 4.929874072565337,
        5.0572858193187376, 5.183856715323707, 5.302187479857232, 5.425233252841499, 5.545537924432401,
        5.665829743176893, 5.786232172779206, 5.905022225954679, 6.0258369227675415,
        6.138636508498319, 6.255698645622665, 6.368823971955043, 6.483252788007763, 6.59763016077933, 6.712631284288037,
        6.821808987572604, 6.930747854698604, 7.042786867284195,
        7.1569283794315215, 7.262409238073819, 7.373585907964206, 7.482135041196106, 7.59007128705397,
        7.7007077186924295, 7.809759218549317, 7.9138204216959975, 8.023202918358328,
        8.127074833291443, 8.233259515944884, 8.342830412649826, 8.450962847401101, 8.552629916299894,
        8.656124052323225, 8.763687606605245, 8.866273209021838, 8.967552235363398,
        9.075691478656855, 9.183898759360103, 9.286054607208035, 9.386649353128306, 9.493095959513699, 9.59256817840535,
        9.695189855326634, 9.798081193077408, 9.904232549402886,
        10.00775435077294, 10.104328978940778, 10.209250799915372, 10.309536048894504, 10.415619200278199,
        10.512914630611792, 10.616925067998318, 10.719545840790364, 10.82569052636717,
        10.920026204996434, 11.019123642207628, 11.123726243878718, 11.223396862538673, 11.322457925709962,
        11.423608717091355, 11.52137264956311, 11.626939696768193, 11.724363636568516,
        11.829078290281942, 11.926756834670329, 12.023866977410192, 12.128754637861284, 12.229005579173206,
        12.327150549599835, 12.42320073379204, 12.527505575937475, 12.624263995550486,
        12.728288529421276, 12.831179430890277, 12.929373577091864, 13.019296117317667, 13.12576615740427,
        13.228316221379316, 13.326752271746773, 13.421767475692928, 13.528376240165276,
        13.621213484612975, 13.720237321024324, 13.820223741511278, 13.916508156762665, 14.014710829324216,
        14.114795925887437, 14.214748909756986, 14.317607141627759, 14.411039427599727,
        14.512167399083273, 14.610786681259539, 14.703982090462606, 14.805431627124811, 14.908040642743243,
        15.002575847531137, 15.101189984138106, 15.198367008875683, 15.29764620343702,
        15.396719987595997, 15.496524643113691, 15.596348694738474, 15.695865032147397, 15.798223622015554,
        15.896844304936113, 15.994636909250309, 16.089005357220906, 16.187523173262452,
        16.286877065339063, 16.38136510301903, 16.48611312927631, 16.580694250205873, 16.676140007828867,
        16.777566155452426, 16.86851472623604, 16.972276170062106, 17.07236604676929,
        17.17002998393461, 17.263518742143038, 17.371960381122054, 17.46679324741868, 17.566527606164218,
        17.657190108223872, 17.76294506487139, 17.85800143424382, 17.95533926398797,
        18.054084208438717, 18.152947633230575, 18.24784882593423, 18.348745139145336, 18.447309008046254,
        18.544926914067986, 18.646181353693144, 18.74249797154525, 18.838536968682273,
        18.938410159329106, 19.039708295100098, 19.12963606837031, 19.23961832162645, 19.328789120506478,
        19.428831868588325, 19.524135208177725, 19.623132432405413, 19.725036388007293,
        19.822908851049924, 19.916681049410915, 20.017173948675897, 20.117788115518543, 20.212385261481966,
        20.309657314200834, 20.410678718805556, 20.50591946714346, 20.6039071615597,
        20.710023214016594, 20.795935078224993, 20.892445919823906, 20.997837703659616, 21.091883072100714,
        21.189986857759582, 21.28474246979953, 21.383595675156343, 21.489635305247493,
        21.58658072126945, 21.67923306964906, 21.777151543542594, 21.877503807332698, 21.9777124345802,
        22.074420570417356, 22.171924206337504, 22.26412444140273, 22.370400437296468,
        22.46440200316392, 22.564445638277135, 22.661225884857526, 22.759392829896036, 22.855759327870224,
        22.961025126269117, 23.04967330445861, 23.15086367387428, 23.24693124118719,
        23.342482534019904, 23.440903850791386, 23.54353641327109, 23.639265951519356, 23.735095878576605,
        23.83628836843447, 23.937952900572512, 24.029500321184912, 24.135181019639223,
        24.232097520335735, 24.326548923665776, 24.420821030333116, 24.524166950199877, 24.62264131885209,
        24.720467087699234, 24.8079699036987, 24.915334091759682, 25.006039372376208,
        25.11092926075667, 25.206329650925014, 25.303373364546196, 25.40449321448272, 25.498342521804513,
        25.59578509653883, 25.691699415485985, 25.792891227265198, 25.889764286973644,
        25.989357277209372, 26.09162357758695, 26.18373528480038, 26.281826407237375, 26.37521869912794,
        26.47482057242853, 26.580942063470875, 26.678337201255403, 26.776286672356854,
        26.864896626788482, 26.965377322197643, 27.069142380346584, 27.163463331840834, 27.257275123072255,
        27.35603533182855, 27.454050023145445, 27.55381560630426, 27.652028171736465,
        27.753068180356422, 27.851517014690277, 27.941954988568334, 28.04625825672801, 28.143590313973302,
        28.239825994178066, 28.33633209342876, 28.431863693763457, 28.53321745058022,
        28.63225613663344, 28.724961183491533, 28.828884530768693, 28.924708511933343, 29.025410664836738,
        29.121282759809816, 29.220270025971008, 29.310135957157353, 29.413877402170062,
        29.510914274473848, 29.616368558872022, 29.706650112894227, 29.80755514838043, 29.90229574136799,
        30.00459742503278, 30.098843898329033, 30.194621572805133, 30.296481970402994,
        30.393556467012353, 30.49110980213624, 30.586373509120993, 30.685775036441395, 30.78946619562423,
        30.883725941744665, 30.97720911177579, 31.07983605573461, 31.176156653596653,
        31.27896288818872, 31.368149662172968, 31.469379909589453, 31.568463567862715, 31.665013490312656
};

constexpr int gauss_laguerre_size = 200;

constexpr double gauss_laguerre_nodes[] = {
        0.007210969203825846, 0.0379946533149587, 0.09337884753136749, 0.17337925254703734, 0.27800228779091074,
        0.4072547401404738, 0.5611446730468262, 0.7396815982422049, 0.9428765238609318, 1.1707419731003792,
        1.4232919942768794, 1.7005421681364932, 2.0025096143712564, 2.3292129980799605, 2.680672536484583,
        3.056910006046025, 3.4579487500517723, 3.8838136867158406, 4.334531317815899, 4.810129737885178,
        5.31063864397255, 5.836089345982861, 6.386514777607974, 6.961949507860055, 7.562429753216415, 8.187993390388847,
        8.838679969727261, 9.51453072927025, 10.215588609455258, 10.94189826850071, 11.693506098474709,
        12.470460242063924, 13.272810610058054, 14.100608899565815, 14.953908612978505, 15.832765077699271,
        16.737235466655537, 17.66737881961392, 18.62325606531771, 19.604930044467626, 20.61246553356731,
        21.64592926965728, 22.705389975960202, 23.790918388463044, 24.90258728346199, 26.0404715060976,
        27.20464799990846, 28.395195837432958, 29.612196251891337, 30.855732669979048, 32.12589074580691,
        33.422758396022445, 34.74642583615081, 36.09698561819301, 37.47453266952271, 38.8791643331242,
        40.31098040921585, 41.77008319830519, 43.25657754572524, 44.770570887701496, 46.31217329900461,
        47.881497542242975, 49.47865911885442, 51.103776321858085, 52.75697029043024, 54.43836506637088,
        56.148087652532084, 57.886268073281144, 59.65303943707597, 61.44853800123349, 63.27290323897623,
        65.1262779088458, 67.00880812657769, 68.92064343953433, 70.86193690379984, 72.8328451640466, 74.83352853628409,
        76.8641510936134, 78.92488075511012, 81.0158893779706, 83.13735285305938, 85.28945120400508, 87.47236868999957,
        89.68629391246309, 91.93141992574655, 94.20794435205381, 96.51606950077273, 98.85600249241861,
        101.22795538740122, 103.63214531984268, 106.06879463668189, 108.53813104231789, 111.04038774905847,
        113.57580363365481, 116.14462340022166, 118.74709774985809, 121.38348355730339, 124.05404405498524,
        126.75904902483619, 129.49877499827804, 132.2735054648016, 135.08353108959335, 137.929149940689,
        140.8106677261691, 143.72839804193876, 146.682662630675, 149.67379165256298, 152.70212396847847,
        155.76800743632836, 158.87179922129988, 162.01386612082908, 165.19458490515078, 168.41434267435565,
        171.6735372329475, 174.97257748296127, 178.31188383678787, 181.6918886509284, 185.1130366819985,
        188.57578556640217, 192.08060632520176, 195.62798389583412, 199.21841769244773, 202.85242219678307,
        206.53052758166962, 210.25328036938683, 214.02124412731897, 217.8350002035438, 221.69514850521543,
        225.60230832285188, 229.5571192039125, 233.5602418793479, 237.61235924714128, 241.71417741722613,
        245.86642682257587, 250.06986340170963, 254.3252698583655, 258.6334570046502, 262.9952651945983,
        267.4115658557703, 271.88326312730106, 276.4112956136785, 280.9966382645146, 285.6403043916765,
        290.3433478363775, 295.10686530023975, 299.9319988559197, 304.81993865468473, 309.7719258503762,
        314.7892557615144, 319.8732812959589, 325.0254166655744, 330.2471414218363, 335.54000484731984,
        340.90563074263855, 346.3457226537445, 351.8620695907119, 357.45655229633655, 363.13115013132364,
        368.88794865270114, 374.7291479737157, 380.6570720071561, 386.67417871028186, 392.7830714688308,
        398.9865117806249, 405.2874334269487, 411.6889583531987, 418.1944145206763, 424.8073560405489,
        431.53158596118175, 438.37118215415387, 445.33052683608435, 452.41434037790884, 459.62772019704227,
        466.9761857097865, 474.4657305533142, 482.10288358484604, 489.8947805527052, 497.84924884106726,
        505.9749083617108, 514.2812925652623, 522.7789947634429, 531.4798466286464, 540.3971380721143,
        549.5458910098196, 558.943204294893, 568.6086941091803, 578.5650646478555, 588.8388601672838, 599.4614751769695,
        610.470541577057, 621.9118827650128, 633.842350678966, 646.3340958894291, 659.4812831369854, 673.4112473253756,
        688.3043771556801, 704.4330522841468, 722.2487752087823, 742.6220869994196, 767.8146922967122
};

constexpr double gauss_laguerre_weights[] = {
        0.018372766795478238, 0.041472857188875674, 0.061656364810778734, 0.0776165403238259, 0.0885554492774702,
        0.09421106624235254, 0.09483302598034032, 0.09109285312127291, 0.08394939562242641, 0.07449647889595658,
        0.06381897014494235, 0.052877536993020155, 0.04243360820337018, 0.03301681769439271, 0.024929571603679004,
        0.018278549978639536, 0.013021240724436198, 0.009016530857203353, 0.006071019745434235, 0.0039760354756828535,
        0.0025334636457925784, 0.0015708920359582537, 0.0009480293423239927, 0.0005569357682463247,
        0.0003185297029969012, 0.00017737846986940766, 9.618225438294225e-5, 5.0788142691503066e-5,
        2.611724967923984e-5, 1.308005359734636e-5, 6.380037628796025e-6, 3.030952664415161e-6, 1.4024395053537512e-6,
        6.32035852452258e-7, 2.774293198971438e-7, 1.1860793112056644e-7, 4.938793958118153e-8, 2.002933614570185e-8,
        7.911176724229093e-9, 3.043216153641342e-9, 1.1400607679656006e-9, 4.159212350729819e-10,
        1.4776265274486515e-10, 5.1117502184990794e-11, 1.721888917364375e-11, 5.64740315088098e-12,
        1.8033359474191282e-12, 5.606130324364673e-13, 1.6966078499466524e-13, 4.998078431106501e-14,
        1.4331713767864812e-14, 3.999773462138043e-15, 1.0863811692285752e-15, 2.8714721566508464e-16,
        7.385283124650752e-17, 1.848132493637661e-17, 4.499501175665134e-18, 1.0656669466627782e-18,
        2.455068831719016e-19, 5.501089815605696e-20, 1.198762522748413e-20, 2.5402241476371243e-21,
        5.2338239521497e-22, 1.0483992560083534e-22, 2.0414829463456282e-23, 3.863899188602396e-24,
        7.107467526357688e-25, 1.2704561733115303e-25, 2.206507310255883e-26, 3.723022122623214e-27,
        6.102005590035238e-28, 9.713527292266692e-29, 1.501577632605412e-29, 2.2538373730605852e-30,
        3.284257546440432e-31, 4.6454300389295945e-32, 6.377071651806957e-33, 8.49482011143389e-34,
        1.0978770691045016e-34, 1.3764081112421736e-35, 1.6736330452190763e-36, 1.9734060043041126e-37,
        2.2559922872802142e-38, 2.5000205428997027e-39, 2.685048906431027e-40, 2.794340361137083e-41,
        2.817335904359463e-42, 2.7513243422471602e-43, 2.6019455939159634e-44, 2.3824021743045994e-45,
        2.1115300771774496e-46, 1.8111203978959266e-47, 1.5030219658608392e-48, 1.2065599369550572e-49,
        9.366839834329689e-51, 7.030566659474735e-52, 5.100706398618385e-53, 3.576020614934395e-54,
        2.422051990335041e-55, 1.5843888761825218e-56, 1.0007211076056123e-57, 6.101157426399236e-59,
        3.5894756384551313e-60, 2.0372192821198905e-61, 1.1150580501722755e-62, 5.883984247057542e-64,
        2.9923857452824502e-65, 1.466187947873666e-66, 6.918910785982279e-68, 3.143459076419942e-69,
        1.3744902845828099e-70, 5.781982540891726e-72, 2.339085415945253e-73, 9.096575028334792e-75,
        3.39934894025185e-76, 1.2201637073340897e-77, 4.2049253465220886e-79, 1.3906704954181126e-80,
        4.411825248647524e-82, 1.3419520966734076e-83, 3.911748095585153e-85, 1.0922049220198121e-86,
        2.9195401858563654e-88, 7.467457479367728e-90, 1.826597423536408e-91, 4.270523610191112e-93,
        9.537539224539778e-95, 2.0335232840006943e-96, 4.1366792524166654e-98, 8.023575500096553e-100,
        1.4828996863391997e-101, 2.609684924402513e-103, 4.37010170097307e-105, 6.958351062348703e-107,
        1.0527023044487484e-108, 1.511996528113801e-110, 2.0601145506921448e-112, 2.6605038391586223e-114,
        3.253804831037912e-116, 3.7651578724207657e-118, 4.118442870206282e-120, 4.254214916005335e-122,
        4.145759447111929e-124, 3.8074245415083625e-126, 3.2917492714187464e-128, 2.676064113105895e-130,
        2.0432706643784e-132, 1.4634508855363796e-134, 9.819584246619774e-137, 6.1643271901145805e-139,
        3.615299949165377e-141, 1.9780135295694693e-143, 1.0080231678226168e-145, 4.7771179578343084e-148,
        2.101740132004247e-150, 8.569118526265068e-153, 3.2316492266214514e-155, 1.1250905145174907e-157,
        3.608479784613956e-160, 1.0638577506267038e-162, 2.876476060344005e-165, 7.11530379360982e-168,
        1.6060401125261061e-170, 3.298800898542803e-173, 6.147871014922091e-176, 1.0363662381700279e-178,
        1.575011268378154e-181, 2.1503067829809905e-184, 2.6273660647637203e-187, 2.8614268550662127e-190,
        2.7656285176601983e-193, 2.361106516022e-196, 1.7715225591626906e-199, 1.1617269853889192e-202,
        6.619120130428007e-206, 3.255495584911086e-209, 1.372381988302243e-212, 4.9203548564633624e-216,
        1.48751899411141e-219, 3.7563075021120645e-223, 7.840142016022902e-227, 1.3367484169799703e-230,
        1.8374261097003344e-234, 2.0060014900132215e-238, 1.7102031672114025e-242, 1.1166357636129291e-246,
        5.459448502775133e-251, 1.9468705545063559e-255, 4.908776409725221e-260, 8.431015556076061e-265,
        9.426927576899392e-270, 6.48688854581158e-275, 2.5581673170610254e-280, 5.268602581343383e-286,
        4.998469790084349e-292, 1.8292375750678093e-298, 1.976265083311595e-305, 4.0513560427e-313, 6.77e-322, 0.0
};

static const boost::math::cubic_b_spline<double> _eta_cs_spline{_log_eta_cs.begin(),
        _log_eta_cs.end(), _log_xmin, _log_xstep};


class DarkSUN : public ThermallyDecoupledModel {
protected:
    double m_lam; // Confinement scale
    unsigned int m_N; // N in SU(N)
    double m_L1; // Low-energy constant for 4-pt eta interaction
    double m_c; // Exponential supression factor for 2eta<->2delta, M~exp(-cN)
    double m_a; // Exponential supression factor for initial delta density
    double m_mu_eta; // Ratio of eta mass to Lambda / sqrt(N)
    double m_mu_delta; // Ratio of delta mass to N * Lambda
    bool m_has_dp; // Flag specifying if a dark photon is present in the model

    ThermodynamicParticle m_eta; // The dark eta'
    ThermodynamicParticle m_delta; // The dark Delta
    ThermodynamicParticle m_dark_photon; // The dark photon

    double m_xi_fo{-1};
    double m_Tsm_fo{-1};

    double m_xi_bbn{-1};
    double m_xi_cmb{-1};

    double m_rd_eta{};
    double m_rd_delta{};

    double m_dneff_cmb{};
    double m_dneff_bbn{};

    double m_eta_si_per_mass{};
    double m_delta_si_per_mass{};

    ODESolution m_solution{};

    static double scaled_cross_section_2eta_4eta(double);

    static double scaled_thermal_cross_section_2eta_4eta(double);

    friend class DarkSUNBoltzmann;

public:
    DarkSUN(double t_lam, unsigned int t_N, double t_L1, double t_c, double t_a,
            double t_mu_eta, double t_mu_delta, double t_xi_inf, bool t_has_dp)
            : m_lam(t_lam), m_N(t_N), m_L1(t_L1), m_c(t_c), m_a(t_a), m_mu_eta(t_mu_eta),
              m_mu_delta(t_mu_delta), m_has_dp(t_has_dp),
              m_eta{t_mu_eta * t_lam / sqrt(t_N), 1.0, 0},
              m_delta{t_mu_delta * t_lam * t_N, double(t_N + 1), t_N},
              m_dark_photon{0.0, t_has_dp ? 2.0 : 0.0, 2} {
        m_xi_inf = t_xi_inf;
        m_hd_inf = 7.0 / 2.0 * m_N + 2 * m_N * m_N - 2.0 + m_dark_photon.get_g();

        if (m_has_dp && m_dark_photon.get_mass() < m_eta.get_mass()) {
            m_gl = m_dark_photon.get_g();
            m_ml = m_dark_photon.get_mass();
        } else {
            m_gl = m_eta.get_g();
            m_ml = m_eta.get_mass();
        }
    }

    // Getters
    double get_lam() const { return m_lam; }

    unsigned int get_N() const { return m_N; }

    double get_L1() const { return m_L1; }

    double get_c() const { return m_c; }

    double get_a() const { return m_a; }

    double get_mu_eta() const { return m_mu_eta; }

    double get_mu_delta() const { return m_mu_delta; }

    double get_xi_inf() const { return m_xi_inf; }

    bool get_has_dp() const { return m_has_dp; }

    double get_xi_fo() const { return m_xi_fo; }

    double get_Tsm_fo() const { return m_Tsm_fo; }

    double get_xi_bbn() const { return m_xi_bbn; }

    double get_xi_cmb() const { return m_xi_cmb; }

    double get_rd_eta() const { return m_rd_eta; }

    double get_rd_delta() const { return m_rd_delta; }

    double get_dneff_cmb() const { return m_dneff_cmb; }

    double get_dneff_bbn() const { return m_dneff_bbn; }

    double get_eta_si_per_mass() const { return m_eta_si_per_mass; }

    double get_delta_si_per_mass() const { return m_delta_si_per_mass; }

    ODESolution get_solution() const { return m_solution; }

    //Setters
    void set_lam(double lam) {
        m_lam = lam;
        m_eta.set_mass(m_mu_eta * m_lam / std::sqrt(double(m_N)));
        m_delta.set_mass(m_mu_delta * m_lam * m_N);
    }

    void set_N(unsigned int N) {
        m_N = N;
        m_eta.set_mass(m_mu_eta * m_lam / std::sqrt(double(m_N)));
        m_delta.set_mass(m_mu_delta * m_lam * m_N);
        m_hd_inf = 7.0 / 2.0 * m_N + 2 * m_N * m_N - 2.0 + m_dark_photon.get_g();
    }

    void set_L1(double L1) { m_L1 = L1; }

    void set_c(double c) { m_c = c; }

    void set_a(double a) { m_a = a; }

    void set_mu_eta(double mu_eta) {
        m_mu_eta = mu_eta;
        m_eta.set_mass(m_mu_eta * m_lam / std::sqrt(double(m_N)));
    }

    void set_mu_delta(double mu_delta) {
        m_mu_delta = mu_delta;
        m_delta.set_mass(m_mu_delta * m_lam * m_N);
    }

    void set_xi_inf(double xi_inf) { m_xi_inf = xi_inf; }

    void set_has_dp(bool has_dp) { m_has_dp = has_dp; }

    void set_xi_fo(double xi_fo) { m_xi_fo = xi_fo; }

    void set_Tsm_fo(double Tsm_fo) { m_Tsm_fo = Tsm_fo; }

    void set_xi_bbn(double xi_bbn) { m_xi_bbn = xi_bbn; }

    void set_xi_cmb(double xi_cmb) { m_xi_cmb = xi_cmb; }

    void set_rd_eta(double val) { m_rd_eta = val; }

    void set_rd_delta(double val) { m_rd_delta = val; }

    void set_dneff_cmb(double val) { m_dneff_cmb = val; }

    void set_dneff_bbn(double val) { m_dneff_bbn = val; }

    void set_eta_si_per_mass(double val) { m_eta_si_per_mass = val; }

    void set_delta_si_per_mass(double val) { m_delta_si_per_mass = val; }

    void set_solution(ODESolution sol) { m_solution = std::move(sol); }

    double dark_heff(const double Td) const override {
        return m_eta.heff(Td) + m_delta.heff(Td) + m_dark_photon.heff(Td);
    }

    double dark_geff(const double Td) const {
        return m_eta.geff(Td) + m_delta.geff(Td) + m_dark_photon.geff(Td);
    }

    double sqrt_gstar(double Tsm, double xi) const {
        double gd = dark_geff(Tsm * xi);
        double gsm = sm_geff(Tsm);
        return sm_sqrt_gstar(Tsm) * sqrt(gsm / (gsm + gd * xi * xi * xi * xi));
    }

    double compute_xi(double Tsm) const {
        /* Conditions where we need to compute xi by root-solve:
         *  1) We have a dark photon
         *  2) Tsm is larger than Tsm_fo (might happen if solver re-samples
         *     past values? or if haven't frozen out, in which case m_Tsm_fo = -1)
         */
        if (m_dark_photon.get_g() > 0.0 || Tsm > m_Tsm_fo) {
            return compute_xi_const_tsm(Tsm);
        }
        /* If eta has frozen out, we compute xi using:
         *  1) if eta freezes out relativistically, xi remains constant,
         *  2) if eta freezes out non-relativistically, xi red-shifts.
         */
        if (m_Tsm_fo * m_xi_fo > m_eta.get_mass()) {
            return m_xi_fo;
        }
        return m_xi_fo * Tsm / m_Tsm_fo;
    }

    double thermal_cross_section_2eta_4eta(double x) const;

    double thermal_cross_section_2eta_2delta(double x) const;

    double cross_section_2eta_2eta() const;

    double cross_section_2delta_2delta() const;

    double delta_n_eff_cmb() const;

    double delta_n_eff_bbn() const;

    void solve_boltzmann(double, double, const std::string &);

};

struct DarkSUNBoltzmann : public ODEFunction {
    std::shared_ptr<DarkSUN> model;

    explicit DarkSUNBoltzmann(std::shared_ptr<DarkSUN> t_model) : model(std::move(t_model)) {}

    void dudt(Vector<double> &dw, const Vector<double> &w, double logx) override {
        double x = exp(logx);
        double meta = model->m_eta.get_mass();
        double Tsm = meta / x;
        double xi = model->compute_xi(Tsm);
        double Td = xi * Tsm;
        double s = sm_entropy_density(Tsm);

        double neq = model->m_eta.neq(Td);
        double weq = log(neq / s);

        // Determine if the eta' has frozen out
        if (w(0) - weq > 0.4 && model->m_xi_fo == -1.0) {
            model->m_xi_fo = xi;
            model->m_Tsm_fo = Tsm;
        }
        // save xi at BBN
        if (Tsm < kT_BBN && model->m_xi_bbn == -1.0) {
            model->m_xi_bbn = xi;
        }

        double sig_ee_eeee = model->thermal_cross_section_2eta_4eta(meta / Td);
        double sig_eeee_ee = sig_ee_eeee / neq / neq;
        sig_eeee_ee = std::isnan(sig_eeee_ee) ? 0.0 : sig_eeee_ee;

        double sig_ee_dd = model->thermal_cross_section_2eta_2delta(meta / Td);

        double pf_e = -sqrt(M_PI / 45) * kM_PLANK * model->sqrt_gstar(Tsm, xi) * Tsm * s * s;
        double pf_d = sqrt(M_PI / 45) * kM_PLANK * model->sqrt_gstar(Tsm, xi) * meta / (x * x);

        // dW_e / dlogx
        dw(0) = pf_e * sig_eeee_ee * exp(w(0)) * (exp(2.0 * w(0)) - exp(2.0 * weq));

        // dY_d / dlogx
        dw(1) = pf_d * sig_ee_dd * exp(2.0 * w(0));
    }

    void dfdu(Matrix<double> &J, const Vector<double> &w, double logx) override {
        double x = exp(logx);
        double meta = model->m_eta.get_mass();
        double Tsm = meta / x;
        double xi = model->compute_xi(Tsm);
        double Td = xi * Tsm;
        double s = sm_entropy_density(Tsm);

        double neq = model->m_eta.neq(Td);
        double weq = log(neq / s);

        double sig_ee_eeee = model->thermal_cross_section_2eta_4eta(meta / Td);
        double sig_eeee_ee = sig_ee_eeee / neq / neq;
        sig_eeee_ee = std::isnan(sig_eeee_ee) ? 0.0 : sig_eeee_ee;

        double sig_ee_dd = model->thermal_cross_section_2eta_2delta(meta / Td);

        double pf_e = -sqrt(M_PI / 45) * kM_PLANK * model->sqrt_gstar(Tsm, xi) * Tsm * s * s;
        double pf_d = sqrt(M_PI / 45) * kM_PLANK * model->sqrt_gstar(Tsm, xi) * meta / (x * x);

        // df_d / dW_e
        J(0, 0) = pf_e * sig_eeee_ee * exp(w(0)) * (3.0 * exp(2.0 * w(0)) - exp(2.0 * weq));

        // df_e / dY_d
        J(0, 1) = 0.0;

        // df_d / dW_e
        J(1, 0) = 2.0 * pf_d * sig_ee_dd * exp(2.0 * w(0));

        // df_d / dY_d
        J(1, 1) = 0.0;
    }
};


/**
 * Cross section for 4eta -> 2eta scaled by model-dependent constants.
 * @param z Center-of-mass scaled by eta' mass: z = cme / meta.
 * @return Scaled cross section.
 */
double DarkSUN::scaled_cross_section_2eta_4eta(double z) {
    double logz = log10(z);
    if (logz <= _log_xmin) {
        return 0.0;
    } else if (logz >= _log_xmax) {
        return pow(z, _eta_cs_slope) * pow(10.0, _eta_cs_intercept);
    }
    return pow(10.0, _eta_cs_spline(logz));
}

/**
 * Thermal cross section for 4eta -> 2eta scaled by model-dependent constants
 * and without a cutoff
 * @param x Mass of the eta' divided by the dark temperature.
 * @return Scaled thermal cross section.
 */
double DarkSUN::scaled_thermal_cross_section_2eta_4eta(double x) {
    using namespace boost::math;
    using namespace lanre::special_functions;
    static const double z_min = 4.0;

    double integral = 0.0;
    for (int i = 0; i < gauss_laguerre_size; i++) {
        double zp = gauss_laguerre_nodes[i] / x + z_min;
        double weight = gauss_laguerre_weights[i];
        integral += weight * zp * sqrt(zp * zp - 4.0) *
                besselk1e(x * zp) *
                scaled_cross_section_2eta_4eta(zp);
    }
    const double bessel2 = besselk2(x);
    const double pf = x / (8.0 * bessel2 * bessel2);
    return pf * integral * exp(-z_min * x) / x;
}


/**
 * Thermal cross section for 4eta -> 2eta without cut-off
 * @param x Mass of the eta' divided by the dark temperature.
 * @return Thermal cross section.
 */
double DarkSUN::thermal_cross_section_2eta_4eta(double x) const {
    using namespace boost::math;
    double meta = m_eta.get_mass();
    double feta = sqrt(double(m_N)) * m_lam / (4.0 * M_PI);
    double pf = m_L1 * m_L1 * pow<7>(meta) / pow<4>(feta * m_lam);
    pf *= pf;

    return pf * scaled_thermal_cross_section_2eta_4eta(x);
}

/**
 * Thermal cross section for 2eta -> 2delta.
 * @param x Mass of the eta' divided by the dark temperature.
 * @return
 */
double DarkSUN::thermal_cross_section_2eta_2delta(double x) const {
    using namespace boost::math;
    using namespace boost::math::quadrature;
    using namespace special_functions;
    /*
    double me = m_eta.get_mass();
    double md = m_delta.get_mass();
    double pf = exp(-2 * m_c * m_N) * m_N * x / (128 * pow<7>(me) * M_PI * pow<2>(besselk2(x)));
    double z_min = 2.0 * md / me;

    double integral = 0.0;
    for (int i = 0; i < gauss_laguerre_size; i++) {
        double zp = gauss_laguerre_nodes[i] / x + z_min;
        double weight = gauss_laguerre_weights[i];
        integral += weight * sqrt(zp * zp - 4.0) * pow(zp * zp * me * me - 4 * md * md, 3.0 / 2.0) *
                besselk1e(x * zp) / pow<2>(zp * zp - 1.0);
    }


    double result = pf * integral * exp(-z_min * x) / x;
    return isnan(result) ? 0.0 : result;
    */

    // The following seems more accurate. This amplitude comes from assuming
    // an effective interaction of deltabar.delta * del[eta,mu]^2
    double me = m_eta.get_mass();
    double md = m_delta.get_mass();
    double feta = (std::sqrt(m_N) * m_lam / (4.0 * M_PI));
    double pf_num = std::exp(-2.0 * m_c * m_N) * pow<4>(me) * x;
    double pf_den = 2.0 * M_PI * pow<2>(4.0 * feta * feta * m_lam * besselk2e(x));
    double pf = pf_num / pf_den;
    double r = md / me;

    auto integrand = [x, r](double z) {
        return exp(x * (2.0 - z)) * sqrt(-4 + z * z) * pow<2>(-2 + z * z) *
                pow(-4 * r * r + z * z, 1.5) * besselk1e(x * z);
    };

    double integral = gauss_kronrod<double, 15>::integrate(
            integrand, 2.0 * r, std::numeric_limits<double>::infinity(), 5, 1e-9
    );

    return integral * pf;
}

/**
 * Compute the zero-temperature self-interaction cross section for
 * eta + eta -> eta + eta.
 * @return Self-interaction cross-section 2eta->2eta.
 */
double DarkSUN::cross_section_2eta_2eta() const {
    using namespace boost::math;

    //compute dark temperature today in order to compute eta velocity
    double Tsm_today = 2.7255 * 8.6173303e-14;
    double xi_today = compute_xi(Tsm_today);
    double Td_today = xi_today * Tsm_today;

    double me = m_eta.get_mass();
    double vrel = sqrt(3.0 * Td_today / me);
    double s = 4 * me * me / (1 - vrel * vrel);

    double fe = m_lam * sqrt(m_N) / (4.0 * M_PI);
    return (8.0 * m_L1 * m_L1 * (376.0 * pow<8>(me) -
            576.0 * pow<6>(me) * s +
            396.0 * pow<4>(me) * s * s -
            136.0 * me * me * pow<3>(s) +
            21 * pow<4>(s))) / (15.0 * pow<4>(fe * m_lam));
}

/**
 * Compute the zero-temperature self-interaction cross section for
 * delta + delta -> delta + delta.
 * @return Self-interaction cross-section delta->delta.
 */
double DarkSUN::cross_section_2delta_2delta() const {
    /* diagrams with sigma exchange
     * ----------------------------
     * p1 -->--.-->-- p3       p1 -->--.-->-- p4
     *         .                       .
     *         .          +            .
     *         .                       .
     * p2 -->--.-->-- p4       p2 -->--.-->-- p3
     */
    using namespace boost::math;
    double g = 1.0; // TODO: not sure what to use for this
    double msigma = m_lam;
    double mdelta = m_delta.get_mass();
    // NOTE: I removed extra factor of N:
    //  Hiren gave a good argument for why there should be no factor of N
    //  from spin sums. The spin-flipping pieces of the amplitude will go
    //  like S.p where S is the spin and p is some momentum scale. At zero
    //  velocity, this term should diappear.
    return 1.5 * pow<4>(g) * mdelta * mdelta / (32.0 * M_PI) / pow<4>(msigma);
}

/**
 * Compute the contribution to delta N_eff from the dark SU(N) sector at CMB.
 * @return Delta N_eff at CMB
 */
double DarkSUN::delta_n_eff_cmb() const {
    using namespace boost::math;
    return 4.0 / 7.0 * pow(11.0 / 4.0, 4.0 / 3.0) *
            dark_geff(kT_CMB * m_xi_cmb) * pow<4>(m_xi_cmb);
}

/**
 * Compute the contribution to delta N_eff from the dark SU(N) sector at BBN.
 * @return Delta N_eff at BBN
 */
double DarkSUN::delta_n_eff_bbn() const {
    using namespace boost::math;
    return 4.0 / 7.0 * dark_geff(kT_BBN * m_xi_bbn) * pow<4>(m_xi_bbn);
}


/**
 * Compute the solution to the Boltzmann equation.
 * @param reltol Relative tolerance
 * @param abstol Absolute tolerance
 * @param t_alg Algorithm to use
 */
void DarkSUN::solve_boltzmann(
        double reltol,
        double abstol,
        const std::string &t_alg = "radau"
) {
    using namespace diffeq;
    using namespace cosmology;

    // Reset parameters
    m_xi_fo = -1;
    m_Tsm_fo = -1;
    m_xi_bbn = -1;
    m_xi_cmb = -1;

    std::shared_ptr<DarkSUN> model = std::make_shared<DarkSUN>(*this);
    DarkSUNBoltzmann boltz{model};

    double Td_init = m_lam / 2.0;
    double xi_init = compute_xi_const_td(Td_init);
    double Tsm_init = Td_init / xi_init;

    double xstart = m_eta.get_mass() / Tsm_init;
    double xfinal = m_eta.get_mass() / kT_CMB;
    double w_eta_init = log(m_eta.neq(Td_init) / sm_entropy_density(Tsm_init));
    double y_delta_init = exp(-m_a * m_N) * m_delta.neq(Td_init) / sm_entropy_density(Tsm_init);

    auto logx_span = std::make_pair(std::log(xstart), std::log(xfinal));

    Vector<double> winit{2};
    winit(0) = w_eta_init;
    winit(1) = y_delta_init;

    ODEProblem problem{std::make_shared<DarkSUNBoltzmann>(boltz), winit, logx_span};

    ODEIntegratorOptions opts{};
    opts.abstol = abstol;
    opts.reltol = reltol;
    if (t_alg == "rodas") {
        Rodas alg{};
        m_solution = solve(problem, alg, opts);
    } else {
        Radau5 alg{};
        m_solution = solve(problem, alg, opts);
    }
    // Exponentiate the first component of the solutions, which are currently
    // W_eta = log(Y_eta).
    for (auto &vec: m_solution.us) {
        vec(0) = exp(vec(0));
    }


    m_xi_fo = model->get_xi_fo();
    m_Tsm_fo = model->get_Tsm_fo();
    m_xi_bbn = model->get_xi_bbn();
    try {
        m_xi_cmb = model->compute_xi(kT_CMB);
    } catch (...) {
        m_xi_cmb = std::nan("");
    }

    double eta_yinf, delta_yinf;
    if (m_solution.retcode == Success) {
        eta_yinf = m_solution.us.back()[0];
        delta_yinf = m_solution.us.back()[1];
        m_rd_eta = m_eta.get_mass() * eta_yinf * kS_TODAY / kRHO_CRIT;
        m_rd_delta = m_delta.get_mass() * delta_yinf * kS_TODAY / kRHO_CRIT;
        m_dneff_cmb = delta_n_eff_cmb();
        m_dneff_bbn = delta_n_eff_bbn();
        m_eta_si_per_mass = cross_section_2eta_2eta() / m_eta.get_mass();
        m_delta_si_per_mass = cross_section_2delta_2delta() / m_delta.get_mass();
    } else {
        m_rd_eta = std::nan("");
        m_rd_delta = std::nan("");
        m_dneff_cmb = std::nan("");
        m_dneff_bbn = std::nan("");
        m_eta_si_per_mass = cross_section_2eta_2eta() / m_eta.get_mass();
        m_delta_si_per_mass = cross_section_2delta_2delta() / m_delta.get_mass();
    }
}


}
}

#endif //LANRE_DM_MODELS_DARKSUN_HPP
