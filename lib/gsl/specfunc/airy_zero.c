/* specfunc/airy_zero.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
#include "stdafx.h"
#include <config.h.in>
#include <gsl/gsl_math.h>
#include <gsl/err/gsl_errno.h>
#include <gsl/specfunc/gsl_sf_airy.h>

#include "error.h"

static const double zero_Ai[] = {
  0,
  -2.338107410459767039,
  -4.087949444130970617,
  -5.520559828095551059,
  -6.786708090071758999,
  -7.944133587120853123,
  -9.022650853340980380,
  -10.04017434155808593,
  -11.00852430373326289,
  -11.93601556323626252,
  -12.82877675286575720,
  -13.69148903521071793,
  -14.52782995177533498,
  -15.34075513597799686,
  -16.13268515694577144,
  -16.90563399742994263,
  -17.661300105697057509,
  -18.401132599207115416,
  -19.126380474246952144,
  -19.838129891721499701,
  -20.537332907677566360,
  -21.224829943642096955,
  -21.901367595585130707,
  -22.567612917496502831,
  -23.224165001121681061,
  -23.871564455535918567,
  -24.510301236589677490,
  -25.140821166148963748,
  -25.763531400982756459,
  -26.378805052137232374,
  -26.986985111606367686,
  -27.588387809882444812,
  -28.183305502632644923,
  -28.772009165237435382,
  -29.354750558766287963,
  -29.931764119086555913,
  -30.503268611418505287,
  -31.069468585183755604,
  -31.63055565801265934,
  -32.18670965295205069,
  -32.73809960900026913,
  -33.28488468190140188,
  -33.82721494950865194,
  -34.36523213386365906,
  -34.89907025034531210,
  -35.42885619274788846,
  -35.95471026189862926,
  -36.47674664437480896,
  -36.99507384699450161,
  -37.50979509200501613,
  -38.02100867725525443,
  -38.52880830509424882,
  -39.03328338327251391,
  -39.53451930072301805,
  -40.03259768075417603,
  -40.52759661388971821,
  -41.01959087233248966,
  -41.50865210780525018,
  -41.99484903432643004,
  -42.47824759730839188,
  -42.95891113021656009,
  -43.43690049989685412,
  -43.91227424156370168,
  -44.38508868433939023,
  -44.85539806814583243,
  -45.32325465267043011,
  -45.78870881905730086,
  -46.25180916491254629,
  -46.71260259315651633,
  -47.17113439520631705,
  -47.62744832892739292,
  -48.08158669175325711,
  -48.53359038933679845,
  -48.98349900006458366,
  -49.43135083573678341,
  -49.87718299868941729,
  -50.32103143561221860,
  -50.76293098829428522,
  -51.20291544151056412,
  -51.64101756824489758,
  -52.07726917242964943,
  -52.51170112936766183,
  -52.94434342398931824,
  -53.37522518708567514,
  -53.80437472964785717,
  -54.23181957543308298,
  -54.65758649186871111,
  -55.08170151939748312,
  -55.50418999935962251,
  -55.92507660050055598,
  -56.34438534418670066,
  -56.76213962840595327,
  -57.17836225062417808,
  -57.59307542956407782,
  -58.00630082596830627,
  -58.41805956240450934,
  -58.82837224216613231,
  -59.23725896731927534,
  -59.64473935594259360,
  -60.05083255860419805,
  -60.45555727411669871
};
static const size_t size_zero_Ai = sizeof(zero_Ai)/sizeof(double);


static const double zero_Bi[] = {
  0,
  -1.173713222709127925,
  -3.271093302836352716,
  -4.830737841662015933,
  -6.169852128310251260,
  -7.376762079367763714,
  -8.491948846509388013,
  -9.538194379346238887,
  -10.52991350670535792,
  -11.47695355127877944,
  -12.38641713858273875,
  -13.26363952294180555,
  -14.11275680906865779,
  -14.93705741215416404,
  -15.739210351190482771,
  -16.521419550634379054,
  -17.285531624581242533,
  -18.033113287225001572,
  -18.765508284480081041,
  -19.483880132989234014,
  -20.189244785396202420,
  -20.882495994193175768,
  -21.564425284712977653,
  -22.235737881803385167,
  -22.897065554219793474,
  -23.548977079642448269,
  -24.191986850649000086,
  -24.826562012152892172,
  -25.453128427085131994,
  -26.072075698466804494,
  -26.683761425120990449,
  -27.288514830076298204,
  -27.886639871735962459,
  -28.478417925678661737,
  -29.064110107777650305,
  -29.643959295918396591,
  -30.218191897047274645,
  -30.787019397921766297,
  -31.350639731255585371,
  -31.90923848358456965,
  -32.46298996683685318,
  -33.01205817205683814,
  -33.55659762084006113,
  -34.09675412765602851,
  -34.63266548426775468,
  -35.16446207582101720,
  -35.69226743681080479,
  -36.21619875398748222,
  -36.73636732230120657,
  -37.25287895916828697,
  -37.76583438165180116,
  -38.27532955056003997,
  -38.78145598496327279,
  -39.28430105019802461,
  -39.78394822205711298,
  -40.28047732954369150,
  -40.77396477829068148,
  -41.26448375650675678,
  -41.75210442510106287,
  -42.23689409345656643,
  -42.71891738216253539,
  -43.19823637387693118,
  -43.67491075336673948,
  -44.14899793766617113,
  -44.62055319719727274,
  -45.08962976861312825,
  -45.55627896004907928,
  -46.02055024940102076,
  -46.48249137619078661,
  -46.94214842752602207,
  -47.39956591861496210,
  -47.85478686825452176,
  -48.30785286967246692,
  -48.75880415707066192,
  -49.20767966818603897,
  -49.65451710315861501,
  -50.09935297997125482,
  -50.54222268670364757,
  -50.98316053082286586,
  -51.42219978571468262,
  -51.85937273464332870,
  -52.29471071231240525,
  -52.72824414418606069,
  -53.16000258371716397,
  -53.59001474761792882,
  -54.01830854929815828,
  -54.44491113058688729,
  -54.86984889184461534,
  -55.29314752056546491,
  -55.71483201856140365,
  -56.13492672781406761,
  -56.55345535507366411,
  -56.97044099527886475,
  -57.38590615386647834,
  -57.79987276803497897,
  -58.21236222702161974,
  -58.62339539144885603,
  -59.03299261179210306,
  -59.44117374601743460,
  -59.84795817643466996,
  -60.25336482580837088
};
static const size_t size_zero_Bi = sizeof(zero_Bi)/sizeof(double);


static const double zero_Aip[] = {
  0,
  -1.018792971647471089,
  -3.248197582179836738,
  -4.820099211178735639,
  -6.163307355639486822,
  -7.372177255047770177,
  -8.488486734019722133,
  -9.535449052433547471,
  -10.52766039695740728,
  -11.47505663348024529,
  -12.384788371845747325,
  -13.262218961665210382,
  -14.111501970462995282,
  -14.935937196720517467,
  -15.738201373692538303,
  -16.520503825433793542,
  -17.284695050216437357,
  -18.032344622504393395,
  -18.764798437665954740,
  -19.483221656567231178,
  -20.188631509463373154,
  -20.881922755516737701,
  -21.563887723198974958,
  -22.235232285348913331,
  -22.896588738874619001,
  -23.548526295928801574,
  -24.191559709526353841,
  -24.826156425921155001,
  -25.452742561777649948,
  -26.071707935173912515,
  -26.683410328322449767,
  -27.288179121523985029,
  -27.886318408768461192,
  -28.478109683102278108,
  -29.063814162638199090,
  -29.643674814632015921,
  -30.217918124468574603,
  -30.786755648012502519,
  -31.350385379083034671,
  -31.90899295843046320,
  -32.46275274623847982,
  -33.01182877663428709,
  -33.55637560978942190,
  -34.09653909480913771,
  -34.63245705463586589,
  -35.16425990255340758,
  -35.69207119851046870,
  -36.21600815233519918,
  -36.73618207994680321,
  -37.25269881785414827,
  -37.76565910053887108,
  -38.27515890473087933,
  -38.78128976408036876,
  -39.28413905729859644,
  -39.78379027246823278,
  -40.28032324990371935,
  -40.77381440566486637,
  -41.26433693758643383,
  -41.75196101547722703,
  -42.23675395695976012,
  -42.71878039026198233,
  -43.19810240513270670,
  -43.67477969292950869,
  -44.14886967681966886,
  -44.62042763293925724,
  -45.08950680327102630,
  -45.55615850092696446,
  -46.02043220845493728,
  -46.48237566972975615,
  -46.94203497593635633,
  -47.39945464610575493,
  -47.85467770262241617,
  -48.30774574208398774,
  -48.75869900186057804,
  -49.20757642267037247,
  -49.65441570746105074,
  -50.09925337686182515,
  -50.54212482144867502,
  -50.98306435104524282,
  -51.42210524126365311,
  -51.85927977747301469,
  -52.29461929636838876,
  -52.72815422529939506,
  -53.15991411950524351,
  -53.58992769739169611,
  -54.01822287397517367,
  -54.44482679260982599,
  -54.86976585510479430,
  -55.29306575033103518,
  -55.71475148140987392,
  -56.13484739156885235,
  -56.55337718874437424,
  -56.97036396900508167,
  -57.38583023886477265,
  -57.79979793654895377,
  -58.21228845227477578,
  -58.62332264760009139,
  -59.03292087389367419,
  -59.44110298997521892,
  -59.84788837897058171,
  -60.25329596442479317
};
static const size_t size_zero_Aip = sizeof(zero_Aip)/sizeof(double);


static const double zero_Bip[] = {
  0,
  -2.294439682614123247,
  -4.073155089071828216,
  -5.512395729663599496,
  -6.781294445990305390,
  -7.940178689168578927,
  -9.019583358794239067,
  -10.037696334908545802,
  -11.006462667712289940,
  -11.934261645014844663,
  -12.827258309177217640,
  -13.690155826835049101,
  -14.526645763485711410,
  -15.339693082242404109,
  -16.131724782385900578,
  -16.904759411889649958,
  -17.660498743114976102,
  -18.400394367181703280,
  -19.125697156412638066,
  -19.837494718415910503,
  -20.536740241453273980,
  -21.224275044889266569,
  -21.900846445139208281,
  -22.567122080497200470,
  -23.223701521208962116,
  -23.871125771677973595,
  -24.509885117016242729,
  -25.140425655367878908,
  -25.763154776913454319,
  -26.378445791146615697,
  -26.986641859775034987,
  -27.588059359225600573,
  -28.182990771292975456,
  -28.771707180886056250,
  -29.354460444612957224,
  -29.931485082026055160,
  -30.502999931936645516,
  -31.069209608721234058,
  -31.63030578754333679,
  -32.18646834257807369,
  -32.73786635840274752,
  -33.28465903151424981,
  -33.82699647630635587,
  -34.36502044767239631,
  -34.89886499060196419,
  -35.42865702564380962,
  -35.95451687785511190,
  -36.47655875580547918,
  -36.99489118631672770,
  -37.50961740986809593,
  -38.02083574095788210
};
static const size_t size_zero_Bip = sizeof(zero_Bip)/sizeof(double);



/* [Abramowitz+Stegun, 10.4.105] */
static double
zero_f(double z)
{
  const double pre = pow(z, 2.0/3.0);
  const double zi2 = 1.0/(z*z);
  const double zi4 = zi2 * zi2;
  const double t1  =  5.0/48.0 * zi2;
  const double t2  = -5.0/36.0 * zi4;
  const double t3  =  77125.0/82944.0 * zi4 * zi2;
  const double t4  = -108056875.0/6967296.0 * zi4 * zi4;
  return pre * (1.0 + t1 + t2 + t3 + t4);
  
}
static double
zero_g(double z)
{
  const double pre = pow(z, 2.0/3.0);
  const double zi2 = 1.0/(z*z);
  const double zi4 = zi2 * zi2;
  const double t1  = -7.0/48.0 * zi2;
  const double t2  =  35.0/288.0 * zi4;
  const double t3  = -181223.0/207360.0 * zi4 * zi2;
  const double t4  =  18683371.0/1244160.0 * zi4 * zi4;
  return pre * (1.0 + t1 + t2 + t3 + t4);
}



int
gsl_sf_airy_zero_Ai_e(unsigned int s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s < 1) {
    DOMAIN_ERROR_MSG("s is less than 1", result);
  }
  else if(s < size_zero_Ai) {
    result->val = zero_Ai[s];
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double z = 3.0*M_PI/8.0 * (4.0*s - 1.0);
    const double f = zero_f(z);
    result->val = -f;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_airy_zero_Bi_e(unsigned int s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s < 1) {
    DOMAIN_ERROR_MSG("s is less than 1", result);
  }
  else if(s < size_zero_Bi) {
    result->val = zero_Bi[s];
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double z = 3.0*M_PI/8.0 * (4.0*s - 3.0);
    const double f = zero_f(z);
    result->val = -f;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_airy_zero_Ai_deriv_e(unsigned int s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s < 1) {
    DOMAIN_ERROR_MSG("s is less than 1", result);
  }
  else if(s < size_zero_Aip) {
    result->val = zero_Aip[s];
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double z = 3.0*M_PI/8.0 * (4.0*s - 3.0);
    const double g = zero_g(z);
    result->val = -g;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_airy_zero_Bi_deriv_e(unsigned int s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s < 1) {
    DOMAIN_ERROR_MSG("s is less than 1", result);
  }
  else if(s < size_zero_Bip) {
    result->val = zero_Bip[s];
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double z = 3.0*M_PI/8.0 * (4.0*s - 1.0);
    const double g = zero_g(z);
    result->val = -g;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_airy_zero_Ai(unsigned int s)
{
  EVAL_RESULT(gsl_sf_airy_zero_Ai_e(s, &result));
}

double gsl_sf_airy_zero_Bi(unsigned int s)
{
  EVAL_RESULT(gsl_sf_airy_zero_Bi_e(s, &result));
}

double gsl_sf_airy_zero_Ai_deriv(unsigned int s)
{
  EVAL_RESULT(gsl_sf_airy_zero_Ai_deriv_e(s, &result));
}

double gsl_sf_airy_zero_Bi_deriv(unsigned int s)
{
  EVAL_RESULT(gsl_sf_airy_zero_Bi_deriv_e(s, &result));
}
