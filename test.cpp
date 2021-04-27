#include "pch.h"
#include "Satellite.h"
#include "CentralBody.h"
#include "GroundStation.h"
#include "PI.h"
#include "datevec.h"
#include <string>
#include <vector>
#include <complex>
#include <cmath>
using namespace std;
typedef complex<double> dcomp;

double tolerance = 1e-04;

TEST(CentralBodyTest, MuTest) {
	string body = "Earth";
	CentralBody cb(body);
	EXPECT_EQ(cb.getMu(), 398600.433);
}

TEST(CentralBodyTest, RTest) {
	string body = "Mercury";
	CentralBody cb(body);
	EXPECT_EQ(cb.getR(), 2440.53);
}

TEST(SatelliteTest, EccVecTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	vector<double> evect = sc.calcEccentricityVector();
	vector<double> expected_evect = { 0.00846604405136776, -0.00073541277437619, -0.00301627159784152 };
	EXPECT_TRUE(abs(evect[0] - expected_evect[0]) < tolerance);
	EXPECT_TRUE(abs(evect[1] - expected_evect[1]) < tolerance);
	EXPECT_TRUE(abs(evect[2] - expected_evect[2]) < tolerance);
}

TEST(SatelliteTest, EccTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	double e = sc.calcEccentricity();
	double expected_e = 0.00901735150586692;
	EXPECT_TRUE(abs(e - expected_e) < tolerance);
}

TEST(SatelliteTest, ThetaTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	double theta = sc.calcTheta();
	double expected_theta = 44.6081237118121 * PI / 180.0;
	EXPECT_TRUE(abs(theta - expected_theta) < tolerance);
}

TEST(SatelliteTest, BetaTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	double beta = sc.calcBeta();
	double expected_beta = 0.360504754845545 * PI / 180.0;
	EXPECT_TRUE(abs(beta - expected_beta) < tolerance);
}

TEST(SatelliteTest, Alpha0Test) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	double alpha0 = sc.calcAlpha0();
	double expected_alpha0 = 0.000146434895578188;
	EXPECT_TRUE(abs(alpha0 - expected_alpha0) < tolerance);
}

TEST(SatelliteTest, EccAnomTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	dcomp E = sc.calcE();
	dcomp expected_E = 44.2464545943028 * PI / 180.0;
	EXPECT_TRUE(abs(E - expected_E) < tolerance);
}

TEST(SatelliteTest, MeanAnomTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	dcomp M = sc.calcM();
	dcomp expected_M = 0.765954492848508;
	EXPECT_TRUE(abs(M - expected_M) < tolerance);
}

TEST(SatelliteTest, HyperAnomTest) {
	string body = "Saturn";
	CentralBody cb(body);
	vector<double> r0 = { -321601.0957, -584995.9962, -78062.5449 };	// km
	vector<double> v0 = { 8.57101142, 7.92783797, 1.90640217 };	// km/s
	Satellite sc(r0, v0, cb);
	double H = sc.calcH();
	double expected_H = -0.77403995924571;
	EXPECT_TRUE(abs(H - expected_H) < tolerance);
}

TEST(SatelliteTest, InitOrbitPropTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	dcomp x = sc.initOrbitProp();
	dcomp expected_x = -4.48166451835941e-06;
	EXPECT_TRUE(abs(x - expected_x) < tolerance);
}

TEST(SatelliteTest, SCInputTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	dcomp x = sc.initOrbitProp();
	dcomp SCInput = sc.calcSCInput(x);
	dcomp expected_SCInput = 2.94119127633455e-15;
	EXPECT_TRUE(abs(SCInput - expected_SCInput) < tolerance);
}

TEST(SatelliteTest, STest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	dcomp x = sc.initOrbitProp();
	dcomp S = sc.calcS(x);
	dcomp expected_S = 0.166666666666667;
	EXPECT_TRUE(abs(S - expected_S) < tolerance);
}

TEST(SatelliteTest, CTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	dcomp x = sc.initOrbitProp();
	dcomp C = sc.calcC(x);
	dcomp expected_C = 0.5;
	EXPECT_TRUE(abs(C - expected_C) < tolerance);
}

TEST(SatelliteTest, PropOrbitTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	double dt = 5 * 3600;
	dcomp x = sc.propagateOrbit(dt);
	dcomp expected_x = 1664.26223079728;
	EXPECT_TRUE(abs(x - expected_x) < tolerance);
}

TEST(SatelliteTest, RtTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	double dt = 5 * 3600;
	vector<double> rt = sc.calcrt(dt);
	vector<double> expected_rt = { -1872.47618067505, 5800.33676900589, 3143.6189000774 };
	EXPECT_TRUE(abs(rt[0] - expected_rt[0]) < tolerance);
	EXPECT_TRUE(abs(rt[1] - expected_rt[1]) < tolerance);
	EXPECT_TRUE(abs(rt[2] - expected_rt[2]) < tolerance);
}

TEST(SatelliteTest, VtTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	double dt = 5 * 3600;
	vector<double> vt = sc.calcvt(dt);
	vector<double> expected_vt = { -7.01875940054339, -2.74976102374641, 1.0249091900684 };
	EXPECT_TRUE(abs(vt[0] - expected_vt[0]) < tolerance);
	EXPECT_TRUE(abs(vt[1] - expected_vt[1]) < tolerance);
	EXPECT_TRUE(abs(vt[2] - expected_vt[2]) < tolerance);
}

TEST(SatelliteTest, Date2JDTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	vector<int> date = { 2021, 4, 23, 15, 23, 0 };
	double jd = sc.date2jd(date);
	double expected_jd = 2459328.14097222;
	EXPECT_TRUE(abs(jd - expected_jd) < tolerance);
}

TEST(SatelliteTest, JD2GMSTTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	vector<int> date = { 2021, 4, 23, 15, 23, 0 };
	double jd = sc.date2jd(date);
	double gmst = sc.jd2gmst(jd);
	double expected_gmst = 82.6430367763902 * PI / 180;
	EXPECT_TRUE(abs(gmst - expected_gmst) < tolerance);
}

TEST(SatelliteTest, JD2GASTTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	vector<int> date = { 2021, 4, 23, 15, 23, 0 };
	double jd = sc.date2jd(date);
	double gast = sc.jd2gast(jd);
	double expected_gast = 82.6394346411274 * PI / 180;
	EXPECT_TRUE(abs(gast - expected_gast) < tolerance);
}

TEST(SatelliteTest, ECI2ECEFTest) {
	string body = "Earth";
	CentralBody cb(body);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> v0 = { -3.931046491, 5.498676921, 3.665980697 };	// km/s
	Satellite sc(r0, v0, cb);
	vector<int> date = { 2021, 4, 23, 15, 23, 0 };
	double jd = sc.date2jd(date);
	vector<double> r0_ecef = sc.eci2ecef(jd, r0);
	vector<double> expected_r0_ecef = { 4654.7683545423, -4936.34140384493, 2.95581 };
	EXPECT_TRUE(abs(r0_ecef[0] - expected_r0_ecef[0]) < tolerance);
	EXPECT_TRUE(abs(r0_ecef[1] - expected_r0_ecef[1]) < tolerance);
	EXPECT_TRUE(abs(r0_ecef[2] - expected_r0_ecef[2]) < tolerance);
}

TEST(GroundStationTest, LLA2ECEFTest) {
	string body = "Earth";
	CentralBody cb(body);
	double lat = 15;
	double lon = 23;
	double alt = 0;
	GroundStation site(lat, lon, alt, cb);
	vector<double> ecef = site.lla2ecef();
	vector<double> expected_ecef = { 5672.32496343096, 2407.75909633351, 1640.10014019589 };
	EXPECT_TRUE(abs(ecef[0] - expected_ecef[0]) < tolerance);
	EXPECT_TRUE(abs(ecef[1] - expected_ecef[1]) < tolerance);
	EXPECT_TRUE(abs(ecef[2] - expected_ecef[2]) < tolerance);
}

TEST(GroundStationTest, AzElTest) {
	string body = "Earth";
	CentralBody cb(body);
	double lat = 15;
	double lon = 23;
	double alt = 0;
	GroundStation site(lat, lon, alt, cb);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	vector<double> azel = site.azel(r0);
	vector<double> expected_azel = { 2.4109374785519, 0.00475632632480929 };
	EXPECT_TRUE(abs(azel[0] - expected_azel[0]) < tolerance);
	EXPECT_TRUE(abs(azel[1] - expected_azel[1]) < tolerance);
}

TEST(GroundStationTest, InVisTrueTest) {
	string body = "Earth";
	CentralBody cb(body);
	double lat = 90;
	double lon = 0;
	double alt = 0;
	GroundStation site(lat, lon, alt, cb);
	vector<double> r0 = { 0, 0, 2.95581 };	// km
	EXPECT_FALSE(site.inVis(r0));
}

TEST(GroundStationTest, InVisFalseTest) {
	string body = "Earth";
	CentralBody cb(body);
	double lat = 15;
	double lon = 23;
	double alt = 0;
	GroundStation site(lat, lon, alt, cb);
	vector<double> r0 = { 5492.00034, 3984.0014, 2.95581 };	// km
	EXPECT_FALSE(site.inVis(r0));
}

TEST(DateVecTest, ValidDateTrueTest) {
	vector<int> date = { 1996, 11, 29, 15, 23, 0 };
	EXPECT_TRUE(valid_date(date));
}

TEST(DateVecTest, ValidDateFalseTest) {
	vector<int> date = { 1995, 12, 15, 24, 60, 60 };
	EXPECT_FALSE(valid_date(date));
}

TEST(DateVecTest, ValidDatePairTrueTest) {
	vector<int> start = { 1995, 12, 15, 15, 23, 0 };
	vector<int> end = { 1996, 11, 29, 15, 23, 0 };
	EXPECT_TRUE(valid_date_pair(start, end));
}

TEST(DateVecTest, ValidDatePairFalseTest) {
	vector<int> start = { 1996, 11, 29, 15, 23, 0 };
	vector<int> end = { 1995, 12, 15, 15, 23, 0 };
	EXPECT_FALSE(valid_date_pair(start, end));
}

TEST(DateVecTest, DateDiffTest) {
	vector<int> start = { 1995, 12, 15, 15, 23, 0 };
	vector<int> end = { 1996, 11, 29, 15, 23, 0 };
	double dt = datediff(start, end);
	double expected_dt = 30240000;
	EXPECT_TRUE(abs(dt - expected_dt) < tolerance);
}