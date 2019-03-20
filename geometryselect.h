#ifndef GEOMETRYSELECT_H
#define GEOMETRYSELECT_H

#include <QObject>
#include <math.h>
#include <QDebug>

class GeometrySelect : public QObject
{
    Q_OBJECT
public:
    explicit GeometrySelect(QObject *parent = 0);
    double JD_p;
    double TS_p;
    double a_p;
    double e_p;
    double i_p;
    double Omega_p;
    double w_p;
    double M_p;
    double f_p;

    double JD_s;
    double TS_s;
    double a_s;
    double e_s;
    double i_s;
    double Omega_s;
    double w_s;
    double M_s;
    double f_s;


    void setMainSatellitePara(double JD_p2,double TS_p2,double a_p2,double e_p2,double i_p2,
                             double Omega_p2,double w_p2 ,double M_p2);

    void setOtherSatellitePara(double JD_s2,double TS_s2,double a_s2,double e_s2,double i_s2,
                             double Omega_s2,double w_s2 ,double M_s2);

    void setMain_F(double f_p2);
    void setOther_F(double f_s2);
    double  cal_K1();
    double  cal_K2();
    double  cal_K3();
    double  cal_SinI_R(double K1,double K2,double K3);
    double  cal_CosI_R(double sinI_R);
    double  cal_CosDelta_p(double sinI_R);
    double  cal_SinDelta_p(double sinI_R);
    double  cal_CosDelta_s(double sinI_R);
    double  cal_SinDelta_s(double sinI_R);

    double  cal_Delta_p(double SinDelta_p,double CosDelta_p);
    double  cal_Delta_s(double SinDelta_s,double CosDelta_s);

    double  cal_u_R_p(double f_p, double Delta_p);
    double  cal_u_R_s(double f_s, double Delta_s);

    double  cal_CosR(double u_R_p, double u_R_s, double CosI_R);
    double  cal_A(double u_R_p,double f_p);
    double  cal_B(double u_R_p,double f_p);
    double  cal_C(double u_R_s,double f_s);
    double  cal_D(double u_R_s,double f_s);
    double  cal_r_p(double f_p);
    double  cal_r_s(double f_s);

    double  cal_CosE_p(double f_p);
    double  cal_CosE_s(double f_s);
    double  cal_F(double r_p,double f_p,double r_s,double u_R_s,double CosI_R,double A,double B);
    double  cal_G(double r_s,double f_s,double r_p,double u_R_p,double CosI_R,double C,double D);

    double  cal_ZF_Zf_p(double r_p,double CosE_p,double r_s,double CosR);
    double  cal_ZF_Zf_s(double r_s,double f_s,double A,double B,double C,double D,double CosI_R);
    double  cal_ZG_Zf_p(double r_p,double f_p,double A,double B,double C,double D,double CosI_R);
    double  cal_ZG_Zf_s(double r_s,double CosE_s,double r_p,double CosR);

    double  cal_Delta_f_p(double F,double G,double ZF_Zf_p,double ZF_Zf_s,double ZG_Zf_p,double ZG_Zf_s );
    double  cal_Delta_f_s(double F,double G,double ZF_Zf_p,double ZF_Zf_s,double ZG_Zf_p,double ZG_Zf_s );

    int    cal_run();
signals:

public slots:
};

#endif // GEOMETRYSELECT_H
