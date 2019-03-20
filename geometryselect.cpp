#include "geometryselect.h"

GeometrySelect::GeometrySelect(QObject *parent) : QObject(parent)
{

}

void GeometrySelect::setMainSatellitePara(double JD_p2,double TS_p2,double a_p2,double e_p2,double i_p2,
                         double Omega_p2,double w_p2 ,double M_p2)
{
    JD_p = JD_p2;
    TS_p = TS_p2;
    a_p = a_p2;
    e_p = e_p2;
    i_p = i_p2;
    Omega_p = Omega_p2;
    w_p = w_p2;
    M_p = M_p2;


}

void GeometrySelect::setOtherSatellitePara(double JD_s2,double TS_s2,double a_s2,double e_s2,double i_s2,
                         double Omega_s2,double w_s2 ,double M_s2)
{
    JD_s = JD_s2;
    TS_s = TS_s2;
    a_s = a_s2;
    e_s = e_s2;
    i_s = i_s2;
    Omega_s = Omega_s2;
    w_s = w_s2;
    M_s = M_s2;


}

void GeometrySelect::setMain_F(double f_p2)
{
    f_p =f_p2;
}

void GeometrySelect::setOther_F(double f_s2)
{
    f_s=f_s2;

}

double  GeometrySelect::cal_K1()
{
    double K1;
    K1 = sin(Omega_s)*sin(i_s)*cos(i_p) - sin(Omega_p)*sin(i_p)*cos(i_s);
    return K1;
}

double  GeometrySelect::cal_K2()
{
    double K2 = cos(Omega_p)*sin(i_p)*cos(i_s) - cos(Omega_s)*sin(i_s)*cos(i_p);
    return K2;
}

double  GeometrySelect::cal_K3()
{
    double K3 = sin(i_s)*sin(i_p)*sin(Omega_p - Omega_s);
    return K3;

}

double  GeometrySelect::cal_SinI_R(double K1, double K2, double K3)
{
    double SinI_R = sqrt(K1*K1+K2*K2+K3*K3);
    return SinI_R;

}

double  GeometrySelect::cal_CosI_R(double sinI_R)
{
    double CosI_R = sqrt(1.0-sinI_R*sinI_R);
    return CosI_R;

}
double  GeometrySelect::cal_CosDelta_p(double sinI_R)
{
    double cosDelta_p = 1.0/sinI_R *(sin(i_p)*cos(i_s) - sin(i_s)*cos(i_p)*cos(Omega_p - Omega_s)) ;
    return cosDelta_p;

}

double  GeometrySelect::cal_SinDelta_p(double sinI_R)
{
    double SinDelta_p = 1.0/sinI_R *(sin(i_s)*sin(Omega_p - Omega_s));
    return SinDelta_p;

}

double GeometrySelect::cal_CosDelta_s(double sinI_R)
{
    double CosDelta_s = 1.0/sinI_R *(sin(i_p)*cos(i_s)*cos(Omega_p - Omega_s) - sin(i_s)*cos(i_p));
    return CosDelta_s;

}

double  GeometrySelect::cal_SinDelta_s(double sinI_R)
{
    double SinDelta_s = 1.0/sinI_R * (sin(i_p)*sin(Omega_p - Omega_s));
    return SinDelta_s;

}

double  GeometrySelect::cal_Delta_p(double SinDelta_p,double CosDelta_p)
{

    double Delta_p = 0.0;

    if( fabs(CosDelta_p) < 0.00001)
    {
        if(SinDelta_p >= 0)
        {
            Delta_p = M_PI / 2.0;
        }else
        {
            Delta_p = 3.0 * M_PI / 2.0;
        }
    }else if(CosDelta_p > 0)
    {
        if(SinDelta_p > 0)
        {
            Delta_p = atan(SinDelta_p/CosDelta_p);
        }else
        {
            Delta_p = 2.0 * M_PI + atan(SinDelta_p/CosDelta_p);
        }

    }else if(CosDelta_p < 0)
    {
         Delta_p = M_PI + atan(SinDelta_p/CosDelta_p);
    }


//        if(dabs(cosdp).lt.1d-5) then
//                 if(sindp.gt.0) then
//                     dp=pi/2.
//                 else
//                     dp=3*pi/2
//                 end if
//             else
//                 if(cosdp.gt.0) then
//                     if(sindp.gt.0) then
//                         dp=atan(sindp/cosdp)
//                     else
//                         dp=2*pi+atan(sindp/cosdp)
//                     end if
//                 else
//                     dp=pi+atan(sindp/cosdp)
//                 end if
//             end if

    return Delta_p;
}

double   GeometrySelect::cal_Delta_s(double SinDelta_s,double CosDelta_s)
{
     double Delta_s = cal_Delta_p(SinDelta_s,CosDelta_s);
     return Delta_s;

}

double  GeometrySelect::cal_u_R_p(double f_p,double Delta_p)
{
    double u_R_p=w_p + f_p - Delta_p ;
    return u_R_p;
}


double   GeometrySelect::cal_u_R_s(double f_s, double Delta_s)
{
    double u_R_s = w_s + f_s -Delta_s;
    return u_R_s;

}

double   GeometrySelect::cal_CosR(double u_R_p,double u_R_s,double CosI_R)
{
    double CosR = cos(u_R_p)*cos(u_R_s) + sin(u_R_p)*sin(u_R_s)*CosI_R;
    return CosR;

}

double  GeometrySelect::cal_A(double u_R_p,double f_p)
{
    double A = sin(u_R_p) + e_p*sin(u_R_p - f_p);
    return A;
}

double GeometrySelect::cal_B(double u_R_p,double f_p)
{
    double B = cos(u_R_p) + e_p*cos(u_R_p - f_p);
    return B;
}

double  GeometrySelect::cal_C(double u_R_s,double f_s)
{
    double C = sin(u_R_s) + e_s*sin(u_R_s - f_s);
    return C;

}

double  GeometrySelect::cal_D(double u_R_s,double f_s)
{
    double D = cos(u_R_s) + e_s*cos(u_R_s - f_s);
    return D;
}

double  GeometrySelect::cal_r_p(double f_p)
{
    double  r_p = a_p*(1.0-e_p*e_p)/(1.0+e_p*cos(f_p));
    return r_p;

}

double  GeometrySelect::cal_r_s(double f_s)
{
    double r_s = a_s*(1.0-e_s*e_s)/(1.0+e_s*cos(f_s));
    return r_s;
}


double  GeometrySelect::cal_CosE_p(double f_p)
{
    double CosE_p = (e_p+cos(f_p))/(1.0+e_p*cos(f_p));
    return CosE_p;
}

double   GeometrySelect::cal_CosE_s(double f_s)
{
    double CosE_s = (e_s + cos(f_s))/(1.0 + e_s*cos(f_s));
    return CosE_s;
}

double  GeometrySelect::cal_F(double r_p,double f_p,double r_s,double u_R_s,double CosI_R,double A,double B)
{
    double F = r_p*e_p*sin(f_p)+r_s*(A*cos(u_R_s)-B*sin(u_R_s)*CosI_R);
    return F;

}

double  GeometrySelect::cal_G(double r_s,double f_s,double r_p,double u_R_p,double CosI_R,double C,double D)
{
    double G = r_s*e_s*sin(f_s)+r_p*(C*cos(u_R_p)-D*sin(u_R_p)*CosI_R);
    return G;

}

double  GeometrySelect::cal_ZF_Zf_p(double r_p,double CosE_p,double r_s,double CosR)
{
    double ZF_Zf_p = r_p*e_p*CosE_p + r_s*CosR;
    return ZF_Zf_p;
}

double  GeometrySelect::cal_ZF_Zf_s(double r_s,double f_s,double A,double B,double C,double D,double CosI_R)
{
    double ZF_Zf_s = -r_s*(A*C + B*D*CosI_R)/(1.0+e_s*cos(f_s)) ;
    return ZF_Zf_s;

}

double GeometrySelect::cal_ZG_Zf_p(double r_p,double f_p,double A,double B,double C,double D,double CosI_R)
{
    double ZG_Zf_p = -r_p*(A*C +B*D*CosI_R)/(1.0+e_p*cos(f_p));
    return ZG_Zf_p;
}


double  GeometrySelect::cal_ZG_Zf_s(double r_s,double CosE_s,double r_p,double CosR)
{
    double ZG_Zf_s = r_s*e_s*CosE_s + r_p * CosR ;
    return ZG_Zf_s;
}

double  GeometrySelect::cal_Delta_f_p(double F,double G,double ZF_Zf_p,double ZF_Zf_s,double ZG_Zf_p,double ZG_Zf_s )
{
    double Delta_f_p = (G*ZF_Zf_s - F*ZG_Zf_s)/(ZF_Zf_p * ZG_Zf_s - ZF_Zf_s*ZG_Zf_p);
    return Delta_f_p;

}

double  GeometrySelect::cal_Delta_f_s(double F,double G,double ZF_Zf_p,double ZF_Zf_s,double ZG_Zf_p,double ZG_Zf_s )
{
    double Delta_f_s = (F*ZG_Zf_p - G*ZF_Zf_p)/(ZF_Zf_p * ZG_Zf_s - ZG_Zf_p*ZF_Zf_s);
    return Delta_f_s;

}

int GeometrySelect::cal_run()
{
    double limitDelta = 0.00000015;
    double limitRange = 10000.0/6378137.0;
    double K1 = cal_K1();
    double K2 = cal_K2();
    double K3 = cal_K3();
    double SinI_R = cal_SinI_R(K1,K2,K3);
    double CosI_R = cal_CosI_R(SinI_R);
    double CosDelta_p = cal_CosDelta_p(SinI_R);
    double SinDelta_p = cal_SinDelta_p(SinI_R);
    double CosDelta_s = cal_CosDelta_s(SinI_R);
    double SinDelta_s = cal_SinDelta_s(SinI_R);
    double Delta_p = cal_Delta_p(SinDelta_p,CosDelta_p);
    double Delta_s = cal_Delta_s(SinDelta_s,CosDelta_s);

    if(CosI_R > cos(0.01*M_PI/180.0))
    {
        return  -1;
    }

    int loopCount =0; // loop counter
    double f_p;
    double f_s;
    double u_R_p;
    double u_R_s;
    double CosR;
    double A,B,C,D;
    double r_p;
    double r_s;
    double CosE_p;
    double CosE_s;
    double F,G;
    double ZF_Zf_p;
    double ZF_Zf_s;
    double ZG_Zf_p;
    double ZG_Zf_s;
    double next_f_p;
    double next_f_s;
    double Delta_f_p;
    double Delta_f_s;


    while(1)
    {
        if(loopCount == 0)
        {
            f_p = Delta_p - w_p;
            f_s = Delta_s - w_s;
            u_R_p = cal_u_R_p(f_p,Delta_p);
            u_R_s = cal_u_R_s(f_s,Delta_s);
        }else
        {
            u_R_p = cal_u_R_p(f_p,Delta_p);
            u_R_s = cal_u_R_s(f_s,Delta_s);
        }

        CosR = cal_CosR(u_R_p,u_R_s,CosI_R);
        A = cal_A(u_R_p,f_p);
        B = cal_B(u_R_p,f_p);
        C = cal_C(u_R_s,f_s);
        D = cal_D(u_R_s,f_s);
        r_p = cal_r_p(f_p);
        r_s = cal_r_s(f_s);
        CosE_p = cal_CosE_p(f_p);
        CosE_s = cal_CosE_s(f_s);
        F = cal_F(r_p,f_p,r_s,u_R_s,CosI_R,A,B);
        G = cal_G(r_s,f_s,r_p,u_R_p,CosI_R,C,D);
        ZF_Zf_p = cal_ZF_Zf_p(r_p,CosE_p,r_s,CosR);
        ZF_Zf_s = cal_ZF_Zf_s(r_s,f_s,A,B,C,D,CosI_R);
        ZG_Zf_p = cal_ZG_Zf_p(r_p,f_p,A,B,C,D,CosI_R);
        ZG_Zf_s = cal_ZG_Zf_s(r_s,CosE_s,r_p,CosR);
        Delta_f_p = cal_Delta_f_p(F,G,ZF_Zf_p,ZF_Zf_s,ZG_Zf_p,ZG_Zf_s);
        Delta_f_s = cal_Delta_f_s(F,G,ZF_Zf_p,ZF_Zf_s,ZG_Zf_p,ZG_Zf_s);

        next_f_p = f_p +Delta_f_p;
        next_f_s = f_s +Delta_f_s;

        f_p = next_f_p;
        f_s = next_f_s;

        if(fabs(Delta_f_p)<limitDelta && fabs(Delta_f_s)<limitDelta)
        {
            qDebug()<<"loopCount = "<<loopCount;
            qDebug()<<"CosR="<<QString::number(CosR,'f',15);
            qDebug()<<"CosI_R="<<QString::number(CosI_R,'f',15);
            qDebug()<<"u_R_p="<<QString::number(u_R_p,'f',15);
            break;

        }


        loopCount = loopCount + 1;
    }

    double r_p_star = a_p*(1.0 - e_p*e_p)/(1.0 + e_p * cos(f_p));
    double r_s_star = a_s*(1.0 - e_s*e_s)/(1.0 + e_s * cos(f_s));
    double r_rel = sqrt(r_p_star * r_p_star + r_s_star * r_s_star -2*r_p_star*r_s_star*CosR);

    double f_p_2_star = f_p + M_PI ;
    double f_s_2_star = f_s + M_PI ;

    double r_p_star2 = a_p*(1.0 - e_p*e_p)/(1.0 + e_p * cos(f_p_2_star));
    double r_s_star2 = a_s*(1.0 - e_s*e_s)/(1.0 + e_s * cos(f_s_2_star));
    double r_rel2 = sqrt(r_p_star2 * r_p_star2 + r_s_star2 * r_s_star2 -2*r_p_star2*r_s_star2*CosR);

    qDebug()<<"f_p,f_s:  "<<f_p*180.0/M_PI<<f_s*180.0/M_PI;

    qDebug()<<"r_rel,r_rel2:"<<r_rel*6378137.0<<r_rel2*6378137.0;

    if(r_rel >limitRange && r_rel2>limitRange)
    {
        return -1;
    }else
    {
        return 1;
    }



}
