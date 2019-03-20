#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include "geometryselect.h"
#include <math.h>
#include "geometryselect.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    double a,b;
    a = cos(0.0*M_PI/180.0);
    b = cos(0.01*M_PI/180.0);
   // qDebug()<<QString::number(a,'f',10)<<QString::number(b,'f',10);

    //25263  46800.000 7007937.883  0.4000 45.0    139.188109  67.794207  30.933685
                     //7007937.883  0.4000 45.0000 319.188109  67.794207  30.933685
    GeometrySelect *select = new GeometrySelect(NULL);
    select->setMainSatellitePara(25263,46800,7007937.883/6378137.0,0.04,45.0/180.0*M_PI,139.188109/180.0*M_PI,67.794207/180.0*M_PI,30.933685/180.0*M_PI);
    select->setOtherSatellitePara(25263,46800,7007937.883/6378137.0,0.04,45.0/180.0*M_PI,319.188109/180.0*M_PI,67.794207/180.0*M_PI,30.933685/180.0*M_PI);
    select->setMain_F(33.3945/180.0*M_PI);
    select->setOther_F(33.3945/180.0*M_PI);

    int aaa  = select->cal_run();
    qDebug()<<aaa;

}
