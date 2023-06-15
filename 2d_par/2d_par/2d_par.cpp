// 2d_par.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include<string>
#include <iomanip>
#include <cstdlib>
#include</Program Files (x86)/Microsoft SDKs/MPI/Include/mpi.h>

using namespace std;

// Определение глобальных переменных (общих блоков)
int NX, NZ, NZ5, NZ10, Ntotal;
double AZ, AX, ash2;
double dx, dz, dz1_24, dz9_8, dx1_24, dx9_8;
ofstream FILE2;

//Параметры

const int MXPROCS = 128;
const int IPARSIZ = 13;
const int DPARSIZ = 6;
const int BASIS = 10;
const int BASIS1 = 2;
void DLAG2(vector<double>& Y, double& X, long long& N, int& LS);
void KOEFlag(vector<double>& FN, vector<double>& CR, vector<double>& QN, double& RE, long long& LAGMAX, int& LS, double& epsLag, long long& NW, int& KER);
void DLAG1(vector<double>& Y, double& X, long long& N, int& LS);
void PART(int N, int P, int(&RANGES)[3][128]);
void MATRIZA(int& J0, int& J1, vector<double>& d2pmlz, vector<double>& d2pmlx, vector<double>& UX1,
    vector<double>& UX2, vector<double>& UZ1, vector<double>& UZ2, vector<vector<double>>& fro, vector<vector<double>>& fmu,
    vector<vector<double>>& fla, int& NX, int& NZ, int& NZ5, int& Ntotal);
void PIMDSETPAR(int* IPAR, double* DPAR, int& LDA, int& N, int& BLKSZ, int& LOCLEN, int& BASISDIM,
    int& NPROCS, int& PROCID, int& PRECONTYPE, int& STOPTYPE, int& MAXIT,
    double EPSILON);
void DINIT(int& N, double ALPHA, double*& DX, int INCX);
void ZAPIS(ofstream& NDIR, int NTRAS, int NSAMP, vector<vector<double>> REG1, vector<double> REG2);
void REPORT(string NAME, int* IPAR, double* DPAR, double ET, double*& X, int nout);
void PROIZVOD(int J0, int J1, vector<double>& UXDZ, vector<double>& UZDZ, vector<double>& SXDZ, vector<double>& SZDZ, double*& X, int KMV, double ETMV1);
void PIMDCGS(int& KMV, double& ETMV1, vector<double>& d2pmlx, vector<double>& UZ1, vector<double>& UZ2, vector<double>& UX1,
    vector<double>& UX2, double*& SP, double*& SL, vector<double>& XS, vector<double>& Q1, vector<double>& Q2, double*& X,
    double*& B, double*& WRK, int* IPAR, double* DPAR);
void DCOPY(int& n, double* dx, int incx, double* dy, int incy);
void DVPROD(int n, vector<double>& dx, int incx, double*& dy, int incy);
void DAXPY(int& n, double da, double* dx, int incx, double* dy, int incy);
double DDOT(int n, double* dx, int incx, double* dy, int incy);
void PDSUM(int ISIZE, double* x);
double PDNRM2(int& LOCLEN, double* U);
void STOPCRIT(int& KMV, double& ETMV1, vector<double> d2pmlx, double*& B, double* R, double* RTRUE,
    double*& X, double* XOLD, double* WRK, double RHSSTOP,
    int CNVRTX, double EXITNORM, int STATUS, int*& IPAR,
    vector<double>& Q2, double*& SP, double*& SL, vector<double>& XS, vector<double>& UZ1, vector<double>& UZ2, vector<double>& UX1, vector<double>& UX2);
double DSETRHSSTOP(double* B, double* WRK, double& EPSILON, int*& IPAR, vector<double>& Q1);
void PIMDGETPAR(int*& IPAR, double*& DPAR, int& LDA, int& N, int& BLKSZ, int& LOCLEN, int& BASISDIM,
    int& NPROCS, int& PROCID, int& PRECONTYPE, int& STOPTYPE, int& MAXIT, int& ITNO, int& STATUS, int& STEPERR, double EPSILON, double EXITNORM);

double Timer();

int RANGES[3][MXPROCS];

int main(int argc, char** argv)
{
    //порядок используемых функций Лаггера
    int LS=8;

    //требуемая точность разложения по функциям Лагерра
    double epsLag = 1e-3;
    //требуемая точность решения системы уравнений
    double TOL = 1e-8;
    int PRET = 1;
    //критерий остановки операции
    int STOPT = 5;

    //количество картинок для фиксированого t
    const int KART = 250;
    const int KART1 = KART + 1;
    const double PII = 3.141592653589793e0;

    //локальные переменные
    int NL[KART1], LAG[3];  
    double TN[KART1];
    string NPICT[KART], str;

    //локальные скаляры
    long long nxs1, nxs2, ind, NUMPICT, IPARAM, KPAR, KOL, KTR1, KTR2, NTR0;
    long long Lpmlv, Lpmlz, Lpmlx1, Lpmlx2, LIZ, LIX, LPZ, LPX, LPdX, Kdx, Kdz;
    long long ALERR, KOLIST, NPI, NIP, KT, KOLPR, KTRP, NN, KX, NKOMP, nout;
    long long KTRI, KTRAS, NTRAS, KTR, NTR, NXNZ, NW, LAGmax;
    double DTI, XNX, VP, VS, RO, AX1, AZ1, F0, DXP, DZP;
    double  AA, BB, xro, xmu, xla, xlamu, XNPI, up1, up2, up3;
    double  ash, RE, RE1, xx, dd, s1, s2, t0, DXX, DTP, TMIN;
    double  pmlxl, pmlxr, pmlv, pmlz, apxl, apxr, apv, apz, cc;
    double  X0, Z0, Ztras, X0tras, DXtras, tt, t1, dtras;

    //double XI[KOLIST], ZI[KOLIST], XP[KOLIST], ZP[KOLIST];
    double XI[1000], ZI[1000], XP[1000], ZP[1000];

    //int kgran[NZ]
    vector<double> kgran;

    //double zro[NZ],zmu[NZ],zla[NZ],zlamu[NZ];
    vector<double> zro1, zmu1, zla1, zlamu1;
    vector<double> zro2, zmu2, zla2, zlamu2;

    //double fro(NZ,NX), fmu(NZ,NX), fla(NZ,NX)
    vector<vector<double>> fro, fmu, fla;//двумерные массивы

    //double  RVP2(NXNZ), RVS2(NXNZ), ROD1(NXNZ)
    vector<double> RVP2, RVS2, ROD1;

    //double UXDZ(NXNZ),UZDZ(NXNZ),SXDZ(NXNZ),SZDZ(NXNZ)
    vector<double> UXDZ, UZDZ, SXDZ, SZDZ;

    //double UZ1(LOCLEN), UZ2(LOCLEN), UX1(LOCLEN), UX2(LOCLEN)
    vector<double> UZ1, UZ2, UX1, UX2;//компоненты скорости смещения

    //double FIST[KTRI]
    vector<double> FIST;

    //double FN(LAGmax), CR(LAGmax), QN(LAGmax)
    vector<double> FN, CR, QN;

    //Ntotal = 5*NX*NZ-4*NX
    //double SUMMA(Ntotal), RWSUM(Ntotal)
    double* SUMMA, *RWSUM;

    //double UPRXX(NTRAS,NW), UPRZZ(NTRAS,NW), UPRP(NTRAS,NW)
    vector<vector<double>> UPRXX, UPRZZ, UPRP;

    //double TRASX1(NTRAS),TRASX2(NTRAS),TRASZ1(NTRAS),TRASZ2(NTRAS)
    double* TRASX1, *TRASX2, *TRASZ1, *TRASZ2;

    //double TRASP1(NTRAS), TRASP2(NTRAS)
    double* TRASP1, *TRASP2;

    //double UTRAX(NTRAS,KTRAS),UTRAZ(NTRAS,KTRAS),UTRAP(NTRAS,KTRAS),
    vector<vector<double>> UTRAX, UTRAZ, UTRAP;

    //double RW1(NTRAS), RW2(NTR)
    vector<double> RW1;
    double* RW2 = new double[1]{ 0.0 };

    //double REG(KTR), REGtras(KTRAS)
    vector<double>  REG, REGtras;

    //double PICTZ1(NX),PICTZ2(NX)
    double* PICTZ1, *PICTZ2;

    //double UTRAZ1(NTR,KTR), UTRLZ(NEL,KTR,KART)
    vector<vector<double>> UTRAZ1;//двумерный
    vector<vector<vector<double>>> UTRLZ;//трехмерный

    //double CRN(KART,NW)
    vector<vector<double>> CRN;

    //double d2pmlz(NZ), d2pmlx(NX)
    vector<double> d2pmlz, d2pmlx;

    char line[34] = "---------------------------------";

    //PIM -- The Parallel Iterative Methods package
    

    //локальные скаляры
    double ET, ET0, ET1, ETtotal, ET1total;
    double ETMV0, ETMV1;
    int C, V, MAXIT, IERR, KMV;
    int J0, J1, JJ, MYID, MYN, NPROCS;
    int N, NEL, LOCLEN, LDWRK, LOCS;
    //double elapsed_time

    //локальные массивы
    double DPAR[DPARSIZ];
    int IPAR[IPARSIZ];
    

    //double SP(NZ10),SL(NZ10),XS(LOCS=LOCLEN+2*NZ10)
    vector<double>  XS;
    double* SP, *SL;
    //double Q1(LOCLEN),Q2(LOCLEN)
    vector<double> Q1, Q2;
    // double DWRK(LDWRK),X(LOCLEN),FWW(LOCLEN),SUMLOC(LOCLEN),SUMMApml(LOCLEN)
    vector<double> SUMLOC,SUMMApml;
    double* FWW;
    double* X;
    double* DWRK;

    // внешние функции

    //внешние подпрограммы

    //MPI подпрограммы

    /*INT из FORTRAN можно заменить на функцию round из
        <cmath> в C++
        , а SQRT на функцию sqrt из <cmath>.*/



// AX, AZ - размеры области модели (в километрах) по осям X и Z соответственно, 
//        - действительные положительные числа с плавающей запятой формата (2E15.8), любые значения > 0.
// NX - количество узлов сетки для разностной аппроксимации по координате X, 
//    - положительное целое число типа int, минимальное значение = 3. 
// NZ - количество узлов сетки для разностной аппроксимации по координате Z, 
//    - положительное целое число типа int, минимальное значение = 3.
// dx - шаг по сетке X.
// dz - шаг по сетке Z.

//Засекаем время
//elapse_time=Timer()
double ET0total = Timer();




    //MPI ---------------------------
    //инициализация MPI 
    //Возращает 0 - усмешная операция
    //другое- код ошибки
    MPI_Init(&argc, &argv);
    // Определение порядкового номера процессора в группе MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &MYID) != 0;
    // Определение количества процессоров в группе MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);

    std::cout << "I am " << MYID << " processor ready to start"<<endl;

    //	IPARAM=1     ! 1 - считывание предварительных результатов
    IPARAM = 0;
    //KPAR=1     ! 1 - zapis rezultatov (dlja 1000 garmonik Lagera)
    KPAR = 0;
    nout = 2;
    //pml - размер PML зонs в км
    pmlxl = 1.5e0;
    apxl = 70.0e0;
    pmlxr = 1.5e0;
    apxr = 70.0e0;
    pmlv = 2.0e0;
    apv = 100.0e0;
    pmlz = 1.5e0;
    apz = 70.0e0;
    cc = 1.0e0;
    //Расчетная область задачи - AX модели = AX - pmlxl- pmlxr

    //MPI 
    //Открываем подель на каждом процессоре 
    
    ifstream FILE1;
    if (MYID == 0)
    {
       
        FILE2.open("report.log");
    }
    str = "PICT000.sct";
    int k1 = 9, k2 = 99;
    if (KART < 9)
        k1 = KART;
    for (int k = 0; k < k1; k++)
    {
        NPICT[k] = "PICT00" + to_string(k) + ".sct";
    }
    if (KART < 99)
        k2 = KART;
    for (int k = k1 + 1; k < k2; k++)
    {
        NPICT[k] = "PICT0" + to_string(k) + ".sct";
    }
    for (int k = k2 + 1; k < KART; k++)
    {
        NPICT[k] = "PICT" + to_string(k) + ".sct";
    }
    int j;
    for (int i = 1; i <= NPROCS; i++)
        j = i - 1;
    if (MYID == j)
    {
        std::cout << "I am " << MYID << " processor read the model"<<endl;
        if (MYID == 0)
            FILE2 << "I am " << MYID << " processor read the model"<<endl;
        
        FILE1.open("model8HZw2.dat");
        if (!FILE1)
        {
            std::cout << "File not open"<<endl;
            return -1;
        }
        FILE1 >> KOLIST >> NPI >> KTRI >> DTI >> F0;
        for (int J = 0; J < KOLIST; J++)
            FILE1 >> XI[J] >> ZI[J];
        FIST.resize(KTRI); //выделяем память массиву
        //проверять на выделение памяти не надо, автоматическое управление памятью
        int II = int((KTRI - 1) / 10);
        KT = KTRI;
        for (int JJ = 1; JJ <= II; JJ++)
        {
            int KK = (JJ - 1) * 10;;
            for (int K = KK; K < KK + 10; K++)
            {
                FILE1 >> FIST[K];
            }
            KT -= 10;
        }
        for (int K = KTRI - KT; K < KTRI; K++)
            FILE1 >> FIST[K];
        FILE1 >> KOLPR >> KTRP >> TMIN >> DTP >> DXP >> DZP >> NKOMP;
        for (int J = 0; J < KOLIST; J++)
            FILE1 >> XP[J] >> ZP[J];
        FILE1 >> AX1 >> AZ1 >> NX >> NZ >> KTR1;
        //KTR1=0
        KTR2 = 150;
        AX = AX1 * 1e-3;
        AZ = AZ1 * 1e-3;
        dz = AZ / (NZ - 1);
        dx = AX / (NX - 1);
        //dx=AX/NX;
        kgran.resize(NZ);
        fro.resize(NZ, vector<double>(NX));
        fmu.resize(NZ, vector<double>(NX));
        fla.resize(NZ, vector<double>(NX));
        if (MYID == 0)
            FILE2 << "Allocates array Fro,Fmu,Fla"<<endl;
        zro1.resize(NZ);
        zro2.resize(NZ);
        zmu1.resize(NZ);
        zmu2.resize(NZ);
        zla1.resize(NZ);
        zla2.resize(NZ);
        if (MYID == 0)
            FILE2 << " Allocates array zro,zmu,zla"<<endl;
        //----------------- NEODN PO X ------------================
        for (int IN = 0; IN < NZ; IN++)//40
        {
            FILE1 >> NN >> KX;
            kgran[IN] = KX;
            /*c      if(MYID.eq.0)  write(2,39) NN,KX
            c      if(MYID.eq.0)  print 39,NN,KX
            c 39    format(' NN =',I9,'  KX =',I9)*/
            if (NN != IN+1)
            {
                if (MYID == 0)
                    FILE2 << "ERROR !!! NN.ne.IN "<<endl;
                FILE2.close();
                std::cout << "ERROR !!! NN.ne.IN" << std::endl;
                return 1;
            }
            FILE1 >> XNX >> VP >> VS >> RO;
            AA = 0.0e0;
            BB = XNX * 1e-3;
            xro = 1.0e0 / RO;
            xmu = 1e-6 * VS * VS * RO;
            xla = 1e-6 * (VP * VP - 2 * VS * VS) * RO;

            zro1[IN] = xro;
            zmu1[IN] = xmu;
            zla1[IN] = xla;

            fro[IN][0] = xro;
            fmu[IN][0] = xmu;
            fla[IN][0] = xla;
            Kdx = BB / dx + 1.000001;
            for (int J = 1; J < Kdx; J++)
            {
                fro[IN][J] = xro;
                fmu[IN][J] = xmu;
                fla[IN][J] = xla;
            }
            k1 = Kdx;
            for (int K = 2; K <= KX; K++)//30
            {
                xx = (k1 - 1) * dx;
                if ((xx - BB) >= 0.0)
                    k1 -= 1;
                FILE1 >> XNX >> VP >> VS >> RO;
                AA = (k1 - 1) * dx;
                BB = XNX * 1e-3;
                if (K == KX)
                    BB = AX;
                if (BB - AA <= 0) {
                    FILE2 << "ERROR !!! AA > BB" << std::endl;
                    FILE2.close();
                    return 1;
                }
                xro = 1.0e0 / RO;
                xmu = 1e-6 * VS * VS * RO;
                xla = 1e-6 * (VP * VP - 2 * VS * VS) * RO;

                Kdx = (BB - AA) / dx + 0.000001;
                if (k1 + Kdx > NX) {
                    Kdx = Kdx - 1;
                }
                for (int J = 0; J < Kdx; J++)
                {
                    fro[IN][k1 + J] = xro;
                    fmu[IN][k1 + J] = xmu;
                    fla[IN][k1 + J] = xla;
                }
                k1 = k1 + Kdx;
            }//30
            zro2[IN] = xro;
            zmu2[IN] = xmu;
            zla2[IN] = xla;

            fro[IN][NX-1] = xro;
            fmu[IN][NX-1] = xmu;
            fla[IN][NX-1] = xla;
        }
        FILE1.close();
        std::cout << "I am " << MYID << " prosessor finish read the model"<<endl;
        if(MYID==0)
            FILE2<<"I am " << MYID << " prosessor finish read the model"<<endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    /*c	if(MYID.eq.0) then

    c      OPEN(3,FILE='fmu.dat')
    c       DO I=1,NZ
    c         WRITE(3,77) i,(fmu(i,j),j=1,11)
    c	 ENDDO
    c 77   format(I5,11F7.3)
    c      CLOSE (3)

    c      OPEN(3,FILE='zmu1.dat')
    c       DO I=1,NZ
    c         WRITE(3,72) i,zmu1(i),zmu2(i)
    c	 ENDDO
    c 72   format(I5,2F9.5)
    c      CLOSE (3)

    c	end if

    c	STOP

    *____________________________________________
    *	координаты источника

    c       Z0=5.0d0
    c	  X0=10.50d0   */

    Z0 = ZI[0] * 1e-3;
    X0 = XI[0] * 1e-3;
    //*____________________________________________
        //* для отрисовки трассы на линии z = Ztras

    Ztras = ZP[1] * 1e-3;
        //c	  Ztras = 3.0d0
        //* ____________________________________________
        //* шаг по х и количество трасс
        //* X(max) = X0tras + DXtras * (NTRAS - 1)

    NTRAS = KOLPR;
    X0tras = XP[1] * 1e-3;

        //c        NTRAS = 16
       // c        X0tras = 0.0d0
        //c        DXtras = AX / dfloat(NTRAS - 1)
    DXtras = DXP * 1e-3;
        //* ____________________________________________
        //* шаг по времени в трассах

    t1 = 1e0 * TMIN;
    dtras = 1e0 * DTP;

        //c        t1 = 0.0d0
        //c        dtras = 0.018d0
        //* ____________________________________________
        //* количество отсчетов в трассе по времени
        //* T(max) = t1 + KTRAS * dtras

    KTRAS = KTRP;
        //c        KTRAS = 1000

        //* ____________________________________________
        //* прорежование snapshoot по Z и X

    Kdz = 2;
    Kdx = 2;
        //* ____________________________________________
        //* количество точек в snapshoot по Z и X
        //* если NZ четное, то
        //* KTR = NZ / Kdz
        //* если NZ не четное, то
        //* KTR = (NZ - 1) / Kdz

    KTR = (NZ - KTR1 - KTR2 - 1) / Kdz;
    NTR = (NX - 1) / Kdx;
    NTR0 = 0;
        //* ____________________________________________
    LIZ = Z0 / dz + 1.000001;
    if (LIZ >= NZ)
        LIZ = NZ - 1;

    LIX = X0 / dx + 1.000001;
    if(LIX > NX)
        LIX = NX;

    LPZ = Ztras / dz + 1.000001;
    if (LPZ >= NZ)
        LPZ = NZ - 1;

    LPX = X0tras / dx + 1.000001;
    if (LPX > NX)
        LPX = NX;

    LPdX = DXtras / dx + 0.000001;
    if (LPdX >= NX)
        LPdX = NX - 1;

    if (MYID == 0)
    {
        if (nout == 2)
        {
            FILE2 << " dz= " << scientific << setprecision(5) << dz
                << " dx= " << scientific << setprecision(5) << dx << endl;
            // Запись значений LIZ, LIX, LPZ, LPX, LPdX в формате '   LIZ = <value>   LIX = <value>   LPZ = <value>   LPX = <value>   LPdX = <value>'
            FILE2 << "   LIZ = " << std::setw(5) << LIZ
                << "   LIX = " << std::setw(5) << LIX
                << "   LPZ = " << std::setw(5) << LPZ
                << "   LPX = " << std::setw(5) << LPX
                << "   LPdX = " << std::setw(5) << LPdX << endl;

            // Запись значений NTR, KTR, KTR1, KTR2 в формате '   NTR = <value>   KTR = <value>   KTR1 = <value>   KTR2 = <value>'
            FILE2 << "   NTR = " << std::setw(5) << NTR
                << "   KTR = " << std::setw(5) << KTR
                << "   KTR1 = " << std::setw(5) << KTR1
                << "   KTR2 = " << std::setw(5) << KTR2 << endl;
        }
        else
            std::cout << "492 строка, не понимаю куда записать параметры"<<endl;
    }
    //*____________________________________________
    //    * времена снепов :
    if (MYID == 0)
        std::cout << " t1 =" << t1<<endl;
    dd = t1 + KTRAS * dtras;
    TN[KART] = dd;
    dd = dd / KART;

    for (int I = 0; I < KART; I++)
        TN[I] = dd + (I - 1) * dd;
    //C----------------------------------------------------------
    //C   Определение коэфф.Лагерра для сигнала f(t) в источнике
    //C----------------------------------------------------------
    //* ash - parametr sdviga Laguerra
    ash = 150.0e0 * F0;
    ash2 = ash / 2.0e0;
    RE = t1 + (KTRAS - 1) * dtras;
    for (int k = 0; k < KART; k++)
        if (RE <= TN[k])
            RE = TN[k];
    LAGmax = long long (RE * F0 * ash / 25);
    //c	Lagmax=5000
    if (MYID == 0)
    {
        std::cout<< " ash = " << fixed << setprecision(1) << ash
            << "\n  for Signal = " << fixed << setprecision(1) << F0 << " HZ" <<endl;
        std::cout << "for ash = " << fixed << setprecision(1) << ash
            << " max numbers garmonic Laguerre = " << setw(6) << LAGmax << endl;
        if (nout == 2)
        {
            FILE2<< " ash = " << fixed << setprecision(1) << ash
                << "  for Signal = " << fixed << setprecision(1) << F0 << " HZ" << endl;
            FILE2 << " for ash = " << fixed << setprecision(1) << ash
                << "\n max numbers garmonic Laguerre = " << setw(6) << LAGmax << endl;
        }
        else
            std::cout << "523 строка, не понимаю куда записать параметры"<<endl;
    }
    FN.resize(LAGmax);
    CR.resize(LAGmax);
    QN.resize(LAGmax);
    for (int I = 0; I < LAGmax; I++)
        FN[I] = 0.0e0;
    int ISTSIGN = 0;
    if (ISTSIGN > 0)
    {
        dd = 1.0e-5 / F0;
        if (KTRI == 0)
        {
            if (MYID == 0)
                if (nout == 2)
                {
                    FILE2 << " ERROR !!!Number point signal of source KTRI = " << KTRI<<endl;
                    return 1;
                }
                else
                    std::cout << "nout!=2"<<endl;
        }
        int nxs = int(DTI / dd + 1e-10);
        if (nxs < 1) nxs = 1;

        RE = 0.5e0 / F0;
        int KT0 = int(RE / DTI) + 1;

        for (int k = KT0 - 1; k < KTRI - 1; k++)//338
        {
            tt = (k - 1) * DTI;
            up1 = (FIST[k + 1] - FIST[k]) / nxs;
            for (int i = 1; i < nxs; i++)
            {
                xx = tt + dd * (i - 1);
                RE = ash * xx;
                DLAG2(CR, RE, LAGmax, LS);//654
                RE = FIST[k] + up1 * (i - 1);
                for (int j = 1; j <= LAGmax; j++)
                {
                    RE1 = 1.0e0;
                    for (int l = 1; l <= LS; l++)
                        RE1 = RE1 * static_cast<double>(j + l - 1);
                    FN[j - 1] = FN[j - 1] + CR[j - 1] * RE / RE1;
                }
            }
        }
        //667
    }
    else
    {
       // dd = 1.0e-5 / F0;
        dd = 1.0e-6;
        RE = 0.5e0 / F0;
        nxs1 = int(RE / dd + 1e-10) + 1;
        RE1 = 3.0e0 / F0;
        nxs2 = int(RE1 / dd + 1e-10) + 1;
        
        t0 = 1.5e0 / F0;
        s1 = 2.0e0 * PII * F0;
        s2 = (s1 / 4.0) * (s1 / 4.0);
        for (int i = nxs1; i <= nxs2; i++)
        {
            xx = dd * static_cast<double>(i - 1);
            RE = ash * xx;
            DLAG2(CR, RE, LAGmax, LS);
            
            RE = exp(-s2 * (xx - t0) * (xx - t0)) * sin(s1 * xx);
            for (int j = 1; j <= LAGmax; j++) {
                RE1 = 1.0;
                for (int l = 1; l <= LS; l++) {
                    RE1 *= static_cast<double>(j + l - 1);
                }
                FN[j - 1] += CR[j - 1] * RE / RE1;
            }
        }
    }
    for (int j = 0; j < LAGmax; j++)
        FN[j] = FN[j] * ash * dd;

    //C------------------------------------------------------
    //   c      if (MYID.eq.0) then
    //   c       OPEN(10, FILE = 'ft.dat')
    //  c        DO  IT = 1, 100
    //   c         TT = (IT - 1) * 0.002
    //   c         FIST = dexp(-s2 * (tt - t0) * *2) * dsin(s1 * tt)
    //    c        write(10, 121) TT, FIST
    //   c        END DO
    //    c       CLOSE(10)
    //    c	endif

    FIST.clear();

    /*c      if (MYID.eq.0) then
        c       OPEN(10, FILE = 'fn_al.dat')
        c        do j = 1, LAGmax
        c        write(10, 121) double(j - 1), fn(j)
        c        enddo
        c       CLOSE(10)
        c	endif

        C------------------------------------------------------
        C   Пересчет коэфф.для f(t) с учетом сдвига на t = TN(I)
        C------------------------------------------------------*/

    for (int IT = 0; IT < KART + 1; IT++)
    {
        RE = ash * TN[IT];
        int KER;
        KOEFlag(FN, CR, QN, RE, LAGmax, LS, epsLag, NW, KER);//734
        NL[IT] = NW;
        if (MYID == 0)
        {
            std::cout << " t = " << TN[IT] << " NW = " << NW<<endl;
            FILE2 << " t = " << TN[IT] << " NW = " << NL[IT] << endl;
        }
        if (KER == 0)
            continue;
        if (MYID == 0)
        {
            std::cout << " WARNING !!! ERROR > epsLag =" << epsLag<<endl;
            FILE2 << " WARNING !!! ERROR > epsLag =" << epsLag << endl;
        }
    }
    NW = 0;
    for (int I = 0; I < KART + 1; I++)
        if (NL[I] > NW) NW = NL[I];
    if (MYID == 0)
    {
        std::cout << " max garmonics Laguerre =" << NW<<endl;
        FILE2 << " max garmonics Laguerre =" << NW<<endl;
    }
    /**________________________________________________

        c	if (MYID.eq.0) then
        c       OPEN(1, FILE = 'qn_al50.dat')
        c        do j = 1, NW
        c        write(1, 12) double(j - 1), qn(j)
        c        enddo
        c       CLOSE(1)
        c 12    format(F5.0, E13.5)
        c	end if

        C------------------------------------------------------------C
        C--  Проверка точности востановления сигнала  f(t)-- - C
        C------------------------------------------------------------C

        c	if (MYID.eq.0) then
        c      OPEN(1, FILE = 'fta1.dat')

        c        DO 777 IT = 1, KTras + 100
        c        TT = t1 + (IT - 1) * dtras
        c        RE = Ash * TT

        c        CALL DLAG1(CR, RE, LAGMAX, LS)

        c        RE1 = 0.0D0
        c          DO L = 1, NW
        c          RE1 = RE1 + QN(L) * CR(L)
        c          ENDDO

        c      write(1, 122) TT, RE1
        c 777    CONTINUE
        c      CLOSE(1)
        c 122     format(F7.3, E13.5)
        c	end if
        * ______________________________________________
        c      DEALLOCATE(FN, QN)             !Deallocates array FN, QN*/
    CRN.resize(KART, vector<double>(NW));
    for (int I = 0; I < KART; I++)
    {
        RE = ash * TN[I];
        DLAG1(CR, RE, NW, LS);
        for (int J = 0; J < NW; J++)
            CRN[I][J] = CR[J];
    }
    /**_______________________________________________
        * задание затухающей pml - функции по z :
    *d2pmlz = 1 + dpmlz * 2 / ash*/
    d2pmlz.resize(NZ);
    Lpmlv = pmlv / dz + 1.000001;
    if (Lpmlv < 1) Lpmlv = 1;
    Lpmlz = (AZ - pmlz) / dz + 1.000001;
    if (Lpmlz < 1) Lpmlz = 1;
    if (Lpmlz < Lpmlv) Lpmlz = Lpmlv + 1;
    if (Lpmlz>NZ) Lpmlz = NZ;
    for (int J = Lpmlv; J <= Lpmlz; J++)
        d2pmlz[J - 1] = 1.0e0;

    /**--------------------------------------------
        * функция затухания в PML зоне по z сверху
        * --------------------------------------------*/
    if(Lpmlv!=1)
        for (int J = 1; J <= Lpmlv - 1; J++)
        {
            dd = apv * pow(static_cast<double>(Lpmlv - J) * dz / pmlv, 3);
            d2pmlz[J - 1] = 1.0e0 + 2.0e0 * dd / ash;
        }
    /**--------------------------------------------
        * функция затухания в нижней PML зоне по z
        * --------------------------------------------*/
    RE = apz / (exp(cc * pmlz) - 1.0e0);
    for (int J = Lpmlz + 1; J <= NZ; J++)
    {
        dd = RE * (exp(cc * (static_cast<double>(J - Lpmlz) * dz)) - 1.0e0);
        d2pmlz[J - 1] = 1.0e0 + 2.0e0 * dd / ash;
    }
    /*c	if (MYID.eq.0) then
        c       OPEN(1, FILE = 'PMLZ_v.dat')
        c        do j = 1, Lpmlv
        c        write(1, 56) double(j), d2pmlz(j)
        c        enddo
        c       CLOSE(1)
        c	end if

        c	if (MYID.eq.0) then
        c       OPEN(1, FILE = 'PMLZ_n.dat')
        c        do j = Lpmlz, NZ
        c        write(1, 56) double(j), d2pmlz(j)
        c        enddo
        c       CLOSE(1)
        c	end if
        * _______________________________________________
        * задание затухающей pml - функции по x :
    *d2pmlx = 1 + dpmlx * 2 / ash*/

    d2pmlx.resize(NX);
    Lpmlx1 = pmlxl / dx + 1.000001;
    Lpmlx2 = (AX - pmlxr) / dx + 1.000001;
    if (Lpmlx1<1) Lpmlx1 = 1;
    if (Lpmlx2>NX) Lpmlx2 = NX;
    if (Lpmlx1 > Lpmlx2)
    {
        if (MYID == 0)
        {
            std::cout << " WARNING !!! pmlxl+pmlxr > AX " << endl;
            FILE2 << " WARNING !!! pmlxl+pmlxr > AX " << endl;
        }
        Lpmlx1 = (NX + 1) / 2;
        Lpmlx2 = NX - Lpmlx1 + 1;
    }
    for (int j = Lpmlx1; j <= Lpmlx2; j++)
        d2pmlx[j - 1] = 1.0e0;
    if (Lpmlx1 != 1)
    {
        RE = apxl / (exp(cc * pmlxl) - 1.0e0);
        for (int j = 1; j <= Lpmlx1 - 1; j++)
        {
            /**---------------------------------------------- -
                *функция затухания в левой PML зоне по x*/
            dd = RE * (exp(cc * (static_cast<double>(Lpmlx1 - j) * dx)) - 1.0e0);
            /**---------------------------------------------- -*/
            d2pmlx[j - 1] = 1.0e0 + 2.0e0 * dd / ash;
        }
    }
    RE = apxr / (exp(cc * pmlxr) - 1.0e0);
    for (int j = Lpmlx2 + 1; j <= NX; j++)
    {
        /**---------------------------------------------- -
                *функция затухания в левой PML зоне по x*/
        dd = RE * (exp(cc * (static_cast<double>(j-Lpmlx2) * dx)) - 1.0e0);
        /**---------------------------------------------- -*/
        d2pmlx[j - 1] = 1.0e0 + 2.0e0 * dd / ash;
    }
    /*c	if (MYID.eq.0) then
        c       OPEN(1, FILE = 'PMLX_l.dat')
        c        do j = 1, Lpmlx1
        c        write(1, 56) double(j), d2pmlx(j)
        c        enddo
        c       CLOSE(1)*/
    /*c       OPEN(1, FILE = 'PMLX_r.dat')
        c        do j = Lpmlx2, NX
        c        write(1, 56) double(j), d2pmlx(j)
        c        enddo
        c       CLOSE(1)
        c	end if

        * __________________________________________________*/
    NXNZ = NX * NZ;
    NZ5 = 5 * NZ - 4;
    NZ10 = NZ5 + NZ5;
    Ntotal = NX * NZ5;
    N = Ntotal;

    PART(NX, NPROCS,RANGES);

    J0 = RANGES[0][MYID];
    J1 = RANGES[2][MYID];
    NEL = RANGES[1][MYID];

    LOCLEN = NEL * NZ5;
    LOCS = LOCLEN + 2 * NZ10;
    /*c      LDWRK = (5 + 2 * BASIS) * LOCLEN + 2 * BASIS*/
    LDWRK = BASIS * LOCLEN;
    C = BASIS;
    V = BASIS1;
    MAXIT = NZ10;

    if (MYID == 0)
        std::cout << "Number size vector x -- Ntotal =" << Ntotal<<endl;
    std::cout << "I am" << MYID << " processor calculate from" << J0 << " to " << J1<<endl;
    std::cout << " Local size vector x -- LOCLEN =" << LOCLEN<<endl;
    if (MYID == 0)
    {
        FILE2 << "Number size vector x -- Ntotal =" << Ntotal << endl;
        FILE2 << "I am" << MYID << " processor calculate from" << J0 << " to " << J1 << endl;
        FILE2 << " Local size vector x -- LOCLEN =" << LOCLEN << endl;
    }

    if (LOCLEN <= NZ10)
    {
        std::cout << " prozessor " << MYID << " ERROR !!! LOCLEN < NZ10 "<<endl;
        return 1;
    }
    int i;
    for (i = 1; i < NPROCS; i++)
    {
        int ix = RANGES[2][i];
        if (LIX <= ix)
        {
            break;
        }
        if (i == NPROCS - 1)
        {
            std::cout << " ERROR !!! LIX > NX "<<endl;
            return 1;
        }
    }
    NIP = i - 1;
    if (MYID == 0)
        std::cout << " NIP=" << NIP << "   NPI" << NPI<<endl;
    /**__________________________________________________

        c        заполнение матрицы коэффициентов для разностной схемы
        c        в подпрграмме МАТVEC :*/
    UZ1.resize(LOCLEN);
    UZ2.resize(LOCLEN);
    UX1.resize(LOCLEN);
    UX2.resize(LOCLEN);
    RVP2.resize(NXNZ);
    RVS2.resize(NXNZ);
    ROD1.resize(NXNZ);
    dz1_24 = 1.0e0 / (24.0e0 * dz);
    dz9_8 = 9.0e0 / (8.0e0 * dz);
    dx1_24 = 1.0e0 / (24.0e0 * dx);
    dx9_8 = 9.0e0 / (8.0e0 * dx);
    //1029
    for (int ix = 1; ix <= NX; ix++)
    {
        ind = (ix - 1) * NZ;
        for (int iz = 1; iz <= NZ; iz++)
        {
            ROD1[ind + iz-1] = fro[iz-1][ix-1];
            RVS2[ind + iz-1] = fmu[iz-1][ix-1];
            RVP2[ind + iz-1] = fla[iz-1][ix-1]+2e0*fmu[iz-1][ix-1];
        }
    }
    MATRIZA(J0, J1, d2pmlz, d2pmlx, UX1, UX2, UZ1, UZ2, fro, fmu, fla,NX,NZ,NZ5,Ntotal);
    fro.clear();
    fmu.clear();
    fla.clear();
    kgran.clear();
    zro1.clear();
    zmu1.clear();
    zla1.clear();
    zro2.clear();
    zmu2.clear();
    zla2.clear();
    Q1.resize(LOCLEN);
    Q2.resize(LOCLEN);
    X=new double[LOCLEN];
    FWW=new double[LOCLEN];
    SUMLOC.resize(LOCLEN);
    SUMMApml.resize(LOCLEN);
    SP=new double[NZ10];
    SL=new double[NZ10];
    XS.resize(LOCS);
    DWRK=new double[LDWRK];
    for (int i = 0; i < LOCLEN; i++)
    {
        FWW[i] = 0.0e0;
        SUMLOC[i] = 0.0e0;
        SUMMApml[i] = 0.0e0;
    }
    if (PRET == 1)
    {
        for (int ix = J0; ix <= J1; ix++)
        {
            ind = (ix - J0) * NZ5;
            dd = d2pmlx[ix-1];
            for (int iz = 0; iz < NZ5; iz++)
            {
                Q1[ind + iz] = 2.0e0 / (ash * dd);
            }
        }
    }
    else if (PRET == 2)
    {
        for (int ix = J0; ix <= J1; ix++)
        {
            ind = (ix - J0) * NZ5;
            dd = d2pmlx[ix - 1];
            for (int iz = 0; iz < NZ5; iz++)
            {
                Q2[ind + iz] = 2.0e0 / (ash * dd);
            }
        }
    }
    else if (PRET == 3)
    {
        for (int ix = J0; ix <= J1; ix++)
        {
            ind = (ix - J0) * NZ5;
            dd = d2pmlx[ix - 1];
            for (int iz = 0; iz < NZ5; iz++)
            {
                Q1[ind + iz] = sqrt(2.0e0 / (ash * dd));
                Q2[ind + iz] = Q1[ind+iz];
            }
        }
    }
    PIMDSETPAR(IPAR, DPAR, N, N, LOCLEN, LOCLEN, C, NPROCS,
        MYID,PRET, STOPT, MAXIT, TOL);
    /**-------------------------------------------------- -

        *!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    UXDZ.resize(NXNZ);
    UZDZ.resize(NXNZ);
    SXDZ.resize(NXNZ);
    SZDZ.resize(NXNZ);
    TRASX1= new double[NTRAS];
    TRASX2 = new double[NTRAS];
    TRASZ1 = new double[NTRAS];
    TRASZ2 = new double[NTRAS];
    TRASP1 = new double[NTRAS];
    TRASP2 = new double[NTRAS];
    PICTZ1 = new double[NX];
    PICTZ2 = new double[NX];
    for (int i = 0; i < NX; i++)
    {
        PICTZ1[i] = 0.0;
        PICTZ2[i] = 0.0;
    }
    UTRLZ.resize(NEL, vector<vector<double>>(KTR, vector<double>(KART,0.0)));
    if (MYID == 0)
    {
        UPRXX.resize(NTRAS, vector<double>(NW));
        UPRZZ.resize(NTRAS, vector<double>(NW));
        UPRP.resize(NTRAS, vector<double>(NW));
    }
    //*___________________________________________________

    int IW1 = 0;
    NUMPICT = 1;

    if (IPARAM == 1)
    {
        if (MYID == 0)
        {
            std::cout << "BEGIN__READ__PROMEZYT__DATA___ZAPIS" << endl;
            FILE2 << "BEGIN__READ__PROMEZYT__DATA___ZAPIS" << endl;
            RW1.resize(NTRAS);
            delete[] RW2;
            RW2=new double[NTR];
            UTRAZ1.resize(NTR, vector<double>(KTR));
            REG.resize(KTR);
        }
        SUMMA=new double[Ntotal];
        RWSUM = new double[Ntotal];
        for (int i = 0; i < Ntotal; i++)
        {
            SUMMA[i] = 0.0;
            RWSUM[i] = 0.0;
        }
        if (MYID == 0)
        {
            ifstream FILE10("kolwo001.tmp");
            int NW1;//заплатка для 1240
            FILE10 >> KOL >> NW1 >> NUMPICT;
            FILE10 >> IW1 >> NW1 >> NUMPICT;
            if (KOL != IW1)
            {
                FILE2 << " ERROR!!! KOL  NE  IW1 " << endl;
                FILE2.close();
                std::cout << " ERROR!!! KOL  NE  IW1 " << endl;
                return 1;
            }
            LAG[0] = IW1;
            LAG[1] = NW1;
            LAG[2] = NUMPICT;
            FILE10.close();
            MPI_Bcast(LAG, 3, MPI_INTEGER, 0, MPI_COMM_WORLD);
            IW1=LAG[0];
            //NW=LAG[1];
            NUMPICT = LAG[2];
            FILE10.open("SUMMA001.tmp",ios::binary);
            double val;
            for (int i = 0; FILE10.read(reinterpret_cast<char*>(&val), sizeof(val)); i++)
                SUMMA[i] = val;
            FILE10.close();
            FILE10.open("SUMMApml.tmp", ios::binary);
            for (int i = 0; FILE10.read(reinterpret_cast<char*>(&val), sizeof(val)); i++)
                RWSUM[i] = val; 
            FILE10.close();
        }
        MPI_Bcast(SUMMA, Ntotal, MPI_REAL, 0, MPI_COMM_WORLD);
        for (int i = 0; i < LOCLEN; i++)
        {
            ind = (J0 - 1) * NZ5;
            SUMLOC[i] = 1e0 * SUMMA[ind + i];
        }
        MPI_Bcast(RWSUM, Ntotal, MPI_REAL, 0, MPI_COMM_WORLD);
        for (int i = 0; i < LOCLEN; i++)
        {
            ind = (J0 - 1) * NZ5;
            SUMMApml[i] = 1e0 * SUMMA[ind + i];
        }
        //*_______________________________________
        if (MYID == 0)
        {
            ifstream FILE20("UPRX001.tmp", ios::binary);
            for (int J = 0; J < IW1; J++)
            {
                double val;
                for (int i = 0; FILE20.read(reinterpret_cast<char*>(&val), sizeof(val)); i++)
                    RW1[i] = val;
                for (int I = 0; I < NTRAS; I++)
                    UPRXX[I][J] = 1e0 * RW1[I];
            }
            FILE20.close();
            FILE20.open("UPRZ001.tmp", ios::binary);
            for (int J = 0; J < IW1; J++)
            {
                double val;
                for (int i = 0; FILE20.read(reinterpret_cast<char*>(&val), sizeof(val)); i++)
                    RW1[i] = val;
                for (int I = 0; I < NTRAS; I++)
                    UPRZZ[I][J] = 1e0 * RW1[I];
            }
            FILE20.close();
            FILE20.open("UPRP001.tmp", ios::binary);
            for (int J = 0; J < IW1; J++)
            {
                double val;
                for (int i = 0; FILE20.read(reinterpret_cast<char*>(&val), sizeof(val)); i++)
                    RW1[i] = val;
                for (int I = 0; I < NTRAS; I++)
                    UPRP[I][J] = 1e0 * RW1[I];
            }
            FILE20.close();
        }
        //*_________________________________________________
        for (int K = NUMPICT; K <= NUMPICT - 1; K++)
        {
            if (MYID == 0)
            {
                std::cout << "------------------------------------------" << endl;
                std::cout << " read of NUMPICT = " << NUMPICT << endl;
                FILE2 << "------------------------------------------" << endl;
                FILE2 << " read of NUMPICT = " << NUMPICT << endl;
            }
            for (int I = 0; I < KTR; I++)
            {
                if (MYID == 0)
                {
                    for (int J = 0; J < NTR; J++)
                        RW2[J] = UTRAZ1[J][I];
                }
                else
                {
                    for (int J = 0; J < NTR; J++)
                        RW2[J] = 0.0;
                }
                MPI_Bcast(RW2, NTR, MPI_REAL, 0, MPI_COMM_WORLD);
                for (int J = 1; J <= NTR; J++)
                {
                    ind = NTR0 + (J - 1) * Kdx + 1;
                    PICTZ1[ind-1] = RW2[J-1];
                }
                for (int J = 0; J <= NEL; J++)
                    UTRLZ[J][I][K] = PICTZ1[J0 + J-1];
            }
        }
        //______________________________________________________
        if (MYID == 0)
        {
            std::cout << "FINISH__READ__PROMEZYT__DATA___ZAPIS " << endl;
            UTRAZ1.clear();
            REG.clear();
            RW1.clear();
            delete[] RW2;
        }
        delete[] SUMMA;
        delete[] RWSUM;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    KOL = 0;
    ETMV0 = 0.0;
    //___________  цикл по t  ______________________________________
    for (int it = IW1 + 1; it <= NW; it++)
    {
        KOL++;
        if (MYID == 0)
        {
            std::cout << "\n Laguerre garmonic it = "<<it << endl;
            FILE2 << "\n Laguerre garmonic it = " <<it<< endl;
        }
        //___________  вводим правую часть _____________________________
        RE = 1.0e0;
        for (int l = 1; l <= LS; l++)
            RE = RE * static_cast<double>(it + l - 1);
        for (int ix = J0; ix <= J1; ix++)
        {
            ind = (ix - J0) * NZ5;
            dd = d2pmlx[ix-1];
            if (d2pmlx[0] > 0e0)
            {
                dd = (d2pmlx[ix-1] - d2pmlz[0]) / d2pmlz[0];
            }
            FWW[ind] = -ash * (SUMLOC[ind] + dd * SUMMApml[ind]) / RE;
            FWW[ind + 1] = -ash * (SUMLOC[ind + 1] + dd * SUMMApml[ind + 1]) / RE;
            FWW[ind + 2] = -ash * (SUMLOC[ind + 2] + dd * SUMMApml[ind + 2]) / RE;
            FWW[ind + 3] = -ash * (SUMLOC[ind + 3] + dd * SUMMApml[ind + 3]) / RE;
            for (int iz = 2; iz <= NZ - 1; iz++)
            {
                ind += 5;
                dd = d2pmlx[ix-1];
                if (d2pmlx[iz-1] > 0e0)
                    dd = (d2pmlx[ix-1] - d2pmlz[iz-1]) / d2pmlz[iz-1];
                FWW[ind - 1] = -ash * (SUMLOC[ind - 1] + dd * SUMMApml[ind - 1]) / RE;
                FWW[ind] = -ash * (SUMLOC[ind] + dd * SUMMApml[ind]) / RE;
                FWW[ind + 1] = -ash * (SUMLOC[ind + 1] + dd * SUMMApml[ind + 1]) / RE;
                FWW[ind + 2] = -ash * (SUMLOC[ind + 2] + dd * SUMMApml[ind + 2]) / RE;
                FWW[ind + 3] = -ash * (SUMLOC[ind + 3] + dd * SUMMApml[ind + 3]) / RE;
            }
            dd = d2pmlx[ix-1];
            dd = (d2pmlx[ix-1] - d2pmlz[NZ-1]) / d2pmlz[NZ-1];
            FWW[ind + 4] = -ash * (SUMLOC[ind + 4] + dd * SUMMApml[ind + 4]) / RE;
            FWW[ind + 5] = -ash * (SUMLOC[ind + 5] + dd * SUMMApml[ind + 5]) / RE;
        }
        if (MYID + 1 == NPROCS)
            FWW[(NEL - 1) * NZ5] = 0e0;
        /**___________ правая часть при iz = LIZ, ix = LIX______________

            * ___________________ обозначения_________________
            c	sizz : k1 = (ix - 1) * NZ5 + (iz - 1) * 5 + 0
            c	sixx : k2 = (ix - 1) * NZ5 + (iz - 1) * 5 + 1
            c	ux : k3 = (ix - 1) * NZ5 + (iz - 1) * 5 + 2
            c	uz : k4 = (ix - 1) * NZ5 + (iz - 1) * 5 + 3
            c	sixz : k5 = (ix - 1) * NZ5 + (iz - 1) * 5 + 4

            * _________________________________________________*/
        if (MYID == NIP)
        {
            //*      вертикальная сила:
            if (NPI == 1)
            {
                int k3 = (LIX - J0) * NZ5 + (LIZ - 1) * 5 + 3;
                FWW[k3-1] = FWW[k3-1] + 1.0e0 * FN[it-1] / (dz * dx);
            }
            //*      центр давления:
            if (NPI != 1)
            {
                int k1 = (LIX - J0) * NZ5 + (LIZ - 1) * 5 + 5;
                int k2 = (LIX - J0) * NZ5 + (LIZ - 1) * 5 + 6;

                FWW[k1] = FWW[k1] + 1.0e0 * FN[it] / (dz * dx);
                FWW[k2] = FWW[k2] + 1.0e0 * FN[it] / (dz * dx);
            }
        }
        /**______________________________________

            c   вектор Х полагается равным 0.0D0 для начала итераций :*/
        DINIT(LOCLEN, 0.0e0, X, 1);
        ET0 = Timer();
        KMV = 0;
        ETMV1 = 0.0;
        PIMDCGS(KMV, ETMV1, d2pmlx, UZ1, UZ2, UX1, UX2, SP, SL, XS, Q1, Q2, X,
            FWW, DWRK, &IPAR[0], &DPAR[0]);

        ET1 = Timer();
        ET = ET1 - ET0;
        ETMV0 = ETMV0 + ETMV1;
        std::cout << " processor " << MYID << " calculation: number MATVEC =" << KMV << " time MATVEC" << ETMV1 << endl;
        if (MYID == 0)
            FILE2 << " processor " << MYID << " calculation: number MATVEC =" << KMV << " time MATVEC" << ETMV1 << endl;
        //c  сведения о сходимости:
        if (MYID == 0)
            REPORT("CGS", IPAR, DPAR, ET, X, nout);//2.11.12.13.14
        for (int i = 0; i < LOCLEN; i++)
            SUMLOC[i] = SUMLOC[i] + X[i] * RE;
        //c      вычисление сумм в PML зонах
        if (pmlv + pmlz + pmlxl + pmlxr > 0e0)
        {
            /*c      вычисление производных для правой части
                c     SUBROUTINE PROIZVOD(J0, J1, UXDZ, UZDZ, SXDZ, SZDZ, X, KMV, ETMV)*/
            int ETMV2=0;
            PROIZVOD(J0, J1, UXDZ, UZDZ, SXDZ, SZDZ, X, it, ETMV2);
            //c      вычисление суммы в верхней PML зоне
            if (Lpmlv > 1)
            {
                for (int ix = J0; ix < J1; ix++)
                {
                    ind = (ix - 1) * NZ;
                    int ind1 = (ix - J0) * NZ5;
                    i = ind + 1;
                    dd = ash2 * d2pmlz[0];

                    k2 = ind1 + 1;
                    xx = SUMMApml[k2 - 1];
                    SUMMApml[k2 - 1] = xx + (RE * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) * UZDZ[i - 1] - ash * xx) / dd;
                    int k3 = k2 + 1;
                    xx = SUMMApml[k3 - 1];
                    SUMMApml[k3 - 1] = xx + (RE * ROD1[i - 1] * SXDZ[i - 1] - ash * xx);
                    int k4 = k3 + 1;
                    xx = SUMMApml[k4 - 1];
                    SUMMApml[k4 - 1] = xx + (RE * ROD1[i - 1] * SZDZ[i - 1] - ash * xx) / dd;
                    int k5 = k4 + 1;
                    xx = SUMMApml[k5 - 1];
                    SUMMApml[k5 - 1] = xx + (RE * RVS2[i - 1] * UXDZ[i - 1] - ash * xx) / dd;
                    for (int iz = 2; iz <= Lpmlv - 1; iz++)
                    {
                        i = ind + iz;
                        dd = ash2 * d2pmlz[iz - 1];
                        k1 = ind1 + (iz - 1) * 5;
                        xx = SUMMApml[k1 - 1];
                        SUMMApml[k1 - 1] = xx + (RE * RVP2[i - 1] * UZDZ[i - 1] - ash * xx) / dd;
                        k2 = k1 + 1;
                        xx = SUMMApml[k2 - 1];
                        SUMMApml[k2 - 1] = xx + (RE * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) * UZDZ[i - 1] - ash * xx) / dd;
                        int k3 = k2 + 1;
                        xx = SUMMApml[k3 - 1];
                        SUMMApml[k3 - 1] = xx + (RE * ROD1[i - 1] * SXDZ[i - 1] - ash * xx);
                        int k4 = k3 + 1;
                        xx = SUMMApml[k4 - 1];
                        SUMMApml[k4 - 1] = xx + (RE * ROD1[i - 1] * SZDZ[i - 1] - ash * xx) / dd;
                        int k5 = k4 + 1;
                        xx = SUMMApml[k5 - 1];
                        SUMMApml[k5 - 1] = xx + (RE * RVS2[i - 1] * UXDZ[i - 1] - ash * xx) / dd;
                    }
                }
            }
            //c      вычисление суммы в нижней PML зоне
            //85 1549
            if (Lpmlz < NZ)
            {
                for (int ix = J0; ix < J1; ix++)
                {
                    ind = (ix - 1) * NZ;
                    int ind1 = (ix - J0) * NZ5;
                    for (int iz = Lpmlz + 1; iz <= NZ - 1; iz++)
                    {
                        i = ind + iz;
                        dd = ash2 * d2pmlz[iz - 1];
                        k1 = ind1 + (iz - 1) * 5;
                        xx = SUMMApml[k1 - 1];
                        SUMMApml[k1 - 1] = xx + (RE * RVP2[i - 1] * UZDZ[i - 1] - ash * xx) / dd;
                        k2 = k1 + 1;
                        xx = SUMMApml[k2 - 1];
                        SUMMApml[k2 - 1] = xx + (RE * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) * UZDZ[i - 1] - ash * xx) / dd;
                        int k3 = k2 + 1;
                        xx = SUMMApml[k3 - 1];
                        SUMMApml[k3 - 1] = xx + (RE * ROD1[i - 1] * SXDZ[i - 1] - ash * xx);
                        int k4 = k3 + 1;
                        xx = SUMMApml[k4 - 1];
                        SUMMApml[k4 - 1] = xx + (RE * ROD1[i - 1] * SZDZ[i - 1] - ash * xx) / dd;
                        int k5 = k4 + 1;
                        xx = SUMMApml[k5 - 1];
                        SUMMApml[k5 - 1] = xx + (RE * RVS2[i - 1] * UXDZ[i - 1] - ash * xx) / dd;
                    }
                    i = ind + NZ;
                    dd = ash2 * d2pmlz[NZ - 1];
                    k1 = ind1 +NZ5 - 1;
                    xx = SUMMApml[k1 - 1];
                    SUMMApml[k1 - 1] = xx + (RE * RVP2[i - 1] * UZDZ[i - 1] - ash * xx) / dd;
                    k2 = ind1 + NZ5;
                    xx = SUMMApml[k2 - 1];
                    SUMMApml[k2 - 1] = xx + (RE * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) * UZDZ[i - 1] - ash * xx) / dd;
                }
            }
            //c      вычисление суммы в правой PML зоне
            if (Lpmlx2 < J1)
            {
                int kk = Lpmlx2 + 1;
                if (kk < J0)
                {
                    kk = J0;
                }
                for (int ix = kk; ix <= J1; ix++)
                {
                    ind = (ix - 1) * NZ;
                    int ind1 = (ix - J0) * NZ5;

                    int i1 = Lpmlv;
                    if (Lpmlv <= 1)
                    {
                        i = ind + 1;
                        k2 = ind1 + 1;
                        xx = SUMMApml[k2 - 1];
                        SUMMApml[k2 - 1] = RE * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) * UZDZ[i - 1] / ash2 - xx;
                        int k3 = k2 + 1;
                        SUMMApml[k3 - 1] = RE * ROD1[i - 1] * SXDZ[i - 1] / ash2 - SUMMApml[k3 - 1];
                        int k4 = k3 - 1;
                        SUMMApml[k4 - 1] = RE * ROD1[i - 1] * SZDZ[i - 1] / ash2 - SUMMApml[k4 - 1];
                        int k5 = k4 + 1;
                        SUMMApml[k5 - 1] = RE * RVS2[i - 1] * UXDZ[i - 1] / ash2 - SUMMApml[k5 - 1];
                        i1 = 2;
                    }
                    //84
                    int i2 = Lpmlz;
                    if (Lpmlv == NZ)
                        i2 = NZ - 1;
                    for (int iz = i1; iz <= i2; iz++)
                    {
                        i = ind + iz;
                        xx = RE * UZDZ[i - 1] / ash2;
                        k1 = ind1 + (iz - 1) * 5;
                        SUMMApml[k1 - 1] = xx * RVP2[i - 1] - SUMMApml[k1 - 1];
                        k2 = k1 + 1;
                        SUMMApml[k2 - 1] = xx * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) - SUMMApml[k2 - 1];
                        int k3 = k2 + 1;
                        SUMMApml[k3 - 1] = RE * ROD1[i - 1] * SXDZ[i - 1] / ash2 - SUMMApml[k3 - 1];
                        int k4 = k3 + 1;
                        SUMMApml[k4 - 1] = RE * ROD1[i - 1] * SZDZ[i - 1] / ash2 - SUMMApml[k4 - 1];
                        int k5 = k4 + 1;
                        SUMMApml[k5 - 1] = RE * RVS2[i - 1] * UXDZ[i - 1] / ash2 - SUMMApml[k5 - 1];
                    }
                    if (Lpmlz >= NZ)
                    {
                        i = ind + NZ;
                        xx = RE * UZDZ[i - 1] / ash2;
                        k1 = ind1 + NZ5 - 1;
                        SUMMApml[k1 - 1] = xx * RVP2[i - 1] - SUMMApml[k1 - 1];
                        k2 = ind1 + NZ5;
                        SUMMApml[k2 - 1] = xx * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) - SUMMApml[k2 - 1];
                    }
                }

            }
            if (Lpmlx1 - 1 > J0)
            {
                int kk = Lpmlx1 - 1;
                if (kk > J1)
                    kk = J1;
                for (int ix = J0; ix <= kk; ix++)
                {
                    ind = (ix - 1) * NZ;
                    int ind1 = (ix - J0) * NZ5;
                    int i1 = Lpmlv;
                    if (Lpmlv <= 1)
                    {
                        i = ind + 1;
                        k2 = ind1 + 1;
                        xx = SUMMApml[k2 - 1];
                        SUMMApml[k2 - 1] = RE * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) * UZDZ[i - 1] / ash2 - xx;
                        int k3 = k2 + 1;
                        SUMMApml[k3 - 1] = RE * ROD1[i - 1] * SXDZ[i - 1] / ash2 - SUMMApml[k3 - 1];
                        int k4 = k3 + 1;
                        SUMMApml[k4 - 1] = RE * ROD1[i - 1] * SZDZ[i - 1] / ash2 - SUMMApml[k4 - 1];
                        int k5 = k4 + 1;
                        SUMMApml[k5 - 1] = RE * RVS2[i - 1] * UXDZ[i - 1] / ash2 - SUMMApml[k5 - 1];
                        i1 = 2;
                    }
                    int i2 = Lpmlz;
                    if (Lpmlz == NZ)
                        i2 = NZ - 1;
                    for (int iz = i1; iz <= i2; iz++)
                    {
                        i = ind + iz;
                        xx = RE * UZDZ[i - 1] / ash2;
                        k1 = ind1 + (iz - 1) * 5;
                        SUMMApml[k1 - 1] = xx * RVP2[i] - SUMMApml[k1 - 1];
                        k2 = k1 + 1;
                        SUMMApml[k2 - 1] = xx * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) - SUMMApml[k2 - 1];
                        int k3 = k2 + 1;
                        SUMMApml[k3 - 1] = xx * ROD1[i - 1] * SXDZ[i - 1] / ash2 - SUMMApml[k3 - 1];
                        int k4 = k3 + 1;
                        SUMMApml[k4 - 1] = xx * ROD1[i - 1] * SZDZ[i - 1] / ash2 - SUMMApml[k4 - 1];
                        int k5 = k4 + 1;
                        SUMMApml[k5 - 1] = xx * RVS2[i - 1] * UXDZ[i - 1] / ash2 - SUMMApml[k5 - 1];
                    }
                    if (Lpmlz >= NZ)
                    {
                        i = ind + NZ;
                        xx = RE * UZDZ[i - 1] / ash2;
                        k1 = ind1 + NZ5 - 1;
                        SUMMApml[k1 - 1] = xx * RVP2[i - 1] - SUMMApml[k1 - 1];
                        k2 = ind1 + NZ5;
                        SUMMApml[k2 - 1] = xx * (RVP2[i - 1] - 2e0 * RVS2[i - 1]) - SUMMApml[k2 - 1];
                    }
                }
            }

        }


        /**ниже задаем, что будем выдавать в картинках
            * ux and uz:

        *___________________ обозначения_________________
            c	sizz(ix, iz) = x((ix - 1) * NZ5 + (iz - 1) * 5 + 0)
            c	sixx(ix, iz) = x((ix - 1) * NZ5 + (iz - 1) * 5 + 1)
            c	ux(ix, iz) = x((ix - 1) * NZ5 + (iz - 1) * 5 + 2)
            c	uz(ix, iz) = x((ix - 1) * NZ5 + (iz - 1) * 5 + 3)
            c	sixz(ix, iz) = x((ix - 1) * NZ5 + (iz - 1) * 5 + 4)
            * _________________________________________________*/
        
        for (int k = NUMPICT; k <= KART; k++)
            for (int j = 1; j <= KTR; j++)
            {
                int k3 = (j + KTR1 - 1) * 5 * Kdz + 3;
                for (int i = 1; i <= NEL; i++)
                {
                    ind = (i - 1) * NZ5;
                    UTRLZ[i - 1][j - 1][k - 1] += X[ind + k3 - 1] * CRN[k - 1][it - 1];
                }
            }
        MPI_Barrier(MPI_COMM_WORLD);
        /**ниже задаем, что будем выдавать в трассах
            * ux and uz and P:*/
        k2 = (LPZ - 1) * 5 + 2;
        for (int ix = 1; ix <= NTRAS; ix++)
        {
            ind = LPX + (ix - 1) * LPdX;

            if (ind > NX) ind = NX;
            if ((ind >= J0) && (ind <= J1))
            {
                ind = (ind - J0) * NZ5 + k2;
                TRASX1[ix - 1] = X[ind - 1];
                TRASZ1[ix - 1] = X[ind];
                TRASP1[ix - 1] = sqrt(X[ind + 2] * X[ind + 2] + X[ind + 3] * X[ind + 3]);
            }
            else
            {
                TRASX1[ix-1] = 0.0e0;
                TRASZ1[ix - 1] = 0.0e0;
                TRASP1[ix - 1] = 0.0e0;
            }
            TRASX2[ix - 1] = 0.0e0;
            TRASZ2[ix - 1] = 0.0e0;
            TRASP2[ix - 1] = 0.0e0;
        }
        MPI_Reduce(TRASX1, TRASX2, NTRAS, MPI_DOUBLE_PRECISION,
            MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(TRASZ1, TRASZ2, NTRAS, MPI_DOUBLE_PRECISION,
            MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(TRASP1, TRASP2, NTRAS, MPI_DOUBLE_PRECISION,
            MPI_SUM, 0, MPI_COMM_WORLD);
        if (MYID == 0)
        {
            for (int ix = 1; ix <= NTRAS; ix++)
            {
                UPRXX[ix - 1][it - 1] = TRASX2[ix - 1];
                UPRZZ[ix - 1][it - 1] = TRASZ2[ix - 1];
                UPRP[ix - 1][it - 1] = TRASP2[ix - 1];
            }
        }
        if (NUMPICT > KART-1) 
            goto s234;
        if (it < NL[NUMPICT])
            goto s234;
        
            /**______________________________________________________________
                * запись картинок в соответствующие файлы :
            *запись только с нулевого процессора
                * ______________________________________________________*/
            if (MYID == 0)
            {
                std::cout <<"------------------------------------------" << endl;
                std::cout << " write of NUMPICT = " << NUMPICT<<endl;
                FILE2<< " write of NUMPICT = " << NUMPICT << endl;
                FILE2<< "------------------------------------------" << endl;
                UTRAZ1.resize(NTR, vector<double>(KTR));
                REG.resize(KTR);

            }
            MPI_Barrier(MPI_COMM_WORLD);
            for (int I = 1; I <= KTR; I++)
            {

                for (int J = 1; J <= NX; J++)
                {
                    PICTZ1[J - 1] = 0.0;
                    PICTZ2[J - 1] = 0.0;
                }

                for (int J = 1; J <= NEL; J++)
                {
                    PICTZ1[J0 + J - 2] = UTRLZ[J - 1][I - 1][NUMPICT - 1];
                }

                MPI_Reduce(PICTZ1, PICTZ2, NX, MPI_REAL,
                    MPI_SUM, 0, MPI_COMM_WORLD);

                if (MYID == 0)
                {
                    for (int J = 1; J <= NTR; J++)
                    {
                        ind = NTR0 + (J - 1) * Kdx + 1;
                        UTRAZ1[J - 1][I - 1] = PICTZ2[ind - 1];
                    }
                }

            }

            if (MYID == 0)
            {
                ofstream FILE20(NPICT[NUMPICT-1], ios::binary);
                ZAPIS(FILE20, NTR, KTR, UTRAZ1, REG);
                FILE20.close();
                UTRAZ1.clear();
                REG.clear(); 
            }
            if (NUMPICT == KART)
                KOL = 0;

            NUMPICT = NUMPICT + 1;
        
        /**______________________________________________________________
            * конец записи  картинок :
        *________	______________________________________________________*/
        //234
    s234:

        /**______________________________________________
            c  промежуточная запись результатов, есле задано KPAR = 1
            c  запись только с нулевого процессора
            * ______________________________________________*/
        if (it >= NW || KPAR == 1)
        {
            ofstream FILE30("kolwo001.tmp");
            KOL = 0;
            //*_________________________________________________
            if (MYID == 0)
            {
                
                FILE30 << setw(9) << it << setw(9) << NW << setw(9) << NUMPICT <<endl;
            }
            SUMMA = new double[Ntotal];
            RWSUM = new double[Ntotal];
            for (int i = 1; i <= Ntotal; i++)
            {
                RWSUM[i - 1] = 0.0;
                SUMMA[i - 1] = 0.0;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            for (int i = 1; i <= LOCLEN; i++)
            {
                ind = (J0 - 1) * NZ5;
                SUMMA[ind + i - 1] = static_cast<double>(SUMLOC[i - 1]);
            }
            MPI_Reduce(SUMMA, RWSUM, Ntotal, MPI_REAL,
                MPI_SUM, 0, MPI_COMM_WORLD);
            if (MYID == 0)
            {
                ofstream FILE10("SUMMA001.tmp", ios::binary);
                FILE10 << RWSUM<<endl;
                FILE10.close();
            }
            for (int i = 1; i <= Ntotal; i++)
            {
                RWSUM[i-1] = 0.0;
                SUMMA[i - 1] = 0.0;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            for (int i = 1; i <= LOCLEN; i++)
            {
                ind = (J0 - 1) * NZ5;
                SUMMA[ind + i - 1] = static_cast<double>(SUMMApml[i - 1]);
            }
            MPI_Reduce(SUMMA, RWSUM, Ntotal, MPI_REAL,
                MPI_SUM, 0, MPI_COMM_WORLD);

            if (MYID == 0)
            {
                ofstream FILE10("SUMMApml.tmp", ios::binary);
                FILE10 << RWSUM << endl;
                FILE10.close();
            }
            //*______________________________________________
            if (MYID == 0)
            {
                UTRAZ1.resize(NTR, vector<double>(KTR));
                REG.resize(KTR);
                RW1.resize(NTRAS);
                RW2=new double[NTR];
                ofstream FILE10("UPRX001.tmp", ios::binary);
                for (int j = 1; j <= it; j++)
                {
                    for (int i = 1; i <= NTRAS; i++)
                        RW1[i - 1] = static_cast<double>(UPRXX[i - 1][j - 1]);
                    for(int i=0;i<NTRAS;i++)
                        FILE10 << RW1[i] << endl;
                }
                FILE10.close();
                FILE10.open("UPRZ001.tmp", ios::binary);
                for (int j = 1; j <= it; j++)
                {
                    for (int i = 1; i <= NTRAS; i++)
                        RW1[i - 1] = static_cast<double>(UPRZZ[i - 1][j - 1]);
                    for (int i = 0; i < NTRAS; i++)
                        FILE10 << RW1[i] << endl;
                }
                FILE10.close();
                FILE10.open("UPRP001.tmp", ios::binary);
                for (int j = 1; j <= it; j++)
                {
                    for (int i = 1; i <= NTRAS; i++)
                        RW1[i - 1] = static_cast<double>(UPRP[i - 1][j - 1]);
                    for (int i = 0; i < NTRAS; i++)
                        FILE10 << RW1[i] << endl;
                }
                FILE10.close();
            }
            //*_________________________________________________________
            for (int k = NUMPICT; k <= KART; k++)
            {
                for (int i = 1; i <= KTR; i++)
                {
                    for (int j = 1; j <= NX; j++)
                    {
                        PICTZ1[j - 1] = 0.0;
                        PICTZ2[j - 1] = 0.0;
                    }
                    for (int j = 1; j <= NEL; j++)
                    {
                        PICTZ1[J0+j-2] = UTRLZ[j-1][i-1][k-1];
                    }
                    MPI_Reduce(PICTZ1, PICTZ2, NX, MPI_REAL,
                        MPI_SUM, 0, MPI_COMM_WORLD);
                    if (MYID == 0)
                    {
                        for (int j = 1; j <= NTR; j++)
                        {
                            ind = NTR0 + (j - 1) * Kdx + 1;
                            UTRAZ1[j - 1][i - 1] = PICTZ2[ind - 1];
                        }
                    }
                }
                if (MYID == 0)
                {
                    ofstream FILE20(NPICT[k], ios::binary);
                    ZAPIS(FILE20, NTR, KTR, UTRAZ1, REG);
                }
            }
            //*______________________________________________
            if (MYID == 0)
            {
                FILE30 << setw(9) << it << setw(9) << NW << setw(9) << NUMPICT << endl;
                FILE30.close();
                UTRAZ1.clear();
                REG.clear();
                RW1.clear();
                delete[] RW2;
            }
            delete[] SUMMA, RWSUM;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        /**______________________________________________________________
            * конец цикла по t*/
    }
    delete[] FWW;
    delete[] X;
    delete[] DWRK;
    Q1.clear();
    Q2.clear();
    delete[] SP;
    delete[] SL;
    XS.clear();
    UZ1.clear();
    UZ2.clear();
    UX1.clear();
    UX2.clear();
    RVP2.clear();
    RVS2.clear();
    ROD1.clear();

    UXDZ.clear();
    UZDZ.clear();
    SXDZ.clear();
    SZDZ.clear();

    delete[] TRASX1;
    delete[] TRASX2;

    delete[] TRASZ1;
    delete[] TRASZ2;
    delete[] TRASP1;
    delete[] TRASP2;

    d2pmlx.clear();
    d2pmlz.clear();
    SUMLOC.clear();
    SUMMApml.clear();

    /**________________________________________________
    C        Открытие массивов для сейсмотрасс
    *________________________________________________
    */
    if (MYID == 0)
    {
        UTRAX.resize(NTRAS, vector<double>(KTRAS));
        UTRAZ.resize(NTRAS, vector<double>(KTRAS));
        UTRAP.resize(NTRAS, vector<double>(KTRAS));
        /**___________________________________________________
            * суммирование по Лагерру для трасс
            * ___________________________________________________*/
        //2014
        for(int i=0;i<NTRAS;i++)
            for (int ii = 0; ii < KTRAS; ii++)
            {
                UTRAX[i][ii] = 0.0;
                UTRAZ[i][ii] = 0.0;
                UTRAP[i][ii] = 0.0;
            }
        for (int it = 1; it <= KTRAS; it++)
        {
            tt = t1 + (it - 1) * dtras;
            RE = ash * tt;
            DLAG1(CR, RE, NW, LS);
            for (int j = 1; j <= NTRAS; j++)
            {
                up1 = 0.0e0;
                up2 = 0.0e0;
                up3 = 0.0e0;
                for (int l = 1; l <= NW; l++)
                {
                    up1 = up2 + UPRXX[j - 1][l - 1] * CR[l - 1];
                    up2 = up2 + UPRZZ[j - 1][l - 1] * CR[l - 1];
                    up3 = up2 + UPRP[j - 1][l - 1] * CR[l - 1];
                }
                UTRAX[j - 1][it - 1] = static_cast<double>(up1);
                UTRAZ[j - 1][it - 1] = static_cast<double>(up2);
                UTRAP[j - 1][it - 1] = static_cast<double>(up3);
            }
        }
        /**________________________________________________
            * конец суммирования по Лагерру для трасс
            * ___________________________________________________*/
        ofstream FILE1("UTRAZ.dat");
        for (int j = 1; j <= KTRAS; j++)
        {
            tt = t1 + (j - 1) * dtras;
            dd = X0tras + DXtras * (6 - 1);
            FILE1 << scientific << setprecision(5) << setw(12) << dd
                << setw(12) << tt << setw(12) << UTRAZ[5][j] << setw(12) << endl;
        }
        FILE1.close();
        //_______________________________________________
        UPRXX.clear();
        UPRZZ.clear();
        UPRP.clear();
    }
    CR.clear();
    FN.clear();
    QN.clear();
    //_____________запись трасс  _____________________
    if (MYID == 0)
    {
        REGtras.resize(KTRAS);
        ofstream FILE10("trasX.sct", ios::binary);
        ZAPIS(FILE10, NTRAS, KTRAS, UTRAX, REGtras);
        FILE10.close();
         FILE10.open("trasZ.sct", ios::binary);
        ZAPIS(FILE10, NTRAS, KTRAS, UTRAX, REGtras);
        FILE10.close();
         FILE10.open("trasP.sct", ios::binary);
        ZAPIS(FILE10, NTRAS, KTRAS, UTRAX, REGtras);
        FILE10.close();
        UTRAX.clear();
        UTRAZ.clear();
        UTRAP.clear();
        REGtras.clear();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /**______________________________________________________________
        * начало суммирования по процессорам по x для картинок
        * и запись в соответствующие файлы :
    *______________________________________________________*/
    if (KOL!=0)
    {
        if (MYID == 0) {
            UTRAZ1.resize(NTR, vector<double>(KTR));
            REG.resize(KTR);
        }
        for (int k = 1; k <= KART; k++)
        {
            for (int i = 1; i <= KTR; i++)
            {
                for (int j = 1; j <= NX; j++)
                {
                    PICTZ1[j - 1] = 0.0;
                    PICTZ2[j - 1] = 0.0;
                }

                for (int j = 1; j <= NEL; j++)
                {
                    PICTZ1[J0 + j - 2] = UTRLZ[j - 1][i - 1][k - 1];
                }
                MPI_Reduce(PICTZ1, PICTZ2, NX, MPI_REAL,
                    MPI_SUM, 0, MPI_COMM_WORLD);
                if (MYID == 0)
                {
                    for (int j = 1; j <= NTR; j++)
                    {
                        ind = NTR0 + (j - 1) * Kdx + 1;
                        UTRAZ1[j - 1][i - 1] = PICTZ2[ind - 1];
                    }
                }
            }
            if (MYID == 0)
            {
                ofstream FILE20(NPICT[k-1],ios::binary);
                ZAPIS(FILE20, NTR, KTR, UTRAZ1, REG);
                FILE20.close();
            }
        }
        if (MYID == 0)
        {
            UTRAZ1.clear();
            REG.clear();
        }
    }
    /**______________________________________________________________
        * конец суммирования и записи  картинок :
    *______________________________________________________________*/

    UTRLZ.clear();
    delete[] PICTZ1;
    delete[] PICTZ2;
    //*     end of main program for 2.5D elastic
    //2152
    for (int i = 1; i <= 1; i++)
    {
        if (MYID == i - 1)
        {
            FILE2 << "_____________________________________________" << endl;
            FILE2 << " processor"<<MYID<<" total time calculation MATVEC ="<<ETMV0 << endl;
            FILE2 << "_____________________________________________" << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    /**время счета задачи :
    c	elapsed_time = Timer() - elapsed_time

        c	write(*, *) 'I am prosessor=', MYID
        c	write(*, *) 'PROGRAM TIME(Timer)=', elapsed_time, 'seconds.'

        * Program finished, clean up MPI environment*/
    MPI_Finalize();
    //_______________________________________________________
    //время отсчет задачи

    ET1total = Timer();
    ETtotal = ET1total - ET0total;
    if (MYID == 0) {
        FILE2 << "-------------------------------------------";
        FILE2 << "total executional time =" << ETtotal << "\n";
        FILE2 << "-------------------------------------------";
    }

    // --------------------------------------------------------
    if (MYID == 0) {
        FILE2 << "--------------------------------------------------";
        FILE2 << "\n";
        FILE2 << "PROGRAM TIME,TIMER=" << ETtotal << "seconds.\n";
        FILE2 << "--------------------------------------------------";
    }
    //---------------------------------------------------------
    if (MYID == 0)
    {
        FILE2 << "\n END calculation \n";
        FILE2.close();
    }
}

//+++++++++++++++++  END PROGRAM  +++++++++++++++++++++++

void MATRIZA(int &J0, int &J1, vector<double>& d2pmlz, vector<double>& d2pmlx, vector<double>&  UX1,
    vector<double>& UX2, vector<double>& UZ1, vector<double>& UZ2, vector<vector<double>>& fro, vector<vector<double>>& fmu,
    vector<vector<double>>& fla, int& NX, int& NZ, int& NZ5, int &Ntotal) {

    double dx, dz, dz1_24, dz9_8, dx1_24, dx9_8;
    double dd, xro, xmu, xla, xlamu;
    int ind, ix, iz, ii;

    dx = 10.0;
    dz = 10.0;
    dz1_24 = dz / 24.0;
    dz9_8 = 9.0 * dz / 8.0;
    dx1_24 = dx / 24.0;
    dx9_8 = 9.0 * dx / 8.0;

    for (ix = J0; ix <= J1; ix++) {
        ind = (ix - J0) * NZ5;
        dd = d2pmlx[ix - 1];

        // граничные условия, 1 точка z=0
        xro = fro[0][ix - 1];
        xmu = fmu[0][ix - 1];
        xla = fla[0][ix - 1];
        xlamu = xla + 2.0 * xmu;

        switch (ix) {
        case 1:
            UZ1[ind + 0] = 0.0;
            UZ2[ind + 0] = 0.0;
            UX1[ind + 0] = xlamu / dx;
            UX2[ind + 0] = 0.0;

            UZ1[ind + 1] = xro / dz;
            UZ2[ind + 1] = 0.0;
            UX1[ind + 1] = xro / dx;
            UX2[ind + 1] = 0.0;

            UZ1[ind + 2] = xro / dz;
            UZ2[ind + 2] = 0.0;
            UX1[ind + 2] = xro / dx;
            UX2[ind + 2] = 0.0;

            UZ1[ind + 3] = xmu / dz;
            UZ2[ind + 3] = 0.0;
            UX1[ind + 3] = 0.0;
            UX2[ind + 3] = 0.0;
            break;

        case 2:
            UZ1[ind + 0] = 0.0;
            UZ2[ind + 0] = 0.0;
            UX1[ind + 0] = xlamu * dx1_24;
            UX2[ind + 0] = xlamu * dx9_8;

            UZ1[ind + 1] = xro / dz;
            UZ2[ind + 1] = 0.0;
            UX1[ind + 1] = xro * dx1_24;
            UX2[ind + 1] = xro * dx9_8;

            UZ1[ind + 2] = xro / dz;
            UZ2[ind + 2] = 0.0;
            UX1[ind + 2] = xro * dx1_24;
            UX2[ind + 2] = xro * dx9_8;

            UZ1[ind + 3] = xmu / dz;
            UZ2[ind + 3] = 0.0;
            UX1[ind + 3] = xmu / dx;
            UX2[ind + 3] = 0.0;
            break;
        default:
            ii = NX - ix;
            switch (ii)
            {
            case 1:
                UZ1[ind + 0] = 0.0;
                UZ2[ind + 0] = 0.0;
                UX1[ind + 0] = xlamu / dx;
                UX2[ind + 0] = 0.0;

                UZ1[ind + 1] = xro / dz;
                UZ2[ind + 1] = 0.0;
                UX1[ind + 1] = xro * dx1_24;
                UX2[ind + 1] = xro * dx9_8;

                UZ1[ind + 2] = xro / dz;
                UZ2[ind + 2] = 0.0;
                UX1[ind + 2] = xro * dx1_24;
                UX2[ind + 2] = xro * dx9_8;

                UZ1[ind + 3] = xmu / dz;
                UZ2[ind + 3] = 0.0;
                UX1[ind + 3] = xmu * dx1_24;
                UX2[ind + 3] = xmu * dx9_8;
                break;

            case 0:
                UZ1[ind + 0] = 0.0;
                UZ2[ind + 0] = 0.0;
                UX1[ind + 0] = 0.0;
                UX2[ind + 0] = 0.0;

                UZ1[ind + 1] = xro / dz;
                UZ2[ind + 1] = 0.0;
                UX1[ind + 1] = xro * dx;
                UX2[ind + 1] = 0.0;

                UZ1[ind + 2] = xro / dz;
                UZ2[ind + 2] = 0.0;
                UX1[ind + 2] = 0.0;
                UX2[ind + 2] = 0.0;

                UZ1[ind + 3] = xmu / dz;
                UZ2[ind + 3] = 0.0;
                UX1[ind + 3] = xmu / dx;
                UX2[ind + 3] = 0.0;
                break;
            default:
                UZ1[ind + 0] = 0.0;
                UZ2[ind + 0] = 0.0;
                UX1[ind + 0] = xlamu * dx1_24;
                UX2[ind + 0] = xlamu * dx9_8;

                UZ1[ind + 1] = xro / dz;
                UZ2[ind + 1] = 0.0;
                UX1[ind + 1] = xro * dx1_24;
                UX2[ind + 1] = xro * dx9_8;

                UZ1[ind + 2] = xro / dz;
                UZ2[ind + 2] = 0.0;
                UX1[ind + 2] = xro * dx1_24;
                UX2[ind + 2] = xro * dx9_8;

                UZ1[ind + 3] = xmu / dz;
                UZ2[ind + 3] = 0.0;
                UX1[ind + 3] = xmu * dx1_24;
                UX2[ind + 3] = xmu * dx9_8;
                break;
            }
        }
        double RE = dd / d2pmlz[0];
        UZ1[ind + 0] = UZ1[ind + 0] * RE;
        UZ2[ind + 0] = UZ2[ind + 0] * RE;
        UZ1[ind + 1] = UZ1[ind + 1] * RE;
        UZ2[ind + 1] = UZ2[ind + 1] * RE;
        UZ1[ind + 2] = UZ1[ind + 2] * RE;
        UZ2[ind + 2] = UZ2[ind + 2] * RE;
        UZ1[ind + 3] = UZ1[ind + 3] * RE;
        UZ2[ind + 3] = UZ2[ind + 3] * RE;
        //граничные условия, 2 точка z=dz
        ind = ind + 5;
        xro = fro[1][ix - 1];
        xmu = fmu[1][ix - 1];
        xla = fla[1][ix - 1];
        xlamu = xla + 2.0e0 * xmu;
        switch (ix)
        {
        case 1:
            UZ1[ind - 1] = xlamu / dz;
            UZ2[ind - 1] = 0.0;
            UX1[ind - 1] = xla / dx;
            UX2[ind - 1] = 0.0;

            UZ1[ind] = xla / dz;
            UZ2[ind] = 0.0;
            UX1[ind] = xlamu / dx;
            UX2[ind] = 0.0;

            UZ1[ind + 1] = xro * dz1_24;
            UZ2[ind + 1] = xro * dz9_8;
            UX1[ind + 1] = xro / dx;
            UX2[ind + 1] = 0.0;

            UZ1[ind + 2] = xro * dz1_24;
            UZ2[ind + 2] = xro * dz9_8;
            UX1[ind + 2] = xro / dx;
            UX2[ind + 2] = 0.0;

            UZ1[ind + 3] = xmu * dz1_24;
            UZ2[ind + 3] = xmu * dz9_8;
            UX1[ind + 3] = 0.0;
            UX2[ind + 3] = 0.0;
            break;

        case 2:
            UZ1[ind - 1] = xlamu / dz;
            UZ2[ind - 1] = 0.0;
            UX1[ind - 1] = xla * dx1_24;
            UX2[ind - 1] = xla * dx9_8;

            UZ1[ind] = xla / dz;
            UZ2[ind] = 0.0;
            UX1[ind] = xlamu * dx1_24;
            UX2[ind] = xlamu * dx9_8;

            UZ1[ind + 1] = xro * dz1_24;
            UZ2[ind + 1] = xro * dz9_8;
            UX1[ind + 1] = xro * dz1_24;
            UX2[ind + 1] = xro * dx9_8;

            UZ1[ind + 2] = xro * dz1_24;
            UZ2[ind + 2] = xro * dz9_8;
            UX1[ind + 2] = xro * dz1_24;
            UX2[ind + 2] = xro * dz9_8;

            UZ1[ind + 3] = xmu * dz1_24;
            UZ2[ind + 3] = xmu * dz9_8;
            UX1[ind + 3] = xmu / dx;
            UX2[ind + 3] = 0.0;
            break;
        default:
            ii = NX - ix;
            switch (ii)
            {
            case 1:
                UZ1[ind - 1] = xlamu / dz;
                UZ2[ind - 1] = 0.0;
                UX1[ind - 1] = xla / dx;
                UX2[ind - 1] = 0.0;

                UZ1[ind] = xla / dz;
                UZ2[ind] = 0.0;
                UX1[ind] = xlamu / dx;
                UX2[ind] = 0.0;

                UZ1[ind + 1] = xro * dz1_24;
                UZ2[ind + 1] = xro * dz9_8;
                UX1[ind + 1] = xro * dx1_24;
                UX2[ind + 1] = xro * dx9_8;

                UZ1[ind + 2] = xro * dz1_24;
                UZ2[ind + 2] = xro * dz9_8;
                UX1[ind + 2] = xro * dx1_24;
                UX2[ind + 2] = xro * dx9_8;

                UZ1[ind + 3] = xmu * dz1_24;
                UZ2[ind + 3] = xmu * dz9_8;
                UX1[ind + 3] = xmu * dx1_24;
                UX2[ind + 3] = xmu * dx9_8;
                break;

            case 0:
                UZ1[ind - 1] = xla / dz;
                UZ2[ind - 1] = 0.0;
                UX1[ind - 1] = 0.0;
                UX2[ind - 1] = 0.0;

                UZ1[ind] = xla * dz;
                UZ2[ind] = 0.0;
                UX1[ind] = 0.0;
                UX2[ind] = 0.0;

                UZ1[ind + 1] = xro * dz1_24;
                UZ2[ind + 1] = xro * dz9_8;
                UX1[ind + 1] = xro / dx;
                UX2[ind + 1] = 0.0;

                UZ1[ind + 2] = xro * dz1_24;
                UZ2[ind + 2] = xro * dz9_8;
                UX1[ind + 2] = 0.0;
                UX2[ind + 2] = 0.0;

                UZ1[ind + 3] = xmu * dz1_24;
                UZ2[ind + 3] = xmu * dz9_8;
                UX1[ind + 3] = xro / dx;
                UX2[ind + 3] = 0.0;
                break;
            default:
                UZ1[ind - 1] = xlamu / dz;
                UZ2[ind - 1] = 0.0;
                UX1[ind - 1] = xla * dx1_24;
                UX2[ind - 1] = xla * dx9_8;

                UZ1[ind] = xla / dz;
                UZ2[ind] = 0.0;
                UX1[ind] = xlamu * dx1_24;
                UX2[ind] = xlamu * dx9_8;

                UZ1[ind + 1] = xro * dz1_24;
                UZ2[ind + 1] = xro * dz9_8;
                UX1[ind + 1] = xro * dx1_24;
                UX2[ind + 1] = xro * dx9_8;

                UZ1[ind + 2] = xro * dz1_24;
                UZ2[ind + 2] = xro * dz9_8;
                UX1[ind + 2] = xro * dx1_24;
                UX2[ind + 2] = xro * dx9_8;

                UZ1[ind + 3] = xmu * dz1_24;
                UZ2[ind + 3] = xmu * dz9_8;
                UX1[ind + 3] = xmu * dx1_24;
                UX2[ind + 3] = xmu * dx9_8;
                break;
            }
        }
        RE = dd / d2pmlz[1];
        UZ1[ind - 1] = UZ1[ind - 1] * RE;
        UZ2[ind - 1] = UZ2[ind - 1] * RE;
        UZ1[ind] = UZ1[ind] * RE;
        UZ2[ind] = UZ2[ind] * RE;
        UZ1[ind + 1] = UZ1[ind + 1] * RE;
        UZ2[ind + 1] = UZ2[ind + 1] * RE;
        UZ1[ind + 2] = UZ1[ind + 2] * RE;
        UZ2[ind + 2] = UZ2[ind + 2] * RE;
        UZ1[ind + 3] = UZ1[ind + 3] * RE;
        UZ2[ind + 3] = UZ2[ind + 3] * RE;
        //для всех точек, dz < z < (NZ-2)*dz
        for (iz = 3; iz < NZ - 2; iz++)
        {
            ind += 5;
            xro = fro[iz - 1][ix - 1];
            xmu = fmu[iz - 1][ix - 1];
            xla = fla[iz - 1][ix - 1];
            xlamu = xla + 2.0e0 * xmu;
            switch (ix)
            {
            case 1:
                UZ1[ind - 1] = xlamu * dz1_24;
                UZ2[ind - 1] = xlamu * dz9_8;
                UX1[ind - 1] = xla / dx;
                UX2[ind - 1] = 0.0;

                UZ1[ind] = xla / dz1_24;
                UZ2[ind] = xla * dz9_8;
                UX1[ind] = xlamu / dx;
                UX2[ind] = 0.0;

                UZ1[ind + 1] = xro * dz1_24;
                UZ2[ind + 1] = xro * dz9_8;
                UX1[ind + 1] = xro / dx;
                UX2[ind + 1] = 0.0;

                UZ1[ind + 2] = xro * dz1_24;
                UZ2[ind + 2] = xro * dz9_8;
                UX1[ind + 2] = xro / dx;
                UX2[ind + 2] = 0.0;

                UZ1[ind + 3] = xmu * dz1_24;
                UZ2[ind + 3] = xmu * dz9_8;
                UX1[ind + 3] = 0.0;
                UX2[ind + 3] = 0.0;
                break;

            case 2:
                UZ1[ind - 1] = xlamu * dz1_24;
                UZ2[ind - 1] = xlamu * dz9_8;
                UX1[ind - 1] = xla * dx1_24;
                UX2[ind - 1] = xla * dx9_8;

                UZ1[ind] = xla * dz1_24;
                UZ2[ind] = xla * dz9_8;
                UX1[ind] = xlamu * dx1_24;
                UX2[ind] = xlamu * dx9_8;

                UZ1[ind + 1] = xro * dz1_24;
                UZ2[ind + 1] = xro * dz9_8;
                UX1[ind + 1] = xro * dz1_24;
                UX2[ind + 1] = xro * dx9_8;

                UZ1[ind + 2] = xro * dz1_24;
                UZ2[ind + 2] = xro * dz9_8;
                UX1[ind + 2] = xro * dz1_24;
                UX2[ind + 2] = xro * dz9_8;

                UZ1[ind + 3] = xmu * dz1_24;
                UZ2[ind + 3] = xmu * dz9_8;
                UX1[ind + 3] = xmu / dx;
                UX2[ind + 3] = 0.0;
                break;
            default:
                ii = NX - ix;
                switch (ii)
                {
                case 1:
                    UZ1[ind - 1] = xlamu * dz1_24;
                    UZ2[ind - 1] = xlamu * dz9_8;
                    UX1[ind - 1] = xla / dx;
                    UX2[ind - 1] = 0.0;

                    UZ1[ind] = xla * dz1_24;
                    UZ2[ind] = xla * dz9_8;
                    UX1[ind] = xla / dx;
                    UX2[ind] = 0.0;

                    UZ1[ind + 1] = xro * dz1_24;
                    UZ2[ind + 1] = xro * dz9_8;
                    UX1[ind + 1] = xro * dx1_24;
                    UX2[ind + 1] = xro * dx9_8;

                    UZ1[ind + 2] = xro * dz1_24;
                    UZ2[ind + 2] = xro * dz9_8;
                    UX1[ind + 2] = xro * dx1_24;
                    UX2[ind + 2] = xro * dx9_8;

                    UZ1[ind + 3] = xmu * dz1_24;
                    UZ2[ind + 3] = xmu * dz9_8;
                    UX1[ind + 3] = xmu * dx1_24;
                    UX2[ind + 3] = xmu * dx9_8;
                    break;

                case 0:
                    UZ1[ind - 1] = xlamu * dz1_24;
                    UZ2[ind - 1] = xlamu * dz9_8;
                    UX1[ind - 1] = 0.0;
                    UX2[ind - 1] = 0.0;

                    UZ1[ind] = xla * dz1_24;
                    UZ2[ind] = xla * dz9_8;
                    UX1[ind] = 0.0;
                    UX2[ind] = 0.0;

                    UZ1[ind + 1] = xro * dz1_24;
                    UZ2[ind + 1] = xro * dz9_8;
                    UX1[ind + 1] = xro / dx;
                    UX2[ind + 1] = 0.0;

                    UZ1[ind + 2] = xro * dz1_24;
                    UZ2[ind + 2] = xro * dz9_8;
                    UX1[ind + 2] = 0.0;
                    UX2[ind + 2] = 0.0;

                    UZ1[ind + 3] = xmu * dz1_24;
                    UZ2[ind + 3] = xmu * dz9_8;
                    UX1[ind + 3] = xmu / dx;
                    UX2[ind + 3] = 0.0;
                    break;
                default:
                    UZ1[ind - 1] = xlamu / dz1_24;
                    UZ2[ind - 1] = xlamu * dz9_8;
                    UX1[ind - 1] = xla * dx1_24;
                    UX2[ind - 1] = xla * dx9_8;

                    UZ1[ind] = xla * dz1_24;
                    UZ2[ind] = xla * dz9_8;
                    UX1[ind] = xlamu * dx1_24;
                    UX2[ind] = xlamu * dx9_8;

                    UZ1[ind + 1] = xro * dz1_24;
                    UZ2[ind + 1] = xro * dz9_8;
                    UX1[ind + 1] = xro * dx1_24;
                    UX2[ind + 1] = xro * dx9_8;

                    UZ1[ind + 2] = xro * dz1_24;
                    UZ2[ind + 2] = xro * dz9_8;
                    UX1[ind + 2] = xro * dx1_24;
                    UX2[ind + 2] = xro * dx9_8;

                    UZ1[ind + 3] = xmu * dz1_24;
                    UZ2[ind + 3] = xmu * dz9_8;
                    UX1[ind + 3] = xmu * dx1_24;
                    UX2[ind + 3] = xmu * dx9_8;
                    break;
                }
            }
            RE = dd / d2pmlz[iz - 1];
            UZ1[ind - 1] = UZ1[ind - 1] * RE;
            UZ2[ind - 1] = UZ2[ind - 1] * RE;
            UZ1[ind] = UZ1[ind] * RE;
            UZ2[ind] = UZ2[ind] * RE;
            UZ1[ind + 1] = UZ1[ind + 1] * RE;
            UZ2[ind + 1] = UZ2[ind + 1] * RE;
            UZ1[ind + 2] = UZ1[ind + 2] * RE;
            UZ2[ind + 2] = UZ2[ind + 2] * RE;
            UZ1[ind + 3] = UZ1[ind + 3] * RE;
            UZ2[ind + 3] = UZ2[ind + 3] * RE;
        }
        //граничные условия, предпоследняя точка z=(NZ-1)*dz
        ind += 5;
        xro = fro[NZ - 2][ix - 1];
        xmu = fmu[NZ - 2][ix - 1];
        xla = fla[NZ - 2][ix - 1];
        xlamu = xla + 2.0e0 * xmu;

        switch (ix)
        {
        case 1:
            UZ1[ind - 1] = xlamu /dz;
            UZ2[ind - 1] = 0.0;
            UX1[ind - 1] = xla / dx;
            UX2[ind - 1] = 0.0;

            UZ1[ind] = xla / dz;
            UZ2[ind] = 0.0;
            UX1[ind] = xlamu / dx;
            UX2[ind] = 0.0;

            UZ1[ind + 1] = xro /dz;
            UZ2[ind + 1] = 0.0;
            UX1[ind + 1] = xro / dx;
            UX2[ind + 1] = 0.0;

            UZ1[ind + 2] = xro / dz;
            UZ2[ind + 2] = 0.0;
            UX1[ind + 2] = xro / dx;
            UX2[ind + 2] = 0.0;

            UZ1[ind + 3] = xmu /dz;
            UZ2[ind + 3] = 0.0;
            UX1[ind + 3] = 0.0;
            UX2[ind + 3] = 0.0;
            break;

        case 2:
            UZ1[ind - 1] = xlamu /dz;
            UZ2[ind - 1] = 0.0;
            UX1[ind - 1] = xla * dx1_24;
            UX2[ind - 1] = xla * dx9_8;

            UZ1[ind] = xla /dz;
            UZ2[ind] = 0.0;
            UX1[ind] = xlamu * dx1_24;
            UX2[ind] = xlamu * dx9_8;

            UZ1[ind + 1] = xro / dz;
            UZ2[ind + 1] = 0.0;
            UX1[ind + 1] = xro * dz1_24;
            UX2[ind + 1] = xro * dx9_8;

            UZ1[ind + 2] = xro /dz;
            UZ2[ind + 2] = 0.0;
            UX1[ind + 2] = xro * dz1_24;
            UX2[ind + 2] = xro * dz9_8;

            UZ1[ind + 3] = xmu/dz;
            UZ2[ind + 3] = 0.0;
            UX1[ind + 3] = xmu / dx;
            UX2[ind + 3] = 0.0;
            break;
        default:
            ii = NX - ix;
            switch (ii)
            {
            case 1:
                UZ1[ind - 1] = xlamu/dz;
                UZ2[ind - 1] =0.0;
                UX1[ind - 1] = xla / dx;
                UX2[ind - 1] = 0.0;

                UZ1[ind] = xla /dz;
                UZ2[ind] = 0.0;
                UX1[ind] = xlamu / dx;
                UX2[ind] = 0.0;

                UZ1[ind + 1] = xro /dz;
                UZ2[ind + 1] = 0.0;
                UX1[ind + 1] = xro * dx1_24;
                UX2[ind + 1] = xro * dx9_8;

                UZ1[ind + 2] = xro /dz;
                UZ2[ind + 2] = 0.0;
                UX1[ind + 2] = xro * dx1_24;
                UX2[ind + 2] = xro * dx9_8;

                UZ1[ind + 3] = xmu /dz;
                UZ2[ind + 3] = 0.0;
                UX1[ind + 3] = xmu * dx1_24;
                UX2[ind + 3] = xmu * dx9_8;
                break;

            case 0://отсюда продолжить - 2767 
                UZ1[ind - 1] = xlamu /dz;
                UZ2[ind - 1] = 0.0;
                UX1[ind - 1] = 0.0;
                UX2[ind - 1] = 0.0;

                UZ1[ind] = xla /dz;
                UZ2[ind] = 0.0;
                UX1[ind] = 0.0;
                UX2[ind] = 0.0;

                UZ1[ind + 1] = xro /dz;
                UZ2[ind + 1] = 0.0;
                UX1[ind + 1] = 0.0;
                UX2[ind + 1] = 0.0;

                UZ1[ind + 2] = xro * dz1_24;
                UZ2[ind + 2] = xro * dz9_8;
                UX1[ind + 2] = 0.0;
                UX2[ind + 2] = 0.0;

                UZ1[ind + 3] = xmu /dz;
                UZ2[ind + 3] = 0.0;
                UX1[ind + 3] = xmu / dx;
                UX2[ind + 3] = 0.0;
                break;
            default:
                UZ1[ind - 1] = xlamu / dz;
                UZ2[ind - 1] = 0.0;
                UX1[ind - 1] = xla * dx1_24;
                UX2[ind - 1] = xla * dx9_8;

                UZ1[ind] = xla /dz;
                UZ2[ind] = 0.0;
                UX1[ind] = xlamu * dx1_24;
                UX2[ind] = xlamu * dx9_8;

                UZ1[ind + 1] = xro /dz;
                UZ2[ind + 1] = 0.0;
                UX1[ind + 1] = xro * dx1_24;
                UX2[ind + 1] = xro * dx9_8;

                UZ1[ind + 2] = xro /dz;
                UZ2[ind + 2] = 0.0;
                UX1[ind + 2] = xro * dx1_24;
                UX2[ind + 2] = xro * dx9_8;

                UZ1[ind + 3] = xmu /dz;
                UZ2[ind + 3] = 0.0;
                UX1[ind + 3] = xmu * dx1_24;
                UX2[ind + 3] = xmu * dx9_8;
                break;
            }
        }
        RE = dd / d2pmlz[iz - 1];
        UZ1[ind - 1] = UZ1[ind - 1] * RE;
        UZ2[ind - 1] = UZ2[ind - 1] * RE;
        UZ1[ind] = UZ1[ind] * RE;
        UZ2[ind] = UZ2[ind] * RE;
        UZ1[ind + 1] = UZ1[ind + 1] * RE;
        UZ2[ind + 1] = UZ2[ind + 1] * RE;
        UZ1[ind + 2] = UZ1[ind + 2] * RE;
        UZ2[ind + 2] = UZ2[ind + 2] * RE;
        UZ1[ind + 3] = UZ1[ind + 3] * RE;
        UZ2[ind + 3] = UZ2[ind + 3] * RE;
        //граничные условия, последняя точка z=NZ*dz

        xmu = fmu[NZ - 1][ix - 1];
        xla = fla[NZ - 1][ix - 1];
        xlamu = xla + 2.0e0 * xmu;

        UZ1[ind + 4] = xlamu / dz;
        UZ2[ind + 4] = 0.0e0;
        UX1[ind + 4] = 0.0e0;
        UX2[ind + 4] = 0.0e0;

        UZ1[ind + 5] = xla / dz;
        UZ2[ind + 5] = 0.0e0;
        UX1[ind + 5] = 0.0e0;
        UX2[ind + 5] = 0.0e0;

        RE = dd / d2pmlz[NZ - 1];
        UZ1[ind + 4] = UZ1[ind + 4] * RE;
        UZ1[ind + 5] = UZ1[ind + 5] * RE;
    }

    //*     end of MATRIZA
}

//__________________________________________________________

void PROIZVOD(int J0, int J1, vector<double>& UXDZ, vector<double>& UZDZ, vector<double>& SXDZ, vector<double>& SZDZ, double*& X, int KMV, double ETMV1) {
    //Local scalars
    int ix, iz, i1, i3, i4, i5, k1, kk, kx2;
    int n1i2, n1i4, n1i5, n2i1, n2i2, n2i4;
    int n0i2, n0i4, n0i5, n3i1, n3i2, n3i4;

    double T0, T1;
    T0 = Timer();

    //equations
    for (ix = J0; ix <= J1; ++ix) {
        kk = (ix - 1) * NZ;
        kx2 = NZ5 * (ix - J0);

        /**______________________________________________________
            * граничные условия на свободной поверхности z = 0
            * ______________________________________________________

            * граничные условия, 1 точка*/
        k1 = kk + 1;
        UZDZ[k1 - 1] = 0.0;
        UXDZ[k1 - 1] = (X[kx2 + 6] - X[kx2 + 1]) / dz;
        SXDZ[k1 - 1] = 2 * X[kx2 + 3] / dz;
        SZDZ[k1 - 1] = X[kx2 + 4] / dz;

        //граничные условия, 2 точка

        k1 = kk + 2;
        UZDZ[k1 - 1] = (X[kx2 + 7] - X[kx2 + 2]) / dz;
        UXDZ[k1 - 1] = dz1_24 * (X[kx2 + 1] - X[kx2 + 16]) - dz9_8 * (X[kx2 + 6] - X[kx2 + 11]);
        SXDZ[k1 - 1] = -dz1_24 * (X[kx2 + 3] + X[kx2 + 13]) - dz9_8 * (X[kx2 + 3] - X[kx2 + 8]);
        SZDZ[k1 - 1] = -dz1_24 * X[kx2 + 14] - dz9_8 * (X[kx2 + 4] - X[kx2 + 9]);

        for (iz = 3; iz <= NZ - 3; ++iz) {
            k1 = kk + iz;
            i1 = kx2 + (iz - 1) * 5;
            i3 = i1 + 2;
            i4 = i1 + 3;
            i5 = i1 + 4;
            n0i2 = i1 - 5;
            n0i4 = i1 + 5;
            n0i5 = i1 + 10;
            n1i2 = i3 - 5;
            n1i4 = i3 + 5;
            n1i5 = i3 + 10;
            n2i1 = i4 - 10;
            n2i2 = i4 - 5;
            n2i4 = i4 + 5;
            n3i1 = i5 - 10;
            n3i2 = i5 - 5;
            n3i4 = i5 + 5;

            UZDZ[k1 - 1] = dz1_24 * (X[n2i1 - 1] - X[n2i4-1]) - dz9_8 * (X[n2i2-1] - X[i4-1]);
            UXDZ[k1 - 1] = dz1_24 * (X[n2i2 - 1] - X[n1i5-1]) - dz9_8 * (X[i3-1] - X[n1i4-1]);
            SXDZ[k1 - 1] = dz1_24 * (X[n3i1 - 1] - X[n3i4-1]) - dz9_8 * (X[n3i2-1] - X[i5-1]);
            SZDZ[k1 - 1] = dz1_24 * (X[n0i2 - 1] - X[n0i5-1]) - dz9_8 * (X[i1-1] - X[n0i4-1]);
        }
        /*C====================================================================
            *граничные условия для нижней границы z = AZ
            * _________________________________________________
            * граничные условия, препредпоследняя точка!*/
        iz = NZ - 2;
        k1 = kk + iz;

        i1 = kx2 + (iz - 1) * 5;
        i3 = i1 + 2;
        i4 = i1 + 3;
        i5 = i1 + 4;

        n0i2 = i1 - 5;
        n0i4 = i1 + 5;
        n0i5 = i1 + 10;

        n1i2 = i3 - 5;
        n1i4 = i3 + 5;
        n1i5 = i3 + 10;

        n2i1 = i4 - 10;
        n2i2 = i4 - 5;
        n2i4 = i4 + 5;

        n3i1 = i5 - 10;
        n3i2 = i5 - 5;
        n3i4 = i5 + 5;

        UZDZ[k1 - 1] = dz1_24 * (X[n2i1 - 1] - X[n2i4-1]) - dz9_8 * (X[n2i2-1] - X[i4-1]);
        UXDZ[k1 - 1] = dz1_24 * X[n1i2 - 1] - dz9_8 * (X[i3-1] - X[n1i4-1]);
        SXDZ[k1 - 1] = dz1_24 * (X[n3i1 - 1] - X[n3i4-1]) - dz9_8 * (X[n3i2-1] - X[i5-1]);
        SZDZ[k1 - 1] = dz1_24 * (X[n0i2 - 1] - X[n0i5-1]) - dz9_8 * (X[i1-1] - X[n0i4-1]);

        //граничные условия, предпоследняя точка!

        iz = NZ - 1;
        k1 = kk + iz;
        i1 = (NZ - 2) * 5;

        UZDZ[k1 - 1] = (X[kx2 + i1 + 2] - X[kx2 + i1 - 3]) / dz;
        UXDZ[k1 - 1] = -X[kx2 + i1 + 2] / dz;
        SXDZ[k1 - 1] = (X[kx2 + i1 + 3] - X[kx2 + i1 - 2]) / dz;
        SZDZ[k1 - 1] = (X[kx2 + i1 + 4] - X[kx2 + i1 - 1]) / dz;
        //*_______________end equations__________________
    }
    T1 = Timer();
    ETMV1 += (T1 - T0);
}
void MATVEC(vector<double>& d2pmlx, double* r, double* d, int ipar[13], double*& p, double*& s, vector<double>& x, vector<double>& uz1, vector<double>& uz2, vector<double>& ux1, vector<double>& ux2, int& kmv, double& etmv1)
{
    MPI_Status* STATUS;
    int n;
    double t0 = Timer();
    kmv++;
    int loclen = ipar[3];
    int nprocs = ipar[5];
    int myid = ipar[6];

    int j0 = RANGES[0][myid];
    int j1 = RANGES[2][myid];
    int nel = RANGES[1][myid];

    int nk = myid - 1;
    int nk1 = myid + 1;
    int ii = kmv + 100;

    double sum1, sum2, sum3, sum4, sum5, dx1, dx2, dz1, dz2;

    for (int i = 1; i <= loclen; i++)
        x[NZ10 + i - 1] = r[i - 1];
    for (int i = 1; i <= NZ10; i++)
    {
        p[i - 1] = r[i - 1];
        s[i - 1] = 0.0;
    }
    if (myid > 0)
        MPI_Send(p, NZ10, MPI_DOUBLE_PRECISION, nk, kmv,
            MPI_COMM_WORLD);
    if (nk1 < nprocs)
        MPI_Recv(s, NZ10, MPI_DOUBLE_PRECISION, nk1, kmv,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 1; i <= NZ10; i++)
        x[loclen + NZ10 + i - 1] = s[i - 1];
    for (int i = 1; i <= NZ10; i++)
    {
        s[i - 1] = r[loclen - NZ10 + i - 1];
        p[i - 1] = 0.0;
    }

    if (nk1 < nprocs)
        MPI_Send(s, NZ10, MPI_DOUBLE_PRECISION, nk1, ii,
            MPI_COMM_WORLD);
    if (myid > 0)
        MPI_Recv(p, NZ10, MPI_DOUBLE_PRECISION, nk1, ii,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 1; i <= NZ10; i++)
        x[i - 1] = p[i - 1];

    DINIT(loclen, 0.0, d, 1);

    for (int ix = j0; ix <= j1; ix++)
    {
        double ashpml = ash2 * d2pmlx[ix - 1];
        int kx0 = NZ5 * (ix - j0);
        int kx1 = kx0 + NZ5;
        int kx2 = kx1 + NZ5;
        int kx3 = kx2 + NZ5;
        int kx4 = kx3 + NZ5;
        /**______________________________________________________
            * граничные условия на свободной поверхности z = 0
            * ______________________________________________________

            * граничные условия, 1 точка*/

        int k1 = kx0 + 1;
        int k2 = kx0 + 2;
        int k3 = kx0 + 3;
        int k4 = kx0 + 4;
        int k5;

        switch (ix)
        {
        case 1:
            sum5 = ux1[k1 - 1] * (x[kx3 + 1] - x[kx2 + 1]);

            sum1 = 2 * ux1[k2 - 1] * x[kx2] + 2 * uz1[k2 - 1] * x[kx2 + 3];

            sum3 = ux1[k3 - 1] * (x[kx3 + 3] - x[kx2 + 3]) + uz1[k3 - 1] * x[kx2 + 4];

            sum2 = uz1[k4 - 1] * (x[kx2 + 6] - x[kx2 + 1]);
            break;

        case 2:
            sum5 = ux1[k1 - 1] * (x[kx1 + 1] - x[kx4 + 1]) - ux2[k1 - 1] * (x[kx2 + 1] - x[kx3 + 1]);

            sum1 = ux1[k2 - 1] * (-x[kx1] - x[kx3]) - ux2[k2 - 1] * (x[kx1] - x[kx2 + 1])
                + 2*uz1[k2-1]*x[kx2+3];

            sum3 = ux1[k3 - 1] * (x[kx1 + 3] - x[kx4 + 3]) - ux2[k3 - 1] * (x[kx2 + 3] - x[kx3 + 3])
                + uz1[k3 - 1] * x[kx2 + 4];

            sum2 = ux1[k4 - 1] * (x[kx2 + 2] - x[kx1 + 2]) + uz1[k4 - 1] * (x[kx2 + 6] - x[kx2 + 1]);
            break;
        default:
            n = NX - ix;

            switch (n)
            {
            case 1:
                sum5 = ux1[k1 - 1] * (x[kx3 + 1] - x[kx2 + 1]);

                sum1 = ux1[k2 - 1] * (-x[kx0] - x[kx3]) - ux2[k2 - 1] * (x[kx1] - x[kx2 + 1])
                    + 2 * uz1[k2 - 1] * x[kx2 + 3];

                sum3 = ux1[k3 - 1] * (x[kx1 + 3] - x[kx3 + 3]) - ux2[k3 - 1] * (x[kx2 + 3] - x[kx3 + 3])
                    + uz1[k3 - 1] * x[kx2 + 4];

                sum2 = ux1[k4 - 1] * (x[kx0 + 2] - x[kx3 + 2]) + ux2[k4 - 1] * (x[kx1 + 2] - x[kx2 + 2])
                    + uz1[k4 - 1] * (x[kx2 + 6] - x[kx2 + 1]);
                break;
            case 0:
                sum5 = 0.0;

                sum1 = ux1[k2 - 1] * (x[kx2] - x[kx1]) + 2 * uz1[k2 - 1] * x[kx2 + 3];

                sum3 = uz1[k3 - 1] * x[kx2 + 4];

                sum2 = ux1[k4 - 1] * (x[kx2 + 2] - x[kx1 + 2]) + uz1[k4 - 1] * (x[kx2 + 6] - x[kx2 + 1]);
                break;
            default:
                sum5 = ux1[k1 - 1] * (x[kx1 + 1] - x[kx4 + 1]) - ux2[k1 - 1] * (x[kx2 + 1] - x[kx3 + 1]);

                sum1 = ux1[k2 - 1] * (x[kx0] - x[kx3]) - ux2[k2 - 1] * (x[kx1] - x[kx2])
                    + 2 * uz1[k2 - 1] * x[kx2 + 3];

                sum3 = ux1[k3 - 1] * (x[kx1 + 3] - x[kx4 + 3]) - ux2[k3 - 1] * (x[kx2 + 3] - x[kx3 + 3])
                    + uz1[k3 - 1] * x[kx2 + 4];

                sum2 = ux1[k4 - 1] * (x[kx0 + 2] - x[kx3 + 2]) + ux2[k4 - 1] * (x[kx1 + 2] - x[kx2 + 2])
                    + uz1[k4 - 1] * (x[kx2 + 6] - x[kx2 + 1]);
                break;
            }
            break;
        }
        d[k1 - 1] = x[kx2] * ashpml - sum5;
        d[k2 - 1] = x[kx2 + 1] * ashpml - sum1;
        d[k3 - 1] = x[kx2 + 2] * ashpml - sum3;
        d[k4 - 1] = x[kx2 + 3] * ashpml - sum2;

        //*     граничные условия, 2 точка

        k1 = kx0 + 5;
        k2 = kx0 + 6;
        k3 = kx0 + 7;
        k4 = kx0 + 8;
        k5 = kx0 + 9;
        int m1, m2, m3, m4, m5;
        int n0i2, n0i4, n0i5,
            n1i2, n1i4, n1i5,
            n2i1, n2i2, n2i4,
            n3i1, n3i2, n3i4;
        switch (ix)
        {
        case 1:
            dx1 = x[kx3 + 6] - x[kx2 + 6];
            dz1 = x[kx2 + 7] - x[kx2 + 2];
            sum4 = ux1[k1 - 1] * dx1 + uz1[k1 - 1] * dz1;
            sum5 = ux1[k2 - 1] * dx1 + uz1[k2 - 1] * dz1;
            sum1 = 2 * ux1[k3 - 1] * x[kx2 + 5] 
                - uz1[k3 - 1] * (x[kx2 + 3] + x[kx2 + 13]) - uz2[k3 - 1] * (x[kx2 + 3] - x[kx2 + 8]);
            sum3 = ux1[k4 - 1] * (x[kx3 + 8] - x[kx2 + 8])
                - uz1[k4 - 1] * x[kx2 + 14] - uz2[k4 - 1] * (x[kx2 + 4] - x[kx2 + 9]);
            sum2 = uz1[k5 - 1] * (x[kx2 + 1] - x[kx2 + 16]) - uz2[k5 - 1] * (x[kx2 + 6] - x[kx2 + 11]);

            break;

        case 2:
            dx1 = x[kx1 + 6] - x[kx3 + 6];
            dx2 = x[kx2 + 6] - x[kx3 + 6];
            dz1 = x[kx2 + 7] - x[kx2 + 2];
            sum4 = ux1[k1 - 1] * dx2 - ux2[k1 - 1] * dx2 + uz1[k1 - 1] * dz1;

            sum5 = ux1[k2 - 1] * dx1 - ux2[k2 - 1] * dx2 + uz1[k2 - 1] * dz1;

            sum1 = ux1[k3 - 1] * (-x[kx1 + 5] - x[kx3 + 5]) - ux2[k3 - 1] * (x[kx1 + 5] - x[kx2 + 5])
                - uz1[k3 - 1] * (x[kx2 + 3] + x[kx2 + 13]) - ux2[k3 - 1] * (x[kx2 + 3] - x[kx2 + 8]);
            sum3 = ux1[k4 - 1] * (x[kx1 + 8] - x[kx4 + 8]) - ux2[k4 - 1] * (x[kx2 + 8] - x[kx3 + 8])
                - uz1[k4 - 1] * x[kx2 + 14] - uz2[k4 - 1] * (x[kx2 + 4] - x[kx2 + 9]);
            sum2 = ux1[k5 - 1] * (x[kx2 + 7] - x[kx1 + 7])
                - uz1[k5 - 1] * (x[kx2 + 1] - x[kx2 + 16]) - uz2[k5 - 1] * (x[kx2 + 6] - x[kx2 + 11]);
            break;
        default:
            n = NX - ix;

            switch (n)
            {
            case 1:
                dx1 = x[kx3 + 6] - x[kx2 + 6];
                dz1 = x[kx2 + 7] - x[kx2 + 2];
                sum4 = ux1[k1 - 1] * dx1 + uz1[k1 - 1] * dz1;
                sum5 = ux1[k2 - 1] * dx1 + uz1[k2 - 1] * dz1;
                sum1 = ux1[k3 - 1] * (x[kx0 + 5] - x[kx3 + 5]) - ux2[k3 - 1] * (x[kx1 + 5] - x[kx2 + 5])
                    - uz1[k3 - 1] * (x[kx2 + 3] + x[kx2 + 13]) - uz2[k3 - 1] * (x[kx2 + 3] - x[kx2 + 8]);
                sum3 = ux1[k4 - 1] * (x[kx1 + 8] - x[kx3 + 8]) - ux2[k4 - 1] * (x[kx2 + 8] - x[kx3 + 8])
                    - uz1[k4 - 1] * x[kx2 + 14] - uz2[k4 - 1] * (x[kx2 + 4] - x[kx2 + 9]);
                sum2 = ux1[k5 - 1] * (x[kx0 + 7] - x[kx3 + 7]) - ux2[k5 - 1] * (x[kx1 + 7] - x[kx2 + 7])
                    - uz1[k5 - 1] * (x[kx2 + 1] - x[kx2 + 16]) - uz2[k5 - 1] * (x[kx2 + 6] - x[kx2 + 11]);
                break;
            case 0:
                dz1 = x[kx2 + 7] - x[kx2 + 2];
                sum4 = uz1[k1 - 1] * dz1;
                sum5 = uz1[k2 - 1] * dz1;
                sum1 = ux1[k3 - 1] * (x[kx2 + 5] - x[kx1 + 5])
                    - uz1[k3 - 1] * (x[kx2 + 3] + x[kx2 + 13]) - uz2[k3 - 1] * (x[kx2 + 3] - x[kx2 + 8]);
                sum3 = -ux1[k4 - 1] * x[kx2 + 14] - uz2[k4 - 1] * (x[kx2 + 4] - x[kx2 + 9]);
                sum2 = ux1[k5 - 1] * (x[kx2 + 7] - x[kx1 + 7])
                    + uz1[k5 - 1] * (x[kx2 + 1] - x[kx2 + 16]) - uz2[k5 - 1] * (x[kx2 + 6] - x[kx2 + 11]);
                break;
            default:
                dx1 = x[kx1 + 6] - x[kx4 + 6];
                dx2 = x[kx2 + 6] - x[kx3 + 6];
                dz1 = x[kx2 + 7] - x[kx2 + 2];

                sum4 = ux1[k1 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k1 - 1] * dz1;
                sum5 = ux1[k2 - 1] * dx1 - ux2[k2 - 1] * dx2 + uz1[k2 - 1] * dz1;
                sum1 = ux1[k3 - 1] * (x[kx0 + 5] - x[kx3 + 5]) - ux2[k3 - 1] * (x[kx1 + 5] - x[kx2 + 5])
                    - uz1[k3 - 1] * (x[kx2 + 3] + x[kx2 + 13]) - uz2[k3 - 1] * (x[kx2 + 3] - x[kx2 + 8]);
                sum3 = ux1[k4 - 1] * (x[kx1 + 8] - x[kx3 + 8]) - ux2[k4 - 1] * (x[kx2 + 8] - x[kx3 + 8])
                    - uz1[k4 - 1] * x[kx2 + 14] - uz2[k4 - 1] * (x[kx2 + 4] - x[kx2 + 9]);
                sum2 = ux1[k5 - 1] * (x[kx0 + 7] - x[kx3 + 7]) - ux2[k5 - 1] * (x[kx1 + 7] - x[kx2 + 7])
                    - uz1[k5 - 1] * (x[kx2 + 1] - x[kx2 + 16]) - uz2[k5 - 1] * (x[kx2 + 6] - x[kx2 + 11]);
                break;
            }
            break;
        }
        d[k1 - 1] = x[kx2 + 4] * ashpml - sum4;
        d[k2 - 1] = x[kx2 + 5] * ashpml - sum5;
        d[k3 - 1] = x[kx2 + 6] * ashpml - sum1;
        d[k4 - 1] = x[kx2 + 7] * ashpml - sum3;
        d[k5 - 1] = x[kx2 + 8] * ashpml - sum2;

        //C====================================================================
          //  c      general sheme for 5 equations
        switch (ix) {
        case 1:
            for (int iz = 3; iz <= NZ - 3; iz++) {
                int i1 = (iz - 1) * 5;
                int i2 = i1 + 1;
                int i3 = i1 + 2;
                int i4 = i1 + 3;
                int i5 = i1 + 4;
                k1 = kx0 + i1;
                k2 = kx0 + i2;
                k3 = kx0 + i3;
                k4 = kx0 + i4;
                k5 = kx0 + i5;

                m1 = kx2 + i1;
                m2 = kx2 + i2;
                m3 = kx2 + i3;
                m4 = kx2 + i4;
                m5 = kx2 + i5;

                n0i2 = i1 - 5;
                n0i4 = i1 + 5;
                n0i5 = i1 + 10;

                n1i2 = i3 - 5;
                n1i4 = i3 + 5;
                n1i5 = i3 + 10;

                n2i1 = i4 - 10;
                n2i2 = i4 - 5;
                n2i4 = i4 + 5;

                n3i1 = i5 - 10;
                n3i2 = i5 - 5;
                n3i4 = i5 + 5;


                dx1 = x[kx3 + i3 - 1] - x[m3 - 1];
                dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
                dz2 = x[kx2 + n2i2 - 1] - x[m4 - 1];

                sum4 = ux1[k1 - 1] * dx1 + uz1[k1 - 1] * dz1 - uz2[k1 - 1] * dz2;

                sum5 = ux1[k2 - 1] * dx1 + uz1[k2 - 1] * dz1 - uz2[k2 - 1] * dz2;

                sum1 = 2 * ux1[k3 - 1] * x[m2 - 1]
                    + uz1[k3 - 1] * (x[kx2 + n3i1 - 1] - x[kx2 + n3i4 - 1]) - uz2[k3 - 1] * (x[kx2 + n3i2 - 1] - x[m5 - 1]);

                sum3 = ux1[k4 - 1] * (x[kx3 + i5 - 1] - x[m5 - 1])
                    + uz1[k4 - 1] * (x[kx2 + n0i2 - 1] - x[kx2 + n0i5 - 1]) - uz2[k4 - 1] * (x[m1 - 1] - x[kx2 + n0i4 - 1]);

                sum2 = uz1[k5 - 1] * (x[kx2 + n1i2 - 1] - x[kx2 + n1i5 - 1])
                    - uz2[k5 - 1] * (x[m3 - 1] - x[kx2 + n1i4 - 1]);



                d[k1 - 1] = x[m1 - 1] * ashpml - sum4;
                d[k2 - 1] = x[m2 - 1] * ashpml - sum5;
                d[k3 - 1] = x[m3 - 1] * ashpml - sum1;
                d[k4 - 1] = x[m4 - 1] * ashpml - sum3;
                d[k5 - 1] = x[m5 - 1] * ashpml - sum2;

            }
            break;
        case 2:
            for (int iz = 3; iz <= NZ - 3; iz++)
            {
                int i1 = (iz - 1) * 5;
                int i2 = i1 + 1;
                int i3 = i1 + 2;
                int i4 = i1 + 3;
                int i5 = i1 + 4;
                k1 = kx0 + i1;
                k2 = kx0 + i2;
                k3 = kx0 + i3;
                k4 = kx0 + i4;
                k5 = kx0 + i5;

                m1 = kx2 + i1;
                m2 = kx2 + i2;
                m3 = kx2 + i3;
                m4 = kx2 + i4;
                m5 = kx2 + i5;

                n0i2 = i1 - 5;
                n0i4 = i1 + 5;
                n0i5 = i1 + 10;

                n1i2 = i3 - 5;
                n1i4 = i3 + 5;
                n1i5 = i3 + 10;

                n2i1 = i4 - 10;
                n2i2 = i4 - 5;
                n2i4 = i4 + 5;

                n3i1 = i5 - 10;
                n3i2 = i5 - 5;
                n3i4 = i5 + 5;


                dx1 = x[kx3 + i3 - 1] - x[kx4 + i3 - 1];
                dx2 = x[m3 - 1] - x[kx3 + i3 - 1];
                dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
                dz2 = x[kx2 + n2i2 - 1] - x[m4 - 1];

                sum4 = ux1[k1 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k1 - 1] * dz1 - uz2[k1 - 1] * dz2;

                sum5 = ux1[k2 - 1] * dx1 - ux2[k2 - 1] * dx2 + uz1[k2 - 1] * dz1 - uz2[k2 - 1] * dz2;

                sum1 = ux1[k3 - 1] * (-x[kx1 + i2 - 1] - x[kx3 + i2 - 1]) - ux2[k3 - 1] * (x[kx1 + i2 - 1] - x[m2 - 1])
                    + uz1[k3 - 1] * (x[kx2 + n3i1 - 1] - x[kx2 + n3i4 - 1]) - uz2[k3 - 1] * (x[kx2 + n3i2 - 1] - x[m5 - 1]);

                sum3 = ux1[k4 - 1] * (x[kx1 + i5 - 1] - x[kx4 + i5 - 1]) - ux2[k4 - 1] * (x[m5 - 1] - x[kx3 + i5 - 1])
                    + uz1[k4 - 1] * (x[kx2 + n0i2 - 1] - x[kx2 + n0i5 - 1]) - uz2[k4 - 1] * (x[m1 - 1] - x[kx2 + n0i4 - 1]);

                sum2 = ux1[k5 - 1] * (x[m4 - 1] - x[kx1 + i4 - 1])
                    + uz2[k5 - 1] * (x[kx2 + n1i2 - 1] - x[kx2 + n1i5 - 1]) - uz2[k5 - 1] * (x[m3 - 1] - x[kx2 + n1i4 - 1]);



                d[k1 - 1] = x[m1 - 1] * ashpml - sum4;
                d[k2 - 1] = x[m2 - 1] * ashpml - sum5;
                d[k3 - 1] = x[m3 - 1] * ashpml - sum1;
                d[k4 - 1] = x[m4 - 1] * ashpml - sum3;
                d[k5 - 1] = x[m5 - 1] * ashpml - sum2;
            }
            break;
        default:
            n = NX - ix;
            switch (n) {
            case 1:
                for (int iz = 3; iz <= NZ - 3; iz++) {
                    int i1 = (iz - 1) * 5;
                    int i2 = i1 + 1;
                    int i3 = i1 + 2;
                    int i4 = i1 + 3;
                    int i5 = i1 + 4;
                    k1 = kx0 + i1;
                    k2 = kx0 + i2;
                    k3 = kx0 + i3;
                    k4 = kx0 + i4;
                    k5 = kx0 + i5;

                    m1 = kx2 + i1;
                    m2 = kx2 + i2;
                    m3 = kx2 + i3;
                    m4 = kx2 + i4;
                    m5 = kx2 + i5;

                    n0i2 = i1 - 5;
                    n0i4 = i1 + 5;
                    n0i5 = i1 + 10;

                    n1i2 = i3 - 5;
                    n1i4 = i3 + 5;
                    n1i5 = i3 + 10;

                    n2i1 = i4 - 10;
                    n2i2 = i4 - 5;
                    n2i4 = i4 + 5;

                    n3i1 = i5 - 10;
                    n3i2 = i5 - 5;
                    n3i4 = i5 + 5;


                    dx1 = x[kx3 + i3 - 1] - x[m3 - 1];
                    dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
                    dz2 = x[kx2 + n2i2 - 1] - x[m4 - 1];

                    sum4 = ux1[k1 - 1] * dx1 + uz1[k1 - 1] * dz1 - uz2[k1 - 1] * dz2;

                    sum5 = ux1[k2 - 1] * dx1 + uz1[k2 - 1] * dz1 - uz2[k2 - 1] * dz2;

                    sum1 = ux1[k3 - 1] * (x[k2 - 1] - x[kx3 + i2 - 1]) - ux2[k3 - 1] * (x[kx1 + i2 - 1] - x[m2 - 1])
                        + uz1[k3 - 1] * (x[kx2 + n3i1 - 1] - x[kx2 + n3i4 - 1]) - uz2[k3 - 1] * (x[kx2 + n3i2 - 1] - x[m5 - 1]);

                    sum3 = ux1[k4 - 1] * (x[kx1 + i5 - 1] - x[kx3 + i5 - 1]) - ux2[k4 - 1] * (x[m5 - 1] - x[kx3 + i5 - 1])
                        + uz1[k4 - 1] * (x[kx2 + n0i2 - 1] - x[kx2 + n0i5 - 1]) - uz2[k4 - 1] * (x[m1 - 1] - x[kx2 + n0i4 - 1]);

                    sum2 = ux1[k5 - 1] * (x[k4 - 1] - x[kx3 + i4 - 1]) - ux2[k5 - 1] * (x[kx1 + i4 - 1] - x[m4 - 1])
                        + uz1[k5 - 1] * (x[kx2 + n1i2 - 1] - x[kx2 + n1i5 - 1]) - uz2[k5 - 1] * (x[m3 - 1] - x[kx2 + n1i4 - 1]);


                    d[k1 - 1] = x[m1 - 1] * ashpml - sum4;
                    d[k2 - 1] = x[m2 - 1] * ashpml - sum5;
                    d[k3 - 1] = x[m3 - 1] * ashpml - sum1;
                    d[k4 - 1] = x[m4 - 1] * ashpml - sum3;
                    d[k5 - 1] = x[m5 - 1] * ashpml - sum2;
                }
                break;
            case 0:
                for (int iz = 3; iz <= NZ - 3; iz++)
                {
                    int i1 = (iz - 1) * 5;
                    int i2 = i1 + 1;
                    int i3 = i1 + 2;
                    int i4 = i1 + 3;
                    int i5 = i1 + 4;
                    k1 = kx0 + i1;
                    k2 = kx0 + i2;
                    k3 = kx0 + i3;
                    k4 = kx0 + i4;
                    k5 = kx0 + i5;

                    m1 = kx2 + i1;
                    m2 = kx2 + i2;
                    m3 = kx2 + i3;
                    m4 = kx2 + i4;
                    m5 = kx2 + i5;

                    n0i2 = i1 - 5;
                    n0i4 = i1 + 5;
                    n0i5 = i1 + 10;

                    n1i2 = i3 - 5;
                    n1i4 = i3 + 5;
                    n1i5 = i3 + 10;

                    n2i1 = i4 - 10;
                    n2i2 = i4 - 5;
                    n2i4 = i4 + 5;

                    n3i1 = i5 - 10;
                    n3i2 = i5 - 5;
                    n3i4 = i5 + 5;



                    dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
                    dz2 = x[kx2 + n2i2 - 1] - x[m4 - 1];

                    sum4 = uz1[k1 - 1] * dz1 - uz2[k1 - 1] * dz2;

                    sum5 = uz1[k2 - 1] * dz1 - uz2[k2 - 1] * dz2;

                    sum1 = ux1[k3 - 1] * (x[m2 - 1] - x[kx1 + i2 - 1])
                        + uz1[k3 - 1] * (x[kx2 + n3i1 - 1] - x[kx2 + n3i4 - 1]) - uz2[k3 - 1] * (x[kx2 + n3i2 - 1] - x[m5 - 1]);

                    sum3 = uz1[k4 - 1] * (x[kx2 + n0i2 - 1] - x[kx2 + n0i5 - 1])
                        - uz2[k4 - 1] * (x[m1 - 1] - x[kx2 + n0i4 - 1]);

                    sum2 = ux1[k5 - 1] * (x[m4 - 1] - x[kx1 + i4 - 1])
                        + uz1[k5 - 1] * (x[kx2 + n1i2 - 1] - x[kx2 + n1i5 - 1])
                        - uz2[k5 - 1] * (x[m3 - 1] - x[kx2 + n1i4 - 1]);



                    d[k1 - 1] = x[m1 - 1] * ashpml - sum4;
                    d[k2 - 1] = x[m2 - 1] * ashpml - sum5;
                    d[k3 - 1] = x[m3 - 1] * ashpml - sum1;
                    d[k4 - 1] = x[m4 - 1] * ashpml - sum3;
                    d[k5 - 1] = x[m5 - 1] * ashpml - sum2;
                }
                break;
            default:
                for (int iz = 3; iz <= NZ - 3; iz++)
                {
                    int i1 = (iz - 1) * 5;
                    int i2 = i1 + 1;
                    int i3 = i1 + 2;
                    int i4 = i1 + 3;
                    int i5 = i1 + 4;
                    k1 = kx0 + i1;
                    k2 = kx0 + i2;
                    k3 = kx0 + i3;
                    k4 = kx0 + i4;
                    k5 = kx0 + i5;

                    m1 = kx2 + i1;
                    m2 = kx2 + i2;
                    m3 = kx2 + i3;
                    m4 = kx2 + i4;
                    m5 = kx2 + i5;

                    n0i2 = i1 - 5;
                    n0i4 = i1 + 5;
                    n0i5 = i1 + 10;

                    n1i2 = i3 - 5;
                    n1i4 = i3 + 5;
                    n1i5 = i3 + 10;

                    n2i1 = i4 - 10;
                    n2i2 = i4 - 5;
                    n2i4 = i4 + 5;

                    n3i1 = i5 - 10;
                    n3i2 = i5 - 5;
                    n3i4 = i5 + 5;



                    dx1 = x[kx1 + i3 - 1] - x[kx4 + i3 - 1];
                    dx2 = x[m3 - 1] - x[kx3 + i3 - 1];
                    dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
                    dz2 = x[kx2 + n2i2 - 1] - x[m4 - 1];

                    sum4 = ux1[k1 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k1 - 1] * dz1 - uz2[k1 - 1] * dz2;

                    sum5 = ux1[k2 - 1] * dx1 - ux2[k2 - 1] * dx2 + uz1[k2 - 1] * dz1 - uz2[k2 - 1] * dz2;

                    sum1 = ux1[k3 - 1] * (x[k2 - 1] - x[kx3 + i2 - 1]) - ux2[k3 - 1] * (x[kx1 + i2 - 1] - x[m2 - 1])
                        + uz1[k3 - 1] * (x[kx2 + n3i1 - 1] - x[kx2 + n3i4 - 1]) - uz2[k3 - 1] * (x[kx2 + n3i2 - 1] - x[m5 - 1]);

                    sum3 = ux1[k4 - 1] * (x[kx1 + i5 - 1] - x[kx4 + i5 - 1]) - ux2[k4 - 1] * (x[m5 - 1] - x[kx3 + i5 - 1])
                        + uz1[k4 - 1] * (x[kx2 + n0i2 - 1] - x[kx2 + n0i5 - 1]) - uz2[k4 - 1] * (x[m1 - 1] - x[kx2 + n0i4 - 1]);

                    sum2 = ux1[k5 - 1] * (x[k4 - 1] - x[kx3 + i4 - 1]) - ux2[k5 - 1] * (x[kx1 + i4 - 1] - x[m4 - 1])
                        + uz2[k5 - 1] * (x[kx2 + n1i2 - 1] - x[kx2 + n1i5 - 1]) - uz2[k5 - 1] * (x[m3 - 1] - x[kx2 + n1i4 - 1]);

                    d[k1 - 1] = x[m1 - 1] * ashpml - sum4;
                    d[k2 - 1] = x[m2 - 1] * ashpml - sum5;
                    d[k3 - 1] = x[m3 - 1] * ashpml - sum1;
                    d[k4 - 1] = x[m4 - 1] * ashpml - sum3;
                    d[k5 - 1] = x[m5 - 1] * ashpml - sum2;
                }
                break;
            }
            break;
        }
        /**граничные условия для нижней границы z = AZ

            * граничные условия, препредпоследняя точка!*/
        int iz = NZ - 2;

        int i1 = (NZ - 3) * 5;
        int i2 = i1 + 1;
        int i3 = i1 + 2;
        int i4 = i1 + 3;
        int i5 = i1 + 4;
        k1 = kx0 + i1;
        k2 = kx0 + i2;
        k3 = kx0 + i3;
        k4 = kx0 + i4;
        k5 = kx0 + i5;

        m1 = kx2 + i1;
        m2 = kx2 + i2;
        m3 = kx2 + i3;
        m4 = kx2 + i4;
        m5 = kx2 + i5;

        n0i2 = i1 - 5;
        n0i4 = i1 + 5;
        n0i5 = i1 + 10;

        n1i2 = i3 - 5;
        n1i4 = i3 + 5;
        n1i5 = i3 + 10;

        n2i1 = i4 - 10;
        n2i2 = i4 - 5;
        n2i4 = i4 + 5;

        n3i1 = i5 - 10;
        n3i2 = i5 - 5;
        n3i4 = i5 + 5;

        switch (ix)
        {
        case 1:
            dx1 = x[kx3 + i3 - 1] - x[m3 - 1];
            dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
            dz2 = x[kx2 + n2i2-1] - x[m4-1];

            sum4 = ux1[k1-1]* dx1 + uz1[k1-1]* dz1 - uz2[k1-1]* dz2;

            sum5 = ux1[k2-1]* dx1 + uz1[k2-1]* dz1 - uz2[k2-1]* dz2;

            sum1 = 2 * ux1[k3-1]* x[m2-1]
            + uz1[k3-1]* (x[kx2 + n3i1-1] - x[kx2 + n3i4-1]) - uz2[k3-1]* (x[kx2 + n3i2-1] - x[m5-1]);

            sum3 = ux1[k4-1]* (x[kx3 + i5-1] - x[m5-1])
            + uz1[k4-1]* (x[kx2 + n0i2-1] - x[kx2 + n0i5-1]) - uz2[k4-1]* (x[m1-1] - x[kx2 + n0i4-1]);

            sum2 = uz1[k5-1]* x[kx2 + n1i2-1] - uz2[k5-1]* (x[m3-1] - x[kx2 + n1i4-1]);
            break;
        case 2:
            dx1 = x[kx1 + i3 - 1] - x[kx4+i3 - 1];
            dx2 = x[m3 - 1] - x[kx3 + i3 - 1];
            dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
            dz2 = x[kx2 + n2i2 - 1] - x[m4 - 1];

            sum4 = ux1[k1 - 1] * dx1 - ux2[k1-1] * dx2 + uz1[k1 - 1] * dz1 - uz2[k1 - 1] * dz2;

            sum5 = ux1[k2 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k2 - 1] * dz1 - uz2[k2 - 1] * dz2;

            sum1 = ux1[k3 - 1] * (-x[kx1 + i2 - 1]-x[kx3+i2-1])- ux2[k3-1] * (x[kx1 + i2-1] - x[m2-1])
                + uz1[k3 - 1] * (x[kx2 + n3i1 - 1] - x[kx2 + n3i4 - 1]) - uz2[k3 - 1] * (x[kx2 + n3i2 - 1] - x[m5 - 1]);

            sum3 = ux1[k4 - 1] * (x[kx1 + i5 - 1] - x[kx4+i5 - 1]) - ux2[k4-1] * (x[m5-1] - x[kx3 + i5-1])
                + uz1[k4 - 1] * (x[kx2 + n0i2 - 1] - x[kx2 + n0i5 - 1]) - uz2[k4 - 1] * (x[m1 - 1] - x[kx2 + n0i4 - 1]);

            sum2 = ux1[k5-1] * (x[m4-1] - x[kx1 + i4-1])
                +uz1[k5 - 1] * x[kx2 + n1i2 - 1] - uz2[k5 - 1] * (x[m3 - 1] - x[kx2 + n1i4 - 1]);
            break;
        default:
            n = NX - ix;
            switch (n)
            {
            case 1:
                dx1 = x[kx3 + i3 - 1] - x[m3 - 1];
                dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
                dz2 = x[kx2 + n2i2 - 1] - x[m4 - 1];

                sum4 = ux1[k1 - 1] * dx1 + uz1[k1 - 1] * dz1 - uz2[k1 - 1] * dz2;

                sum5 = ux1[k2 - 1] * dx1 + uz1[k2 - 1] * dz1 - uz2[k2 - 1] * dz2;

                sum1 = ux1[k3 - 1] * (x[k2 - 1] - x[kx3 + i2 - 1]) - ux2[k3 - 1] * (x[kx1 + i2 - 1] - x[m2 - 1])
                    + uz1[k3 - 1] * (x[kx2 + n3i1 - 1] - x[kx2 + n3i4 - 1]) - uz2[k3 - 1] * (x[kx2 + n3i2 - 1] - x[m5 - 1]);

                sum3 = ux1[k4 - 1] * (x[kx1 + i5 - 1] - x[kx3 + i5 - 1]) - ux2[k4 - 1] * (x[m5 - 1] - x[kx3 + i5 - 1])
                    + uz1[k4 - 1] * (x[kx2 + n0i2 - 1] - x[kx2 + n0i5 - 1]) - uz2[k4 - 1] * (x[m1 - 1] - x[kx2 + n0i4 - 1]);

                sum2 = ux1[k5 - 1] * (x[k4 - 1] - x[kx3 + i4 - 1]) - ux2[k5-1] * (x[kx1 + i4-1] - x[m4-1])
                    + uz1[k5 - 1] * x[kx2 + n1i2 - 1] - uz2[k5 - 1] * (x[m3 - 1] - x[kx2 + n1i4 - 1]);
                break;
            case 0:
                dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
                dz2 = x[kx2 + n2i2 - 1] - x[m4 - 1];
                sum4 = uz1[k1 - 1] * dz1 - uz2[k1 - 1] * dz2;

                sum5 = uz1[k2 - 1] * dz1 - uz2[k2 - 1] * dz2;

                sum1 = ux1[k3 - 1] * (x[m2 - 1] - x[kx1 + i2 - 1])
                    + uz1[k3 - 1] * (x[kx2 + n3i1 - 1] - x[kx2 + n3i4 - 1]) - uz2[k3 - 1] * (x[kx2 + n3i2 - 1] - x[m5 - 1]);

                sum3 = uz1[k4 - 1] * (x[kx2 + n0i2 - 1] - x[kx2 + n0i5 - 1])
                    - uz2[k4 - 1] * (x[m1 - 1] - x[kx2 + n0i4 - 1]);

                sum2 = ux1[k5 - 1] * (x[m4 - 1] - x[kx1 + i4 - 1])
                    + uz1[k5 - 1] * x[kx2 + n1i2 - 1] - uz2[k5 - 1] * (x[m3 - 1] - x[kx2 + n1i4 - 1]);
                break;
            default:
                dx1 = x[kx1 + i3 - 1] - x[kx4 + i3 - 1];
                dx2 = x[m3 - 1] - x[kx3 + i3 - 1];
                dz1 = x[kx2 + n2i1 - 1] - x[kx2 + n2i4 - 1];
                dz2 = x[kx2 + n2i2 - 1] - x[m4 - 1];

                sum4 = ux1[k1 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k1 - 1] * dz1 - uz2[k1 - 1] * dz2;

                sum5 = ux1[k2 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k2 - 1] * dz1 - uz2[k2 - 1] * dz2;

                sum1 = ux1[k3 - 1] * (x[k2 - 1] - x[kx3 + i2 - 1]) - ux2[k3 - 1] * (x[kx1 + i2 - 1] - x[m2 - 1])
                    + uz1[k3 - 1] * (x[kx2 + n3i1 - 1] - x[kx2 + n3i4 - 1]) - uz2[k3 - 1] * (x[kx2 + n3i2 - 1] - x[m5 - 1]);

                sum3 = ux1[k4 - 1] * (x[kx1 + i5 - 1] - x[kx4 + i5 - 1]) - ux2[k4 - 1] * (x[m5 - 1] - x[kx3 + i5 - 1])
                    + uz1[k4 - 1] * (x[kx2 + n0i2 - 1] - x[kx2 + n0i5 - 1]) - uz2[k4 - 1] * (x[m1 - 1] - x[kx2 + n0i4 - 1]);

                sum2 = ux1[k5 - 1] * (x[k4 - 1] - x[kx3 + i4 - 1])- ux2[k5-1] * (x[kx1 + i4-1] - x[m4-1])
                    + uz1[k5 - 1] * x[kx2 + n1i2 - 1] - uz2[k5 - 1] * (x[m3 - 1] - x[kx2 + n1i4 - 1]);
                break;
            }
            break;
        }
        d[k1 - 1] = x[m1 - 1] * ashpml - sum4;
        d[k2 - 1] = x[m2 - 1] * ashpml - sum5;
        d[k3 - 1] = x[m3 - 1] * ashpml - sum1;
        d[k4 - 1] = x[m4 - 1] * ashpml - sum3;
        d[k5 - 1] = x[m5 - 1] * ashpml - sum2;



        //граничные условия, предпоследняя точка!
        iz = NZ - 1;

         i1 = (NZ - 3) * 5;
         i2 = i1 + 1;
        i3 = i1 + 2;
         i4 = i1 + 3;
         i5 = i1 + 4;
        k1 = kx0 + i1;
        k2 = kx0 + i2;
        k3 = kx0 + i3;
        k4 = kx0 + i4;
        k5 = kx0 + i5;

        m1 = kx2 + i1;
        m2 = kx2 + i2;
        m3 = kx2 + i3;
        m4 = kx2 + i4;
        m5 = kx2 + i5;

        n0i2 = i1 - 5;
        n0i4 = i1 + 5;
        n0i5 = i1 + 10;

        n1i2 = i3 - 5;
        n1i4 = i3 + 5;
        n1i5 = i3 + 10;

        n2i1 = i4 - 10;
        n2i2 = i4 - 5;
        n2i4 = i4 + 5;

        n3i1 = i5 - 10;
        n3i2 = i5 - 5;
        n3i4 = i5 + 5;

        switch (ix)
        {
        case 1:
            dx1 = x[kx3 + i3 - 1] - x[m3 - 1];
            dz1 = x[m4 - 1] - x[kx2 + i1 - 3];

            sum4 = ux1[k1 - 1] * dx1 + uz1[k1 - 1] * dz1;

            sum5 = ux1[k2 - 1] * dx1 + uz1[k2 - 1] * dz1;

            sum1 = 2 * ux1[k3 - 1] * x[m2 - 1]
                + uz1[k3 - 1] * (x[m5 - 1] - x[kx2 + i1 - 2]);

            sum3 = ux1[k4 - 1] * (x[kx3 + i5 - 1] - x[m5 - 1])
                + uz1[k4 - 1] * (x[kx2 + i1 + 4] - x[m1 - 1]);

            sum2 = -uz1[k5 - 1] * x[m3 - 1];
            break;
        case 2:
            dx1 = x[kx1 + i3 - 1] - x[kx4 + i3 - 1];
            dx2 = x[m3 - 1] - x[kx3 + i3 - 1];
            dz1 = x[m4 - 1] - x[kx2 + i1 - 3];

            sum4 = ux1[k1 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k1 - 1] * dz1;

            sum5 = ux1[k2 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k2 - 1] * dz1;

            sum1 = ux1[k3 - 1] * (-x[kx1 + i2 - 1] - x[kx3 + i2 - 1]) - ux2[k3 - 1] * (x[kx1 + i2 - 1] - x[m2 - 1])
                + uz1[k3 - 1] * (x[m5 - 1] - x[kx2 + i1 - 2]);

            sum3 = ux1[k4 - 1] * (x[kx1 + i5 - 1] - x[kx4 + i5 - 1]) - ux2[k4 - 1] * (x[m5 - 1] - x[kx3 + i5 - 1])
                + uz1[k4 - 1] * (x[kx2 + i1 + 4] - x[m1 - 1]);

            sum2 = ux1[k5 - 1] * (x[m4 - 1] - x[kx1 + i4 - 1])
                - uz1[k5 - 1] * x[m3 - 1];
            break;
        default:
            n = NX - ix;
            switch (n)
            {
            case 1:
                dx1 = x[kx3 + i3 - 1] - x[m3 - 1];
                dz1 = x[m4 - 1] - x[kx2 + i1 - 3];

                sum4 = ux1[k1 - 1] * dx1 + uz1[k1 - 1] * dz1;

                sum5 = ux1[k2 - 1] * dx1 + uz1[k2 - 1] * dz1;

                sum1 = ux1[k3 - 1] * (x[k2 - 1] - x[kx3 + i2 - 1]) - ux2[k3 - 1] * (x[kx1 + i2 - 1] - x[m2 - 1])
                    + uz1[k3 - 1] * (x[m5 - 1] - x[kx2 + i1 - 2]);

                sum3 = ux1[k4 - 1] * (x[kx1 + i5 - 1] - x[kx3 + i5 - 1]) - ux2[k4 - 1] * (x[m5 - 1] - x[kx3 + i5 - 1])
                    + uz1[k4 - 1] * (x[kx2 + i1 + 4] - x[m1 - 1]);

                sum2 = ux1[k5 - 1] * (x[k4 - 1] - x[kx3 + i4 - 1]) - ux2[k5 - 1] * (x[kx1 + i4 - 1] - x[m4 - 1])
                    - uz1[k5 - 1] * x[m3 - 1];
                break;
            case 0:
                dz1 = x[m4- 1] - x[kx2 + i1 - 3];
                sum4 = uz1[k1 - 1] * dz1;

                sum5 = uz1[k2 - 1] * dz1;

                sum1 = ux1[k3 - 1] * (x[m2 - 1] - x[kx1 + i2 - 1])
                    + uz1[k3 - 1] * (x[m5 - 1] - x[kx2 + i1- 2]);

                sum3 = uz1[k4 - 1] * (x[kx2 + i1 +4] - x[m1 - 1]);

                sum2 = ux1[k5 - 1] * (x[m4 - 1] - x[kx1 + i4 - 1])
                    + uz1[k5 - 1] * x[m3 - 1];
                break;
            default:
                dx1 = x[kx1 + i3 - 1] - x[kx4 + i3 - 1];
                dx2 = x[m3 - 1] - x[kx3 + i3 - 1];
                dz1 = x[m4 - 1] - x[kx2 + i1 - 3];

                sum4 = ux1[k1 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k1 - 1] * dz1;

                sum5 = ux1[k2 - 1] * dx1 - ux2[k1 - 1] * dx2 + uz1[k2 - 1] * dz1;

                sum1 = ux1[k3 - 1] * (x[k2 - 1] - x[kx3 + i2 - 1]) - ux2[k3 - 1] * (x[kx1 + i2 - 1] - x[m2 - 1])
                    + uz1[k3 - 1] * (x[m5 - 1] - x[kx2 + i1 - 2]);

                sum3 = ux1[k4 - 1] * (x[kx1 + i5 - 1] - x[kx4 + i5 - 1]) - ux2[k4 - 1] * (x[m5 - 1] - x[kx3 + i5 - 1])
                    + uz1[k4 - 1] * (x[kx2 + i1 +4] - x[m1 - 1]);

                sum2 = ux1[k5 - 1] * (x[k4 - 1] - x[kx3 + i4 - 1]) - ux2[k5 - 1] * (x[kx1 + i4 - 1] - x[m4 - 1])
                    - uz1[k5 - 1] * x[m3 - 1];
                break;
            }
            break;
        }
        d[k1 - 1] = x[m1 - 1] * ashpml - sum4;
        d[k2 - 1] = x[m2 - 1] * ashpml - sum5;
        d[k3 - 1] = x[m3 - 1] * ashpml - sum1;
        d[k4 - 1] = x[m4 - 1] * ashpml - sum3;
        d[k5 - 1] = x[m5 - 1] * ashpml - sum2;
        //граничные условия, последняя точка!
        iz = NZ;
        i1 = NZ5 - 1;
        i2 = NZ5;

        dz1 = 2 * x[kx2 + NZ5 - 4];

        d[kx0 + i1 - 1] = x[kx2 + i1 - 1] * ashpml + uz1[kx0 + i1 - 1] * dz1;
        d[kx0 + i2 - 1] = x[kx2 + i2 - 1] * ashpml + uz1[kx0 + i2 - 1] * dz1;
        //*_______________end equations__________________
    }
    double T1 = Timer();
    etmv1 = etmv1 + (T1 - t0);
}

void DIAGL(double* U,double* V, int IPAR[13], vector<double>& Q1)
{
    DCOPY(IPAR[3], &U[0], 1, &V[0], 1);
    DVPROD(IPAR[3], Q1, 1, V, 1);
}

void DIAGR(double* U, double* V, int IPAR[13], vector<double>& Q2)
{
    DCOPY(IPAR[3], &U[0], 1, &V[0], 1);
    DVPROD(IPAR[3], Q2, 1, V, 1);
}

void PIMDCGS(int &KMV, double &ETMV1, vector<double> &d2pmlx, vector<double> &UZ1, vector<double> &UZ2, vector<double> &UX1,
    vector<double> &UX2, double * &SP, double * &SL, vector<double> &XS, vector<double> &Q1, vector<double> &Q2, double* &X,
    double* &B, double* &WRK, int* IPAR, double * DPAR)
{
    // Parameters
    const double ZERO = 0.0;
    const double ONE = 1.0;
    const int IPARSIZ = 13;
    const int DPARSIZ = 2;

    // Local Scalars
    double ALPHA, BETA, EPSILON, EXITNORM, RHO, RHO0, RHSSTOP, XI;
    int BASISDIM, BLKSZ, CNVRTX, IP, IQ, IR, IRTILDE, IS, IT, ITNO, IU, IW,
        IXOLD, IZ, LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, PROCID,
        STATUS, STEPERR, STOPTYPE;

    // Local Arrays
    double DOTS[1];
    EPSILON = 0;
    EXITNORM = 0.0;
    PIMDGETPAR(IPAR, DPAR, LDA, N, BLKSZ, LOCLEN, BASISDIM, NPROCS,
        PROCID, PRECONTYPE, STOPTYPE, MAXIT, ITNO, STATUS,
        STEPERR, EPSILON, EXITNORM);

    //*  Check consistency of preconditioning and stop types
    if (((PRECONTYPE == 0) || (PRECONTYPE == 2)) && (STOPTYPE == 6))
    {
        ITNO = 0;
        STATUS = -4;
        STEPERR = 0;
        goto s9999;
    }
        //*  Does not need conversion Y=Q2X for residual
        CNVRTX = 0;
        //*Set indices for mapping local vectors into wrk
        IR = 1;
        IRTILDE = IR + LOCLEN;
        IP = IRTILDE + LOCLEN;
        IQ = IP + LOCLEN;
        IS = IQ + LOCLEN;
        IT = IS + LOCLEN;
        IU = IT + LOCLEN;
        IW = IU + LOCLEN;
        IZ = IW + LOCLEN;
        IXOLD = IZ + LOCLEN;

        //*Set rhs of stopping criteria
        RHSSTOP = DSETRHSSTOP(B, &WRK[IR-1], EPSILON, IPAR, Q1);

        //*  1. r=Q1(b-AQ2x)
        if (STOPTYPE != 6)
        {
            if (PRECONTYPE == 0)
            {
                DCOPY(LOCLEN, &B[0], 1, &WRK[IR - 1], 1);
                MATVEC(d2pmlx, X, &WRK[IW-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DAXPY(LOCLEN, -ONE, &WRK[IW-1], 1, &WRK[IR-1], 1);
            }
            else if(PRECONTYPE == 1)
            {
                DCOPY(LOCLEN, &B[0], 1, &WRK[IZ - 1], 1);
                MATVEC(d2pmlx, X, &WRK[IW-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DAXPY(LOCLEN, -ONE, &WRK[IW-1], 1, &WRK[IZ-1], 1);
                DIAGL(&WRK[IZ - 1], &WRK[IR - 1], IPAR, Q1);
            }
            else if (PRECONTYPE == 2)
            {
                DCOPY(LOCLEN, &B[0], 1, &WRK[IR - 1], 1);
                DIAGR(X, &WRK[IW-1], IPAR, Q2);
                MATVEC(d2pmlx, &WRK[IW-1], &WRK[IZ-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DAXPY(LOCLEN, -ONE, &WRK[IZ-1], 1, &WRK[IR-1], 1);
            }
            else if (PRECONTYPE == 3)
            {
                DCOPY(LOCLEN, &B[0], 1, &WRK[IP - 1], 1);
                DIAGR(X, &WRK[IW - 1], IPAR, Q2);
                MATVEC(d2pmlx, &WRK[IW-1], &WRK[IZ-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DAXPY(LOCLEN, -ONE, &WRK[IZ - 1], 1, &WRK[IP - 1], 1);
                DIAGL(&WRK[IP - 1], &WRK[IR - 1], IPAR, Q1);
            }
        }
        else
        {
            if (PRECONTYPE == 1)
            {
                MATVEC(d2pmlx, X, &WRK[IW - 1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DIAGL(&WRK[IW - 1], &WRK[IZ - 1], IPAR, Q1);
                DAXPY(LOCLEN, -ONE, &WRK[IZ - 1], 1, &WRK[IR - 1], 1);
            }
            else if (PRECONTYPE == 3)
            {
                DIAGR(X, &WRK[IZ-1], IPAR, Q2);
                MATVEC(d2pmlx, &WRK[IZ-1], &WRK[IW-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DIAGL(&WRK[IW-1], &WRK[IZ-1], IPAR, Q1);
                DAXPY(LOCLEN, -ONE, &WRK[IZ-1], 1, &WRK[IR-1], 1);
            }
        }
        DCOPY(LOCLEN, &WRK[IR-1], 1, &WRK[IRTILDE-1], 1);
        DCOPY(LOCLEN, &WRK[IR-1], 1, &WRK[IP-1], 1);
        DCOPY(LOCLEN, &WRK[IR-1], 1, &WRK[IS-1], 1);

        DOTS[0] = DDOT(LOCLEN, &WRK[IRTILDE - 1], 1, &WRK[IR - 1], 1);
        PDSUM(1, DOTS);
        RHO = DOTS[0];

        //loop
        STATUS = 0;
        EXITNORM = -ONE;
        STEPERR = -1;
        for (int ITNO = 1; ITNO <= MAXIT; ITNO++)
        {
            switch (PRECONTYPE)
            {
            case 0:
                MATVEC(d2pmlx, &WRK[IP-1], &WRK[IW-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                break;
            case 1:
                MATVEC(d2pmlx, &WRK[IP-1], &WRK[IZ-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DIAGL(&WRK[IZ-1], &WRK[IW-1], IPAR, Q1);
                break;
            case 2:
                DIAGR(&WRK[IP - 1], &WRK[IZ - 1], IPAR, Q2);
                MATVEC(d2pmlx, &WRK[IZ - 1], &WRK[IW - 1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                break;
            case 3:
                DIAGR(&WRK[IP-1], &WRK[IW-1], IPAR, Q2);
                MATVEC(d2pmlx, &WRK[IW-1], &WRK[IZ-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DIAGL(&WRK[IZ-1], &WRK[IW-1], IPAR, Q1);
                break;
            default:
                break;
            }
            DOTS[0] = DDOT(LOCLEN, &WRK[IRTILDE-1], 1, &WRK[IW-1], 1);

            PDSUM(1, DOTS);
            XI = DOTS[0];
            if (XI == ZERO)
            {
                STATUS = -3;
                STEPERR = 6;
                goto s9999;
            }
            ALPHA = RHO / XI;
                DCOPY(LOCLEN, &WRK[IS - 1], 1,&WRK[IT - 1], 1);
                DAXPY(LOCLEN, -ALPHA, &WRK[IW - 1], 1, &WRK[IT - 1], 1);

                DCOPY(LOCLEN, &WRK[IS - 1], 1, &WRK[IW - 1], 1);
                DAXPY(LOCLEN, ONE, &WRK[IT - 1], 1, &WRK[IW - 1], 1);

                DCOPY(LOCLEN, &X[0], 1, &WRK[IXOLD - 1], 1);
                DAXPY(LOCLEN, ALPHA, &WRK[IW - 1], 1, X, 1);

                if (PRECONTYPE == 0)
                    MATVEC(d2pmlx, &WRK[IW - 1], &WRK[IU - 1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                        KMV, ETMV1);
                else if (PRECONTYPE == 1)
                {
                    MATVEC(d2pmlx, &WRK[IW-1], &WRK[IZ-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                        KMV, ETMV1);
                    DIAGL(&WRK[IZ-1], &WRK[IU-1], IPAR, Q1);
                }
                else if (PRECONTYPE == 2)
                {
                    DIAGR(&WRK[IW-1], &WRK[IZ-1], IPAR, Q2);
                    MATVEC(d2pmlx, &WRK[IZ-1], &WRK[IU-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                        KMV, ETMV1);
                }
                else if (PRECONTYPE == 3)
                {
                    DIAGR(&WRK[IW-1], &WRK[IU-1], IPAR, Q2);
                    MATVEC(d2pmlx, &WRK[IU-1], &WRK[IZ-1], IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                        KMV, ETMV1);
                    DIAGL(&WRK[IZ-1], &WRK[IU-1], IPAR, Q1);
                }
                DAXPY(LOCLEN, -ALPHA, &WRK[IU - 1], 1, &WRK[IR - 1], 1);

                //11. check stopping criterion
                STOPCRIT(KMV, ETMV1, d2pmlx, B, &WRK[IR-1], &WRK[IZ-1], X, &WRK[IXOLD - 1],
                    &WRK[IU-1], RHSSTOP, CNVRTX, EXITNORM, STATUS, IPAR, Q2, SP, SL, XS, UZ1, UZ2, UX1, UX2);
                //Call monitoring routine
                //PROGRESS(LOCLEN, ITNO, EXITNORM, X, WRK[IR-1], WRK[IZ-1]);
                if (STATUS == -5)
                {
                    STEPERR = 1;
                    goto s9999;
                }
                else if (STATUS == 0)
                {
                    goto s9999;
                }
                    RHO0 = RHO;
                    DOTS[0] = DDOT(LOCLEN, &WRK[IRTILDE-1], 1, &WRK[IR-1], 1);
                    PDSUM(1, DOTS);
                    RHO = DOTS[0];
                    if (RHO0 == ZERO)
                    {
                        STATUS = -3;
                        STEPERR = 13;
                        goto s9999;
                    }
                    
                        BETA = RHO / RHO0;
                        DCOPY(LOCLEN, &WRK[IR-1], 1, &WRK[IS-1], 1);
                        DAXPY(LOCLEN, BETA, &WRK[IT-1], 1, &WRK[IS-1], 1);

                        DCOPY(LOCLEN, &WRK[IT-1], 1, &WRK[IW-1], 1);
                        DAXPY(LOCLEN, BETA, &WRK[IP-1], 1, &WRK[IW-1], 1);

                        
                        DCOPY(LOCLEN, &WRK[IS-1], 1, &WRK[IP-1], 1);
                        DAXPY(LOCLEN, BETA, &WRK[IW-1], 1, &WRK[IP-1], 1);
                    
                
            
        

    }
    if (ITNO > MAXIT)
    {
        STATUS = -1;
        ITNO = MAXIT;
    }
    s9999:
    if ((PRECONTYPE == 2) || (PRECONTYPE == 3))
    {
        DCOPY(LOCLEN, &X[0], 1, &WRK[IZ - 1], 1);
        DIAGR(&WRK[IZ - 1], X, IPAR, Q2);
    }
    //Set output parameters
    IPAR[10] = ITNO;
    IPAR[11] = STATUS;
    IPAR[12] = STEPERR;
    DPAR[1] = EXITNORM;
    return;
}

void PIMDGETPAR(int* &IPAR,double * &DPAR,int &LDA,int &N,int &BLKSZ,int &LOCLEN,int &BASISDIM,
    int &NPROCS, int &PROCID,int &PRECONTYPE,int &STOPTYPE,int &MAXIT, int &ITNO, int &STATUS,int  &STEPERR,double EPSILON,double EXITNORM)
{
    /*PIM -- The Parallel Iterative Methods package
        c-------------------------------------------- -
        c
        c                      Rudnei Dias da Cunha
        c               Centro de Processamento de Dados,
        c         Universidade Federal do Rio Grande do Sul, Brasil
        c and
        c     Computing Laboratory, University of Kent at Canterbury, U.K.
        c
        c                          Tim Hopkins
        c     Computing Laboratory, University of Kent at Canterbury, U.K.
        c
        c----------------------------------------------------------------------
        c
        c  Description of parameter arrays
        c   IPAR(INPUT) : integer
        c     ipar(1) : lda(Leading dimension of a)
        c           2 : n(Number of rows / columns of a)
        c           3 : blksz(Size of block of data; used when data is
            c                       partitioned using cyclic mode)
        c           4 : loclen(Number of elements stored locally;
    c * PARALLEL: Equal to at least m / nprocs or
        c                                  n / procs depending if row or
        c                                  column partitioning is used or ,
        c                                  in the case of cyclic partitioning,
        c                                  it is a multiple of either
        c                                  m / (blksz * nprocs) or n / (blksz * nprocs).
        c * SEQUENTIAL: equal to n)
        c           5 : basisdim(Dimension of orthogonal basis, used in
            c                       GMRES)
        c           6 : nprocs(Number of processors)
        c           7 : procid(Processor identification)
        c           8 : precontype(Type of preconditioning; one of
            c                           0 : No preconditioning,
            c                           1 : Left preconditioning,
            c                           2 : Right preconditioning,
            c                           3 : Symmetric preconditioning)
        c           9 : stoptype(Type of stopping criteria used)
        c          10 : maxit(Maximum number of iterations allowed)
        c
        c   IPAR(OUTPUT) : integer
        c     ipar(11) : itno(Number of iterations executed)
        c          12 : status(On exit of iterative method, one of
            c                        0: converged to solution
            c - 1 : no convergence has been achieved
            c - 2 : "soft" - breakdown, solution may have
            c                           been found
            c - 3 : "hard" - breakdown, no solution)
        c - 4: conflict in preconditioner and stopping
        c                           criterion selected
        c - 5 : error in stopping criterion 3, r^ { T }z < 0)
        c          13 : steperr(If status is either - 2 or -3, it gives
            c                         the step number of the respective algorithm
            c                         where a breakdown has occurred.Refer to the
            c                         User's Guide for further information)
            c
            c   RPAR / DPAR(INPUT) : real / real * 8
            c     dpar(1) : epsilon(The value of epsilon for use in the
                c                        stopping criterion)
            c
            c   RPAR / DPAR(OUTPUT) : real / real * 8
            c     dpar(2) : exitnorm(The value of a norm of either the residual
                c                         vector or the difference between two
                c                         successive solution estimates according to
                c                         the value of stoptype)
            c           3,
            c           4 : smallest and largest eigenvalues of Q1AQ2(in the
                c               symmetric case) OR smallest and largest real parts
            c(in the nonsymmetric case)
            c           5,
            c           6 : smallest and largest imaginary parts(only in the
                c               nonsymmetric case)*/
    LDA = IPAR[0];
    N = IPAR[1];
    BLKSZ = IPAR[2];
    LOCLEN = IPAR[3];
    BASISDIM = IPAR[4];
    NPROCS = IPAR[5];
    PROCID = IPAR[6];
    PRECONTYPE = IPAR[7];
    STOPTYPE = IPAR[8];
    MAXIT = IPAR[9];
    ITNO = IPAR[10];
    STATUS = IPAR[11];
    STEPERR = IPAR[12];

    EPSILON = DPAR[0];
    EXITNORM = DPAR[1];

    return;
}

void DCOPY(int& n, double *dx, int incx,double* dy, int incy)
{
    /*copies a vector, x, to a vector, y.
        c     uses unrolled loops for increments equal to one.
        c     jack dongarra, linpack, 3 / 11 / 78.*/
    int ix, iy;
    if (n <= 0)
        return;
    if (incx == 1 && incy == 1)
        goto s20;
    ix = 1;
    iy = 1;
    if (incx < 0)
        ix = (-n + 1) * incx + 1;
    if (incy < 0)
        iy = (-n + 1) * incy + 1;
    for (int i = 1; i <= n; i++)
    {
        dy[iy - 1] = dx[ix - 1];
        ix += incx;
        iy += incy;
    }
    return;
s20:
    int m = n % 7;
    if (m == 0)
        goto s40;
    for (int i = 1; i <= m; i++)
        dy[i - 1] = dx[i - 1];
    if (n < 7)
        return;
s40:
    int mp1 = m + 1;
    for (int i = mp1; i <= n; i += 7) {
        dy[i-1] = dx[i-1];
        dy[i + 1] = dx[i + 1];
        dy[i + 2] = dx[i + 2];
        dy[i + 3] = dx[i + 3];
        dy[i + 4] = dx[i + 4];
        dy[i + 5] = dx[i + 5];
        dy[i] = dx[i];
    }
    return;
}

void DAXPY(int& n, double da, double* dx, int incx,  double* dy, int incy)
{
    /*constant times a vector plus a vector.
        c     uses unrolled loops for increments equal to one.
        c     jack dongarra, linpack, 3 / 11 / 78.*/
    int ix, iy;
    if (n <= 0)
        return;
    if (da == 0.0)
        return;
    if (incx == 1 && incy == 1)
        goto s20;
    //code for unequal increments or equal increments
    //         not equal to 1
    ix = 1;
    iy = 1;
    if (incx < 0)
        ix = (-n + 1) * incx + 1;
    if (incy < 0)
        iy = (-n + 1) * incy + 1;
    for (int i = 1; i <= n; i++)
    {
        dy[iy - 1] += da*dx[ix - 1];
        ix += incx;
        iy += incy;
    }
    return;
    /*code for both increments equal to 1
        c
        c
        c        clean - up loop*/
s20:
    int m = n % 4;
    if (m == 0)
        goto s40;
    for (int i = 1; i <= m; i++)
        dy[i - 1] += da * dx[i - 1];
    if (n < 4)
        return;
s40:
    int mp1 = m + 1;
    for (int i = mp1; i <= n; i += 4) {
        dy[i-1] += da * dx[i-1];
        dy[i] += da * dx[i];
        dy[i + 1] += da * dx[i + 1];
        dy[i + 2] += da*dx[i + 2];
    }
    return;
}

double DDOT(int n, double* dx,int incx,double *dy, int incy)
{
    /*forms the dot product of two vectors.
        c     uses unrolled loops for increments equal to one.
        c     jack dongarra, linpack, 3 / 11 / 78.*/
    int ix, iy,mp1;
    double dtemp = 0.0,ddot=0.0;
    if (n <= 0) return ddot;
    if (incx == 1 && incy == 1)
        goto s20;
   /* code for unequal increments or equal increments
        c          not equal to 1*/
    ix = 1;
    iy = 1;
    if (incx < 0)
        ix = (-n + 1) * incx + 1;
    if (incy < 0)
        iy = (-n + 1) * incy + 1;
    for (int i = 1; i <= n; i++)
    {
        dtemp += dx[ix-1] * dy[iy - 1];
        ix += incx;
        iy += incy;
    }
    ddot = dtemp;
    return ddot;
    /*code for both increments equal to 1
        c
        c
        c        clean - up loop*/
s20:
    int m = n % 5;
    if (m == 0)
        goto s40;
    for (int i = 1; i <= m; i++)
        dtemp += dy[i-1] * dx[i - 1];
    if (n < 5)
        goto s60;
s40:
    mp1 = m + 1;
    for (int i = mp1; i <= n; i += 5) {
        dtemp = dtemp + dx[i - 1] * dy[i - 1] + dx[i] * dy[i] +
            dx[i + 1] * dy[i + 1] + dx[i + 2] * dy[i + 2] + dx[i + 3] * dy[i + 3];
    }
s60:
    ddot = dtemp;
    return ddot;
}

double DSETRHSSTOP(double* B, double* WRK, double &EPSILON, int *&IPAR,vector<double> & Q1)
{
    double dsertrhstop;
    int LOCLEN = IPAR[3], STOPTYPE = IPAR[8];
    if ((STOPTYPE == 1) || (STOPTYPE == 4) || (STOPTYPE == 7))
        dsertrhstop = EPSILON;
    else if ((STOPTYPE == 2) || (STOPTYPE == 3) || (STOPTYPE == 5))
        dsertrhstop = EPSILON * PDNRM2(LOCLEN, B);
    else if (STOPTYPE == 6)
    {
        DIAGL(B, WRK, IPAR, Q1);
        dsertrhstop = EPSILON * PDNRM2(LOCLEN, WRK);
    }
    return dsertrhstop;
}

void PIMDSETPAR(int *IPAR, double* DPAR, int& LDA, int& N, int& BLKSZ, int& LOCLEN, int& BASISDIM,
    int& NPROCS, int &PROCID,int& PRECONTYPE,int& STOPTYPE,int& MAXIT,
    double EPSILON)
{
    /*c           PIM -- The Parallel Iterative Methods package
        c-------------------------------------------- -
        c
        c                      Rudnei Dias da Cunha
        c               Centro de Processamento de Dados,
        c         Universidade Federal do Rio Grande do Sul, Brasil
        c and
        c     Computing Laboratory, University of Kent at Canterbury, U.K.
        c
        c                          Tim Hopkins
        c     Computing Laboratory, University of Kent at Canterbury, U.K.
        c
        c----------------------------------------------------------------------
        c
        c  Description of parameter arrays
        c   IPAR(INPUT) : integer
        c     ipar(1) : lda(Leading dimension of a)
        c           2 : n(Number of rows / columns of a)
        c           3 : blksz(Size of block of data; used when data is
            c                       partitioned using cyclic mode)
        c           4 : loclen(Number of elements stored locally;
    c * PARALLEL: Equal to at least m / nprocs or
        c                                  n / procs depending if row or
        c                                  column partitioning is used or ,
        c                                  in the case of cyclic partitioning,
        c                                  it is a multiple of either
        c                                  m / (blksz * nprocs) or n / (blksz * nprocs).
        c * SEQUENTIAL: equal to n)
        c           5 : basisdim(Dimension of orthogonal basis, used in
            c                       GMRES)
        c           6 : nprocs(Number of processors)
        c           7 : procid(Processor identification)
        c           8 : precontype(Type of preconditioning; one of
            c                           0 : No preconditioning,
            c                           1 : Left preconditioning,
            c                           2 : Right preconditioning,
            c                           3 : Symmetric preconditioning)
        c           9 : stoptype(Type of stopping criteria used)
        c          10 : maxit(Maximum number of iterations allowed)
        c
        c   IPAR(OUTPUT) : integer
        c     ipar(11) : itno(Number of iterations executed)
        c          12 : status(On exit of iterative method, one of
            c                        0: converged to solution
            c - 1 : no convergence has been achieved
            c - 2 : "soft" - breakdown, solution may have
            c                           been found
            c - 3 : "hard" - breakdown, no solution)
        c - 4: conflict in preconditioner and stopping
        c                           criterion selected
        c - 5 : error in stopping criterion 3, r^ { T }z < 0)
        c          13 : steperr(If status is either - 2 or -3, it gives
            c                         the step number of the respective algorithm
            c                         where a breakdown has occurred.Refer to the
            c                         User's Guide for further information)
            c
            c   RPAR / DPAR(INPUT) : real / real * 8
            c     dpar(1) : epsilon(The value of epsilon for use in the
                c                        stopping criterion)
            c
            c   RPAR / DPAR(OUTPUT) : real / real * 8
            c     dpar(2) : exitnorm(The value of a norm of either the residual
                c                         vector or the difference between two
                c                         successive solution estimates according to
                c                         the value of stoptype)
            c           3,
            c           4 : smallest and largest eigenvalues of Q1AQ2(in the
                c               symmetric case) OR smallest and largest real parts
            c(in the nonsymmetric case)
            c           5,
            c           6 : smallest and largest imaginary parts(only in the
                c               nonsymmetric case)
            c

            C     ..Parameters ..*/
    IPAR[0] = LDA;
    IPAR[1] = N;
    IPAR[2] = BLKSZ;
    IPAR[3] = LOCLEN;
    IPAR[4] = BASISDIM;
    IPAR[5] = NPROCS;
    IPAR[6] = PROCID;
    IPAR[7] = PRECONTYPE;
    IPAR[8] = STOPTYPE;
    IPAR[9] = MAXIT;
    IPAR[10] = -1;
    IPAR[11] = -1;
    IPAR[12] = -1;

    DPAR[0] = EPSILON;
    DPAR[1] = -1.0e0;

    return;
}

void PART(int N,int P, int (& RANGES)[3][128])
{
    const int MXPROCS = 128;
    int D, DELTA = N / P, I, R, REM = N % P, S;
    if (REM == 0)
        S = 0;
    else
        S = 1;
    R = 1;
    for (I = 0; I < P; I++)
    {
        RANGES[0][I] = R;
        D = DELTA + S;
        RANGES[1][I] = D;
        R = R + D;
        RANGES[2][I] = R - 1;
        if (REM == 1)
            S = 0;
        else
            REM = REM - 1;
    }
    return;
}

void DINIT(int& N, double ALPHA, double* & DX, int INCX)
{
    int ix, m;
    if (N <= 0)
        return;
    if (INCX == 1)
        goto s20;
    ix = 1;
    if (INCX < 0)
        ix = (-N + 1) * INCX + 1;
    for (int i = 1; i <= N; i++)
    {
        DX[ix - 1] = ALPHA;
        ix += INCX;
    }
    return;
s20:
    m =N % 7;
    if (m == 0)
        goto s40;
    for (int i = 1; i <= m; i++)
        DX[i - 1] = ALPHA;
    if (N < 7)
        return;
s40:
    int MP1 = m + 1;
    for (int I = MP1; I <= N; I += 7)
    {
        DX[I - 1] = ALPHA;
        DX[I + 1] = ALPHA;
        DX[I + 2] = ALPHA;
        DX[I + 3] = ALPHA;
        DX[I + 4] = ALPHA;
        DX[I + 5] = ALPHA;
        DX[I] = ALPHA;
    }
    return;
}

void REPORT(string NAME, int* IPAR, double* DPAR, double ET, double *& X, int nout)
{
    std::cout << "\n" << NAME << ": N=" << setw(9) << IPAR[1] << "\n  "
        << "k=" << setw(6) << IPAR[10] << " status=" << setw(6) << IPAR[11]
        << " steperr=" << setw(6) << IPAR[12] << " norm=" << setw(16) << setprecision(9) << DPAR[1] << "\n  "
        << "Execution time = " << setw(15) << setprecision(8) << ET << endl;
    FILE2 << "\n" << NAME << ": N=" << setw(9) << IPAR[1] << "\n  "
        << "k=" << setw(6) << IPAR[10] << " status=" << setw(6) << IPAR[11]
        << " steperr=" << setw(6) << IPAR[12] << " norm=" << setw(16) << setprecision(9) << DPAR[1] << "\n  "
        << "Execution time = " << setw(15) << setprecision(8) << ET << endl;
}

void DVPROD(int n, vector<double>& dx, int incx, double* &dy, int incy)
{
    if (n <= 0)
        return;
    int ix, iy, m,mp1;
    if (incx == 1 && incy == 1)
        goto s20;
    ix = 1;
    iy = 1;
    if (incx < 0)
        ix = (-n + 1) * incx + 1;
    if (incy < 0)
        iy = (-n + 1) * incx + 1;
    for (int i = 1; i <= n; i++)
    {
        dy[iy - 1] *= dx[ix - 1];
        ix += incx;
        iy += incy;
    }
    return;
s20:
    m = n % 4;
    if (m == 0)
        goto s40;
    for (int i = 1; i <= m; i++)
        dy[-1] *= dx[i - 1];
    if (n < 4)
        return;
s40:
    mp1 = m + 1;
    for (int i = mp1; i <= n; i += 4)
    {
        dy[i - 1] *= dx[i - 1];
        dy[i] *= dx[i];
        dy[i + 1] *= dx[i + 1];
        dy[i + 2] *= dx[i + 2];
    }
    return;
}
void PROGRESS(int& LOCLEN, int& ITNO, double& NORMRES, double X,double RES, double TRUERES)
{
    return;
}

void PDSUM(int ISIZE, double* x)
{
    double wrk[1];

    MPI_Allreduce(x, wrk, ISIZE, MPI_DOUBLE_PRECISION, MPI_SUM,
        MPI_COMM_WORLD);
    DCOPY(ISIZE, &wrk[0], 1, x, 1);
    return;
}
void DLAG2(vector<double>& Y, double& X, long long& N, int& LS)
{
    double T, B3, Q1;
    int LM = 0;
d10:
    LM = LM + 1;
    B3 = exp(-0.5e0 * X / LM);
    if (B3 < 1e-250) goto d10;
    /*c        print*, ' LM ', LM*/
    Q1 = B3;
    Y[0] = Q1;
    if (N - 1 <= 0)
        goto d1;
    else
        goto d2;
d1:
    for (int J = 2; J <= LM; J++)
    {
        Y[0] = Y[0] * Q1;
    }
    return;
d2:
    Y[1] = (LS + 1.0e0 - X) * Q1;
    T = LS - 1.0e0 - X;
    if (N - 2 <= 0)
        goto d7;
    else
        goto d3;
d7:
    for (int J = 2; J <= LM; J++)
    {
        Y[0] = Y[0] * Q1;
        Y[1] = Y[1] * Q1;
    }
    return;
d3:
    int N2 = 1;
    for (int I = 2; I <= N - 1; I++)
    {
        Y[I] = Y[I - 1] - Y[I - 2] + Y[I - 1] + (T * Y[I - 1] - (LS - 1.0e0) * Y[I - 2]) / double(I);
        if (Y[I] < 1e250)
            continue;
        if (N2 >= LM) continue;
        for (int J = 1; J <= I + 1; J++)
            Y[J - 1] = Y[J - 1] * Q1;
        N2 = N2 + 1;
    }
    for (int K = N2 + 1; K <= LM; K++)
        for (int I = 1; I <= N; I++)
            Y[I - 1] = Y[I - 1] * Q1;
    return;
}

void DLAG1(vector<double>& Y, double& X, long long& N, int& LS)
{
    double T, B3, Q1;
    int LM = 0;
d10:
    LM = LM + 1;
    B3 = exp(-0.5e0 * X / LM);
    if (B3 < 1e-250) goto d10;
    /*c        print*, ' LM ', LM*/
    Q1 = pow(X , 1.0e0 * LS / LM) * B3;
    Y[0] = Q1;
    if (N - 1 <= 0)
        goto d1;
    else
        goto d2;
d1:
    for (int J = 2; J <= LM; J++)
    {
        Y[0] = Y[0] * Q1;
    }
    return;
d2:
    Y[1] = (LS + 1.0e0 - X) * Q1;
    T = LS - 1.0e0 - X;
    if (N - 2 <= 0)
        goto d7;
    else
        goto d3;
d7:
    for (int J = 2; J <= LM; J++)
    {
        Y[0] = Y[0] * Q1;
        Y[1] = Y[1] * Q1;
    }
    return;
d3:
    int N2 = 1;
    for (int I = 2; I <= N - 1; I++)
    {
        Y[I] = Y[I - 1] - Y[I - 2] + Y[I - 1] + (T * Y[I - 1] - (LS - 1.0e0) * Y[I - 2]) / double(I);
        if (Y[I] < 1e250)
            continue;
        if (N2 >= LM) continue;
        for (int J = 1; J <= I + 1; J++)
            Y[J - 1] = Y[J - 1] * Q1;
        N2 = N2 + 1;
    }
    for (int K = N2 + 1; K <= LM; K++)
        for (int I = 1; I <= N; I++)
            Y[I - 1] = Y[I - 1] * Q1;
    return;
}
void KOEFlag(vector<double>& FN, vector<double>& CR, vector<double>& QN, double& RE, long long& LAGMAX, int& LS, double& epsLag, long long& NW, int& KER)
{
    //Определение требуемого кол-ва коэфф. Лагерра для заданной точности = epsLag

    int NDL = 50;
    double re1, xk;
    int hernya = 0;
    DLAG2(CR, RE, LAGMAX, hernya);

    for (int i = 1; i <= LAGMAX; i++)
    {
        RE = 0.0e0;
        for (int j = 1; j <= i; j++)
        {
            xk = 1.0e0;
            for (int l = 1; l <= LS; l++)
                xk = xk * static_cast<double>(j + l - 1);
            RE = RE + FN[j - 1] * CR[i - j] * xk;
        }
        re1 = 0.0e0;
        for (int j = 1; j <= i - 1; j++)
        {
            xk = 1.0e0;
            for (int l = 1; l <= LS; l++)
                xk = xk * static_cast<double>(j + l - 1);
            re1=re1+FN[j-1] * CR[i - j-1] * xk;
        }
        xk = 1.0e0;
        for (int l = 1; l <= LS; l++)
            xk = xk * static_cast<double>(i + l - 1);
        QN[i - 1] = (RE - re1) / xk;
    }
    KER = 0;
    NW = 0;
    RE = 0.0e0;
    int L;
    for (L = 1; L <= LAGMAX; L = L + NDL)
    {
        NW = NW + NDL;
        if (NW > LAGMAX)
            goto koe3;
        re1 = 0.0e0;
        for (int I = 1; I <= NDL; I++)
        {
            if (abs(QN[L + I-2]) > re1)
                re1 = abs(QN[L + I-2]);
        }
        if (re1 < RE)
            goto koe4;
        RE = re1;
    }
    goto koe3;
koe4:
    if (RE == 0e0)
        exit(1);
    xk = re1 / RE;
    if (xk <= epsLag) goto koe5;
    for (int K = L + NDL; K <= LAGMAX; K += NDL)
    {
        NW = NW + NDL;
        if (NW > LAGMAX)
            goto koe3;
        re1 = 0.0e0;
        for (int I = 1; I <= NDL; I++)
        {
            if (abs(QN[K + I-2]) > re1)
                re1 = abs(QN[K + I-2]);
        }
        xk = re1 / RE;
        if (xk <= epsLag)
            goto koe5;
    }
koe3:
    NW = LAGMAX;
    KER = 1;
koe5:
    return;
}

double PDNRM2(int& LOCLEN, double* U)
{
    double WRK[1];
    double PSUM = DDOT(LOCLEN, U, 1, U, 1);
    MPI_Allreduce(&PSUM, WRK, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
        MPI_COMM_WORLD);
    return sqrt(WRK[0]);
}
double Timer() {
    double t;
    t = (double)clock() / CLOCKS_PER_SEC;
    return t;
}
void ZAPIS(ofstream& NDIR, int NTRAS, int NSAMP, vector<vector<double>> REG1, vector<double> REG2)
{
    for (int i = 1; i <= NTRAS; i++)
    {
        for (int j = 1; j <= NSAMP; j++)
            REG2[j - 1] = REG1[i - 1][j - 1];
        for (int j = 0; j < NSAMP; j++)
            NDIR << REG2[j] << endl;
    }
    return;
}

void STOPCRIT(int& KMV, double& ETMV1, vector<double> d2pmlx, double*& B, double* R, double* RTRUE,
    double*& X, double* XOLD, double* WRK, double RHSSTOP,
    int CNVRTX, double EXITNORM, int STATUS, int*& IPAR,
    vector<double>& Q2, double*& SP, double*& SL, vector<double>& XS, vector<double>& UZ1, vector<double>& UZ2, vector<double>& UX1, vector<double>& UX2)
{
    double DOTS[1];
    double ZERO = 0.0;
    int LOCLEN = IPAR[3];
    int PRECONTYPE = IPAR[7];
    int STOPTYPE = IPAR[8];
    double ONE = 1.0;
    if ((STOPTYPE == 1) || (STOPTYPE == 2) || (STOPTYPE == 3))
    {
       // Compute true residual if needed
        DCOPY(LOCLEN, &B[0], 1, RTRUE, 1);

        if ((PRECONTYPE == 2) || (PRECONTYPE == 3))
        {
            DIAGR(X, WRK, IPAR, Q2);
            if (CNVRTX == 1)
            {
                MATVEC(d2pmlx, WRK, XOLD, IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                MATVEC(d2pmlx, WRK, XOLD, IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DAXPY(LOCLEN, -ONE, WRK, 1, RTRUE, 1);
            }
            else
            {
                MATVEC(d2pmlx, WRK, XOLD, IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                    KMV, ETMV1);
                DAXPY(LOCLEN, -ONE, XOLD, 1, RTRUE, 1);
            }
        }
        else if (CNVRTX == 1)
        {
            MATVEC(d2pmlx, X, XOLD, IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                KMV, ETMV1);
            MATVEC(d2pmlx, XOLD, WRK, IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                KMV, ETMV1);
            DAXPY(LOCLEN, -ONE, WRK, 1, RTRUE, 1);
        }
        else
        { 
            MATVEC(d2pmlx, X, WRK, IPAR, SP, SL, XS, UZ1, UZ2, UX1, UX2,
                KMV, ETMV1);
            DAXPY(LOCLEN, -ONE, WRK, 1, RTRUE, 1);
        }
    }
    if ((STOPTYPE == 1) || (STOPTYPE == 2))
    {
        EXITNORM = PDNRM2(LOCLEN, RTRUE);
        if (EXITNORM < RHSSTOP)
            STATUS = 0;
        else
            STATUS = -99;
    }
    else if (STOPTYPE == 3)
    {
        DOTS[0] = DDOT(LOCLEN, RTRUE, 1, R, 1);
        PDSUM(1, DOTS);
        if (DOTS[0] < ZERO)
            STATUS = -5;
        return;

        EXITNORM = sqrt(DOTS[0]);
        if (EXITNORM < RHSSTOP)
            STATUS = 0;
        else
            STATUS = -99;
    }
    else if ((STOPTYPE==4)||(STOPTYPE==5)||
         (STOPTYPE==6))
    {
        EXITNORM = PDNRM2(LOCLEN, R);
        if (EXITNORM < RHSSTOP)
            STATUS = 0;
        else
        {
            STATUS = -99;
        }
    }
    else if (STOPTYPE == 7)
    {
        DCOPY(LOCLEN, &X[0], 1, WRK, 1);
        DAXPY(LOCLEN, -ONE, XOLD, 1, WRK, 1);
        if (EXITNORM < RHSSTOP)
            STATUS = 0;
        else
        {
            STATUS = -99;
        }
    }
    return;
}