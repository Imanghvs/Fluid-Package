#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

class complex_number
{
public:
    double real;
    double img;
    complex_number(){
        real = 0;
        img = 0;
    }
    complex_number(double real1,double img1){
        real = real1;
        img = img1;
    }
    friend complex_number summ(complex_number, complex_number);
    friend complex_number mult(complex_number, complex_number);
    friend complex_number poww(complex_number,double);
};
complex_number summ(complex_number c1, complex_number c2) {
    complex_number summed;
    summed.real = c1.real + c2.real;
    summed.img = c1.img + c2.img;
    return summed;
}
complex_number mult(complex_number c1, complex_number c2) {
    complex_number multiplied;
    multiplied.real = c1.real * c2.real - c1.img * c2.img ;
    multiplied.img = c1.real * c2.img + c1.img * c2.real ;
    return multiplied;
}
complex_number poww(complex_number c1,double i) {
    double pi = 3.1415926535897;
    complex_number powed;
    double r = sqrt(c1.real*c1.real + c1.img * c1.img);
    double theta = atan2(c1.img,c1.real);
    if (theta < 0) {
        theta = theta + 2 * pi;
    }
    double newr = pow(r,i);
    double newtheta = theta*i;
    powed.real = newr * cos(newtheta);
    powed.img = newr * sin(newtheta);
    return powed;
}
class component
{
private:
    double Gm_Res;
    double Hm_Res;
    double Sm_Res;
    double Gm_ig;
    double Hm_ig;
    double Sm_ig;
    double quality;
    double P_sat;
    double f_omega;
    double omega;
    double Tc;
    double Pc;
    double w = -1;
    double u = 2;
    double a;
    double T;
    double P;
    double R = 8.3145;
    double Vm;
    double Vm_l;
    double Vm_g;
    double V;
    double Gm;
    double Gm_l;
    double Gm_g;
    double G;
    double Hm;
    double Hm_l;
    double Hm_g;
    double Hm_fg;
    double H;
    double Z;
    double Um;
    double Um_l;
    double Um_g;
    double U;
    double Sm;
    double Sm_l;
    double Sm_g;
    double S;
    double m_dot; //mole/sec
    double b;
    double Z_l;
    double Z_g;
    double fi;
    double fi_l;
    double fi_g;
    double Cpm_ig = 2.5 * R;
    double Cvm_ig = Cpm_ig - R;
    double Cpm;
    double Cvm;
    double Cp;
    double Cv;
    double T_ref = 273.15;
    double P_ref = 101325.1;
    std::string phase;
public:
    void Setu(int u1) {u = u1;}
    void Setw(int w1) {w = w1;}
    void SetPhase(std::string phase1) {phase = phase1;}
    void SetOmega(double omega1) {omega = omega1;}
    void SetTc(double Tc1) {Tc = Tc1;}
    void SetPc(double Pc1) {Pc = Pc1;}
    void Setf_omega(double f_omega1) {f_omega=f_omega1;}
    void SetQuality(double Quality) {quality = Quality;}
    std::string GetPhase() {return phase;}
    double GetOmega() {return omega;}
    double GetTc() {return Tc;}
    double GetPc() {return Pc;}
    double Getf_omega() {return f_omega;}
    double GetQuality() {return quality;}
    void Seta(double a1) {a = a1;}
    void Setb(double b1) {b = b1;}
    void SetT(double T1) {T = T1;}
    void SetP(double P1) {P = P1;}
    void SetGm(double Gm1) {Gm = Gm1;}
    void SetGm_l(double Gm_l1) {Gm_l = Gm_l1;}
    void SetGm_g(double Gm_g1) {Gm_g = Gm_g1;}
    void SetG(double G1) {G = G1;}
    void SetUm(double Um1) {Um = Um1;}
    void SetUm_l(double Um_l1) {Um_l = Um_l1;}
    void SetUm_g(double Um_g1) {Um_g = Um_g1;}
    void SetU(double U1) {U = U1;}
    void SetVm(double Vm1) {Vm = Vm1;}
    void SetVm_l(double Vm_l1) {Vm_l = Vm_l1;}
    void SetVm_g(double Vm_g1) {Vm_g = Vm_g1;}
    void SetV(double V1) {Vm = V1;}
    void SetHm(double Hm1) {Hm = Hm1;}
    void SetHm_g(double Hm_g1) {Hm_g = Hm_g1;}
    void SetHm_l(double Hm_l1) {Hm_l = Hm_l1;}
    void SetCpm_ig(double Cpm_ig1) {Cpm_ig = Cpm_ig1;}
    void Setmdot(double mdot) {m_dot = mdot;}
    double SetCpm(double Cpm1) {Cpm = Cpm1;}
    void SetCvm(double Cvm1) {Cvm = Cvm1;}
    void SetCp(double Cp1) {Cp = Cp1;}
    void SetCv(double Cv1) {Cv = Cv1;}
    double Geta() {return a;}
    double GetT() {return T;}
    double GetP() {return P;}
    double GetG() {return G;}
    double GetH() {return H;}
    double GetU() {return U;}
    double GetV() {return V;}
    double GetS() {return S;}
    double GetCpm() {return Cpm;}
    double GetCvm() {return Cvm;}
    double GetCp() {return Cp;}
    double GetCv() {return Cv;}
    double Getm_dot() {return m_dot;}
    double GetGm() {return Gm;}
    double GetGm_l() {return Gm_l;}
    double GetGm_g() {return Gm_g;}
    double GetUm() {return Um;}
    double GetUm_l() {return Um_l;}
    double GetUm_g() {return Um_g;}
    double GetVm() {return Vm;}
    double GetVm_l() {return Vm_l;}
    double GetVm_g() {return Vm_g;}
    double GetHm() {return Hm;}
    double GetHm_l() {return Hm_l;}
    double GetHm_g() {return Hm_g;}
    double GetHm_fg() {return Hm_fg;}
    double GetSm() {return Sm;}
    double GetSm_l() {return Sm_l;}
    double GetSm_g() {return Sm_g;}
    double GetCpm_ig() {return Cpm_ig;}
    double GetP_sat() {return P_sat;}
    double GetZ() {return Z;}
    double GetZ_l() {return Z_l;}
    double GetZ_g() {return Z_g;}
    double Getfi_l() {return fi_l;}
    double Getfi_g() {return fi_g;}
    std::vector<double> Z_vector;
    std::vector<double> fi_vector;
    std::vector<double> SolveZ()
    {
        std::vector<double> Zs;
        Setf_omega((.37464 + 1.54226 * omega - 0.26992 * (pow(omega, 2))));
        double Tr = T/Tc;
        Seta((0.45724*R*R*Tc*Tc/Pc*pow((1+f_omega*(1-pow(Tr,0.5))),2)));
        Setb(0.07780*R*Tc/Pc);
        double A_star = a*P/(R*R*T*T);
        double B_star = b*P/(R*T);

        // daraje 3

        double aa = 1;
        double bb = - (1 + B_star - u*B_star);
        double cc = A_star + w * B_star * B_star - u * B_star - u * B_star * B_star;
        double dd = - A_star * B_star - w * B_star * B_star - w * pow(B_star,3);
        double delta_0 = bb*bb - 3*aa*cc;
        double delta_1 = 2 * pow(bb,3) - 9*aa*bb*cc + 27* aa*aa*dd;
        double zir_radical = delta_1 * delta_1 - 4 * pow(delta_0,3);
        complex_number kesi;
        kesi.real = -0.5;
        kesi.img = sqrt(3.)/2;
        complex_number com_aa;
        com_aa.real = aa;
        com_aa.img = 0;
        complex_number com_bb;
        com_bb.real = bb;
        com_bb.img = 0;
        complex_number com_delta_0;
        com_delta_0.real = delta_0;
        com_delta_0.img = 0;
        complex_number com_delta_1;
        com_delta_1.real = delta_1;
        com_delta_1.img = 0;
        complex_number com_zir_radical = complex_number(zir_radical,0);
        complex_number com_C = poww(((mult(summ(poww(com_zir_radical,0.5),com_delta_1),poww(complex_number(2,0),-1)))),0.3333333);
        // -1/(3a) -> aval
        complex_number manfi_se;
        manfi_se.real = -3;
        manfi_se.img = 0;
        complex_number aval = poww(mult(manfi_se,com_aa),-1);
        // kesi * C
        complex_number dovom = mult(kesi,com_C);
        // delta_sefr / folan
        complex_number sevvom = mult(com_delta_0,poww(dovom,-1));
        // cheharom -> delta_0/C
        complex_number cheharom = mult(com_delta_0,poww(com_C,-1));
        complex_number Z_1 = mult(aval,summ(com_bb,summ(com_C,cheharom)));
        complex_number Z_2 = mult(summ(com_bb, summ(sevvom,dovom)),aval);
        complex_number Z_3 = mult(summ(com_bb, summ(mult(sevvom,poww(kesi,-1)),mult(dovom,kesi))),aval);
        if (std::abs(Z_1.img)<(0.001 * std::abs(Z_1.real)))
        {
            Zs.push_back(Z_1.real);
        }
        if (std::abs(Z_2.img)<(0.001 * std::abs(Z_2.real)))
        {
            Zs.push_back(Z_2.real);
        }
        if (std::abs(Z_3.img)<(0.001 * std::abs(Z_3.real)))
        {
            Zs.push_back(Z_3.real);
        }
        sort(Zs.begin(),Zs.end());
        if (Zs.size()>1)
        {
            Zs.erase(Zs.begin()+1);
        }
        Z_vector = Zs;
        return Zs;
    }
    //calculate fugacity factors
    std::vector <double> FugFactors() {
        std::vector<double> Z_ha = SolveZ();
        std::vector<double> FugFacs;
        Setf_omega((.37464 + 1.54226 * omega - 0.26992 * (pow(omega, 2))));
        double Tr = T/Tc;
        Seta((0.45724*R*R*Tc*Tc/Pc*pow((1+f_omega*(1-pow(Tr,0.5))),2)));
        Setb(0.07780*R*Tc/Pc);
        double A_star = a*P/(R*R*T*T);
        double B_star = b*P/(R*T);
        for (double Z : Z_ha)
        {
            FugFacs.push_back(exp(Z - 1 - log(Z - B_star) - A_star / (sqrt(8.) * B_star) * log((Z + (1 + sqrt(2.)) * B_star) / (Z + (1 - sqrt(2.)) * B_star))));
        }
        if (Z_ha.size()>1 && (((FugFacs[0]-FugFacs[1])/FugFacs[0])>0.001))
        {
            SetPhase("Gas");
        }
        else if (Z_ha.size()>1 && (((FugFacs[1]-FugFacs[0])/FugFacs[1])>0.001))
        {
            SetPhase("Liquid");
        }
        else if (Z_ha.size()>1)
        {
            SetPhase("TwoPhases");
        }
        else
        {
            SetPhase("Gas");
        }
        fi_vector = FugFacs;
        return FugFacs;
    }
    void IntensivePropertyCalculator() {
        fi_vector = FugFactors();
//        Z_vector = SolveZ();
        Setf_omega((.37464 + 1.54226 * omega - 0.26992 * (pow(omega, 2))));
        double Tr = T / Tc;
        Seta((0.45724 * R * R * Tc * Tc / Pc * pow((1 + f_omega * (1 - pow(Tr, 0.5))), 2)));
        Setb(0.07780 * R * Tc / Pc);
        double A_star = a * P / (R * R * T * T);
        double B_star = b * P / (R * T);
        if (GetPhase() == "TwoPhases")
        {
            Z_l = Z_vector[0];
            Z_g = Z_vector[1];
            fi_l = fi_vector[0];
            fi_g = fi_vector[1];
            P_sat = P;
            double Gm_Res_l = R * T * log(fi_l);
            double Gm_Res_g = R * T * log(fi_g);
            double ln_yechizi_l = log((Z_l + (1 + sqrt(2.)) * B_star) / (Z_l + (1 - sqrt(2.)) * B_star));
            double ln_yechizi_g = log((Z_g + (1 + sqrt(2.)) * B_star) / (Z_g + (1 - sqrt(2.)) * B_star));
            double moshtagh_a = -a * f_omega / sqrt(pow((1 + f_omega * (1 - sqrt(Tr))), 2) * T * Tc);
            double Sm_Res_l = R * log(Z_l - B_star) - 1 / (sqrt(8) * b * R * T) * moshtagh_a * ln_yechizi_l;
            double Sm_Res_g = R * log(Z_g - B_star) - 1 / (sqrt(8) * b * R * T) * moshtagh_a * ln_yechizi_g;
            double Hm_Res_l = Gm_Res_l + T * Sm_Res_l;
            double Hm_Res_g = Gm_Res_g + T * Sm_Res_g;
            Hm_ig = Cpm_ig * (T - T_ref);
            Sm_ig = Cpm_ig * log(T / T_ref) - R * log(P / P_ref);
            Gm_ig = Hm_ig - T * Sm_ig;
            Gm_l = Gm_Res_l + Gm_ig;
            Gm_g = Gm_Res_g + Gm_ig;
            Gm = Gm_l * (1 - quality) + Gm_g * quality;
            Hm_l = Hm_Res_l + Hm_ig;
            Hm_g = Hm_Res_g + Hm_ig;
            Hm_fg = Hm_g - Hm_l;
            Hm = Hm_l * (1 - quality) + Hm_g * quality;
            Sm_l = Sm_Res_l + Sm_ig;
            Sm_g = Sm_Res_g + Sm_ig;
            Sm = Sm_l * (1 - quality) + Sm_g * quality;
            Vm_l = Z_l * R * T / P;
            Vm_g = Z_g * R * T / P;
            Vm = Vm_l * (1 - quality) + Vm_g * quality;
            Um_l = Hm_l - P*Vm;
            Um_g = Hm_g - P*Vm;
            Um = Um_l * (1 - quality) + Um_g * quality;
        }
        else if (GetPhase() == "Gas")
        {
            if (Z_vector.size() > 1)
            {
                Z_g = Z_vector[1];
                Z = Z_g;
                fi_g = fi_vector[1];
                fi = fi_g;
            }
            else if (Z_vector.size() == 1)
            {
                Z_g = Z_vector[0];
                Z = Z_g;
                fi_g = fi_vector[0];
                fi = fi_g;
            }
            Gm_Res = R * T * log(fi);
            double ln_yechizi = log((Z + (1 + sqrt(2.)) * B_star) / (Z + (1 - sqrt(2.)) * B_star));
            double moshtagh_a = - a * f_omega / sqrt(pow((1 + f_omega * (1 - sqrt(Tr))), 2) * T * Tc);
            Sm_Res = R * log(Z - B_star) - 1 / (sqrt(8) * b * R * T) * moshtagh_a * ln_yechizi;
            Hm_Res = Gm_Res + T * Sm_Res;
            Hm_ig = Cpm_ig * (T - T_ref);
            Sm_ig = Cpm_ig * log(T / T_ref) - R * log(P / P_ref);
            Gm_ig = Hm_ig - T * Sm_ig;
            Gm = Gm_Res + Gm_ig;
            Hm = Hm_Res + Hm_ig;
            Sm = Sm_Res + Sm_ig;
            Vm = Z * R * T / P;
            Um = Hm - P*Vm;
        }
        else if (GetPhase() == "Liquid")
        {
            Z_l = Z_vector[0];
            Z = Z_l;
            fi_l = fi_vector[0];
            fi = fi_l;
            Gm_Res = R * T * log(fi);
            double ln_yechizi = log((Z + (1 + sqrt(2.)) * B_star) / (Z + (1 - sqrt(2.)) * B_star));
            double moshtagh_a = -a * f_omega / sqrt(pow((1 + f_omega * (1 - sqrt(Tr))), 2) * T * Tc);
            Sm_Res = R * log(Z - B_star) - 1 / (sqrt(8) * b * R * T) * moshtagh_a * ln_yechizi;
            Hm_Res = Gm_Res + T * Sm_Res;
            Hm_ig = Cpm_ig * (T - T_ref);
            Sm_ig = Cpm_ig * log(T / T_ref) - R * log(P / P_ref);
            Gm_ig = Hm_ig - T * Sm_ig;
            Gm = Gm_Res + Gm_ig;
            Hm = Hm_Res + Hm_ig;
            Sm = Sm_Res + Sm_ig;
            Vm = Z * R * T / P;
        }
    }
    void PropertyCalculator()
    {
        IntensivePropertyCalculator();
        V = Vm * m_dot;
        U = Um * m_dot;
        H = Hm * m_dot;
        G = Gm * m_dot;
        S = Sm * m_dot;
    }
};
component PrCalculator(component Component1, std::string parametrha, double parametr1, double parametr2)
{
    if (parametrha == "TP")
    {
        Component1.SetT(parametr1);
        Component1.SetP(parametr2);
        Component1.PropertyCalculator();
        component TempComponent = Component1;
        TempComponent.SetT(Component1.GetT()+0.1);
        if (Component1.GetPhase() == "TwoPhases" || TempComponent.GetPhase() == "TwoPhases")
        {
           Component1.SetCpm(99999999999999999.);
           Component1.SetCvm(99999999999999999.);
           Component1.SetCp(99999999999999999. * Component1.Getm_dot());
           Component1.SetCv(99999999999999999. * Component1.Getm_dot());
        }
        else
        {
           TempComponent.PropertyCalculator();
           Component1.SetCpm((TempComponent.GetHm()-Component1.GetHm())/(TempComponent.GetT()-Component1.GetT()));
        }
    }
    if (parametrha == "TV") {
        //newton
        Component1.SetT(parametr1);
        double Vm = parametr2;
        double Vmi;
        std::vector<double> hads;
        std::vector<double> Vms;
        for (int i = -8; i < 2; i++) {
            hads.push_back(pow(3, i) * 8.3145 * Component1.GetT() / Vm);
            Component1.SetP(pow(3, i) * 8.3145 * Component1.GetT() / Vm);
            Component1.PropertyCalculator();
            Vms.push_back(Component1.GetVm());
        }
        int indexing = 0;
        for (int i = 0; i < 10; i++) {
            int points = 0;
            for (int j = 0; j < 10; j++) {
                if (i == j) {
                    continue;
                } else if (std::abs(Vms[i] - Vm) / Vms[i] <= std::abs(Vms[j] - Vm) / Vms[j]) {
                    points++;
                }
            }
            if (points == 9) {
                indexing = i;
            }
        }
        Component1.SetP(hads[indexing]);
        Component1.PropertyCalculator();
        double Pi = hads[indexing];
        component TempComp = Component1;
        double moshtagh;
        int t = 0;
        while (true) {
            Component1.SetP(Pi);
            Component1.PropertyCalculator();
            Vmi = Component1.GetVm();
            TempComp.SetP(Pi + 10);
            TempComp.PropertyCalculator();
            moshtagh = (TempComp.GetVm() - Vmi) / 10;
            Pi = (Vm - Vmi) / moshtagh + Pi;
            if ((std::abs(Vm - Vmi) / Vm < 0.01)) {
                break;
            }
            if (Pi < 10 || Pi > pow(10.0, 9)) {
                Pi = (rand() % 100 + .0000001) * 8.3145 * Component1.GetT() / Vm;
            }
        };
    }
    if (parametrha == "TH")
    {
        //newton
        Component1.SetT(parametr1);
        double Hm = parametr2;
        double Hmi;
        std::vector<double> hads;
        std::vector<double> Hms;
        for (int i = 1; i < 11; i++) {
            hads.push_back(pow(10,i));
            Component1.SetP(pow(10,i));
            Component1.PropertyCalculator();
            Hms.push_back(Component1.GetHm());
        }
        int indexing = 0;
        for (int i = 0; i < 10; i++) {
            int points = 0;
            for (int j = 0; j < 10; j++) {
                if (i == j) {
                    continue;
                } else if ((std::abs((Hms[i] - Hm) / Hms[i])) < (std::abs((Hms[j] - Hm) / Hms[j])))
                {
                    points++;
                }
            }
            if (points == 9) {
                indexing = i;
            }
        }
        Component1.SetP(hads[indexing]);
        Component1.PropertyCalculator();
        double Pi = hads[indexing];
        component TempComp = Component1;
        double moshtagh;
        int t = 0;
        while (true) {
            Component1.SetP(Pi);
            Component1.PropertyCalculator();
            if (Component1.GetPhase()=="TwoPhases")
            {
                if (Hm > Component1.GetHm_g())
                {
                    Hmi = Hmi + 2*Component1.GetHm_fg();
                }
                else if (Hm < Component1.GetHm_l())
                {
                    Hmi = Hmi - 2* Component1.GetHm_fg();
                }
                else
                {
                    double Quality = (Hm - Component1.GetHm_l())/Component1.GetHm_fg();
                    Component1.SetQuality(Quality);
                }
            }
            Hmi = Component1.GetHm();
            TempComp.SetP(Pi + 10);
            TempComp.PropertyCalculator();
            moshtagh = (TempComp.GetHm() - Hmi) / 10;
            Pi = (Hm - Hmi) / moshtagh + Pi;
            if ((std::abs(Hm - Hmi) / Hm < 0.01)) {
                break;
            }
            if (Pi < 10 || Pi > pow(10.0, 9))
            {
                Pi = (rand() % 1000000000 + .0000001);
            }
        }
        };
    return Component1;

};
int main() {
    /*component Water;
    Water.SetOmega(0.344);
    Water.SetTc(647.3);
    Water.SetPc(22120000);
    Water.Setmdot(100);
    Water.SetQuality(.5);
    Water = CalculateIntensives(Water,"TH",300,-1600);
    std::cout<<"Z_l = "<<Water.GetZ_l()<<", Z_g = "<<Water.GetZ_g()<<std::endl;
    std::cout<<"Z = "<<Water.GetZ()<<std::endl;
    std::cout<<"Hm = "<<Water.GetHm()<<std::endl;
    std::cout<<"P = "<<Water.GetP()<<std::endl;
    std::cout<<"Phase: "<<Water.GetPhase()<<std::endl;*/
    std::cout<<"Welcome!"<<std::endl;
    std::cout<<"Please enter the value of Omega: "<<std::endl;
    double Omega;
    double Tc;
    double Pc;
    double mdot;
    double Quality;
    std::cin>>Omega;
    component NewComponent;
    NewComponent.SetOmega(Omega);
    std::cout<<"Please enter the critical temperature (Kelvins) : "<<std::endl;
    std::cin>>Tc;
    std::cout<<"Please enter the critical pressure (Pa): "<<std::endl;
    std::cin>>Pc;
    std::cout<<"Please enter the molar flowrate (mole/s . It's optional. But if you want to know extensive properties, it would be necessary. you can enter a random number.): "<<std::endl;
    std::cin>>mdot;
    std::cout<<"Please enter Cp_ig: (J/(K.mole))"<<std::endl;
    NewComponent.SetTc(Tc);
    NewComponent.SetPc(Pc);
    NewComponent.Setmdot(mdot);
    double T;
    double P;
    double Vm;
    double Hm;
    while (true)
    {
        std::cout << "Which two variables do you have?" << std::endl;
        std::cout << "1: T and P          2: T and Vm         3: T and Hm" << std::endl;
        int which;
        std::string which_two;
        std::cin >> which;
        if (which == 1)
        {
            which_two = "TP";
            std::cout<<"Please enter the value of T (Kelvins) : "<<std::endl;
            std::cin>>T;
            std::cout<<"Please enter the value of P (Pa) : "<<std::endl;
            std::cin>>P;
            NewComponent = PrCalculator(NewComponent,which_two,T,P);
            if (NewComponent.GetPhase()=="TwoPhases")
            {
                std::cout<<"It has two phases. We recommend entering the quality (anyway you will receive the right number for liquid and gas phases separately."<<std::endl;
                NewComponent.SetQuality(Quality);
                NewComponent = PrCalculator(NewComponent,which_two,T,P);
            }
            break;
        }
        else if (which == 2) {
            which_two = "TV";
            std::cout<<"Please enter the value of T (Kelvins) : "<<std::endl;
            std::cin>>T;
            std::cout<<"Please enter the value of Vm: (Cubic meters per mole)"<<std::endl;
            std::cin>>Vm;
            NewComponent = PrCalculator(NewComponent,which_two,T,Vm);
            if (NewComponent.GetPhase()=="TwoPhases")
            {
                std::cout<<"It has two phases. We recommend entering the quality (anyway you will receive the right number for liquid and gas phases separately."<<std::endl;
                std::cin>>Quality;
                NewComponent.SetQuality(Quality);
                NewComponent = PrCalculator(NewComponent,which_two,T,P);
            }
            break;
        }
        else if (which == 3) {
            which_two = "TH";
            std::cout<<"Please enter the value of T (Kelvins) : "<<std::endl;
            std::cin>>T;
            std::cout<<"Please enter the value of Hm (J/mole) : "<<std::endl;
            std::cin>>Hm;
            NewComponent = PrCalculator(NewComponent,which_two,T,Hm);
            break;
        }
        else
        {
            std::cout<<"Wrong number! Try again!"<<std::endl;
        }
    };
    std::cout<<"T = "<<NewComponent.GetT()<<std::endl;
    std::cout<<"P = "<<NewComponent.GetP()<<std::endl;
    std::cout<<"Vm = "<<NewComponent.GetVm()<<std::endl;
    std::cout<<"Hm = "<<NewComponent.GetHm()<<std::endl;
    std::cout<<"Um = "<<NewComponent.GetUm()<<std::endl;
    std::cout<<"Gm = "<<NewComponent.GetGm()<<std::endl;
    std::cout<<"V = "<<NewComponent.GetV()<<std::endl;
    std::cout<<"H = "<<NewComponent.GetH()<<std::endl;
    std::cout<<"U = "<<NewComponent.GetU()<<std::endl;
    std::cout<<"Cpm = "<<NewComponent.GetCpm()<<std::endl;

    std::cout<<"G = "<<NewComponent.GetG()<<std::endl;
    if (NewComponent.GetPhase()=="TwoPhases")
    {
        std::cout<<"Vm_l = "<<NewComponent.GetVm_l()<<std::endl;
        std::cout<<"Vm_g = "<<NewComponent.GetVm_g()<<std::endl;
        std::cout<<"Hm_l = "<<NewComponent.GetHm_l()<<std::endl;
        std::cout<<"Hm_g = "<<NewComponent.GetHm_g()<<std::endl;
        std::cout<<"Um_l = "<<NewComponent.GetUm_l()<<std::endl;
        std::cout<<"Um_g = "<<NewComponent.GetUm_g()<<std::endl;
        std::cout<<"Gm_l = "<<NewComponent.GetGm_l()<<std::endl;
        std::cout<<"Gm_g = "<<NewComponent.GetGm_g()<<std::endl;
        std::cout<<"Quality = "<<NewComponent.GetQuality()<<std::endl;
    }
}