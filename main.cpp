#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

namespace alg1 {
    double Variant(double r, double k);
    double func(double x, double y, double z);
    double Rnk(double x, double y);
    double Qnk(double x, double y);
    double Qqn(double x, double y, double z);
    double U(double x);
    double T(double x);
    double RText(double x, double y, double z, std::string text);
    double CText(double max, std::string text);
    double Max(double x, double y, double z, double u);
    double GText(std::string text);
    double Rrr(double f, double r);
    double Trr(double f, double r);
    double Yrr(double f, double r);
    double Y(double x);
}

namespace alg2 {
    double Rnk(double x, double y);
    double Qnk1(double x, double y);
    double Qqn1(double x, double y, double z);
    double U1(double x);
    double T1(double x);
}

namespace alg3 {
    double funk(double x, double y, double z);
    double Qnk2(double x, double y);
    double Qqn2(double x, double y, double z);
    double U1(double x);
    double T1(double x);
}

//implementing
namespace alg1 {
    double Variant(double r, double k) {
        cout << "Called alg1::Variant(" << r << "," << k << ")" << endl;
        return 0.8973*r+0.1027*k;
    }
    double func(double x, double y, double z) {
        cout << "Called alg1::func(" << x << "," << y << ", " << z << ")" << endl;
        try {
        return alg1::Rnk(x, y) + alg1::Rnk(y, z)*alg1::Rnk(x, y);
        } catch(int err) {
            if(err == 64){
                return alg3::funk(x, y, z);
            }
        }
    }
    double Rnk(double x, double y) {
        cout << "Called alg1::Rnk(" << x << "," << y << ")" << endl;
        try {
            return x*alg1::Qnk(x, y)+y*alg1::Qnk(y,x);
        } catch(int err) {
            if(err == 2) { //calculate Rnk by alg2
                return alg2::Rnk(x, y);
            } else throw err;
        }
    }
    double Qnk(double x, double y) {
        cout << "Called alg1::Qnk(" << x << "," << y << ")" << endl;
        return alg1::Qqn(x, y, x+1) - alg1::Qqn(y, x, x-1);
    }
    double Qqn(double x, double y, double z) {
        cout << "Called alg1::Qqn(" << x << "," << y << ", " << z << ")" << endl;
        return x/alg1::U(x) + y * alg1::T(y) - alg1::U(z) * alg1::T(z);
    }
    double U(double x) {
        cout << "Called alg1::U(" << x << ")" << endl;
        ifstream fin;
        fin.open("/home/ob3r0n/Desktop/dat1.dat");
        if(!fin.is_open()) { //6.1
            throw 2;
        }

        if(abs(x) <= 5) { //6.2
            throw 2;
        }

        map<double, double> table;
        while (!fin.eof()) {
            double x_, u;
            fin >> x_ >> u;
            table.insert(pair<double, double>(x_, u));
        }
        fin.close();
        vector<double> x_values;
        for (auto i : table)
            if (x == i.first)
                return i.second;
            else
                x_values.push_back(i.first);

        double x_i = 0;
        for (unsigned int i = 0; i < x_values.size() - 1; i++) //6.3
            if (x_values[i] < x && x < x_values[i + 1]) {
                x_i = x_values[i];
                break;
            }
        return alg1::U(x_i);
    }
    double T(double x) {
        cout << "Called alg1::T(" << x << ")" << endl;
        ifstream fin;
        fin.open("/home/ob3r0n/Desktop/dat2.dat");
        if(!fin.is_open()) { //6.4
            throw 64;
        }

        if(abs(x) <= 5) { //6.5
            throw 2;
        }

        map<double, double> table;
        while (!fin.eof()) {
            double x_, u;
            fin >> x_ >> u;
            table.insert(pair<double, double>(x_, u));
        }
        fin.close();
        vector<double> x_values;
        for (auto i : table)
            if (x == i.first)
                return i.second;
            else
                x_values.push_back(i.first);

        double x_i = 0;
        for (unsigned int i = 0; i < x_values.size() - 1; i++) //6.6
            if (x_values[i] < x && x < x_values[i + 1]) {
                x_i = x_values[i];
                break;
            }
        return alg1::T(x_i);
    }
    double RText(double x, double y, double z, std::string text) {
        cout << "Called alg1::Rnk(" << x << "," << y << ", " << z << ", " << text << ")" << endl;
        return alg1::CText(alg1::Max(x, y, x+z, y+z), text);
    }
    double CText(double x, std::string text) {
        cout << "Called alg1::Rnk(" << x << "," << text << ")" << endl;
        if(x > 0) return alg1::GText(text) + x;
        else if(text == "") return alg1::GText("set")+alg1::GText("get")-x;
        else if (x <= 0) return alg1::GText("set")+alg1::GText(text);
    }

    double Max(double x, double y, double z, double u) {
        cout << "Called alg1::Max(" << x << "," << y << ", " << z << ", " << u << ")" << endl;
        return max(max(x, y), max(z, u));
    }
    double GText(std::string text) {
        cout << "Called alg1::GText(" << text << ")" << endl;
        ifstream fin;
        fin.open("/home/ob3r0n/Desktop/dat3.dat");
        if (!fin.is_open()) { //10.1
            cout << "Can't open dat3.dat";
            exit(0);
        }
        map<string, double> table;
        while (!fin.eof()) {
            string x;
            double u;
            fin >> x >> u;
            table.insert(pair<string, double>(x, u));
        }
        fin.close();
        for (auto i = table.begin(); i != table.end(); ++i) {
            if (text == (*i).first) //10.2
                return (*i).second; //10.3
        }
        return 0; //10.4

    }
    double Rrr(double f, double r) {
        cout << "Called alg1::Rrr(" << f << "," << r << ")" << endl;
        return f*alg1::Trr(f, r)+r*alg1::Trr(f, 2*r);
    }
    double Trr(double f, double r) {
        cout << "Called alg1::Trr(" << f << "," << r << ")" << endl;
        if(4*f*f - r < 0) throw 121;
        return sqrt(4*f*f-r) + 0.5*alg1::Yrr(r, f);
    }
    double Yrr(double f, double r) {
        cout << "Called alg1::Yrr(" << f << "," << r << ")" << endl;
        return alg1::Y(f)*r+0.5*alg1::Y(r);
    }
    double Y(double x) {
        cout << "Called alg1::Y(" << x << ")" << endl;
        if(100-x*x < 0) throw 121; //14.1
        if(x*sqrt(100-x*x) < 1) throw 121; //14.2
        return log(x*sqrt(100-x*x));
    }
}

namespace alg2 {
    double Rnk(double x, double y) {
        cout << "Called alg2::Rnk(" << x << "," << y << ")" << endl;
        return x*alg2::Qnk1(x, y) + y*alg2::Qnk1(y, x) - 0.03 * alg2::Qnk1(x,y) * alg2::Qnk1(y, x);
    }
    double Qnk1(double x, double y) {
        cout << "Called alg2::Qnk1(" << x << "," << y << ")" << endl;
        return 1.1 * alg2::Qqn1(x, y, x+y) - 0.9 * alg2::Qqn1(y, x, x-y);
    }
    double Qqn1(double x, double y, double z) {
        cout << "Called alg2::Qqn1(" << x << "," << y << ", " << z << ")" << endl;
        return x/alg2::U1(x) + y*alg2::T1(y)-alg2::U1(z)*alg2::T1(z);
    }
    double arccot(double x) {
        const double pi = 3.14159265358979323846;
        return pi/2 - atan(x);
    };
    double U1(double x) {
        cout << "Called alg2::U1(" << x << ")" << endl;
        return arccot(asin(sin(3*x)));
    }
    double T1(double x) {
        cout << "Called alg2::T1(" << x << ")" << endl;
        return arccot(acos(sin(2*x)));
    }
}

namespace alg3 {
    double funk(double x, double y, double z) {
        cout << "Called alg3::funk(" << x << "," << y << ", " << z << ")" << endl;
        return 1.75*x*alg3::Qnk2(x, y) + 1.25*y*alg3::Qnk2(y, x) - 1.5*alg3::Qnk2(x, y)*alg3::Qnk2(y, x);
    }
    double Qnk2(double x, double y) {
        cout << "Called alg3::Qnk2(" << x << "," << y << ")" << endl;
        return 1.3 * alg3::Qqn2(x, y, x) - 0.7*alg3::Qqn2(y, x, x);
    }
    double Qqn2(double x, double y, double z) {
        cout << "Called alg3::Qqn2(" << x << "," << y << ", " << z << ")" << endl;
        return x/alg3::U1(x)+y*alg3::T1(x)-0.9*alg3::U1(z)*alg3::T1(z);
    }
    double arccot(double x) {
        const double pi = 3.14159265358979323846;
        return pi/2 - atan(x);
    };
    double U1(double x) {
        cout << "Called alg3::U1(" << x << ")" << endl;
        return arccot(asin(sin(3*x)));
    }
    double T1(double x) {
        cout << "Called alg3::T1(" << x << ")" << endl;
        return arccot(acos(sin(2*x)));
    }
}

int main() {
    int x, y, z;
    cout << "Input x, y, z: ";
    cin >> x >> y >> z;
    string text;
    cout << "Input text: ";
    cin >> text;
    double r = alg1::func(x, y, z), k;
    try {
        k = alg1::RText(x, y, z, text);
    } catch(int err) {
        if(err == 121) {
            k = 0; //12.1
        }
    }
    double result = alg1::Variant(r, k);
    std::cout << "Result is: " << result << std::endl;
    return 0;
}
