#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
#include<map>
#include<unordered_set>
#include<unordered_map>
#include<map>
#include<set>
#include<random>
#include<algorithm>
#include<queue>
#include<stack>
#include<deque>
#include<bitset>
#include<cstdio>
#include<cassert>
#include<sstream>
#include<set>

#define int long long int

using namespace std;

class Rational {
public:
    
    int a, b;
    
    Rational() : a(0), b(1) {}
    
    Rational(int _a, int _b = 1) : a(_a), b(_b) {
        int g = gcd(abs(a), abs(b));
        a /= g; b /= g;
        if (b < 0) a *= -1;
        b = abs(b);
    }
    
    
    Rational operator- () const {
        return Rational(-a, b);
    }
    
    Rational operator+ (const Rational& other) const {
        int na = a * other.b + other.a * b;
        int nb = b * other.b;
        int g = gcd(abs(na), abs(nb));
        na /= g; nb /= g;
        return Rational(na, nb);
    }
    
    Rational operator- (const Rational& other) const {
        return *this + (-other);
    }
    
    Rational operator* (const Rational& other) const {
        int na = a * other.a;
        int nb = b * other.b;
        int g = gcd(abs(na), abs(nb));
        na /= g; nb /= g;
        return Rational(na, nb);
    }
    
    Rational operator/ (const Rational& other) const {
        int na = a * other.b;
        int nb = b * other.a;
        int g = gcd(abs(na), abs(nb));
        na /= g; nb /= g;
        return Rational(na, nb);
    }
    
    void operator+= (const Rational& other) {
        *this = *this + other;
    }
    
    void operator-= (const Rational& other) {
        *this = *this - other;
    }
    
    void operator*= (const Rational& other) {
        *this = *this * other;
    }
    
    void operator/= (const Rational& other) {
        *this = *this / other;
    }
    
    string to_string() const {
        if (b == 1) return std::to_string(a);
        else return std::to_string(a) + "/" + std::to_string(b);
    }
    
    friend ostream& operator<< (ostream& out, const Rational& a) {
        return out << a.to_string();
    }
    
    friend istream& operator>> (istream& in, Rational& a) {
        int f, s = 1; char c = ' ';
        in >> f >> noskipws >> c >> skipws;
        if (c == '/') in >> s;
        a = Rational(f, s);
        return in;
    }
    
    bool operator >= (const Rational& other) const {
        return a * other.b >= other.a * b;
    }
    
    bool operator > (const Rational& other) const {
        return a * other.b > other.a * b;
    }
    
    bool operator <= (const Rational& other) const {
        return a * other.b <= other.a * b;
    }
    
    bool operator < (const Rational& other) const {
        return a * other.b < other.a * b;
    }
    
    bool operator== (const Rational& other) const {
        return a == other.a && b == other.b;
    }
    
    bool operator!= (const Rational& other) const {
        return !(*this == other);
    }
    
};

Rational fabs(const Rational& a) {
    return Rational(abs(a.a), abs(a.b));
}

const string mask = "xyzwqhrt";

bool isDigit(char c) {
    return '0' <= c && c <= '9';
}

int getInd(char c) {
    return find(mask.begin(), mask.end(), c) - mask.begin();
}

int readNum(const string& s, int& i) {
    int ans = 0, mn = 1;
    if (s[i] == '-') {
        mn = -1;
        ++i;
    }
    while (i < s.size() && isDigit(s[i])) {
        ans = ans * 10 + (s[i] - '0');
        ++i;
    }
    return ans * mn;
}

class Monom {
public:
    
    Rational a;
    vector<int> x;
    
    Monom (Rational a = 1, vector<int> x = {}) : a(a), x(x) {}
    
    Monom (const Monom& other) : a(other.a), x(other.x) {}
    
    bool operator== (const Monom& other) const {
        int len = max(size(), other.size());
        for (int i = 0; i < len; ++i) {
            if (at(i) != other.at(i)) return 0;
        }
        return 1;
    }
    
    int getSign() const {
        if (a >= 0) return 1;
        else return -1;
    }
    
    bool zero() const {
        return a == 0;
    }
    
    Monom operator+ (const Monom& other) const {
        return Monom(a + other.a, x);
    }
    
    Monom operator+= (const Monom& other) {
        *this = *this + other;
        return *this;
    }
    
    Monom operator- (const Monom& other) const {
        return Monom(a - other.a, x);
    }
    
    Monom operator-= (const Monom& other) {
        *this = *this - other;
        return *this;
    }
    
    Monom operator* (const Monom& other) const {
        Monom newM(*this);
        newM.a *= other.a;
        int len = max(newM.size(), other.size());
        newM.x.resize(len);
        for (int i = 0; i < other.size(); ++i) newM[i] += other[i];
        return newM;
    }
    
    Monom operator*= (const Monom& other) {
        *this = *this * other;
        return *this;
    }
    
    Monom operator/ (const Monom& other) const {
        Monom newM(*this);
        newM.a /= other.a;
        int len = max(newM.size(), other.size());
        newM.x.resize(len);
        for (int i = 0; i < other.size(); ++i) newM[i] -= other[i];
        return newM;
    }
    
    Monom operator/= (const Monom& other) {
        *this = *this / other;
        return *this;
    }
    
    int at(int ind) const {
        if (ind >= size()) return 0;
        else return x[ind];
    }
    
    int& operator[] (int ind) {
        return x[ind];
    }
    
    int operator[] (int ind) const {
        return x[ind];
    }
    
    int size() const {
        return x.size();
    }
    
    void clear() {
        a = 1;
        x.clear();
    }
    
    friend ostream& operator<< (ostream& out, const Monom& a) {
        if (a.getSign() == -1) out << " - ";
        if (a.x.size() > 0) {
            if (fabs(a.a) != 1) out << fabs(a.a);
            for (int i = 0; i < a.size(); ++i) {
                if (a[i] == 0) continue;
                out << mask[i];
                if (a[i] > 1) {
                    if (a[i] >= 10) out << "^{" << a[i] << "}";
                    else out << "^" << a[i];
                }
            }
        } else out << fabs(a.a);
        return out;
    }
    
    friend istream& operator>> (istream& in, Monom& a) {
        a.clear();
        string s; in >> s;
        int ptr = 0;
        if (s[ptr] == '-' || isDigit(s[ptr])) {
            a.a = readNum(s, ptr);
            if (a.a == 0) a.a = -1;
        }
        while (ptr < s.size()) {
            int ind = getInd(s[ptr++]);
            int deg = 1;
            if (ptr < s.size() && s[ptr] == '^') {
                ++ptr;
                deg = readNum(s, ptr);
            }
            a.x.resize(max(a.size(), ind + 1));
            a[ind] = deg;
        }
        return in;
    }
};

bool operator< (const Monom& a, const Monom& b) {
    int len = max(a.size(), b.size());
    for (int i = 0; i < len; ++i) {
        if (a.at(i) == b.at(i)) continue;
        if (a.at(i) < b.at(i)) return 1;
        else return 0;
    }
    return 0;
}

class Polynom {
private:
    vector<Monom> a;
    
public:
    Polynom(vector<Monom> a = {}) : a(a) {
        fix();
    }
    
    Polynom(const Polynom& other) : a(other.a) {}
    
    Polynom operator+ (const Monom& q) const {
        Polynom newP(*this);
        for (int i = 0; i < newP.size(); ++i) {
            if (a[i] == q) {
                newP[i] += q;
                newP.fix();
                return newP;
            }
        }
        newP.a.push_back(q);
        newP.fix();
        return newP;
    }
    
    Polynom operator+= (const Monom& q) {
        *this = *this + q;
        return *this;
    }
    
    Polynom operator- (const Monom& q) const {
        Polynom newP(*this);
        for (int i = 0; i < newP.size(); ++i) {
            if (a[i] == q) {
                newP[i] -= q;
                newP.fix();
                return newP;
            }
        }
        newP.a.push_back(q * Rational(-1));
        newP.fix();
        return newP;
    }
    
    Polynom operator-= (const Monom& q) {
        *this = *this - q;
        return *this;
    }
    
    Polynom operator* (const Monom& q) const {
        Polynom newP(*this);
        for (int i = 0; i < newP.size(); ++i) newP[i] *= q;
        newP.fix();
        return newP;
    }
    
    Polynom operator*= (const Monom& q) {
        *this = *this * q;
        return *this;
    }
    
    Polynom operator+ (const Polynom& other) const {
        Polynom newP(*this);
        for (int i = 0; i < other.size(); ++i) newP += other[i];
        newP.fix();
        return newP;
    }
    
    Polynom operator+= (const Polynom& other) {
        *this = *this + other;
        return *this;
    }
    
    Polynom operator- (const Polynom& other) const {
        Polynom newP(*this);
        for (int i = 0; i < other.size(); ++i) newP -= other[i];
        newP.fix();
        return newP;
    }
    
    Polynom operator-= (const Polynom& other) {
        *this = *this - other;
        return *this;
    }
    
    Polynom operator* (const Polynom& other) const {
        Polynom newP;
        for (int i = 0; i < other.size(); ++i) {
            for (int j = 0; j < size(); ++j) newP += other[i] * a[j];
        }
        newP.fix();
        return newP;
    }
    
    Polynom operator*= (const Polynom& other) {
        *this = *this * other;
        return *this;
    }
    
    Monom& operator[] (int ind) {
        return a[ind];
    }
    
    Monom operator[] (int ind) const {
        return a[ind];
    }
    
    Monom getL() const {
        return a[0];
    }
    
    void fix() {
        for (int i = size() - 1; i >= 0; --i) {
            if (a[i].zero()) a.erase(a.begin() + i);
        }
        sort(a.rbegin(), a.rend(), [](const Monom& x, const Monom& y) {
            return x < y;
        });
    }
    
    int size() const {
        return a.size();
    }
    
    bool empty() const {
        return size() == 0;
    }
    
    friend istream& operator>> (istream& in, Polynom& p) {
        string line;
        while (line == "" || line == "\n") getline(in, line);
        stringstream ss(line);
        Monom q;
        int mn = 1;
        while (ss >> q) {
            p.a.push_back(q * Rational(mn));
            char c; ss >> c;
            if (c == '-') mn = -1;
            else mn = 1;
        }
        p.fix();
        return in;
    }
    
    friend ostream& operator<< (ostream& out, const Polynom& p) {
        if (p.empty()) return out << 0;
        for (int i = 0; i < p.size(); ++i) {
            out << p[i];
            if (i < p.size() - 1) {
                if (p[i + 1].getSign() == 1) out << " + ";
            }
        }
        return out;
    }
};

Monom lcm(const Monom& a, const Monom& b) {
    int len = max(a.size(), b.size());
    Monom ans(a.a * b.a, vector<int>(len));
    for (int i = 0; i < len; ++i) ans[i] = max(a.at(i), b.at(i));
    return ans;
}

Polynom S(const Polynom& a, const Polynom& b) {
    auto La = a.getL(), Lb = b.getL();
    auto m = lcm(La, Lb);
    auto m1 = m / La;
    auto m2 = m / Lb;
    return a * m1 - b * m2;
}

bool check(const Monom& a, const Monom& b) {
    int len = max(a.size(), b.size());
    for (int i = 0; i < len; ++i) {
        if (a.at(i) < b.at(i)) return 0;
    }
    return 1;
}

bool oneReduction(Polynom& a, const Polynom& b, int ind = 0) {
    for (int i = 0; i < a.size(); ++i) {
        if (check(a[i], b.getL())) {
            auto m = a[i] / b.getL();
            a -= b * m;
            if (ind > 0) cout << "\\xrightarrow[" << m << "]{f_" + to_string(ind) + "} " << a << " ";
            return 1;
        }
    }
    return 0;
}

bool reduction(Polynom& a, const Polynom& b, int ind = 0) {
    int cnt = 0;
    while (oneReduction(a, b, ind)) ++cnt;
    return cnt > 0;
}

void reduction(Polynom& a, vector<Polynom> F) {
    cout << a << " ";
    while (true) {
        bool cont = 0;
        for (int i = 0; i < F.size(); ++i) {
            if (reduction(a, F[i], i + 1)) {
                cont = 1;
                break;
            }
        }
        if (cont) continue;
        else break;
    }
    cout << "$";
}

void findBasis(vector<Polynom>& a) {
    int i = 0;
    while (i < a.size()) {
        int j = i + 1;
        while (j < a.size()) {
            Polynom s = S(a[i], a[j]);
            cout << "\\vspace{\\baselineskip}\n$S(f_" + to_string(i + 1) + ", f_" + to_string(j + 1) + ") = ";
            reduction(s, a);
            if (!s.empty()) {
                a.push_back(s);
                cout << "$\\Rightarrow f = \\{";
                for (int w = 0; w < a.size(); ++w) {
                    cout << a[w];
                    if (w + 1 < a.size()) cout << ", ";
                }
                cout << "\\}$";
            }
            cout << "\n\n";
            ++j;
        }
        ++i;
    }
}

Polynom Pow(const Polynom& a, int k, const Polynom& mod) {
    Polynom ans = a;
    for (int i = 1; i < k; ++i) ans *= a;
    reduction(ans, mod);
    return ans;
}

signed main() {
    Polynom mod, a; cin >> mod >> a;
    cout << mod << " | " << a << "\n\n";
    for (int i = 1; i <= 8; ++i) {
        cout << Pow(a, i, mod);
        if (i < 8) cout << " \\to ";
        else cout << "\n";
    }
    return 0;
}
